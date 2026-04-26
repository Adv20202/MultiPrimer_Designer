"""
Specificity checker for primers using sequence alignment.
Uses BWA or other aligners to check for off-target binding sites.
"""

import os
import subprocess
import tempfile
import shutil
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass
from pathlib import Path

from ..core.models import Primer, PrimerPair, GenomicPosition, Strand
from ..utils.logger import LoggerMixin
from ..utils.config import get_config


@dataclass
class AlignmentHit:
    """Represents an alignment hit for a primer."""
    chromosome: str
    position: int
    strand: str
    mismatches: int
    mapping_quality: int


@dataclass
class SpecificityResult:
    """Result of specificity check for a primer pair."""
    is_specific: bool
    forward_hits: List[AlignmentHit]
    reverse_hits: List[AlignmentHit]
    potential_amplicons: List[Dict[str, Any]]
    warnings: List[str]
    error: str = ""


class SpecificityChecker(LoggerMixin):
    """
    Checks primer specificity by mapping primers to the reference genome.
    Uses BWA-MEM for short-read alignment.
    """

    # Maximum acceptable amplicon size for off-target check
    MAX_AMPLICON_SIZE = 20000

    # Maximum mismatches for primer binding
    MAX_MISMATCHES = 3

    def __init__(
        self,
        bwa_path: str = None,
        genome_index_path: str = None,
        assembly: str = "GRCh38"
    ):
        """
        Initialize specificity checker.

        Args:
            bwa_path: Path to BWA executable
            genome_index_path: Path to BWA genome index
            assembly: Reference assembly
        """
        config = get_config()
        self.bwa_path = bwa_path or config.bwa.bwa_path
        self.genome_index_path = genome_index_path or config.bwa.genome_index_path
        self.assembly = assembly
        self.threads = config.bwa.threads

        # Check if BWA is available
        self._bwa_available = self._check_bwa_available()

        # Temporary directory for intermediate files
        self._temp_dir = None

    def _check_bwa_available(self) -> bool:
        """Check if BWA is installed and accessible."""
        try:
            result = subprocess.run(
                [self.bwa_path, 'version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                self.log_info(f"BWA found: {result.stdout.strip()}")
                return True
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

        self.log_warning("BWA not found. Specificity checking will be limited.")
        return False

    def is_available(self) -> bool:
        """Check if specificity checking is available."""
        return self._bwa_available and bool(self.genome_index_path)

    def check_specificity(
        self,
        primer_pair: PrimerPair,
        expected_chromosome: str = None,
        expected_start: int = None,
        expected_end: int = None,
        progress_callback=None
    ) -> SpecificityResult:
        """
        Check specificity of a primer pair.

        Args:
            primer_pair: PrimerPair to check
            expected_chromosome: Expected target chromosome
            expected_start: Expected amplicon start
            expected_end: Expected amplicon end
            progress_callback: Optional progress callback

        Returns:
            SpecificityResult with specificity information
        """
        if not self.is_available():
            return self._fallback_specificity_check(primer_pair)

        try:
            # Create temp directory
            self._temp_dir = tempfile.mkdtemp(prefix='primer_spec_')

            # Map forward primer
            if progress_callback:
                progress_callback(0.2, "Mapping forward primer...")
            forward_hits = self._map_primer(
                primer_pair.forward.sequence,
                'forward'
            )

            # Map reverse primer
            if progress_callback:
                progress_callback(0.4, "Mapping reverse primer...")
            reverse_hits = self._map_primer(
                primer_pair.reverse.sequence,
                'reverse'
            )

            # Find potential amplicons
            if progress_callback:
                progress_callback(0.6, "Analyzing potential amplicons...")
            potential_amplicons = self._find_potential_amplicons(
                forward_hits, reverse_hits
            )

            # Analyze specificity
            if progress_callback:
                progress_callback(0.8, "Evaluating specificity...")
            is_specific, warnings = self._evaluate_specificity(
                potential_amplicons,
                expected_chromosome,
                expected_start,
                expected_end
            )

            return SpecificityResult(
                is_specific=is_specific,
                forward_hits=forward_hits,
                reverse_hits=reverse_hits,
                potential_amplicons=potential_amplicons,
                warnings=warnings
            )

        except Exception as e:
            self.log_exception(f"Specificity check error: {e}")
            return SpecificityResult(
                is_specific=False,
                forward_hits=[],
                reverse_hits=[],
                potential_amplicons=[],
                warnings=[],
                error=str(e)
            )

        finally:
            # Cleanup temp directory
            if self._temp_dir and os.path.exists(self._temp_dir):
                shutil.rmtree(self._temp_dir, ignore_errors=True)

    def _map_primer(self, sequence: str, name: str) -> List[AlignmentHit]:
        """
        Map a primer sequence to the genome using BWA.

        Args:
            sequence: Primer sequence
            name: Primer name for temporary files

        Returns:
            List of alignment hits
        """
        if not self._temp_dir:
            return []

        # Create FASTA file for primer
        fasta_path = os.path.join(self._temp_dir, f'{name}.fa')
        with open(fasta_path, 'w') as f:
            f.write(f'>{name}\n{sequence}\n')

        # Run BWA MEM
        sam_path = os.path.join(self._temp_dir, f'{name}.sam')

        cmd = [
            self.bwa_path, 'mem',
            '-t', str(self.threads),
            '-k', '10',  # Minimum seed length
            '-T', '20',  # Minimum alignment score
            '-a',        # Output all alignments
            self.genome_index_path,
            fasta_path
        ]

        try:
            with open(sam_path, 'w') as sam_file:
                result = subprocess.run(
                    cmd,
                    stdout=sam_file,
                    stderr=subprocess.PIPE,
                    timeout=60
                )

            if result.returncode != 0:
                self.log_warning(f"BWA error: {result.stderr.decode()}")
                return []

            # Parse SAM file
            return self._parse_sam(sam_path)

        except subprocess.TimeoutExpired:
            self.log_warning("BWA timed out")
            return []
        except Exception as e:
            self.log_error(f"BWA execution error: {e}")
            return []

    def _parse_sam(self, sam_path: str) -> List[AlignmentHit]:
        """Parse SAM file to extract alignment hits."""
        hits = []

        with open(sam_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                flag = int(fields[1])
                if flag & 4:  # Unmapped
                    continue

                chromosome = fields[2]
                position = int(fields[3])
                mapq = int(fields[4])

                # Determine strand
                strand = '-' if flag & 16 else '+'

                # Count mismatches from NM tag
                mismatches = 0
                for field in fields[11:]:
                    if field.startswith('NM:i:'):
                        mismatches = int(field[5:])
                        break

                # Filter by mismatches
                if mismatches <= self.MAX_MISMATCHES:
                    hits.append(AlignmentHit(
                        chromosome=chromosome,
                        position=position,
                        strand=strand,
                        mismatches=mismatches,
                        mapping_quality=mapq
                    ))

        return hits

    def _find_potential_amplicons(
        self,
        forward_hits: List[AlignmentHit],
        reverse_hits: List[AlignmentHit]
    ) -> List[Dict[str, Any]]:
        """
        Find potential amplicons from primer mapping hits.

        Args:
            forward_hits: Forward primer alignments
            reverse_hits: Reverse primer alignments

        Returns:
            List of potential amplicon information
        """
        amplicons = []

        for fwd in forward_hits:
            for rev in reverse_hits:
                # Must be on same chromosome
                if fwd.chromosome != rev.chromosome:
                    continue

                # Check orientation (forward on + strand, reverse on - strand)
                # OR (forward on - strand, reverse on + strand)
                valid_orientation = (
                    (fwd.strand == '+' and rev.strand == '-' and rev.position > fwd.position) or
                    (fwd.strand == '-' and rev.strand == '+' and fwd.position > rev.position)
                )

                if not valid_orientation:
                    continue

                # Calculate amplicon size
                if fwd.strand == '+':
                    size = rev.position - fwd.position
                else:
                    size = fwd.position - rev.position

                # Check size limit
                if 0 < size <= self.MAX_AMPLICON_SIZE:
                    amplicons.append({
                        'chromosome': fwd.chromosome,
                        'start': min(fwd.position, rev.position),
                        'end': max(fwd.position, rev.position),
                        'size': size,
                        'forward_mismatches': fwd.mismatches,
                        'reverse_mismatches': rev.mismatches,
                        'total_mismatches': fwd.mismatches + rev.mismatches
                    })

        # Sort by total mismatches
        amplicons.sort(key=lambda x: x['total_mismatches'])

        return amplicons

    def _evaluate_specificity(
        self,
        potential_amplicons: List[Dict[str, Any]],
        expected_chromosome: str,
        expected_start: int,
        expected_end: int
    ) -> Tuple[bool, List[str]]:
        """
        Evaluate specificity based on potential amplicons.

        Args:
            potential_amplicons: List of potential amplicons
            expected_chromosome: Expected chromosome
            expected_start: Expected start position
            expected_end: Expected end position

        Returns:
            Tuple of (is_specific, warnings)
        """
        warnings = []

        if not potential_amplicons:
            warnings.append("No alignments found - primers may be too specific or have design issues")
            return False, warnings

        # Find the expected amplicon
        expected_found = False
        off_targets = []

        for amp in potential_amplicons:
            # Check if this matches expected target
            if expected_chromosome and expected_start and expected_end:
                is_expected = (
                    amp['chromosome'] == expected_chromosome and
                    abs(amp['start'] - expected_start) < 100 and
                    abs(amp['end'] - expected_end) < 100
                )
            else:
                # If no expected position, first perfect match is assumed to be target
                is_expected = (
                    amp['total_mismatches'] == 0 and not expected_found
                )

            if is_expected:
                expected_found = True
            else:
                off_targets.append(amp)

        if not expected_found:
            warnings.append("Expected target amplicon not found among alignments")

        # Check off-targets
        perfect_off_targets = [a for a in off_targets if a['total_mismatches'] == 0]
        near_off_targets = [a for a in off_targets if 0 < a['total_mismatches'] <= 2]

        if perfect_off_targets:
            warnings.append(
                f"Found {len(perfect_off_targets)} perfect off-target site(s): " +
                ", ".join(f"{a['chromosome']}:{a['start']}" for a in perfect_off_targets[:3])
            )

            # Check for pseudogenes (same chromosome, different location)
            for amp in perfect_off_targets:
                if amp['chromosome'] == expected_chromosome:
                    warnings.append(
                        f"Potential pseudogene or duplicated region on {amp['chromosome']}:{amp['start']}-{amp['end']}"
                    )

        if near_off_targets:
            warnings.append(
                f"Found {len(near_off_targets)} near-match off-target site(s) with 1-2 mismatches"
            )

        # Determine specificity
        is_specific = (
            expected_found and
            len(perfect_off_targets) == 0
        )

        return is_specific, warnings

    def _fallback_specificity_check(self, primer_pair: PrimerPair) -> SpecificityResult:
        """
        Perform basic specificity checks when BWA is not available.

        Args:
            primer_pair: PrimerPair to check

        Returns:
            SpecificityResult with limited information
        """
        warnings = [
            "Full specificity check not available (BWA not configured). "
            "Performing basic sequence analysis only."
        ]

        # Check for poly-X runs
        for primer in [primer_pair.forward, primer_pair.reverse]:
            seq = primer.sequence
            for base in 'ACGT':
                if base * 5 in seq:
                    warnings.append(
                        f"{'Forward' if primer.is_forward else 'Reverse'} primer "
                        f"contains poly-{base} run which may reduce specificity"
                    )

        # Check for low complexity
        for primer in [primer_pair.forward, primer_pair.reverse]:
            seq = primer.sequence
            unique_bases = len(set(seq))
            if unique_bases < 3:
                warnings.append(
                    f"{'Forward' if primer.is_forward else 'Reverse'} primer "
                    "has low sequence complexity"
                )

        return SpecificityResult(
            is_specific=True,  # Assume specific without full check
            forward_hits=[],
            reverse_hits=[],
            potential_amplicons=[],
            warnings=warnings
        )

    def check_primer_sequence(self, sequence: str) -> Dict[str, Any]:
        """
        Perform basic sequence quality checks on a primer.

        Args:
            sequence: Primer sequence

        Returns:
            Dictionary with check results
        """
        seq_len = len(sequence)
        results = {
            'length': seq_len,
            'gc_content': (sequence.count('G') + sequence.count('C')) / seq_len * 100 if seq_len > 0 else 0,
            'issues': []
        }

        # Check length
        if len(sequence) < 15:
            results['issues'].append("Primer is very short (<15 bp)")
        elif len(sequence) > 30:
            results['issues'].append("Primer is very long (>30 bp)")

        # Check GC content
        if results['gc_content'] < 30:
            results['issues'].append("Low GC content (<30%)")
        elif results['gc_content'] > 70:
            results['issues'].append("High GC content (>70%)")

        # Check for poly-X runs
        for base in 'ACGT':
            if base * 4 in sequence:
                results['issues'].append(f"Contains poly-{base} run")

        # Check 3' end stability
        last_5 = sequence[-5:]
        gc_3prime = (last_5.count('G') + last_5.count('C')) / 5 * 100
        if gc_3prime > 80:
            results['issues'].append("3' end may be too GC-rich")
        elif gc_3prime < 20:
            results['issues'].append("3' end may be too AT-rich")

        # Check for repeats
        for repeat_len in range(3, 6):
            for i in range(len(sequence) - repeat_len * 2 + 1):
                repeat = sequence[i:i + repeat_len]
                if repeat * 2 in sequence:
                    results['issues'].append(f"Contains repeat: {repeat}")
                    break

        results['is_acceptable'] = len(results['issues']) == 0

        return results

    def setup_genome_index(self, genome_fasta: str, output_dir: str) -> bool:
        """
        Set up BWA genome index.

        Args:
            genome_fasta: Path to genome FASTA file
            output_dir: Directory for index files

        Returns:
            True if successful
        """
        if not self._bwa_available:
            self.log_error("BWA not available for index creation")
            return False

        try:
            os.makedirs(output_dir, exist_ok=True)
            index_base = os.path.join(output_dir, 'genome')

            # Copy FASTA to index directory
            import shutil
            shutil.copy(genome_fasta, f"{index_base}.fa")

            # Run BWA index
            cmd = [self.bwa_path, 'index', '-p', index_base, f"{index_base}.fa"]

            self.log_info(f"Creating BWA index: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout for indexing
            )

            if result.returncode == 0:
                self.genome_index_path = index_base
                self.log_info(f"BWA index created at {index_base}")
                return True
            else:
                self.log_error(f"BWA index failed: {result.stderr}")
                return False

        except Exception as e:
            self.log_exception(f"Index creation error: {e}")
            return False

"""
Probe designer for qPCR applications.
Designs TaqMan-style probes targeting variant positions.
"""

from typing import List, Optional, Dict, Any
from dataclasses import dataclass

from ..core.models import Variant, Probe, DesignParameters, GenomicPosition
from ..api.ensembl_client import EnsemblClient
from ..utils.logger import LoggerMixin

try:
    import primer3
    PRIMER3_AVAILABLE = True
except ImportError:
    PRIMER3_AVAILABLE = False


@dataclass
class ProbeDesignResult:
    """Result of probe design."""
    success: bool
    probes: List[Probe]
    message: str = ""
    warnings: List[str] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


class ProbeDesigner(LoggerMixin):
    """
    Designs fluorescent probes (TaqMan-style) for qPCR variant detection.
    Positions the probe to span the variant site.
    """

    def __init__(
        self,
        ensembl_client: EnsemblClient = None,
        parameters: DesignParameters = None
    ):
        """
        Initialize probe designer.

        Args:
            ensembl_client: For fetching sequences
            parameters: Design parameters
        """
        if not PRIMER3_AVAILABLE:
            raise ImportError(
                "primer3-py library is required for probe design. "
                "Install it with: pip install primer3-py"
            )

        self.ensembl_client = ensembl_client or EnsemblClient()
        self.parameters = parameters or DesignParameters()

    def design_probe(
        self,
        variant: Variant,
        amplicon_sequence: str = None,
        variant_position_in_amplicon: int = None
    ) -> ProbeDesignResult:
        """
        Design a probe for variant detection.

        Args:
            variant: Variant to design probe for
            amplicon_sequence: Sequence of the amplicon
            variant_position_in_amplicon: Position of variant in amplicon sequence

        Returns:
            ProbeDesignResult with designed probes
        """
        if not variant.genomic_position and not amplicon_sequence:
            return ProbeDesignResult(
                success=False,
                probes=[],
                message="No sequence available for probe design"
            )

        try:
            # Get sequence if not provided
            if not amplicon_sequence:
                sequence, var_pos = self._get_variant_region_sequence(variant)
                if not sequence:
                    return ProbeDesignResult(
                        success=False,
                        probes=[],
                        message="Could not fetch sequence for probe design"
                    )
            else:
                sequence = amplicon_sequence
                var_pos = variant_position_in_amplicon

            # Design probes
            probes = self._design_probes_for_variant(sequence, var_pos, variant)

            if not probes:
                return ProbeDesignResult(
                    success=False,
                    probes=[],
                    message="No suitable probes found",
                    warnings=self._generate_probe_suggestions()
                )

            return ProbeDesignResult(
                success=True,
                probes=probes,
                message=f"Designed {len(probes)} probe(s)"
            )

        except Exception as e:
            self.log_exception(f"Probe design error: {e}")
            return ProbeDesignResult(
                success=False,
                probes=[],
                message=f"Probe design error: {str(e)}"
            )

    def _get_variant_region_sequence(
        self,
        variant: Variant,
        flank_size: int = 100
    ) -> tuple:
        """
        Get sequence around the variant for probe design.

        Args:
            variant: Variant with genomic position
            flank_size: Flanking sequence size

        Returns:
            Tuple of (sequence, variant_position_in_sequence)
        """
        if not variant.genomic_position:
            return None, None

        pos = variant.genomic_position
        start = max(1, pos.start - flank_size)
        end = pos.end + flank_size

        sequence = self.ensembl_client.get_genomic_sequence(
            chromosome=pos.chromosome,
            start=start,
            end=end
        )

        if sequence:
            var_position = pos.start - start
            return sequence, var_position

        return None, None

    def _design_probes_for_variant(
        self,
        sequence: str,
        variant_position: int,
        variant: Variant
    ) -> List[Probe]:
        """
        Design probes that span the variant position.

        Args:
            sequence: DNA sequence
            variant_position: Position of variant in sequence
            variant: Variant information

        Returns:
            List of Probe objects
        """
        probes = []
        p = self.parameters

        # Probe should ideally have the variant in the middle
        # Try different positions with variant in the center
        min_size = p.probe_min_size
        max_size = p.probe_max_size

        for probe_length in range(min_size, max_size + 1, 2):
            # Try centering variant at different positions within the probe
            for var_offset in range(probe_length // 3, 2 * probe_length // 3):
                probe_start = variant_position - var_offset
                probe_end = probe_start + probe_length

                if probe_start < 0 or probe_end > len(sequence):
                    continue

                probe_seq = sequence[probe_start:probe_end]

                # Calculate Tm
                tm = self._calculate_tm(probe_seq)

                if not (p.probe_min_tm <= tm <= p.probe_max_tm):
                    continue

                # Calculate GC content
                gc = (probe_seq.count('G') + probe_seq.count('C')) / len(probe_seq) * 100

                if not (30 <= gc <= 70):
                    continue

                # Check for problems
                if self._has_probe_issues(probe_seq):
                    continue

                # Create probe
                probe = Probe(
                    sequence=probe_seq,
                    start=probe_start,
                    end=probe_end - 1,
                    tm=tm,
                    gc_content=gc,
                    variant_position_in_probe=var_offset
                )
                probes.append(probe)

                # Limit number of probes
                if len(probes) >= 5:
                    break

            if len(probes) >= 5:
                break

        # Sort by Tm closest to optimal
        opt_tm = (p.probe_min_tm + p.probe_max_tm) / 2
        probes.sort(key=lambda x: abs(x.tm - opt_tm))

        return probes[:5]

    def _calculate_tm(self, sequence: str) -> float:
        """Calculate melting temperature for probe."""
        return primer3.calc_tm(
            sequence,
            mv_conc=self.parameters.mv_conc,
            dv_conc=self.parameters.dv_conc,
            dntp_conc=self.parameters.dntp_conc,
            dna_conc=self.parameters.dna_conc
        )

    def _has_probe_issues(self, sequence: str) -> bool:
        """Check if probe sequence has problems."""
        # Check for G at 5' end (quenches fluorescence)
        if sequence.startswith('G'):
            return True

        # Check for runs of same base
        for base in 'ACGT':
            if base * 4 in sequence:
                return True

        # Check for too many Gs (can form G-quadruplexes)
        if len(sequence) > 0 and sequence.count('G') / len(sequence) > 0.35:
            return True

        # Check for palindromes (hairpin formation)
        for length in range(6, min(len(sequence) // 2, 12)):
            for i in range(len(sequence) - length * 2):
                subseq = sequence[i:i + length]
                rev_comp = self._reverse_complement(subseq)
                if rev_comp in sequence[i + length:]:
                    return True

        return False

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(b, b) for b in reversed(sequence))

    def _generate_probe_suggestions(self) -> List[str]:
        """Generate suggestions when no probes found."""
        return [
            "Try relaxing Tm constraints for probe",
            "The variant position may have unfavorable sequence context",
            "Consider using a different variant detection method",
            "Try adjusting probe size parameters"
        ]

    def design_allele_specific_probes(
        self,
        variant: Variant,
        amplicon_sequence: str = None,
        variant_position_in_amplicon: int = None
    ) -> Dict[str, ProbeDesignResult]:
        """
        Design probes for both reference and alternate alleles.

        Args:
            variant: Variant with allele information
            amplicon_sequence: Amplicon sequence
            variant_position_in_amplicon: Position of variant in amplicon (required if amplicon_sequence provided)

        Returns:
            Dictionary with 'ref' and 'alt' ProbeDesignResults
        """
        results = {}

        # Get sequence
        if not amplicon_sequence:
            sequence, var_pos = self._get_variant_region_sequence(variant)
            if not sequence:
                error_result = ProbeDesignResult(
                    success=False,
                    probes=[],
                    message="Could not fetch sequence"
                )
                return {'ref': error_result, 'alt': error_result}
        else:
            sequence = amplicon_sequence
            var_pos = variant_position_in_amplicon

        # Design probe for reference allele (using original sequence)
        results['ref'] = self.design_probe(variant, sequence, var_pos)

        # Create alternate sequence - only if we have valid position and alleles
        if variant.ref_allele and variant.alt_allele and var_pos is not None and var_pos >= 0:
            # Validate position is within sequence bounds
            if var_pos + len(variant.ref_allele) <= len(sequence):
                alt_sequence = (
                    sequence[:var_pos] +
                    variant.alt_allele +
                    sequence[var_pos + len(variant.ref_allele):]
                )
                results['alt'] = self.design_probe(variant, alt_sequence, var_pos)
            else:
                results['alt'] = ProbeDesignResult(
                    success=False,
                    probes=[],
                    message="Variant position exceeds sequence length"
                )
        else:
            results['alt'] = ProbeDesignResult(
                success=False,
                probes=[],
                message="Missing variant position or allele information for alternate sequence"
            )

        return results

    def check_probe_specificity(
        self,
        probe: Probe,
        target_sequence: str
    ) -> Dict[str, Any]:
        """
        Check probe binding specificity.

        Args:
            probe: Probe to check
            target_sequence: Full target sequence

        Returns:
            Dictionary with specificity metrics
        """
        results = {
            'unique_binding': True,
            'binding_sites': [],
            'warnings': []
        }

        probe_seq = probe.sequence

        # Find all binding sites in target
        for i in range(len(target_sequence) - len(probe_seq) + 1):
            site = target_sequence[i:i + len(probe_seq)]
            mismatches = sum(1 for a, b in zip(probe_seq, site) if a != b)

            if mismatches <= 3:  # Allow some mismatches
                results['binding_sites'].append({
                    'position': i,
                    'mismatches': mismatches
                })

        if len(results['binding_sites']) > 1:
            results['unique_binding'] = False
            results['warnings'].append(
                f"Probe has {len(results['binding_sites'])} potential binding sites"
            )

        # Check hairpin potential
        hairpin = primer3.calc_hairpin(
            probe_seq,
            mv_conc=self.parameters.mv_conc,
            dv_conc=self.parameters.dv_conc
        )

        if hairpin.tm > 50:
            results['warnings'].append(
                f"Probe may form hairpin (Tm: {hairpin.tm:.1f}°C)"
            )

        return results

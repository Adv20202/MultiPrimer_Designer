"""
Sequence homology analysis -- maps the variant flanking sequence
to the human reference genome using BLAST+ (blastn) to identify
pseudogenes, duplications and homologous regions.

Uses blastn with -task blastn (NOT megablast) for maximum sensitivity
to divergent homologs (85-100% identity).
"""

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Callable, Tuple

from .blast_setup import BlastSetupManager
from ..utils.logger import LoggerMixin
from ..utils.config import get_config


# ---------------------------------------------------------------------------
#  Data models
# ---------------------------------------------------------------------------

@dataclass
class HomologyHit:
    """Single alignment hit from BLAST+."""
    chromosome: str
    position: int                # Leftmost position on reference (1-based)
    strand: str                  # '+' or '-'
    mapq: int                    # 0 for BLAST (no MAPQ concept)
    cigar: str                   # Synthetic CIGAR built from BLAST alignment
    alignment_score: int         # Bit score (int)
    mismatches: int              # Edit distance / mismatch count
    query_sequence: str          # Query sequence (ungapped)
    md_tag: str                  # Synthetic MD tag built from BLAST alignment
    aligned_length: int          # Alignment length on reference
    percent_identity: float      # From BLAST pident
    is_primary: bool             # True for the on-target hit
    is_supplementary: bool       # Always False for BLAST
    alignment_start: int         # Query start (0-based)
    alignment_end: int           # Query end (0-based, exclusive)
    # BLAST-specific (optional, backward-compatible)
    evalue: float = 0.0
    bit_score: float = 0.0


@dataclass
class HomologyResult:
    """Full homology analysis result for one sequence."""
    query_name: str
    query_length: int
    query_chromosome: str
    query_position: int
    hits: List[HomologyHit] = field(default_factory=list)
    primary_hit: Optional[HomologyHit] = None
    error: str = ""
    blast_stderr: str = ""      # BLAST stderr output
    query_genomic_start: int = 0  # Genomic start of query sequence (for coord translation)


# ---------------------------------------------------------------------------
#  Main analyzer
# ---------------------------------------------------------------------------

class HomologyAnalyzer(LoggerMixin):
    """
    Performs sequence homology analysis by aligning the flanking
    sequence (~500 bp) against the human reference genome using
    BLAST+ (blastn) to find all homologous regions with >85% identity.
    """

    # BLAST output format fields (tab-separated)
    BLAST_OUTFMT = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen "
        "qseq sseq sstrand"
    )
    BLAST_FIELD_COUNT = 17

    def __init__(self):
        self.config = get_config()
        self._blastn_path: str = ""
        self._db_path: str = ""
        self._is_ready = False

    # ------------------------------------------------------------------
    #  Readiness check
    # ------------------------------------------------------------------

    def ensure_blast_ready(
        self, progress_callback: Optional[Callable] = None
    ) -> Tuple[bool, str]:
        """
        Ensure BLAST+ is installed and genome database is built.
        Returns (ready, error_message).
        """
        if self._is_ready:
            return True, ""

        setup = BlastSetupManager()

        # 1. Detect / install blastn
        if progress_callback:
            progress_callback("Checking BLAST+ installation...")

        ok, blastn_path = setup.detect_blastn()
        if not ok:
            if progress_callback:
                progress_callback("Installing BLAST+...")
            ok, result = setup.install_blast(progress_callback)
            if not ok:
                return False, f"BLAST+ not available: {result}"
            blastn_path = result

        # 2. Detect / build BLAST database
        if progress_callback:
            progress_callback("Checking BLAST database...")

        db_ok, db_path = setup.detect_blast_db()
        if not db_ok:
            if progress_callback:
                progress_callback("Building BLAST database (first time only)...")
            db_ok, result = setup.setup_blast_db(progress_callback)
            if not db_ok:
                return False, f"BLAST database not available: {result}"
            db_path = result

        self._blastn_path = blastn_path
        self._db_path = db_path
        self._is_ready = True
        self.log_info(
            f"[HomologyAnalyzer] BLAST ready: {blastn_path}, db: {db_path}"
        )
        return True, ""

    # ------------------------------------------------------------------
    #  Main analysis
    # ------------------------------------------------------------------

    def analyze(
        self,
        sequence: str,
        query_name: str,
        source_chromosome: str,
        source_position: int,
        progress_callback: Optional[Callable] = None,
    ) -> HomologyResult:
        """
        Run homology analysis for the given flanking sequence.

        Args:
            sequence:           Flanking sequence (~501 bp)
            query_name:         Identifier (e.g. "TP53_c.215C>G")
            source_chromosome:  Chromosome of the source variant
            source_position:    Genomic position of the source variant
            progress_callback:  Optional status callback

        Returns:
            HomologyResult with all hits passing the filter criteria.
        """
        if not self._is_ready:
            return HomologyResult(
                query_name=query_name,
                query_length=len(sequence),
                query_chromosome=source_chromosome,
                query_position=source_position,
                error="BLAST is not configured. Call ensure_blast_ready() first.",
            )

        temp_dir = None
        try:
            temp_dir = tempfile.mkdtemp(prefix="homology_")

            # Write query FASTA
            fasta_path = os.path.join(temp_dir, "query.fa")
            with open(fasta_path, "w") as f:
                f.write(f">{query_name}\n{sequence}\n")

            if progress_callback:
                progress_callback("Running BLAST+ homology search...")

            # Run BLAST
            stdout_text, stderr_text = self._run_blast(fasta_path)

            if progress_callback:
                progress_callback("Parsing BLAST results...")

            # Parse results
            hits = self._parse_blast_tabular(
                stdout_text, len(sequence),
                source_chromosome, source_position,
                full_query_sequence=sequence,
            )

            # Identify primary hit
            primary = next((h for h in hits if h.is_primary), None)

            if hits:
                primary_info = (
                    f"{primary.chromosome}:{primary.position}"
                    if primary else "none"
                )
                self.log_info(
                    f"[HomologyAnalyzer] Found {len(hits)} hits "
                    f"(primary: {primary_info})"
                )

            return HomologyResult(
                query_name=query_name,
                query_length=len(sequence),
                query_chromosome=source_chromosome,
                query_position=source_position,
                hits=hits,
                primary_hit=primary,
                blast_stderr=stderr_text,
            )

        except subprocess.TimeoutExpired:
            self.log_warning("[HomologyAnalyzer] BLAST timed out (300 s)")
            return HomologyResult(
                query_name=query_name,
                query_length=len(sequence),
                query_chromosome=source_chromosome,
                query_position=source_position,
                error="BLAST exceeded the 300 s time limit.",
            )
        except Exception as e:
            self.log_exception(f"[HomologyAnalyzer] Error: {e}")
            return HomologyResult(
                query_name=query_name,
                query_length=len(sequence),
                query_chromosome=source_chromosome,
                query_position=source_position,
                error=str(e),
            )
        finally:
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    #  BLAST execution
    # ------------------------------------------------------------------

    def _run_blast(self, fasta_path: str) -> Tuple[str, str]:
        """
        Run blastn and return (stdout, stderr).
        """
        cfg = self.config.blast
        cmd = [
            self._blastn_path,
            '-task', 'blastn',           # Standard blastn (NOT megablast)
            '-query', fasta_path,
            '-db', self._db_path,
            '-outfmt', self.BLAST_OUTFMT,
            '-evalue', str(cfg.evalue),
            '-max_target_seqs', str(cfg.max_target_seqs),
            '-max_hsps', str(getattr(cfg, 'max_hsps', 3)),
            '-num_threads', str(cfg.threads),
            '-dust', 'no',               # Don't filter low-complexity
            '-soft_masking', 'false',
            '-word_size', '11',          # Standard blastn word size
        ]

        self.log_info(f"[HomologyAnalyzer] Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True, text=True,
            timeout=300,
        )

        stderr_text = result.stderr[:2000] if result.stderr else ""

        if result.returncode != 0:
            self.log_warning(
                f"[HomologyAnalyzer] BLAST exited with code "
                f"{result.returncode}: {stderr_text}"
            )

        return result.stdout, stderr_text

    # ------------------------------------------------------------------
    #  BLAST output parsing
    # ------------------------------------------------------------------

    def _parse_blast_tabular(
        self,
        output: str,
        query_length: int,
        source_chromosome: str,
        source_position: int,
        full_query_sequence: str = "",
    ) -> List[HomologyHit]:
        """
        Parse BLAST tabular output (-outfmt 6) into HomologyHit objects.
        """
        cfg = self.config.blast
        hits: List[HomologyHit] = []

        for line in output.strip().split('\n'):
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < self.BLAST_FIELD_COUNT:
                continue

            try:
                # Parse fields: qseqid sseqid pident length mismatch gapopen
                #                qstart qend sstart send evalue bitscore
                #                qlen slen qseq sseq sstrand
                sseqid = fields[1]
                pident = float(fields[2])
                length = int(fields[3])
                mismatch = int(fields[4])
                gapopen = int(fields[5])
                qstart = int(fields[6])
                qend = int(fields[7])
                sstart = int(fields[8])
                send = int(fields[9])
                evalue = float(fields[10])
                bitscore = float(fields[11])
                qlen = int(fields[12])
                # slen = int(fields[13])  # not needed
                qseq = fields[14]
                sseq = fields[15]
                sstrand = fields[16] if len(fields) > 16 else "plus"

                # Filter by thresholds
                if pident < cfg.min_percent_identity:
                    continue
                if length < cfg.min_aligned_length:
                    continue

                # Determine strand and position
                strand = "-" if sstrand == "minus" else "+"
                ref_position = min(sstart, send)

                # Build synthetic CIGAR and MD from aligned sequences
                core_cigar, md_tag = self._build_cigar_and_md(qseq, sseq)

                # Add soft-clips for unaligned query portions
                cigar = self._add_soft_clips(
                    core_cigar, qstart, qend, qlen
                )

                # Use full query sequence (needed by build_alignment_lines
                # which accounts for soft-clips)
                query_for_hit = full_query_sequence if full_query_sequence else qseq.replace('-', '')

                hits.append(HomologyHit(
                    chromosome=sseqid,
                    position=ref_position,
                    strand=strand,
                    mapq=0,
                    cigar=cigar,
                    alignment_score=int(bitscore),
                    mismatches=mismatch + gapopen,
                    query_sequence=query_for_hit,
                    md_tag=md_tag,
                    aligned_length=length,
                    percent_identity=round(pident, 1),
                    is_primary=False,       # Will be set by _identify_primary
                    is_supplementary=False,
                    alignment_start=qstart - 1,  # Convert to 0-based
                    alignment_end=qend,          # Already exclusive in 0-based
                    evalue=evalue,
                    bit_score=bitscore,
                ))

            except (ValueError, IndexError) as e:
                self.log_warning(
                    f"[HomologyAnalyzer] Skipping BLAST line: {e}"
                )
                continue

        # Identify primary hit
        self._identify_primary_hit(hits, source_chromosome, source_position)

        # Sort: primary first, then by alignment_score descending
        hits.sort(key=lambda h: (-int(h.is_primary), -h.alignment_score))

        return hits

    @staticmethod
    def _add_soft_clips(
        core_cigar: str, qstart: int, qend: int, qlen: int
    ) -> str:
        """Add soft-clip operations to CIGAR for unaligned query portions."""
        parts = []
        leading = qstart - 1  # BLAST qstart is 1-based
        trailing = qlen - qend

        if leading > 0:
            parts.append(f"{leading}S")
        parts.append(core_cigar)
        if trailing > 0:
            parts.append(f"{trailing}S")

        return ''.join(parts)

    def _identify_primary_hit(
        self,
        hits: List[HomologyHit],
        source_chromosome: str,
        source_position: int,
    ) -> None:
        """
        Identify and mark the primary (on-target) hit in place.
        The primary hit is the one on the source chromosome closest
        to the source variant position.
        """
        if not hits:
            return

        norm_src = self.normalize_chr(source_chromosome)

        # Find hits on the source chromosome
        same_chr = [
            h for h in hits
            if self.normalize_chr(h.chromosome) == norm_src
        ]

        if same_chr:
            # Pick the one closest to source_position
            best = min(
                same_chr,
                key=lambda h: abs(h.position - source_position)
            )
            best.is_primary = True
        else:
            # No hit on source chromosome -- mark highest scorer
            hits[0].is_primary = True

    # ------------------------------------------------------------------
    #  Synthetic CIGAR / MD builder from BLAST aligned sequences
    # ------------------------------------------------------------------

    @staticmethod
    def _build_cigar_and_md(
        qseq: str, sseq: str
    ) -> Tuple[str, str]:
        """
        Build synthetic CIGAR string and MD tag from BLAST aligned
        query and subject sequences.

        qseq and sseq are equal-length strings where '-' marks gaps.
        """
        if not qseq or not sseq:
            return "", ""

        # Build flat CIGAR ops and MD parts
        cigar_ops: List[str] = []   # 'M', 'I', 'D' per position
        md_parts: List[str] = []
        match_run = 0

        for q, s in zip(qseq, sseq):
            if q == '-':
                # Deletion in query (subject has base, query doesn't)
                cigar_ops.append('D')
                if match_run > 0:
                    md_parts.append(str(match_run))
                    match_run = 0
                # Accumulate deletion bases for MD ^XYZ
                if md_parts and isinstance(md_parts[-1], str) and md_parts[-1].startswith('^'):
                    md_parts[-1] += s.upper()
                else:
                    md_parts.append(f'^{s.upper()}')

            elif s == '-':
                # Insertion in query (query has base, subject doesn't)
                cigar_ops.append('I')
                # Insertions are NOT represented in MD tag

            elif q.upper() == s.upper():
                # Match
                cigar_ops.append('M')
                match_run += 1

            else:
                # Mismatch
                cigar_ops.append('M')
                if match_run > 0:
                    md_parts.append(str(match_run))
                    match_run = 0
                md_parts.append(s.upper())  # MD records ref base at mismatch

        if match_run > 0:
            md_parts.append(str(match_run))

        # Compress CIGAR: consecutive same ops -> count+op
        cigar_str = _compress_cigar(cigar_ops)
        md_str = ''.join(md_parts)

        return cigar_str, md_str

    # ------------------------------------------------------------------
    #  Alignment reconstruction from MD tag
    # ------------------------------------------------------------------

    _CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
    _MD_RE = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")

    @classmethod
    def build_alignment_lines(
        cls, hit: HomologyHit
    ) -> Tuple[str, str, str]:
        """
        Build (query_line, midline, subject_line) triple from the
        MD tag and query sequence for BLAST-like display.

        MD tag format:
          - digits N  -> N matching bases
          - letter X  -> mismatch (ref base = X)
          - ^XYZ      -> deletion in query (ref has XYZ)
        """
        query_line: List[str] = []
        subject_line: List[str] = []
        mid_line: List[str] = []

        md_parts = cls._MD_RE.findall(hit.md_tag) if hit.md_tag else []
        query_seq = hit.query_sequence

        # We also need to handle insertions from CIGAR (I operations),
        # which are NOT represented in MD tag.
        # Walk CIGAR and MD in parallel for full accuracy.
        cigar_ops = cls._CIGAR_RE.findall(hit.cigar)

        # Build a flat list of operations from CIGAR.
        # Only the LEADING soft-clip affects the query start offset;
        # trailing soft-clips must NOT shift q_idx.
        cigar_flat: List[str] = []
        leading_clip = 0
        seen_non_clip = False
        for count_str, op in cigar_ops:
            count = int(count_str)
            if op == "S":
                if not seen_non_clip:
                    leading_clip = count   # only leading soft-clip
                # trailing S is ignored (doesn't affect q_idx)
                continue
            elif op == "H":
                # Hard-clip: nothing in SEQ
                continue
            else:
                seen_non_clip = True
                cigar_flat.extend([op] * count)

        q_idx = leading_clip  # start reading query after leading soft-clip

        # Now walk CIGAR flat ops + MD together
        md_iter = iter(md_parts)
        md_match_remaining = 0
        current_md_part = None

        cigar_idx = 0
        while cigar_idx < len(cigar_flat):
            op = cigar_flat[cigar_idx]

            if op == "I":
                # Insertion in query (not in ref, not in MD)
                if q_idx < len(query_seq):
                    query_line.append(query_seq[q_idx])
                    subject_line.append("-")
                    mid_line.append(" ")
                    q_idx += 1
                cigar_idx += 1

            elif op in ("M", "=", "X"):
                # Match/mismatch -- consumes both query and ref
                # Use MD to determine match vs mismatch
                if md_match_remaining > 0:
                    # Still in a matching stretch
                    if q_idx < len(query_seq):
                        base = query_seq[q_idx]
                        query_line.append(base)
                        subject_line.append(base)
                        mid_line.append("|")
                        q_idx += 1
                    md_match_remaining -= 1
                    cigar_idx += 1
                else:
                    # Need next MD part
                    try:
                        current_md_part = next(md_iter)
                    except StopIteration:
                        # MD exhausted -- treat as match
                        if q_idx < len(query_seq):
                            base = query_seq[q_idx]
                            query_line.append(base)
                            subject_line.append(base)
                            mid_line.append("|")
                            q_idx += 1
                        cigar_idx += 1
                        continue

                    if current_md_part.isdigit():
                        md_match_remaining = int(current_md_part)
                        # Don't advance cigar_idx -- re-enter loop
                        continue
                    elif current_md_part.startswith("^"):
                        # Deletion -- but we're in M op?
                        # This shouldn't happen, but handle gracefully
                        for ref_base in current_md_part[1:]:
                            query_line.append("-")
                            subject_line.append(ref_base)
                            mid_line.append(" ")
                        continue
                    else:
                        # Mismatch -- each letter is a ref base
                        for ref_base in current_md_part:
                            if q_idx < len(query_seq):
                                query_line.append(query_seq[q_idx])
                                subject_line.append(ref_base)
                                mid_line.append(".")
                                q_idx += 1
                            cigar_idx += 1
                        continue

            elif op == "D":
                # Deletion in query -- MD should have ^XYZ
                if md_match_remaining > 0:
                    # Shouldn't happen -- but be defensive
                    cigar_idx += 1
                    continue

                try:
                    current_md_part = next(md_iter)
                except StopIteration:
                    query_line.append("-")
                    subject_line.append("?")
                    mid_line.append(" ")
                    cigar_idx += 1
                    continue

                if current_md_part.startswith("^"):
                    del_bases = current_md_part[1:]
                    for ref_base in del_bases:
                        query_line.append("-")
                        subject_line.append(ref_base)
                        mid_line.append(" ")
                        cigar_idx += 1
                elif current_md_part.isdigit():
                    # Unexpected -- push back and handle as gap
                    md_match_remaining = int(current_md_part)
                    query_line.append("-")
                    subject_line.append("N")
                    mid_line.append(" ")
                    cigar_idx += 1
                else:
                    query_line.append("-")
                    subject_line.append(current_md_part)
                    mid_line.append(" ")
                    cigar_idx += 1

            elif op == "N":
                # Skipped region (intron) -- just advance cigar
                cigar_idx += 1
            else:
                cigar_idx += 1

        return "".join(query_line), "".join(mid_line), "".join(subject_line)

    # ------------------------------------------------------------------
    #  CIGAR helpers
    # ------------------------------------------------------------------

    @classmethod
    def _cigar_ref_length(cls, cigar: str) -> int:
        """Compute alignment length on the reference from CIGAR."""
        length = 0
        for m in cls._CIGAR_RE.finditer(cigar):
            count = int(m.group(1))
            op = m.group(2)
            if op in ("M", "D", "N", "=", "X"):  # ref-consuming ops
                length += count
        return length

    @classmethod
    def _cigar_query_range(
        cls, cigar: str, query_length: int
    ) -> Tuple[int, int]:
        """
        Compute query range (0-based) from CIGAR accounting for
        leading/trailing soft-clips.
        """
        ops = cls._CIGAR_RE.findall(cigar)
        start = 0
        end = query_length

        if ops and ops[0][1] == "S":
            start = int(ops[0][0])
        if len(ops) >= 2 and ops[-1][1] == "S":
            end = query_length - int(ops[-1][0])

        return start, end

    # ------------------------------------------------------------------
    #  Discriminating positions extraction (for primer design)
    # ------------------------------------------------------------------

    @classmethod
    def extract_discriminating_positions(
        cls, result: HomologyResult
    ) -> Dict[int, int]:
        """
        Extract **genomic** positions where our sequence differs from
        homologous regions (pseudogenes, duplications).

        Walks the MD tag + CIGAR of each NON-primary hit to find mismatch
        and indel positions in query coordinates (0-based), then converts
        them to absolute genomic coordinates using
        ``result.query_genomic_start``.

        Returns:
            Dict mapping **genomic position** to the number of secondary
            hits that have a mismatch/indel at that position.
            Higher count = position discriminates from more homologs.
        """
        disc_positions: Dict[int, int] = {}

        secondary_hits = [
            h for h in result.hits
            if not h.is_primary and not h.is_supplementary
        ]

        if not secondary_hits:
            return disc_positions

        genomic_offset = result.query_genomic_start  # 0 if not set

        for hit in secondary_hits:
            mismatch_query_positions = cls._get_mismatch_query_positions(hit)
            for qpos in mismatch_query_positions:
                genomic_pos = genomic_offset + qpos
                disc_positions[genomic_pos] = disc_positions.get(genomic_pos, 0) + 1

        return disc_positions

    @classmethod
    def _get_mismatch_query_positions(cls, hit: HomologyHit) -> List[int]:
        """
        Extract 0-based query positions where mismatches or indels occur
        relative to the homologous hit.

        Uses CIGAR + MD tag to walk through the alignment and identify
        positions where query differs from the hit's reference sequence.
        """
        positions: List[int] = []

        md_parts = cls._MD_RE.findall(hit.md_tag) if hit.md_tag else []
        if not md_parts:
            return positions

        cigar_ops = cls._CIGAR_RE.findall(hit.cigar)

        # Build flat CIGAR ops, track leading soft-clip
        cigar_flat: List[str] = []
        leading_clip = 0
        seen_non_clip = False
        for count_str, op in cigar_ops:
            count = int(count_str)
            if op == "S":
                if not seen_non_clip:
                    leading_clip = count
                continue
            elif op == "H":
                continue
            else:
                seen_non_clip = True
                cigar_flat.extend([op] * count)

        q_idx = leading_clip
        md_iter = iter(md_parts)
        md_match_remaining = 0
        cigar_idx = 0

        while cigar_idx < len(cigar_flat):
            op = cigar_flat[cigar_idx]

            if op == "I":
                # Insertion in query -- query has extra bases not in homolog
                # This position in query is discriminating
                positions.append(q_idx)
                q_idx += 1
                cigar_idx += 1

            elif op in ("M", "=", "X"):
                if md_match_remaining > 0:
                    # Matching base -- not discriminating
                    q_idx += 1
                    md_match_remaining -= 1
                    cigar_idx += 1
                else:
                    try:
                        current_md_part = next(md_iter)
                    except StopIteration:
                        q_idx += 1
                        cigar_idx += 1
                        continue

                    if current_md_part.isdigit():
                        md_match_remaining = int(current_md_part)
                        continue
                    elif current_md_part.startswith("^"):
                        # Deletion -- handled in D ops normally,
                        # but if seen in M context, skip gracefully
                        continue
                    else:
                        # Mismatch letters -- each is a discriminating position
                        for _ in current_md_part:
                            positions.append(q_idx)
                            q_idx += 1
                            cigar_idx += 1
                        continue

            elif op == "D":
                # Deletion in query -- homolog has extra bases
                # The adjacent query positions are discriminating context,
                # but the deletion itself doesn't map to a query position.
                # We mark the current q_idx as "near deletion" discriminating.
                if q_idx not in positions:
                    positions.append(q_idx)

                if md_match_remaining > 0:
                    cigar_idx += 1
                    continue

                try:
                    current_md_part = next(md_iter)
                except StopIteration:
                    cigar_idx += 1
                    continue

                if current_md_part.startswith("^"):
                    del_len = len(current_md_part) - 1  # minus ^
                    # Advance cigar_idx for all D ops in this deletion
                    cigar_idx += del_len
                elif current_md_part.isdigit():
                    md_match_remaining = int(current_md_part)
                    cigar_idx += 1
                else:
                    cigar_idx += 1

            elif op == "N":
                cigar_idx += 1
            else:
                cigar_idx += 1

        return positions

    @classmethod
    def score_primer_discrimination(
        cls,
        primer_start: int,
        primer_length: int,
        is_forward: bool,
        discriminating_positions: Dict[int, int],
    ) -> Tuple[float, int]:
        """
        Score how well a primer discriminates against homologous regions.

        Mismatches at the 3' end of a primer are weighted more heavily
        because they most effectively prevent polymerase extension on
        the non-target (pseudogene) template.

        Args:
            primer_start: 0-based start position of primer in query sequence
            primer_length: Length of primer
            is_forward: True for forward primer, False for reverse
            discriminating_positions: Dict from extract_discriminating_positions()

        Returns:
            Tuple of (weighted_score, count_of_discriminating_positions)
        """
        score = 0.0
        count = 0
        primer_end = primer_start + primer_length  # exclusive

        for pos in range(primer_start, primer_end):
            if pos not in discriminating_positions:
                continue

            count += 1
            n_homologs = discriminating_positions[pos]

            # Calculate distance from the 3' end of the primer
            if is_forward:
                # Forward primer: 3' end is at primer_end - 1
                dist_from_3prime = (primer_end - 1) - pos
            else:
                # Reverse primer: 3' end is at primer_start
                dist_from_3prime = pos - primer_start

            # Weight: 3' end positions are most valuable
            if dist_from_3prime < 5:
                weight = 3.0
            elif dist_from_3prime < 10:
                weight = 2.0
            else:
                weight = 1.0

            score += weight * n_homologs

        return score, count

    # ------------------------------------------------------------------
    #  Chromosome normalisation
    # ------------------------------------------------------------------

    @staticmethod
    def normalize_chr(chrom: str) -> str:
        """Normalize chromosome name (strip 'chr' prefix, uppercase)."""
        return str(chrom).replace("chr", "").upper()

    @staticmethod
    def format_chr(chrom: str) -> str:
        """Ensure chromosome has 'chr' prefix (add if missing, don't duplicate)."""
        s = str(chrom)
        return s if s.startswith("chr") else f"chr{s}"


# ---------------------------------------------------------------------------
#  Module-level helper
# ---------------------------------------------------------------------------

def _compress_cigar(ops: List[str]) -> str:
    """Compress a flat list of CIGAR ops into a standard CIGAR string."""
    if not ops:
        return ""

    parts: List[str] = []
    current_op = ops[0]
    count = 1

    for op in ops[1:]:
        if op == current_op:
            count += 1
        else:
            parts.append(f"{count}{current_op}")
            current_op = op
            count = 1

    parts.append(f"{count}{current_op}")
    return ''.join(parts)

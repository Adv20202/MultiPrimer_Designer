"""
Orchestrates specificity checking across two independent tools:
  1. NCBI Primer-BLAST    (gold-standard, online)
  2. UCSC In-Silico PCR   (quick, online)

Combines results into a single verdict per primer pair.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import List, Optional, Callable
from enum import Enum

from ..core.models import PrimerPair
from ..utils.logger import LoggerMixin
from ..api.primer_blast_client import PrimerBlastClient
from ..api.ucsc_ispcr_client import UCSCisPCRClient


# ======================================================================
#  Data classes for combined results
# ======================================================================

class SpecificityVerdict(Enum):
    """Overall specificity verdict shown to the user."""
    SPECIFIC = "SPECIFIC"
    LIKELY_SPECIFIC = "LIKELY SPECIFIC"
    OFF_TARGET_DETECTED = "OFF-TARGET DETECTED"
    PSEUDOGENE_RISK = "PSEUDOGENE RISK"
    INCONCLUSIVE = "INCONCLUSIVE"


@dataclass
class ToolResult:
    """Result from a single specificity tool."""
    tool_name: str                          # "Primer-BLAST", "UCSC isPCR"
    available: bool                         # Was the tool available/reachable?
    ran_successfully: bool                  # Completed without errors?
    is_specific: bool = True
    off_target_count: int = 0
    warnings: List[str] = field(default_factory=list)
    error: str = ""
    off_target_locations: List[str] = field(default_factory=list)


@dataclass
class CombinedSpecificityResult:
    """Merged result from specificity tools."""
    verdict: SpecificityVerdict
    is_specific: bool
    tool_results: List[ToolResult]
    off_target_summary: List[str]           # All off-target locations
    all_warnings: List[str]
    pseudogene_warnings: List[str]
    tools_used: int
    tools_available: int


# ======================================================================
#  Orchestrator
# ======================================================================

class SpecificityOrchestrator(LoggerMixin):
    """
    Coordinates Primer-BLAST and UCSC isPCR for comprehensive
    primer specificity checking.
    """

    def __init__(self):
        # Remote tools (always available if internet works)
        self.primer_blast = PrimerBlastClient()
        self.ucsc_ispcr = UCSCisPCRClient()

    # ------------------------------------------------------------------
    #  Main entry point
    # ------------------------------------------------------------------

    def check_all(
        self,
        primer_pair: PrimerPair,
        expected_chromosome: Optional[str] = None,
        expected_start: Optional[int] = None,
        expected_end: Optional[int] = None,
        assembly: str = "GRCh38",
        progress_callback: Optional[Callable] = None,
        use_primer_blast: bool = True,
        use_ucsc: bool = True,
    ) -> CombinedSpecificityResult:
        """
        Run selected specificity tools **in parallel** and merge results.
        Only tools with use_* = True are submitted; others get a "skipped" result.
        """
        tool_results: List[ToolResult] = [None, None]  # Preserve order

        # Create "skipped" placeholders for disabled tools
        if not use_ucsc:
            tool_results[0] = ToolResult(
                tool_name="UCSC isPCR", available=False,
                ran_successfully=False, error="Skipped by user"
            )
        if not use_primer_blast:
            tool_results[1] = ToolResult(
                tool_name="Primer-BLAST", available=False,
                ran_successfully=False, error="Skipped by user"
            )

        # Only submit enabled tools
        futures = {}
        active_tools = sum([use_primer_blast, use_ucsc])
        if active_tools == 0:
            return self._combine(tool_results, expected_chromosome)

        with ThreadPoolExecutor(max_workers=active_tools) as executor:
            if use_ucsc:
                futures[executor.submit(
                    self._run_ucsc, primer_pair, assembly,
                    expected_chromosome, expected_start,
                    expected_end, progress_callback
                )] = 0
            if use_primer_blast:
                futures[executor.submit(
                    self._run_primer_blast, primer_pair, assembly,
                    expected_chromosome, expected_start,
                    expected_end, progress_callback
                )] = 1

            for future in as_completed(futures):
                idx = futures[future]
                try:
                    tool_results[idx] = future.result()
                except Exception as e:
                    names = ["UCSC isPCR", "Primer-BLAST"]
                    tool_results[idx] = ToolResult(
                        tool_name=names[idx], available=True,
                        ran_successfully=False, error=str(e)
                    )

        return self._combine(tool_results, expected_chromosome)

    # ------------------------------------------------------------------
    #  Individual tool runners
    # ------------------------------------------------------------------

    def _run_primer_blast(self, pp, assembly, exp_chr, exp_start, exp_end, cb) -> ToolResult:
        try:
            result = self.primer_blast.check_specificity(
                pp.forward.sequence, pp.reverse.sequence,
                expected_chromosome=exp_chr,
                expected_start=exp_start,
                expected_end=exp_end,
                progress_callback=cb,
                assembly=assembly,
            )
            # Build off-target location strings using chromosome info
            # Format: "chr17:NC_000017.11 (description, Nbp, F:Xmm R:Ymm)"
            # so the chromosome is parseable by _combine for pseudogene detection
            off_locs = []
            for h in result.hits:
                if h.is_intended_target:
                    continue
                fwd = '?' if h.forward_mismatches < 0 else str(h.forward_mismatches)
                rev = '?' if h.reverse_mismatches < 0 else str(h.reverse_mismatches)
                if h.chromosome:
                    # Put chromosome in parseable format: "chr<N>:<accession> (...)"
                    loc = (f"chr{h.chromosome}:{h.accession} "
                           f"({h.description[:50]}, {h.product_size}bp, "
                           f"F:{fwd}mm R:{rev}mm)")
                else:
                    loc = (f"{h.accession} ({h.description[:50]}, "
                           f"{h.product_size}bp, F:{fwd}mm R:{rev}mm)")
                off_locs.append(loc)

            return ToolResult(
                tool_name="Primer-BLAST", available=True,
                ran_successfully=result.success,
                is_specific=result.is_specific,
                off_target_count=result.off_target_count,
                warnings=result.warnings,
                error=result.error,
                off_target_locations=off_locs,
            )
        except Exception as e:
            self.log_warning(f"Primer-BLAST check failed: {e}")
            return ToolResult(
                tool_name="Primer-BLAST", available=True,
                ran_successfully=False, error=str(e)
            )

    def _run_ucsc(self, pp, assembly, exp_chr, exp_start, exp_end, cb) -> ToolResult:
        genome = "hg38" if assembly == "GRCh38" else "hg19"
        try:
            result = self.ucsc_ispcr.check_specificity(
                pp.forward.sequence, pp.reverse.sequence,
                genome=genome,
                expected_chromosome=exp_chr,
                expected_start=exp_start,
                expected_end=exp_end,
                progress_callback=cb,
            )
            off_locs = [
                f"{h.chromosome}:{h.start}-{h.end} ({h.product_size}bp)"
                for h in result.hits if not h.is_intended_target
            ]
            return ToolResult(
                tool_name="UCSC isPCR", available=True,
                ran_successfully=result.success,
                is_specific=result.is_specific,
                off_target_count=result.off_target_count,
                warnings=result.warnings,
                error=result.error,
                off_target_locations=off_locs,
            )
        except Exception as e:
            self.log_warning(f"UCSC isPCR check failed: {e}")
            return ToolResult(
                tool_name="UCSC isPCR", available=True,
                ran_successfully=False, error=str(e)
            )

    # ------------------------------------------------------------------
    #  Merge results into a single verdict
    # ------------------------------------------------------------------

    def _combine(
        self,
        tool_results: List[ToolResult],
        expected_chromosome: Optional[str] = None,
    ) -> CombinedSpecificityResult:

        tools_available = sum(1 for t in tool_results if t.available)
        tools_ran = sum(1 for t in tool_results if t.ran_successfully)

        all_warnings: List[str] = []
        all_off_targets: List[str] = []
        pseudogene_warnings: List[str] = []
        any_off_target = False

        for tr in tool_results:
            for w in tr.warnings:
                tag = f"[{tr.tool_name}] "
                all_warnings.append(w if w.startswith(tag) else tag + w)

            if tr.ran_successfully and not tr.is_specific:
                any_off_target = True

            for loc in tr.off_target_locations:
                prefixed = f"[{tr.tool_name}] {loc}"
                all_off_targets.append(prefixed)

                # Pseudogene detection: off-target on SAME chromosome
                if expected_chromosome:
                    loc_chr = self._extract_chromosome_from_location(loc)
                    if loc_chr:
                        norm_exp = str(expected_chromosome).replace('chr', '')
                        norm_loc = loc_chr.replace('chr', '')
                        if norm_loc == norm_exp:
                            pseudogene_warnings.append(
                                f"[{tr.tool_name}] Possible pseudogene/duplication: {loc}"
                            )

        # Determine verdict
        if tools_ran == 0:
            verdict = SpecificityVerdict.INCONCLUSIVE
            # INCONCLUSIVE should NOT count as specific — we couldn't verify,
            # so in continuous mode this pair should not satisfy the "specific
            # pairs found" counter.
            is_specific = False
        elif pseudogene_warnings:
            verdict = SpecificityVerdict.PSEUDOGENE_RISK
            is_specific = False
        elif any_off_target:
            verdict = SpecificityVerdict.OFF_TARGET_DETECTED
            is_specific = False
        elif tools_ran < tools_available:
            verdict = SpecificityVerdict.LIKELY_SPECIFIC
            is_specific = True
        else:
            verdict = SpecificityVerdict.SPECIFIC
            is_specific = True

        return CombinedSpecificityResult(
            verdict=verdict,
            is_specific=is_specific,
            tool_results=tool_results,
            off_target_summary=all_off_targets,
            all_warnings=all_warnings,
            pseudogene_warnings=pseudogene_warnings,
            tools_used=tools_ran,
            tools_available=tools_available,
        )

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_chromosome_from_location(loc: str) -> str:
        """
        Extract a normalised chromosome name from an off-target location string.

        Handles formats from both tools:
          UCSC isPCR:   "chr17:43044295-43044575 (281bp)"
          Primer-BLAST: "chr17:NC_000017.11 (Homo sapiens..., 358bp, ...)"
                        "NC_000017.11 (Homo sapiens..., 358bp, ...)"

        Returns the chromosome part (e.g. "17", "X") or empty string.
        """
        loc = loc.strip()

        # Pattern 1: starts with "chr<N>:" (both tools can produce this)
        #   "chr17:43044295-43044575 ..."  or  "chr17:NC_000017.11 ..."
        import re
        m = re.match(r'chr(\d+|X|Y|MT?)\b', loc, re.IGNORECASE)
        if m:
            return m.group(1)

        # Pattern 2: bare NC_ accession — parse from accession number
        # "NC_000017.11 (Homo sapiens chromosome 17...)"
        m = re.match(r'NC_0000(\d{2})', loc)
        if m:
            num = int(m.group(1))
            if 1 <= num <= 22:
                return str(num)
            elif num == 23:
                return 'X'
            elif num == 24:
                return 'Y'

        # Pattern 3: parse "chromosome N" from the description part
        m = re.search(r'chromosome\s+(\d+|X|Y|MT)', loc, re.IGNORECASE)
        if m:
            return m.group(1)

        return ''

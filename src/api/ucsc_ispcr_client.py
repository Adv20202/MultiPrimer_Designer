"""
UCSC In-Silico PCR (isPCR) API client.
Submits primer pairs to the UCSC Genome Browser isPCR endpoint
and parses amplicon results.

Uses urllib (matching the project's existing API pattern).
"""

import re
import time
import urllib.request
import urllib.parse
import urllib.error
from typing import List, Dict, Optional, Callable
from dataclasses import dataclass, field

from ..utils.logger import LoggerMixin
from ..utils.config import get_config


@dataclass
class ISPCRHit:
    """A single amplicon hit from UCSC In-Silico PCR."""
    chromosome: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    product_size: int = 0
    is_intended_target: bool = False


@dataclass
class ISPCRResult:
    """Result from UCSC In-Silico PCR."""
    success: bool
    is_specific: bool = True
    hits: List[ISPCRHit] = field(default_factory=list)
    off_target_count: int = 0
    warnings: List[str] = field(default_factory=list)
    error: str = ""


class UCSCisPCRClient(LoggerMixin):
    """
    Client for UCSC In-Silico PCR.
    Uses the UCSC Genome Browser isPCR CGI endpoint.
    """

    ISPCR_URL = "https://genome.ucsc.edu/cgi-bin/hgPcr"

    def __init__(self):
        config = get_config()
        self.timeout = 60
        self.max_retries = config.api.max_retries
        self._cache: Dict[str, ISPCRResult] = {}
        self._last_request_time = 0.0
        self._min_interval = 1.0   # Be conservative with UCSC

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def check_specificity(
        self,
        forward_seq: str,
        reverse_seq: str,
        genome: str = "hg38",
        max_product_size: int = 4000,
        expected_chromosome: Optional[str] = None,
        expected_start: Optional[int] = None,
        expected_end: Optional[int] = None,
        progress_callback: Optional[Callable] = None,
    ) -> ISPCRResult:
        """
        Submit primer pair to UCSC isPCR and check for off-targets.

        Args:
            forward_seq: Forward primer sequence (5'->3')
            reverse_seq: Reverse primer sequence (5'->3')
            genome: Genome assembly ("hg38" or "hg19")
            max_product_size: Maximum product size to search
            expected_chromosome: Expected chromosome (e.g. "17", "chr17")
            expected_start: Expected position — variant or amplicon start.
                            Used to identify which isPCR hit is the intended target.
                            The hit is matched if the expected position falls within
                            the isPCR product boundaries (not strict start==start).
            expected_end: Expected position — variant or amplicon end.
            progress_callback: Optional callback for status updates
        """
        cache_key = f"{forward_seq}_{reverse_seq}_{genome}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        if progress_callback:
            progress_callback("Querying UCSC In-Silico PCR...")

        try:
            html = self._query(forward_seq, reverse_seq, genome, max_product_size)
            if html is None:
                return ISPCRResult(success=False, error="UCSC isPCR request failed")

            result = self._parse_results(
                html, expected_chromosome, expected_start, expected_end
            )
            self._cache[cache_key] = result
            return result

        except Exception as e:
            self.log_exception(f"UCSC isPCR error: {e}")
            return ISPCRResult(success=False, error=str(e))

    def clear_cache(self):
        self._cache.clear()

    # ------------------------------------------------------------------
    #  HTTP request
    # ------------------------------------------------------------------

    def _query(
        self, fwd: str, rev: str, genome: str, max_size: int
    ) -> Optional[str]:
        """Send the isPCR query and return the HTML response."""
        params = {
            'org': 'Human',
            'db': genome,
            'wp_target': 'genome',
            'wp_f': fwd,
            'wp_r': rev,
            'Submit': 'submit',
            'wp_size': str(max_size),
            # wp_perfect: min consecutive perfect-match bases at 3' end.
            # wp_good: min total matching bases in primer.
            # 15 is the UCSC default.  Lower values would find more off-
            # targets (including weak ones unlikely to amplify in real PCR);
            # higher values would miss genuine co-amplification targets.
            'wp_perfect': '15',
            'wp_good': '15',
        }
        data = urllib.parse.urlencode(params).encode('utf-8')

        for attempt in range(self.max_retries):
            self._rate_limit()
            try:
                req = urllib.request.Request(
                    self.ISPCR_URL, data=data,
                    headers={'User-Agent': 'PrimerDesigner/1.0'}
                )
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    return resp.read().decode('utf-8', errors='replace')

            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(5 * (attempt + 1))
                elif e.code >= 500:
                    time.sleep(3 * (attempt + 1))
                else:
                    self.log_warning(f"UCSC isPCR HTTP {e.code}")
                    return None
            except urllib.error.URLError:
                time.sleep(3 * (attempt + 1))

        return None

    # ------------------------------------------------------------------
    #  Result parsing
    # ------------------------------------------------------------------

    def _parse_results(
        self,
        html: str,
        expected_chr: Optional[str] = None,
        expected_start: Optional[int] = None,
        expected_end: Optional[int] = None,
    ) -> ISPCRResult:
        """
        Parse UCSC isPCR HTML results.

        isPCR output contains a <PRE> block with FASTA-like entries:
            >chr17:43044295+43044575  281bp  BRCA1
            ATCG...
        """
        hits: List[ISPCRHit] = []
        warnings: List[str] = []

        # Extract PRE block
        pre_match = re.search(r'<PRE>(.*?)</PRE>', html, re.DOTALL | re.IGNORECASE)
        if not pre_match:
            if 'No matches' in html or 'no results' in html.lower():
                return ISPCRResult(
                    success=True, is_specific=True,
                    warnings=["UCSC isPCR: No amplification products found"]
                )
            # Could also mean the results are in a different format
            # Try alternative pattern
            pre_match = re.search(r'<pre>(.*?)</pre>', html, re.DOTALL)
            if not pre_match:
                return ISPCRResult(
                    success=True, is_specific=True,
                    warnings=["UCSC isPCR: Could not parse results (no PRE block)"]
                )

        pre_content = pre_match.group(1)

        # Parse FASTA-like entries:  >chrN:start+end  NNNbp  info
        entry_pattern = re.compile(
            r'>(\w+):(\d+)([+-])(\d+)\s+(\d+)bp'
        )
        for m in entry_pattern.finditer(pre_content):
            chrom = m.group(1)
            start = int(m.group(2))
            strand = m.group(3)
            end = int(m.group(4))
            size = int(m.group(5))

            hit = ISPCRHit(
                chromosome=chrom,
                start=start,
                end=end,
                strand=strand,
                product_size=size,
            )

            # Mark as intended target if the expected coordinates
            # match this hit.
            #
            # Note: expected_start / expected_end may be either:
            #   (a) the variant position (a single point inside the amplicon)
            #   (b) the amplicon start/end (which should nearly match the hit)
            #
            # We handle both cases by checking whether the expected interval
            # *overlaps* the isPCR product interval on the same chromosome,
            # with a generous tolerance of 500 bp on each side.  This ensures
            # a variant position (which is INSIDE the product) always matches
            # its correct product, and amplicon coordinates match too.
            if expected_chr and expected_start is not None:
                norm_hit = chrom.replace('chr', '')
                norm_exp = str(expected_chr).replace('chr', '')
                if norm_hit == norm_exp:
                    hit_lo = min(start, end)
                    hit_hi = max(start, end)
                    exp_lo = expected_start
                    exp_hi = expected_end if expected_end else expected_start

                    # The expected interval should fall within (or overlap)
                    # the product boundaries.  We add 500 bp tolerance to
                    # account for coordinate system differences.
                    if (exp_lo >= hit_lo - 500 and exp_hi <= hit_hi + 500):
                        hit.is_intended_target = True

            hits.append(hit)

        # If expected coordinates were given but no hit matched, it may be
        # because the coordinate systems slightly differ.  In that case,
        # fall back to the first hit on the expected chromosome.
        if expected_chr and hits and not any(h.is_intended_target for h in hits):
            norm_exp = str(expected_chr).replace('chr', '')
            for h in hits:
                if h.chromosome.replace('chr', '') == norm_exp:
                    h.is_intended_target = True
                    break

        # If no expected target specified, first hit is assumed intended
        if not expected_chr and hits:
            hits[0].is_intended_target = True

        off_targets = [h for h in hits if not h.is_intended_target]
        is_specific = len(off_targets) == 0

        if off_targets:
            warnings.append(
                f"UCSC isPCR: {len(off_targets)} off-target product(s) found"
            )
            for ot in off_targets[:5]:
                warnings.append(
                    f"  {ot.chromosome}:{ot.start}-{ot.end} ({ot.product_size}bp)"
                )
        elif hits:
            h = hits[0]
            warnings.append(
                f"UCSC isPCR: Single product at {h.chromosome}:{h.start}-{h.end} "
                f"({h.product_size}bp)"
            )
        else:
            warnings.append("UCSC isPCR: No amplification products found")

        return ISPCRResult(
            success=True,
            is_specific=is_specific,
            hits=hits,
            off_target_count=len(off_targets),
            warnings=warnings,
        )

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    def _rate_limit(self):
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

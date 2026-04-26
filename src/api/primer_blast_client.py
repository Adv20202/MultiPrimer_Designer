"""
NCBI Primer-BLAST API client.
Submits primer pairs to Primer-BLAST and parses results for off-target
amplification products.

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
class PrimerBlastHit:
    """A single off-target hit from Primer-BLAST."""
    accession: str = ""
    description: str = ""
    chromosome: str = ""                # Parsed from accession/description
    product_size: int = 0
    forward_mismatches: int = -1        # -1 = unknown (parsing failed)
    reverse_mismatches: int = -1        # -1 = unknown (parsing failed)
    total_mismatches: int = -1          # -1 = unknown
    is_intended_target: bool = False


@dataclass
class PrimerBlastResult:
    """Result from Primer-BLAST specificity check."""
    success: bool
    is_specific: bool = True
    hits: List[PrimerBlastHit] = field(default_factory=list)
    off_target_count: int = 0
    warnings: List[str] = field(default_factory=list)
    error: str = ""


# Mapping of RefSeq chromosome accessions (GRCh38) to chromosome names
_REFSEQ_CHR_MAP = {
    'NC_000001': '1',  'NC_000002': '2',  'NC_000003': '3',
    'NC_000004': '4',  'NC_000005': '5',  'NC_000006': '6',
    'NC_000007': '7',  'NC_000008': '8',  'NC_000009': '9',
    'NC_000010': '10', 'NC_000011': '11', 'NC_000012': '12',
    'NC_000013': '13', 'NC_000014': '14', 'NC_000015': '15',
    'NC_000016': '16', 'NC_000017': '17', 'NC_000018': '18',
    'NC_000019': '19', 'NC_000020': '20', 'NC_000021': '21',
    'NC_000022': '22', 'NC_000023': 'X',  'NC_000024': 'Y',
    # GRCh37
    'NC_000001.10': '1',  'NC_000002.11': '2',  'NC_000003.11': '3',
    'NC_000004.11': '4',  'NC_000005.9': '5',   'NC_000006.11': '6',
    'NC_000007.13': '7',  'NC_000008.10': '8',  'NC_000009.11': '9',
    'NC_000010.10': '10', 'NC_000011.9': '11',  'NC_000012.11': '12',
    'NC_000013.10': '13', 'NC_000014.8': '14',  'NC_000015.9': '15',
    'NC_000016.9': '16',  'NC_000017.10': '17', 'NC_000018.9': '18',
    'NC_000019.9': '19',  'NC_000020.10': '20', 'NC_000021.8': '21',
    'NC_000022.10': '22', 'NC_000023.10': 'X',  'NC_000024.9': 'Y',
}


class PrimerBlastClient(LoggerMixin):
    """
    Client for NCBI Primer-BLAST.
    Submits primer pairs and checks for off-target amplification.
    """

    SUBMIT_URL = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"

    # Polling configuration — optimized for speed
    MAX_POLL_ATTEMPTS = 30       # up to ~3 min with 5-s interval
    POLL_INTERVAL = 5            # seconds between polls
    INITIAL_WAIT = 3             # seconds before first poll

    def __init__(self):
        config = get_config()
        self.timeout = 60
        self.max_retries = config.api.max_retries
        self.ncbi_api_key = getattr(config.api, 'ncbi_api_key', '')
        # Rate limit: 3 req/s without API key, 10 req/s with key
        self._min_interval = 0.1 if self.ncbi_api_key else 0.34
        self._last_request_time = 0.0
        self._cache: Dict[str, PrimerBlastResult] = {}

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    # Assembly accessions for Primer-BLAST Custom database mode.
    # These are NCBI RefSeq accessions for human genome assemblies.
    _ASSEMBLY_ACCESSIONS = {
        'GRCh38': 'GCF_000001405.40',   # GRCh38.p14
        'GRCh37': 'GCF_000001405.25',   # GRCh37.p13
    }

    def check_specificity(
        self,
        forward_seq: str,
        reverse_seq: str,
        max_product_size: int = 4000,
        expected_chromosome: Optional[str] = None,
        expected_start: Optional[int] = None,
        expected_end: Optional[int] = None,
        progress_callback: Optional[Callable] = None,
        assembly: str = "GRCh38",
    ) -> PrimerBlastResult:
        """
        Submit a primer pair to Primer-BLAST and return specificity info.

        Args:
            forward_seq: Forward primer sequence (5'->3')
            reverse_seq: Reverse primer sequence (5'->3')
            max_product_size: Maximum product size to search
            expected_chromosome: Expected chromosome (e.g. "17", "chr17")
            expected_start: Expected amplicon start (genomic coordinate)
            expected_end: Expected amplicon end (genomic coordinate)
            progress_callback: Optional callback for status updates
            assembly: Reference genome assembly ("GRCh38" or "GRCh37")
        """
        cache_key = f"{forward_seq}_{reverse_seq}_{assembly}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        try:
            # Step 1: Submit the job
            if progress_callback:
                progress_callback("Submitting to Primer-BLAST...")

            html_or_url = self._submit_job(
                forward_seq, reverse_seq, max_product_size, assembly
            )
            if html_or_url is None:
                return PrimerBlastResult(
                    success=False,
                    error="Failed to submit Primer-BLAST job"
                )

            # If we got HTML directly (server-side cache hit), parse it
            if html_or_url.startswith('<') or 'Primer pair specificity' in html_or_url:
                result = self._parse_results(
                    html_or_url, expected_chromosome, expected_start, expected_end
                )
                self._cache[cache_key] = result
                return result

            # Step 2: Poll for results
            if progress_callback:
                progress_callback("Waiting for Primer-BLAST results...")

            time.sleep(self.INITIAL_WAIT)
            html = self._poll_for_results(html_or_url, progress_callback)

            if not html:
                return PrimerBlastResult(
                    success=False,
                    error="Primer-BLAST timed out or returned no results"
                )

            result = self._parse_results(
                html, expected_chromosome, expected_start, expected_end
            )
            self._cache[cache_key] = result
            return result

        except Exception as e:
            self.log_exception(f"Primer-BLAST error: {e}")
            return PrimerBlastResult(success=False, error=str(e))

    def clear_cache(self):
        self._cache.clear()

    # ------------------------------------------------------------------
    #  Job submission
    # ------------------------------------------------------------------

    def _submit_job(
        self, fwd: str, rev: str, max_product: int, assembly: str = "GRCh38"
    ) -> Optional[str]:
        """
        Submit a Primer-BLAST job.
        Returns HTML (if results are immediate) or the poll URL string.
        Returns None on failure.
        """
        # Choose database based on assembly.
        # We use the NCBI "Custom" database mode with the specific RefSeq
        # genome assembly accession so that Primer-BLAST checks specificity
        # against the EXACT genome build the user selected (GRCh38 or GRCh37).
        # This requires two params:
        #   PRIMER_SPECIFICITY_DATABASE = "Custom"
        #   CUSTOM_DB = "<assembly_accession>"
        # When ORGANISM is set with Custom DB, NCBI ignores it (per docs),
        # so specificity is checked purely against the chosen assembly.
        assembly_acc = self._ASSEMBLY_ACCESSIONS.get(assembly)

        params = {
            'PRIMER_LEFT_INPUT': fwd,
            'PRIMER_RIGHT_INPUT': rev,
            'SEARCH_SPECIFIC_PRIMER': 'on',
            'PRIMER_PRODUCT_MAX': str(max_product),
            'TOTAL_MISMATCH_IGNORE': '6',
            'MAX_TARGET_SIZE': '4000',
            'HITSIZE': '50000',
            'EVALUE': '30000',
            'WORD_SIZE': '7',
            'NUM_TARGETS_WITH_PRIMERS': '1000',
            'MAX_TARGET_PER_TEMPLATE': '100',
        }

        if assembly_acc:
            # Custom database mode — use the exact genome assembly
            params['PRIMER_SPECIFICITY_DATABASE'] = 'Custom'
            params['CUSTOM_DB'] = assembly_acc
        else:
            # Fallback for unknown assemblies
            params['PRIMER_SPECIFICITY_DATABASE'] = 'refseq_representative_genomes'
            params['ORGANISM'] = '9606'
        if self.ncbi_api_key:
            params['api_key'] = self.ncbi_api_key

        data = urllib.parse.urlencode(params).encode('utf-8')

        for attempt in range(self.max_retries):
            self._rate_limit()
            try:
                req = urllib.request.Request(
                    self.SUBMIT_URL, data=data,
                    headers={'User-Agent': 'PrimerDesigner/1.0'}
                )
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    result_url = resp.geturl()
                    html = resp.read().decode('utf-8', errors='replace')

                    # Immediate results?
                    if 'Primer pair specificity checking results' in html:
                        return html

                    # Extract job key for polling
                    m = re.search(r'job_key=([A-Za-z0-9_-]+)', html)
                    if m:
                        return f"{self.SUBMIT_URL}?job_key={m.group(1)}"

                    # Fallback: use redirect URL
                    return result_url

            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(5 * (attempt + 1))
                elif e.code >= 500:
                    time.sleep(3 * (attempt + 1))
                else:
                    self.log_warning(f"Primer-BLAST HTTP {e.code}")
                    return None
            except urllib.error.URLError:
                time.sleep(3 * (attempt + 1))

        return None

    # ------------------------------------------------------------------
    #  Polling
    # ------------------------------------------------------------------

    def _poll_for_results(
        self, url: str, progress_callback: Optional[Callable] = None
    ) -> Optional[str]:
        """Poll until results ready. Returns HTML or None."""
        for attempt in range(self.MAX_POLL_ATTEMPTS):
            self._rate_limit()
            try:
                req = urllib.request.Request(
                    url, headers={'User-Agent': 'PrimerDesigner/1.0'}
                )
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    html = resp.read().decode('utf-8', errors='replace')

                    if 'Primer pair specificity checking results' in html:
                        return html

                    if 'No substantial similarity' in html:
                        return html

                    if 'No target templates were found' in html:
                        return html

                    # Check for real Primer-BLAST errors — use specific phrases
                    # that only appear in actual error pages, not in generic HTML
                    if any(phrase in html for phrase in (
                        'Cannot process the request',
                        'Primer-BLAST Error',
                        'Unable to perform specificity',
                        'Search of the input primer pair',
                    )):
                        return html

                    if progress_callback:
                        elapsed = (attempt + 1) * self.POLL_INTERVAL
                        max_time = self.MAX_POLL_ATTEMPTS * self.POLL_INTERVAL
                        progress_callback(
                            f"Primer-BLAST: waiting for NCBI server... "
                            f"({elapsed}s / max {max_time}s)"
                        )
                    time.sleep(self.POLL_INTERVAL)

            except Exception as e:
                self.log_warning(f"Primer-BLAST poll #{attempt + 1} failed: {e}")
                time.sleep(self.POLL_INTERVAL)

        return None

    # ------------------------------------------------------------------
    #  Result parsing
    # ------------------------------------------------------------------

    def _parse_results(
        self,
        html: str,
        expected_chromosome: Optional[str] = None,
        expected_start: Optional[int] = None,
        expected_end: Optional[int] = None,
    ) -> PrimerBlastResult:
        """Parse Primer-BLAST HTML and extract off-target products."""
        # Check for NCBI error pages first.
        # Note: "No target templates were found" is NOT an error — it means
        # no products were amplified, which is a valid (specific) result.
        # It's handled below as "no hits".
        _error_phrases = (
            'Cannot process the request',
            'Primer-BLAST Error',
            'Unable to perform specificity',
            'Search of the input primer pair',
        )
        for phrase in _error_phrases:
            if phrase in html:
                # Extract a more specific error message if possible
                error_detail = phrase
                m = re.search(r'<div class="error[^"]*"[^>]*>([^<]+)', html)
                if m:
                    error_detail = m.group(1).strip()
                self.log_warning(f"Primer-BLAST returned error: {error_detail}")
                return PrimerBlastResult(
                    success=False,
                    is_specific=False,
                    error=f"Primer-BLAST: {error_detail}",
                )

        # No hits — highly specific.
        # "No substantial similarity" = standard BLAST no-hit message.
        # "No target templates were found" = Custom DB mode, primers don't
        # amplify any product in the assembly (also means no off-targets).
        if ('No substantial similarity' in html
                or 'No target templates were found' in html):
            return PrimerBlastResult(
                success=True, is_specific=True,
                warnings=["Primer-BLAST: No off-target products found"]
            )

        hits: List[PrimerBlastHit] = []
        warnings: List[str] = []

        # ---- parse product blocks ----
        # Primer-BLAST lists products as:
        #   >NC_000017.11 Homo sapiens chromosome 17, GRCh38...
        #   product length = 358
        #   Forward primer  1   ATCGATCG...  20
        #   Template        ... ..........   ...
        #   Reverse primer  1   ATCGATCG...  20
        #   Template        ... ..........   ...
        product_pattern = re.compile(
            r'>([\w.]+)\s+([^\n]+?)\s*\n'           # accession + description
            r'.*?product length\s*=\s*(\d+)',        # product length
            re.DOTALL
        )

        for m in product_pattern.finditer(html):
            accession = m.group(1)
            description = m.group(2).strip()[:120]
            size = int(m.group(3))

            # Parse chromosome from accession and description
            chromosome = self._parse_chromosome(accession, description)

            # Use a larger block to ensure we capture the full alignment
            block_end = min(m.start() + 1500, len(html))
            block = html[m.start():block_end]
            fwd_mm = self._count_mismatches(block, 'Forward')
            rev_mm = self._count_mismatches(block, 'Reverse')

            # Total mismatches: sum if both known, else unknown
            total_mm = -1
            if fwd_mm >= 0 and rev_mm >= 0:
                total_mm = fwd_mm + rev_mm

            hits.append(PrimerBlastHit(
                accession=accession,
                description=description,
                chromosome=chromosome,
                product_size=size,
                forward_mismatches=fwd_mm,
                reverse_mismatches=rev_mm,
                total_mismatches=total_mm,
            ))

        # ---- Mark intended target ----
        # Strategy: use expected coordinates when available, fall back to
        # first perfect-match hit
        intended_found = False

        if expected_chromosome and hits:
            norm_exp = str(expected_chromosome).replace('chr', '')
            for h in hits:
                if not h.chromosome:
                    continue
                norm_hit = h.chromosome.replace('chr', '')
                if norm_hit == norm_exp:
                    # Chromosome matches — check if mismatch count is consistent
                    # with being the intended target (0 mismatches or unknown)
                    if h.total_mismatches <= 0:  # 0 = perfect, -1 = unknown
                        h.is_intended_target = True
                        intended_found = True
                        break

        # Fallback: first hit with zero mismatches in both primers
        if not intended_found:
            for h in hits:
                # Only treat as intended if both are genuinely 0, not unknown (-1)
                if h.forward_mismatches == 0 and h.reverse_mismatches == 0:
                    h.is_intended_target = True
                    intended_found = True
                    break

        # Last resort: if we have hits but none matched above (e.g. all mismatches
        # were unknown=-1), mark the first hit as intended
        if not intended_found and hits:
            hits[0].is_intended_target = True

        off_targets = [h for h in hits if not h.is_intended_target]

        if off_targets:
            warnings.append(
                f"Primer-BLAST: {len(off_targets)} off-target product(s) found"
            )
            for ot in off_targets[:5]:
                mm_info = self._format_mismatch_info(ot)
                warnings.append(
                    f"  {ot.accession} ({ot.description[:60]}), "
                    f"{ot.product_size}bp{mm_info}"
                )
        else:
            warnings.append("Primer-BLAST: No off-target products found")

        return PrimerBlastResult(
            success=True,
            is_specific=len(off_targets) == 0,
            hits=hits,
            off_target_count=len(off_targets),
            warnings=warnings,
        )

    # ------------------------------------------------------------------
    #  Mismatch counting
    # ------------------------------------------------------------------

    @staticmethod
    def _count_mismatches(block: str, direction: str) -> int:
        """
        Extract mismatch count for a primer direction from a product block.

        Returns:
            Number of mismatches (0 = perfect match), or -1 if parsing failed.
            The distinction matters: 0 means confirmed perfect match, -1 means
            we couldn't determine the mismatch count.
        """
        # Strategy 1: Look for explicit mismatch count in text
        #   e.g. "Forward primer  1 total mismatch" or "2 total mismatches"
        pattern = rf'{direction}\s+primer.*?(\d+)\s+total\s+mismatch'
        m = re.search(pattern, block, re.IGNORECASE | re.DOTALL)
        if m:
            return int(m.group(1))

        # Strategy 2: Look for the alignment block between primer and template.
        # Primer-BLAST shows alignments like:
        #   Forward primer  1   ATCGATCG  20
        #   Template        58  ........  78
        #
        # Perfect match: all dots/pipes, no spaces/letters between
        # Mismatch: spaces or different chars in alignment line
        #
        # Look for the primer line followed by the template line
        align_pattern = rf'{direction}\s+primer\s+\d+\s+([A-Za-z]+)\s+\d+\s*\n\s+(\S+)'
        am = re.search(align_pattern, block, re.IGNORECASE)
        if am:
            primer_seq = am.group(1)
            alignment_or_template = am.group(2)
            # If next line is shorter — it might be the alignment indicator
            # Count non-matching positions
            mismatches = 0
            for i, c in enumerate(alignment_or_template):
                if i < len(primer_seq) and c != '|' and c != '.':
                    # A space or different letter in alignment = mismatch
                    if c == ' ' or c != primer_seq[i]:
                        mismatches += 1
            if mismatches > 0:
                return mismatches

        # Strategy 3: If we see the direction keyword and then "on template"
        # with no mismatch info → probably 0 mismatches (perfect match display)
        perfect_pattern = rf'{direction}\s+primer.*?on\s+template'
        if re.search(perfect_pattern, block, re.IGNORECASE | re.DOTALL):
            return 0

        # Cannot determine — return -1 (unknown) rather than 0
        # to avoid falsely marking a hit as "perfect match"
        return -1

    # ------------------------------------------------------------------
    #  Chromosome parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_chromosome(accession: str, description: str) -> str:
        """
        Parse chromosome name from a Primer-BLAST hit.

        Uses two strategies:
        1. Map known RefSeq chromosome accessions (NC_000001..NC_000024)
        2. Parse "chromosome N" from the description text

        Returns:
            Chromosome name (e.g. "1", "17", "X") or empty string if unknown.
        """
        # Strategy 1: RefSeq accession lookup (most reliable)
        # Strip version suffix for lookup: "NC_000017.11" -> "NC_000017"
        acc_base = accession.rsplit('.', 1)[0] if '.' in accession else accession
        if acc_base in _REFSEQ_CHR_MAP:
            return _REFSEQ_CHR_MAP[acc_base]
        # Try with full accession (for GRCh37 versioned entries)
        if accession in _REFSEQ_CHR_MAP:
            return _REFSEQ_CHR_MAP[accession]

        # Strategy 2: Parse from description
        # e.g. "Homo sapiens chromosome 17, GRCh38.p14..."
        m = re.search(r'chromosome\s+(\d+|X|Y|MT)', description, re.IGNORECASE)
        if m:
            return m.group(1)

        # Strategy 3: NW_ / NT_ scaffold accessions — try description
        m = re.search(r'chr(\d+|X|Y|M)', description, re.IGNORECASE)
        if m:
            return m.group(1)

        return ''

    # ------------------------------------------------------------------
    #  Formatting helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _format_mismatch_info(hit: PrimerBlastHit) -> str:
        """Format mismatch info for display, handling unknown (-1) values."""
        parts = []
        fwd = '?' if hit.forward_mismatches < 0 else str(hit.forward_mismatches)
        rev = '?' if hit.reverse_mismatches < 0 else str(hit.reverse_mismatches)
        return f", F:{fwd}mm R:{rev}mm"

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    def _rate_limit(self):
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

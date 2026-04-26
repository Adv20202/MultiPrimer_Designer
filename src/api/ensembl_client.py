"""
Ensembl REST API client for fetching transcript, variant, and population data.
Provides access to transcript information, genomic coordinates, and population frequencies.
"""

import json
import time
import urllib.request
import urllib.parse
import urllib.error
from typing import Optional, Dict, List, Any
from dataclasses import dataclass

from ..core.models import (
    Transcript, Exon, Strand, GenomicPosition, PopulationVariant, MANEType
)
from ..utils.logger import LoggerMixin
from ..utils.config import get_config
from ..utils.thread_safe import ThreadSafeCache, get_ensembl_rate_limiter


@dataclass
class EnsemblResponse:
    """Response from Ensembl API."""
    success: bool
    data: Any = None
    error: str = ""


class EnsemblClient(LoggerMixin):
    """
    Client for accessing Ensembl REST API.
    Fetches transcript info, genomic coordinates, and population variant data.
    """

    # Ensembl REST API endpoints
    BASE_URL = "https://rest.ensembl.org"
    GRCH37_URL = "https://grch37.rest.ensembl.org"  # For HG19

    def __init__(self, assembly: str = "GRCh38"):
        """
        Initialize Ensembl client.

        Args:
            assembly: Reference assembly (GRCh38 or GRCh37)
        """
        self.assembly = assembly
        self.base_url = self.BASE_URL if assembly == "GRCh38" else self.GRCH37_URL

        config = get_config()
        self.timeout = config.api.request_timeout
        self.max_retries = config.api.max_retries
        self.retry_delay = config.api.retry_delay

        # Thread-safe cache for API responses
        self._cache = ThreadSafeCache()

        # Thread-safe rate limiter shared with VariantDBClient
        # (both hit rest.ensembl.org, 15 req/s combined limit)
        self._rate_limiter = get_ensembl_rate_limiter()

    def set_assembly(self, assembly: str):
        """
        Change the reference assembly.

        Args:
            assembly: GRCh38 or GRCh37
        """
        self.assembly = assembly
        self.base_url = self.BASE_URL if assembly == "GRCh38" else self.GRCH37_URL
        self._cache.clear()
        self.log_info(f"Assembly set to {assembly}")

    def _make_request(self, endpoint: str, params: Dict[str, str] = None,
                      content_type: str = "application/json") -> EnsemblResponse:
        """
        Make an HTTP request to Ensembl API.

        Args:
            endpoint: API endpoint path
            params: Query parameters
            content_type: Expected content type

        Returns:
            EnsemblResponse with data or error
        """
        self._rate_limiter.wait()

        url = f"{self.base_url}{endpoint}"
        if params:
            url = f"{url}?{urllib.parse.urlencode(params)}"

        for attempt in range(self.max_retries):
            try:
                self.log_debug(f"Ensembl request: {url}")

                request = urllib.request.Request(
                    url,
                    headers={
                        'Content-Type': content_type,
                        'Accept': content_type
                    }
                )

                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    data = response.read().decode('utf-8')

                    if content_type == "application/json":
                        return EnsemblResponse(success=True, data=json.loads(data))
                    else:
                        return EnsemblResponse(success=True, data=data)

            except urllib.error.HTTPError as e:
                self.log_warning(f"Ensembl HTTP error (attempt {attempt + 1}): {e.code}")
                if e.code == 429:  # Rate limited
                    retry_after = int(e.headers.get('Retry-After', self.retry_delay * 5))
                    time.sleep(retry_after)
                elif e.code == 400:  # Bad request
                    return EnsemblResponse(success=False, error=f"Bad request: {e.reason}")
                elif e.code == 404:  # Not found
                    return EnsemblResponse(success=False, error="Not found")
                elif e.code >= 500:  # Server error
                    time.sleep(self.retry_delay * (attempt + 1))
                else:
                    return EnsemblResponse(success=False, error=f"HTTP {e.code}: {e.reason}")

            except urllib.error.URLError as e:
                self.log_warning(f"Ensembl URL error (attempt {attempt + 1}): {e.reason}")
                time.sleep(self.retry_delay * (attempt + 1))

            except Exception as e:
                self.log_exception(f"Ensembl request error: {e}")
                return EnsemblResponse(success=False, error=str(e))

        return EnsemblResponse(success=False, error="Max retries exceeded")

    def get_transcript_info(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get transcript information from Ensembl using RefSeq accession.

        Args:
            accession: RefSeq transcript accession (e.g., NM_002485.4)

        Returns:
            Dictionary with transcript info or None
        """
        cache_key = f"transcript_{accession}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        # Remove version for lookup if needed
        base_accession = accession.split('.')[0]

        # Use xrefs/id endpoint to find Ensembl ID from RefSeq accession
        # This is the correct endpoint for external database IDs like RefSeq
        response = self._make_request(f"/xrefs/id/{base_accession}",
                                      params={'all_levels': '1'})

        ensembl_id = None

        if response.success and response.data:
            # Look for Ensembl transcript ID in xrefs
            for entry in response.data:
                if entry.get('type') == 'transcript' and entry.get('id', '').startswith('ENST'):
                    ensembl_id = entry.get('id')
                    break

        if not ensembl_id:
            # Try with full accession including version
            response = self._make_request(f"/xrefs/id/{accession}",
                                          params={'all_levels': '1'})
            if response.success and response.data:
                for entry in response.data:
                    if entry.get('type') == 'transcript' and entry.get('id', '').startswith('ENST'):
                        ensembl_id = entry.get('id')
                        break

        if not ensembl_id:
            # Try using xrefs/name endpoint as alternative
            response = self._make_request(f"/xrefs/name/human/{base_accession}")
            if response.success and response.data:
                for entry in response.data:
                    if entry.get('type') == 'transcript':
                        ensembl_id = entry.get('id')
                        break

        if not ensembl_id:
            # Try /xrefs/symbol/ — works well for RefSeq accessions
            response = self._make_request(
                f"/xrefs/symbol/homo_sapiens/{base_accession}"
            )
            if response.success and response.data:
                for entry in response.data:
                    eid = entry.get('id', '')
                    if eid.startswith('ENST'):
                        ensembl_id = eid
                        break

        if not ensembl_id:
            self.log_warning(f"No Ensembl transcript found for {accession}")
            return None

        # Get full transcript info
        response = self._make_request(f"/lookup/id/{ensembl_id}",
                                      params={'expand': '1', 'utr': '1'})

        if not response.success:
            return None

        data = response.data
        info = {
            'ensembl_id': ensembl_id,
            'accession': accession,
            'gene_symbol': data.get('display_name', '').split('-')[0] if data.get('display_name') else '',
            'gene_id': data.get('Parent', ''),
            'chromosome': data.get('seq_region_name', ''),
            'start': data.get('start', 0),
            'end': data.get('end', 0),
            'strand': '+' if data.get('strand', 1) == 1 else '-',
            'biotype': data.get('biotype', ''),
            'is_canonical': data.get('is_canonical', False)
        }

        self._cache[cache_key] = info
        return info

    def get_transcript_sequence(self, ensembl_id: str = None,
                                 accession: str = None) -> Optional[str]:
        """
        Get transcript CDS sequence.

        Args:
            ensembl_id: Ensembl transcript ID
            accession: RefSeq accession (alternative)

        Returns:
            CDS sequence string or None
        """
        if not ensembl_id and accession:
            info = self.get_transcript_info(accession)
            if info:
                ensembl_id = info.get('ensembl_id')

        if not ensembl_id:
            return None

        response = self._make_request(
            f"/sequence/id/{ensembl_id}",
            params={'type': 'cds'},
            content_type="text/plain"
        )

        if response.success:
            return response.data.strip().upper()
        return None

    def get_transcript_exons(self, ensembl_id: str = None,
                             accession: str = None) -> List[Exon]:
        """
        Get exon coordinates for a transcript.

        Args:
            ensembl_id: Ensembl transcript ID
            accession: RefSeq accession (alternative)

        Returns:
            List of Exon objects
        """
        if not ensembl_id and accession:
            info = self.get_transcript_info(accession)
            if info:
                ensembl_id = info.get('ensembl_id')

        if not ensembl_id:
            return []

        response = self._make_request(f"/lookup/id/{ensembl_id}",
                                      params={'expand': '1'})

        if not response.success or not response.data:
            return []

        exons = []
        exon_data = response.data.get('Exon', [])

        for i, ex in enumerate(exon_data, 1):
            exon = Exon(
                number=i,
                genomic_start=ex.get('start', 0),
                genomic_end=ex.get('end', 0)
            )
            exons.append(exon)

        return exons

    def get_genomic_sequence(self, chromosome: str, start: int, end: int,
                             strand: int = 1) -> Optional[str]:
        """
        Get genomic sequence for a region.

        Args:
            chromosome: Chromosome name (e.g., '1', 'X')
            start: Start position
            end: End position
            strand: 1 for forward, -1 for reverse

        Returns:
            Sequence string or None
        """
        # Validate coordinates
        if start > end:
            start, end = end, start

        response = self._make_request(
            f"/sequence/region/human/{chromosome}:{start}..{end}:{strand}",
            content_type="text/plain"
        )

        if response.success:
            return response.data.strip().upper()
        return None

    def map_cds_to_genomic(self, accession: str, cds_start: int,
                           cds_end: int) -> Optional[GenomicPosition]:
        """
        Map CDS coordinates to genomic coordinates.

        Args:
            accession: RefSeq transcript accession
            cds_start: CDS start position (1-based)
            cds_end: CDS end position (1-based)

        Returns:
            GenomicPosition or None
        """
        info = self.get_transcript_info(accession)
        if not info or not info.get('ensembl_id'):
            return None

        ensembl_id = info['ensembl_id']

        # Use map endpoint for coordinate translation
        response = self._make_request(
            f"/map/cds/{ensembl_id}/{cds_start}..{cds_end}"
        )

        if not response.success or not response.data:
            return None

        mappings = response.data.get('mappings', [])
        if not mappings:
            return None

        mapping = mappings[0]
        return GenomicPosition(
            chromosome=mapping.get('seq_region_name', ''),
            start=mapping.get('start', 0),
            end=mapping.get('end', 0),
            strand=Strand.POSITIVE if mapping.get('strand', 1) == 1 else Strand.NEGATIVE,
            assembly=self.assembly
        )

    def get_population_variants(self, chromosome: str, start: int, end: int,
                                populations: List[str] = None) -> List[PopulationVariant]:
        """
        Get population variant data (MAF) for a genomic region.

        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
            populations: List of population codes to filter (optional)

        Returns:
            List of PopulationVariant objects
        """
        # Ensembl variant endpoint
        response = self._make_request(
            f"/overlap/region/human/{chromosome}:{start}-{end}",
            params={'feature': 'variation'}
        )

        if not response.success or not response.data:
            return []

        variants = []
        for var in response.data:
            # Get variant details including MAF
            var_id = var.get('id')
            if not var_id:
                continue

            # Get population frequencies
            pop_response = self._make_request(f"/variation/human/{var_id}",
                                             params={'pops': '1'})

            maf_global = 0.0
            maf_by_pop = {}

            if pop_response.success and pop_response.data:
                pop_data = pop_response.data.get('populations', [])
                for pop in pop_data:
                    pop_name = pop.get('population', '')
                    frequency = pop.get('frequency', 0.0)

                    # Map to standard population codes
                    if 'gnomAD' in pop_name or '1000GENOMES' in pop_name:
                        if 'global' in pop_name.lower() or pop_name == '1000GENOMES:phase_3:ALL':
                            maf_global = max(maf_global, frequency)

                        # Extract population code
                        for code in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'NFE', 'FIN', 'ASJ']:
                            if code in pop_name:
                                maf_by_pop[code.lower()] = max(
                                    maf_by_pop.get(code.lower(), 0.0),
                                    frequency
                                )

            pop_variant = PopulationVariant(
                rsid=var_id,
                chromosome=chromosome,
                position=var.get('start', 0),
                ref=var.get('ref', ''),
                alt=var.get('alt', ''),
                maf_global=maf_global,
                maf_by_population=maf_by_pop
            )
            variants.append(pop_variant)

        return variants

    def get_gene_info(self, gene_symbol: str) -> Optional[Dict[str, Any]]:
        """
        Get gene information from Ensembl.

        Args:
            gene_symbol: Gene symbol (e.g., BRCA1)

        Returns:
            Dictionary with gene info or None
        """
        response = self._make_request(f"/lookup/symbol/homo_sapiens/{gene_symbol}",
                                      params={'expand': '1'})

        if not response.success:
            return None

        data = response.data
        return {
            'ensembl_id': data.get('id', ''),
            'symbol': data.get('display_name', ''),
            'description': data.get('description', ''),
            'chromosome': data.get('seq_region_name', ''),
            'start': data.get('start', 0),
            'end': data.get('end', 0),
            'strand': '+' if data.get('strand', 1) == 1 else '-',
            'biotype': data.get('biotype', '')
        }

    def clear_cache(self):
        """Clear the response cache."""
        self._cache.clear()
        self.log_info("Ensembl client cache cleared")

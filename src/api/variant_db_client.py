"""
Client for variant databases (gnomAD, Ensembl) to fetch population allele frequencies.
Used for filtering primer positions based on known polymorphisms.
"""

import json
import time
import urllib.request
import urllib.parse
import urllib.error
from typing import Optional, Dict, List, Any
from dataclasses import dataclass

from ..core.models import PopulationVariant
from ..utils.logger import LoggerMixin
from ..utils.config import get_config
from ..utils.thread_safe import ThreadSafeCache, get_ensembl_rate_limiter


class VariantDBClient(LoggerMixin):
    """
    Client for fetching variant population frequency data.
    Primarily uses Ensembl REST API which includes gnomAD and 1000 Genomes data.
    """

    ENSEMBL_URL = "https://rest.ensembl.org"
    ENSEMBL_GRCH37_URL = "https://grch37.rest.ensembl.org"

    # Population code mappings
    POPULATION_CODES = {
        'global': 'Global (all populations)',
        'afr': 'African/African American',
        'ami': 'Amish',
        'amr': 'Latino/Admixed American',
        'asj': 'Ashkenazi Jewish',
        'eas': 'East Asian',
        'fin': 'Finnish',
        'nfe': 'Non-Finnish European',
        'sas': 'South Asian',
        'oth': 'Other',
        'mid': 'Middle Eastern',
        'remaining': 'Remaining individuals'
    }

    def __init__(self, assembly: str = "GRCh38"):
        """
        Initialize variant database client.

        Args:
            assembly: Reference assembly (GRCh38 or GRCh37)
        """
        self.assembly = assembly
        self.base_url = self.ENSEMBL_URL if assembly == "GRCh38" else self.ENSEMBL_GRCH37_URL

        config = get_config()
        self.timeout = config.api.request_timeout
        self.max_retries = config.api.max_retries
        self.retry_delay = config.api.retry_delay

        # Thread-safe cache for API responses
        self._cache = ThreadSafeCache()

        # Thread-safe rate limiter shared with EnsemblClient
        # (both hit rest.ensembl.org, 15 req/s combined limit)
        self._rate_limiter = get_ensembl_rate_limiter()

    def set_assembly(self, assembly: str):
        """Change reference assembly."""
        self.assembly = assembly
        self.base_url = self.ENSEMBL_URL if assembly == "GRCh38" else self.ENSEMBL_GRCH37_URL
        self._cache.clear()

    def _make_request(self, url: str, content_type: str = "application/json") -> Optional[Any]:
        """
        Make HTTP request with rate limiting and retries.
        """
        self._rate_limiter.wait()

        for attempt in range(self.max_retries):
            try:
                self.log_debug(f"VariantDB request: {url}")

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
                        return json.loads(data)
                    return data

            except urllib.error.HTTPError as e:
                if e.code == 429:
                    retry_after = int(e.headers.get('Retry-After', 5))
                    time.sleep(retry_after)
                elif e.code == 400 or e.code == 404:
                    return None
                elif e.code >= 500:
                    time.sleep(self.retry_delay * (attempt + 1))
                else:
                    self.log_warning(f"HTTP error {e.code}: {e.reason}")
                    return None

            except Exception as e:
                self.log_warning(f"Request error (attempt {attempt + 1}): {e}")
                time.sleep(self.retry_delay * (attempt + 1))

        return None

    def get_variants_in_region(
        self,
        chromosome: str,
        start: int,
        end: int,
        populations: List[str] = None,
        maf_threshold: float = 0.0
    ) -> List[PopulationVariant]:
        """
        Get all known variants in a genomic region with their population frequencies.

        Args:
            chromosome: Chromosome (e.g., '1', 'X')
            start: Start position
            end: End position
            populations: List of population codes to include (None = all)
            maf_threshold: Minimum MAF threshold to include

        Returns:
            List of PopulationVariant objects
        """
        # Validate region size (Ensembl has limits)
        region_size = end - start
        if region_size > 5000000:  # 5Mb limit
            self.log_warning(f"Region too large ({region_size}bp), splitting query")
            return self._get_variants_chunked(chromosome, start, end, populations, maf_threshold)

        cache_key = f"variants_{chromosome}:{start}-{end}"
        if cache_key in self._cache:
            cached = self._cache[cache_key]
            return self._filter_variants(cached, populations, maf_threshold)

        # Query Ensembl overlap endpoint for variations
        url = f"{self.base_url}/overlap/region/human/{chromosome}:{start}-{end}?feature=variation"
        data = self._make_request(url)

        if not data:
            return []

        # Collect all rsIDs for batch fetching
        rsids = []
        other_variants = []

        for var in data:
            var_id = var.get('id', '')

            if var_id and var_id.startswith('rs'):
                rsids.append(var_id)
            elif var.get('start'):
                # Include other variants with basic info (no population data)
                allele_string = var.get('allele_string', '/')
                alleles = allele_string.split('/') if allele_string else ['', '']
                basic_variant = PopulationVariant(
                    rsid=var_id or f"var_{var.get('start', 0)}",
                    chromosome=chromosome,
                    position=var.get('start', 0),
                    ref=alleles[0] if alleles else '',
                    alt=alleles[-1] if len(alleles) > 1 else '',
                    maf_global=0.0,
                    maf_by_population={}
                )
                other_variants.append(basic_variant)

        # Fetch all rsIDs in batch (much faster!)
        variants = other_variants
        if rsids:
            self.log_debug(f"Fetching {len(rsids)} variants in batch...")
            batch_variants = self._get_variants_batch(rsids)
            variants.extend(batch_variants)

        self._cache[cache_key] = variants
        return self._filter_variants(variants, populations, maf_threshold)

    def _get_variants_chunked(
        self,
        chromosome: str,
        start: int,
        end: int,
        populations: List[str],
        maf_threshold: float,
        chunk_size: int = 100000
    ) -> List[PopulationVariant]:
        """Get variants in large region by chunking."""
        all_variants = []
        current_start = start

        while current_start < end:
            current_end = min(current_start + chunk_size, end)
            # Call internal method to avoid recursion check
            variants = self._get_variants_in_region_internal(
                chromosome, current_start, current_end, populations, maf_threshold
            )
            all_variants.extend(variants)
            current_start = current_end + 1

        return all_variants

    def _get_variants_in_region_internal(
        self,
        chromosome: str,
        start: int,
        end: int,
        populations: List[str],
        maf_threshold: float
    ) -> List[PopulationVariant]:
        """Internal method that fetches variants without size check (to avoid recursion)."""
        cache_key = f"variants_{chromosome}:{start}-{end}"
        if cache_key in self._cache:
            cached = self._cache[cache_key]
            return self._filter_variants(cached, populations, maf_threshold)

        url = f"{self.base_url}/overlap/region/human/{chromosome}:{start}-{end}?feature=variation"
        data = self._make_request(url)

        if not data:
            return []

        # Collect all rsIDs for batch fetching
        rsids = []
        other_variants = []

        for var in data:
            var_id = var.get('id', '')

            if var_id and var_id.startswith('rs'):
                rsids.append(var_id)
            elif var.get('start'):
                allele_string = var.get('allele_string', '/')
                alleles = allele_string.split('/') if allele_string else []
                basic_variant = PopulationVariant(
                    rsid=var_id or f"var_{var.get('start', 0)}",
                    chromosome=chromosome,
                    position=var.get('start', 0),
                    ref=alleles[0] if len(alleles) > 0 else '',
                    alt=alleles[1] if len(alleles) > 1 else '',
                    maf_global=0.0,
                    maf_by_population={}
                )
                other_variants.append(basic_variant)

        # Fetch all rsIDs in batch (much faster!)
        variants = other_variants
        if rsids:
            self.log_debug(f"Fetching {len(rsids)} variants in batch...")
            batch_variants = self._get_variants_batch(rsids)
            variants.extend(batch_variants)

        self._cache[cache_key] = variants
        return self._filter_variants(variants, populations, maf_threshold)

    def _get_variants_batch(self, rsids: List[str], batch_size: int = 200) -> List[PopulationVariant]:
        """
        Fetch multiple variants in batch using POST endpoint.
        Much faster than individual requests.

        Args:
            rsids: List of rsIDs to fetch
            batch_size: Maximum number of IDs per batch (Ensembl limit is 1000)

        Returns:
            List of PopulationVariant objects
        """
        all_variants = []

        # Process in batches
        for i in range(0, len(rsids), batch_size):
            batch = rsids[i:i + batch_size]

            # Check cache first
            uncached_ids = []
            for rsid in batch:
                cache_key = f"variant_{rsid}"
                if cache_key in self._cache:
                    all_variants.append(self._cache[cache_key])
                else:
                    uncached_ids.append(rsid)

            if not uncached_ids:
                continue

            # POST batch request using urllib
            # Note: pops parameter must be in URL, not in POST body
            url = f"{self.base_url}/variation/human?pops=1"
            try:
                # Prepare POST data - only ids in the body
                post_data = json.dumps({"ids": uncached_ids}).encode('utf-8')

                request = urllib.request.Request(
                    url,
                    data=post_data,
                    headers={
                        'Content-Type': 'application/json',
                        'Accept': 'application/json'
                    },
                    method='POST'
                )

                with urllib.request.urlopen(request, timeout=60) as response:
                    response_data = response.read().decode('utf-8')
                    data = json.loads(response_data)

                for rsid, var_data in data.items():
                    if not var_data:
                        continue

                    # Parse variant (same logic as _get_variant_details)
                    mappings = var_data.get('mappings', [])
                    if not mappings:
                        continue

                    mapping = mappings[0]
                    chromosome = mapping.get('seq_region_name', '')
                    position = mapping.get('start', 0)
                    allele_string = mapping.get('allele_string', '/')

                    alleles = allele_string.split('/')
                    ref = alleles[0] if alleles else ''
                    alt = alleles[1] if len(alleles) > 1 else ''

                    # Parse population frequencies
                    # IMPORTANT: Ensembl returns frequency for EACH allele separately
                    # We need the MINOR allele frequency (the less common one)
                    # Frequency close to 1.0 means reference allele - we want the alternate!
                    maf_global = 0.0
                    maf_by_pop = {}

                    populations = var_data.get('populations', [])
                    self.log_debug(f"Variant {rsid}: found {len(populations)} population entries")

                    for pop in populations:
                        pop_name = pop.get('population', '')
                        frequency = pop.get('frequency', 0.0)
                        allele = pop.get('allele', '')  # The allele this frequency is for

                        if not frequency:
                            continue

                        # CRITICAL: Skip if this is the reference allele frequency
                        # We want the ALTERNATE allele frequency (MAF)
                        # If frequency > 0.5, this is likely the reference allele
                        # Also skip if allele matches ref
                        if allele and allele.upper() == ref.upper():
                            continue  # Skip reference allele frequencies

                        # If no allele specified but frequency is very high (>0.5),
                        # it's likely the major (reference) allele
                        if not allele and frequency > 0.5:
                            continue

                        pop_lower = pop_name.lower()

                        # Global frequencies - check multiple patterns
                        if any(x in pop_lower for x in ['global', ':all', 'all populations', 'total']):
                            maf_global = max(maf_global, frequency)

                        # gnomAD populations (v3, v4)
                        if 'gnomad' in pop_lower or 'gnomad_' in pop_lower:
                            if ':afr' in pop_lower or '_afr' in pop_lower:
                                maf_by_pop['afr'] = max(maf_by_pop.get('afr', 0.0), frequency)
                            elif ':amr' in pop_lower or '_amr' in pop_lower:
                                maf_by_pop['amr'] = max(maf_by_pop.get('amr', 0.0), frequency)
                            elif ':asj' in pop_lower or '_asj' in pop_lower:
                                maf_by_pop['asj'] = max(maf_by_pop.get('asj', 0.0), frequency)
                            elif ':eas' in pop_lower or '_eas' in pop_lower:
                                maf_by_pop['eas'] = max(maf_by_pop.get('eas', 0.0), frequency)
                            elif ':fin' in pop_lower or '_fin' in pop_lower:
                                maf_by_pop['fin'] = max(maf_by_pop.get('fin', 0.0), frequency)
                            elif ':nfe' in pop_lower or '_nfe' in pop_lower:
                                maf_by_pop['nfe'] = max(maf_by_pop.get('nfe', 0.0), frequency)
                            elif ':sas' in pop_lower or '_sas' in pop_lower:
                                maf_by_pop['sas'] = max(maf_by_pop.get('sas', 0.0), frequency)

                        # 1000 Genomes populations
                        elif '1000genomes' in pop_lower:
                            if ':afr' in pop_lower:
                                maf_by_pop['afr'] = max(maf_by_pop.get('afr', 0.0), frequency)
                            elif ':amr' in pop_lower:
                                maf_by_pop['amr'] = max(maf_by_pop.get('amr', 0.0), frequency)
                            elif ':eas' in pop_lower:
                                maf_by_pop['eas'] = max(maf_by_pop.get('eas', 0.0), frequency)
                            elif ':eur' in pop_lower:
                                maf_by_pop['nfe'] = max(maf_by_pop.get('nfe', 0.0), frequency)
                            elif ':sas' in pop_lower:
                                maf_by_pop['sas'] = max(maf_by_pop.get('sas', 0.0), frequency)
                            elif ':all' in pop_lower:
                                maf_global = max(maf_global, frequency)

                    # If no global MAF was found, use the max of population MAFs
                    if maf_global == 0.0 and maf_by_pop:
                        maf_global = max(maf_by_pop.values())

                    variant = PopulationVariant(
                        rsid=rsid,
                        chromosome=chromosome,
                        position=position,
                        ref=ref,
                        alt=alt,
                        maf_global=maf_global,
                        maf_by_population=maf_by_pop
                    )

                    # Cache it
                    cache_key = f"variant_{rsid}"
                    self._cache[cache_key] = variant
                    all_variants.append(variant)

            except Exception as e:
                self.log_warning(f"Batch request error: {e}, falling back to individual requests")
                # Fallback to individual requests
                for rsid in uncached_ids:
                    var = self._get_variant_details(rsid)
                    if var:
                        all_variants.append(var)

        return all_variants

    def _get_variant_details(self, rsid: str) -> Optional[PopulationVariant]:
        """
        Get detailed variant information including population frequencies.

        Args:
            rsid: rs ID (e.g., rs12345)

        Returns:
            PopulationVariant or None
        """
        cache_key = f"variant_{rsid}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        url = f"{self.base_url}/variation/human/{rsid}?pops=1"
        data = self._make_request(url)

        if not data:
            return None

        # Parse mappings for position
        mappings = data.get('mappings', [])
        if not mappings:
            return None

        # Use first mapping (usually the primary)
        mapping = mappings[0]
        chromosome = mapping.get('seq_region_name', '')
        position = mapping.get('start', 0)
        allele_string = mapping.get('allele_string', '/')

        alleles = allele_string.split('/')
        ref = alleles[0] if alleles else ''
        alt = alleles[1] if len(alleles) > 1 else ''

        # Parse population frequencies
        # IMPORTANT: Ensembl returns frequency for EACH allele separately
        # We need the MINOR allele frequency (the less common one)
        maf_global = 0.0
        maf_by_pop = {}

        populations = data.get('populations', [])
        self.log_debug(f"Variant {rsid}: found {len(populations)} population entries")

        for pop in populations:
            pop_name = pop.get('population', '')
            frequency = pop.get('frequency', 0.0)
            allele = pop.get('allele', '')  # The allele this frequency is for

            if not frequency:
                continue

            # CRITICAL: Skip if this is the reference allele frequency
            # We want the ALTERNATE allele frequency (MAF)
            if allele and allele.upper() == ref.upper():
                continue  # Skip reference allele frequencies

            # If no allele specified but frequency is very high (>0.5),
            # it's likely the major (reference) allele
            if not allele and frequency > 0.5:
                continue

            # Identify population type
            pop_lower = pop_name.lower()

            # Global frequencies - check multiple patterns
            if any(x in pop_lower for x in ['global', ':all', 'all populations', 'total']):
                maf_global = max(maf_global, frequency)

            # gnomAD populations (v3, v4)
            if 'gnomad' in pop_lower or 'gnomad_' in pop_lower:
                if ':afr' in pop_lower or '_afr' in pop_lower or 'afr' in pop_lower:
                    maf_by_pop['afr'] = max(maf_by_pop.get('afr', 0.0), frequency)
                elif ':amr' in pop_lower or '_amr' in pop_lower:
                    maf_by_pop['amr'] = max(maf_by_pop.get('amr', 0.0), frequency)
                elif ':asj' in pop_lower or '_asj' in pop_lower:
                    maf_by_pop['asj'] = max(maf_by_pop.get('asj', 0.0), frequency)
                elif ':eas' in pop_lower or '_eas' in pop_lower:
                    maf_by_pop['eas'] = max(maf_by_pop.get('eas', 0.0), frequency)
                elif ':fin' in pop_lower or '_fin' in pop_lower:
                    maf_by_pop['fin'] = max(maf_by_pop.get('fin', 0.0), frequency)
                elif ':nfe' in pop_lower or '_nfe' in pop_lower or 'nf_european' in pop_lower:
                    maf_by_pop['nfe'] = max(maf_by_pop.get('nfe', 0.0), frequency)
                elif ':sas' in pop_lower or '_sas' in pop_lower or 'south_asian' in pop_lower:
                    maf_by_pop['sas'] = max(maf_by_pop.get('sas', 0.0), frequency)

            # 1000 Genomes populations
            elif '1000genomes' in pop_lower:
                if ':afr' in pop_lower:
                    maf_by_pop['afr'] = max(maf_by_pop.get('afr', 0.0), frequency)
                elif ':amr' in pop_lower:
                    maf_by_pop['amr'] = max(maf_by_pop.get('amr', 0.0), frequency)
                elif ':eas' in pop_lower:
                    maf_by_pop['eas'] = max(maf_by_pop.get('eas', 0.0), frequency)
                elif ':eur' in pop_lower:
                    maf_by_pop['nfe'] = max(maf_by_pop.get('nfe', 0.0), frequency)
                elif ':sas' in pop_lower:
                    maf_by_pop['sas'] = max(maf_by_pop.get('sas', 0.0), frequency)
                elif ':all' in pop_lower:
                    maf_global = max(maf_global, frequency)

        # If no global MAF was found, use the max of population MAFs
        if maf_global == 0.0 and maf_by_pop:
            maf_global = max(maf_by_pop.values())

        variant = PopulationVariant(
            rsid=rsid,
            chromosome=chromosome,
            position=position,
            ref=ref,
            alt=alt,
            maf_global=maf_global,
            maf_by_population=maf_by_pop
        )

        self._cache[cache_key] = variant
        return variant

    def _filter_variants(
        self,
        variants: List[PopulationVariant],
        populations: List[str] = None,
        maf_threshold: float = 0.0
    ) -> List[PopulationVariant]:
        """
        Filter variants by population and MAF threshold.

        A variant passes the filter if **any** of the following is true:
        - maf_global >= maf_threshold  (always checked — a globally common
          SNP must never be ignored regardless of the selected population)
        - maf for any selected non-global population >= maf_threshold

        Args:
            variants: List of variants to filter
            populations: Population codes to check (None = global only)
            maf_threshold: Minimum MAF threshold

        Returns:
            Filtered list of variants
        """
        if maf_threshold <= 0 and not populations:
            return variants

        filtered = []
        for var in variants:
            # Always check global MAF — a globally common variant must be
            # masked regardless of which population the user selected
            if var.maf_global >= maf_threshold:
                filtered.append(var)
                continue

            # Check specific populations (if any selected besides global)
            if populations:
                for pop in populations:
                    if pop == 'global':
                        continue
                    maf = var.maf_by_population.get(pop, 0.0)
                    if maf >= maf_threshold:
                        filtered.append(var)
                        break

        return filtered

    def get_variant_at_position(
        self,
        chromosome: str,
        position: int
    ) -> List[PopulationVariant]:
        """
        Get variants at a specific genomic position.

        Args:
            chromosome: Chromosome
            position: Position

        Returns:
            List of variants at that position
        """
        return self.get_variants_in_region(chromosome, position, position)

    def check_position_polymorphic(
        self,
        chromosome: str,
        position: int,
        populations: List[str] = None,
        maf_threshold: float = 0.01
    ) -> bool:
        """
        Check if a position has a common polymorphism.

        Args:
            chromosome: Chromosome
            position: Position
            populations: Populations to check
            maf_threshold: MAF threshold for "common"

        Returns:
            True if position has common variant
        """
        variants = self.get_variants_in_region(
            chromosome, position, position, populations, maf_threshold
        )
        return len(variants) > 0

    def get_population_labels(self) -> Dict[str, str]:
        """Get human-readable labels for population codes."""
        return self.POPULATION_CODES.copy()

    def clear_cache(self):
        """Clear the response cache."""
        self._cache.clear()
        self.log_info("VariantDB client cache cleared")

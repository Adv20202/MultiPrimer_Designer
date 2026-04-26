"""
NCBI API client for fetching transcript and sequence data.
Uses E-utilities for accessing NCBI databases.
"""

import re
import time
import xml.etree.ElementTree as ET
from typing import Optional, Dict, List, Any
from dataclasses import dataclass
import urllib.request
import urllib.parse
import urllib.error
import json

from ..core.models import Transcript, Exon, Strand, GenomicPosition
from ..utils.logger import LoggerMixin
from ..utils.config import get_config
from ..utils.thread_safe import ThreadSafeCache, ThreadSafeRateLimiter


@dataclass
class NCBIResponse:
    """Response from NCBI API."""
    success: bool
    data: Any = None
    error: str = ""


class NCBIClient(LoggerMixin):
    """
    Client for accessing NCBI E-utilities API.
    Fetches transcript information, sequences, and genomic coordinates.
    """

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    EFETCH_URL = f"{BASE_URL}/efetch.fcgi"
    ESEARCH_URL = f"{BASE_URL}/esearch.fcgi"
    EINFO_URL = f"{BASE_URL}/einfo.fcgi"

    def __init__(self, api_key: str = "", email: str = "", assembly: str = "GRCh38"):
        """
        Initialize NCBI client.

        Args:
            api_key: NCBI API key (optional, but recommended)
            email: Email for NCBI tracking (required by NCBI policy)
            assembly: Reference assembly (GRCh38 or GRCh37)
        """
        config = get_config()
        self.api_key = api_key or config.api.ncbi_api_key
        self.email = email or config.api.ncbi_email
        self.timeout = config.api.request_timeout
        self.max_retries = config.api.max_retries
        self.retry_delay = config.api.retry_delay
        self.assembly = assembly

        # Thread-safe cache for API responses
        self._cache = ThreadSafeCache()

        # Thread-safe rate limiter for NCBI
        self._rate_limiter = ThreadSafeRateLimiter(
            0.34 if not self.api_key else 0.1
        )

    def set_assembly(self, assembly: str):
        """Change reference assembly and clear cache."""
        self.assembly = assembly
        self._cache.clear()

    def _make_request(self, url: str, params: Dict[str, str]) -> NCBIResponse:
        """
        Make an HTTP request to NCBI API with rate limiting and retries.

        Args:
            url: API endpoint URL
            params: Query parameters

        Returns:
            NCBIResponse with data or error
        """
        # Add API key and email to params
        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email

        self._rate_limiter.wait()

        full_url = f"{url}?{urllib.parse.urlencode(params)}"

        for attempt in range(self.max_retries):
            try:
                self.log_debug(f"NCBI request: {full_url}")

                request = urllib.request.Request(
                    full_url,
                    headers={'User-Agent': 'PrimerDesigner/1.0'}
                )

                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    data = response.read().decode('utf-8')
                    return NCBIResponse(success=True, data=data)

            except urllib.error.HTTPError as e:
                self.log_warning(f"NCBI HTTP error (attempt {attempt + 1}): {e.code} - {e.reason}")
                if e.code == 429:  # Rate limited
                    time.sleep(self.retry_delay * (attempt + 1) * 2)
                elif e.code >= 500:  # Server error, retry
                    time.sleep(self.retry_delay * (attempt + 1))
                else:
                    return NCBIResponse(success=False, error=f"HTTP {e.code}: {e.reason}")

            except urllib.error.URLError as e:
                self.log_warning(f"NCBI URL error (attempt {attempt + 1}): {e.reason}")
                time.sleep(self.retry_delay * (attempt + 1))

            except Exception as e:
                self.log_exception(f"NCBI request error: {e}")
                return NCBIResponse(success=False, error=str(e))

        return NCBIResponse(success=False, error="Max retries exceeded")

    def get_transcript_info(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get transcript information from NCBI.

        Args:
            accession: RefSeq transcript accession (e.g., NM_002485.4)

        Returns:
            Dictionary with transcript info or None
        """
        # Check cache
        cache_key = f"transcript_info_{accession}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        # Search for the accession in nucleotide database
        search_params = {
            'db': 'nucleotide',
            'term': accession,
            'retmode': 'json'
        }

        response = self._make_request(self.ESEARCH_URL, search_params)
        if not response.success:
            self.log_error(f"Failed to search for transcript {accession}: {response.error}")
            return None

        try:
            search_result = json.loads(response.data)
            id_list = search_result.get('esearchresult', {}).get('idlist', [])

            if not id_list:
                self.log_warning(f"No results found for transcript {accession}")
                return None

            # Fetch the record
            fetch_params = {
                'db': 'nucleotide',
                'id': id_list[0],
                'rettype': 'gb',
                'retmode': 'xml'
            }

            fetch_response = self._make_request(self.EFETCH_URL, fetch_params)
            if not fetch_response.success:
                return None

            # Parse the GenBank XML
            info = self._parse_genbank_xml(fetch_response.data)
            if info:
                self._cache[cache_key] = info

            return info

        except Exception as e:
            self.log_exception(f"Error parsing NCBI response for {accession}: {e}")
            return None

    def _parse_genbank_xml(self, xml_data: str) -> Optional[Dict[str, Any]]:
        """
        Parse GenBank XML response.

        Args:
            xml_data: XML string from NCBI

        Returns:
            Dictionary with parsed transcript information
        """
        try:
            root = ET.fromstring(xml_data)

            # Navigate to GBSeq element
            gbseq = root.find('.//GBSeq')
            if gbseq is None:
                return None

            info = {
                'accession': '',
                'gene_symbol': '',
                'gene_id': '',
                'organism': '',
                'sequence': '',
                'length': 0,
                'chromosome': '',
                'strand': '+',
                'exons': [],
                'cds_start': 0,
                'cds_end': 0
            }

            # Basic info
            accession_elem = gbseq.find('GBSeq_accession-version')
            if accession_elem is not None and accession_elem.text:
                info['accession'] = accession_elem.text

            organism_elem = gbseq.find('GBSeq_organism')
            if organism_elem is not None and organism_elem.text:
                info['organism'] = organism_elem.text

            sequence_elem = gbseq.find('GBSeq_sequence')
            if sequence_elem is not None and sequence_elem.text:
                info['sequence'] = sequence_elem.text.upper()
                info['length'] = len(info['sequence'])

            # Parse features for gene info and CDS
            features = gbseq.findall('.//GBFeature')
            for feature in features:
                feature_key = feature.find('GBFeature_key')
                if feature_key is None:
                    continue

                key = feature_key.text

                if key == 'gene':
                    quals = feature.findall('.//GBQualifier')
                    for qual in quals:
                        name = qual.find('GBQualifier_name')
                        value = qual.find('GBQualifier_value')
                        if name is not None and value is not None:
                            if name.text == 'gene':
                                info['gene_symbol'] = value.text
                            elif name.text == 'db_xref' and 'GeneID:' in value.text:
                                info['gene_id'] = value.text.replace('GeneID:', '')

                elif key == 'CDS':
                    location = feature.find('GBFeature_location')
                    if location is not None and location.text:
                        # Parse CDS location (e.g., "1..1500" or "join(1..100,200..300)")
                        loc_text = location.text

                        # Handle join() syntax for multi-exon CDS
                        if loc_text.startswith('join('):
                            # Extract all ranges from join()
                            ranges_text = loc_text[5:-1]  # Remove "join(" and ")"
                            ranges = re.findall(r'(\d+)\.\.(\d+)', ranges_text)
                            if ranges:
                                # CDS start is from first range, end is from last range
                                info['cds_start'] = int(ranges[0][0])
                                info['cds_end'] = int(ranges[-1][1])
                        else:
                            # Simple case: single range
                            match = re.match(r'(\d+)\.\.(\d+)', loc_text)
                            if match:
                                info['cds_start'] = int(match.group(1))
                                info['cds_end'] = int(match.group(2))

            return info

        except ET.ParseError as e:
            self.log_error(f"XML parse error: {e}")
            return None

    def get_sequence(self, accession: str, start: int = None, end: int = None) -> Optional[str]:
        """
        Get nucleotide sequence from NCBI.

        Args:
            accession: Accession number
            start: Start position (1-based, optional)
            end: End position (1-based, optional)

        Returns:
            Sequence string or None
        """
        params = {
            'db': 'nucleotide',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text'
        }

        if start is not None and end is not None:
            params['seq_start'] = str(start)
            params['seq_stop'] = str(end)

        response = self._make_request(self.EFETCH_URL, params)
        if not response.success:
            return None

        # Parse FASTA format
        lines = response.data.strip().split('\n')
        sequence_lines = [line for line in lines if not line.startswith('>')]
        sequence = ''.join(sequence_lines).upper()

        return sequence if sequence else None

    def get_genomic_coordinates(self, transcript_accession: str) -> Optional[GenomicPosition]:
        """
        Get genomic coordinates for a transcript.

        Args:
            transcript_accession: RefSeq transcript accession

        Returns:
            GenomicPosition or None
        """
        # Use gene database to get location info
        info = self.get_transcript_info(transcript_accession)
        if not info:
            return None

        # If we have gene ID, try to get more detailed location
        if info.get('gene_id'):
            gene_params = {
                'db': 'gene',
                'id': info['gene_id'],
                'retmode': 'xml'
            }

            response = self._make_request(self.EFETCH_URL, gene_params)
            if response.success:
                try:
                    root = ET.fromstring(response.data)

                    # Extract chromosome from Gene-ref
                    chromosome = ''
                    maploc = root.find('.//Gene-ref_maploc')
                    if maploc is not None and maploc.text:
                        # Format like "17q21.31" - extract just the chromosome number/letter
                        chromosome = re.match(r'^(\d+|[XY])', maploc.text)
                        if chromosome:
                            chromosome = chromosome.group(1)

                    # Find genomic location matching the selected assembly
                    for loc in root.findall('.//Gene-commentary'):
                        heading = loc.find('Gene-commentary_heading')
                        if heading is not None and self.assembly in (heading.text or ''):
                            # Extract coordinates
                            interval_from = loc.find('.//Seq-interval_from')
                            interval_to = loc.find('.//Seq-interval_to')
                            strand_elem = loc.find('.//Seq-interval_strand/Na-strand')

                            if interval_from is not None and interval_to is not None:
                                return GenomicPosition(
                                    chromosome=chromosome or 'unknown',
                                    start=int(interval_from.text) + 1,  # NCBI uses 0-based, convert to 1-based
                                    end=int(interval_to.text) + 1,
                                    strand=Strand.NEGATIVE if strand_elem is not None and strand_elem.get('value') == 'minus' else Strand.POSITIVE,
                                    assembly=self.assembly
                                )
                except Exception as e:
                    self.log_debug(f"Could not parse gene location: {e}")

        return None

    def map_cds_to_genomic(self, accession: str, cds_start: int,
                           cds_end: int) -> Optional[GenomicPosition]:
        """
        Map CDS coordinates to genomic coordinates using NCBI data.

        This method converts coding sequence (CDS) positions to genomic coordinates
        by fetching the transcript structure from NCBI and calculating the genomic
        position based on CDS offset.

        Args:
            accession: RefSeq transcript accession (e.g., NM_002485.4)
            cds_start: CDS start position (1-based)
            cds_end: CDS end position (1-based)

        Returns:
            GenomicPosition or None if mapping fails
        """
        # Get transcript info with CDS boundaries
        info = self.get_transcript_info(accession)
        if not info:
            self.log_warning(f"Could not fetch transcript info for {accession}")
            return None

        # Get genomic coordinates for the whole transcript
        genomic_coords = self.get_genomic_coordinates(accession)
        if not genomic_coords:
            self.log_warning(f"Could not fetch genomic coordinates for {accession}")
            return None

        # Extract CDS info from transcript
        transcript_cds_start = info.get('cds_start', 0)
        transcript_cds_end = info.get('cds_end', 0)

        if not transcript_cds_start or not transcript_cds_end:
            self.log_warning(f"No CDS information found for {accession}")
            return None

        # Calculate genomic position based on CDS offset.
        #
        # IMPORTANT LIMITATION: This method uses a linear offset from the
        # transcript start/end, which does NOT account for introns.  It is
        # only accurate for single-exon transcripts.  For multi-exon genes
        # the result will be wrong.  Prefer Ensembl /map/cds/ endpoint or
        # manual exon-based mapping instead.

        strand = genomic_coords.strand
        genomic_span = genomic_coords.end - genomic_coords.start + 1
        cds_length = transcript_cds_end - transcript_cds_start + 1

        # If the genomic span is much larger than the CDS region, the gene
        # contains introns and linear offset will be inaccurate.  Return
        # None so the caller can try a better mapping method.
        if genomic_span > cds_length * 1.5:
            self.log_debug(
                f"NCBI linear mapping skipped for {accession}: "
                f"genomic span ({genomic_span:,} bp) >> CDS length "
                f"({cds_length:,} bp), likely multi-exon gene"
            )
            return None

        offset_start = cds_start - 1  # Convert to 0-based
        offset_end = cds_end - 1

        if strand == Strand.POSITIVE:
            genomic_variant_start = genomic_coords.start + offset_start
            genomic_variant_end = genomic_coords.start + offset_end
        else:
            genomic_variant_start = genomic_coords.end - offset_end
            genomic_variant_end = genomic_coords.end - offset_start

        # Ensure start <= end
        if genomic_variant_start > genomic_variant_end:
            genomic_variant_start, genomic_variant_end = genomic_variant_end, genomic_variant_start

        return GenomicPosition(
            chromosome=genomic_coords.chromosome,
            start=genomic_variant_start,
            end=genomic_variant_end,
            strand=strand,
            assembly=genomic_coords.assembly
        )

    def get_transcript_exons(self, transcript_accession: str) -> List[Exon]:
        """
        Get exon information for a transcript.

        Args:
            transcript_accession: RefSeq transcript accession

        Returns:
            List of Exon objects
        """
        exons = []

        # This would typically require mapping through the alignment API
        # For now, return empty list - will be populated by EnsemblClient
        self.log_debug(f"Exon fetching via NCBI not fully implemented for {transcript_accession}")

        return exons

    def clear_cache(self):
        """Clear the response cache."""
        self._cache.clear()
        self.log_info("NCBI client cache cleared")

"""
Coordinate translator for converting between CDS (c.) and genomic coordinates.
Handles strand orientation and exon boundaries.
"""

from typing import Optional, List, Tuple, Dict, Any
from dataclasses import dataclass

from .models import (
    Variant, Transcript, Exon, GenomicPosition, Strand
)
from ..api.ensembl_client import EnsemblClient
from ..api.ncbi_client import NCBIClient
from ..utils.logger import LoggerMixin
from ..utils.thread_safe import ThreadSafeCache


@dataclass
class MappingResult:
    """Result of coordinate mapping."""
    success: bool
    genomic_position: Optional[GenomicPosition] = None
    exon: Optional[Exon] = None
    exon_number: int = 0
    distance_to_exon_start: int = 0
    distance_to_exon_end: int = 0
    error: str = ""


class CoordinateTranslator(LoggerMixin):
    """
    Translates CDS coordinates to genomic coordinates.
    Handles exon mapping, strand orientation, and assembly conversion.
    """

    def __init__(
        self,
        ensembl_client: EnsemblClient = None,
        ncbi_client: NCBIClient = None,
        assembly: str = "GRCh38"
    ):
        """
        Initialize coordinate translator.

        Args:
            ensembl_client: EnsemblClient for coordinate mapping
            ncbi_client: NCBIClient for transcript data
            assembly: Reference assembly (GRCh38 or GRCh37)
        """
        self.ensembl_client = ensembl_client or EnsemblClient(assembly)
        self.ncbi_client = ncbi_client or NCBIClient()
        self.assembly = assembly

        # Thread-safe cache for transcript structures (supports deduplication
        # via get_or_compute — multiple variants sharing a transcript will
        # trigger only one set of API calls)
        self._transcript_cache = ThreadSafeCache()

    def set_assembly(self, assembly: str):
        """
        Change reference assembly.

        Args:
            assembly: GRCh38 or GRCh37
        """
        self.assembly = assembly
        self.ensembl_client.set_assembly(assembly)
        self.ncbi_client.set_assembly(assembly)
        self._transcript_cache.clear()
        self.log_info(f"Assembly changed to {assembly}")

    def map_cds_to_genomic(
        self,
        variant: Variant,
        transcript: Transcript = None
    ) -> MappingResult:
        """
        Map CDS position to genomic coordinates.

        Args:
            variant: Variant with CDS position information
            transcript: Optional pre-loaded transcript data

        Returns:
            MappingResult with genomic position
        """
        try:
            # Get transcript data if not provided
            if not transcript:
                transcript = self._get_transcript_data(variant.transcript_accession)

            if not transcript:
                return MappingResult(
                    success=False,
                    error=f"Could not fetch transcript data for {variant.transcript_accession}"
                )

            # Use Ensembl /map/cds/ endpoint first — correctly handles exon
            # structure and introns (unlike NCBI linear offset fallback)
            genomic_pos = self.ensembl_client.map_cds_to_genomic(
                variant.transcript_accession,
                variant.cds_start,
                variant.cds_end
            )

            if not genomic_pos:
                # Try manual mapping using transcript exon structure
                self.log_debug(
                    f"Ensembl CDS mapping failed for {variant.transcript_accession}, "
                    f"trying manual exon-based mapping..."
                )
                genomic_pos = self._manual_cds_to_genomic(
                    variant.cds_start,
                    variant.cds_end,
                    transcript
                )

            if not genomic_pos:
                # Last resort: NCBI linear offset (may be inaccurate for
                # multi-exon genes because it ignores introns)
                self.log_warning(
                    f"Using NCBI linear-offset fallback for "
                    f"{variant.transcript_accession} — result may be "
                    f"inaccurate for multi-exon genes"
                )
                genomic_pos = self.ncbi_client.map_cds_to_genomic(
                    variant.transcript_accession,
                    variant.cds_start,
                    variant.cds_end
                )

            if not genomic_pos:
                return MappingResult(
                    success=False,
                    error="Could not map CDS position to genomic coordinates"
                )

            # Find the exon containing this variant by GENOMIC position
            # (more reliable than CDS position which doesn't account for UTRs)
            exon, exon_number = self._find_exon_for_genomic_position(
                genomic_pos.start, transcript
            )

            # Calculate distances to exon boundaries
            dist_start, dist_end = 0, 0
            if exon:
                dist_start = variant.cds_start - exon.cds_start if exon.cds_start else 0
                dist_end = exon.cds_end - variant.cds_end if exon.cds_end else 0

            # Copy Ensembl exons to the variant's transcript so that
            # downstream code (intron-depth constraint, design region) can
            # use full exon structure
            if variant.transcript and transcript.exons:
                variant.transcript.exons = transcript.exons
                variant.transcript.chromosome = transcript.chromosome or variant.transcript.chromosome
                variant.transcript.strand = transcript.strand

            return MappingResult(
                success=True,
                genomic_position=genomic_pos,
                exon=exon,
                exon_number=exon_number,
                distance_to_exon_start=dist_start,
                distance_to_exon_end=dist_end
            )

        except Exception as e:
            self.log_exception(f"Error mapping coordinates: {e}")
            return MappingResult(
                success=False,
                error=str(e)
            )

    def _get_transcript_data(self, accession: str) -> Optional[Transcript]:
        """
        Get or fetch transcript data with exon structure.

        Uses thread-safe ``get_or_compute`` so that multiple variants
        sharing the same transcript trigger only one set of API calls.

        Args:
            accession: RefSeq transcript accession

        Returns:
            Transcript object or None
        """
        return self._transcript_cache.get_or_compute(
            accession,
            lambda: self._fetch_transcript_data(accession)
        )

    def _fetch_transcript_data(self, accession: str) -> Optional[Transcript]:
        """
        Actually fetch transcript data from NCBI / Ensembl.

        Called at most once per accession thanks to ``get_or_compute``.
        """
        # Try NCBI first (better for RefSeq accessions)
        ncbi_info = self.ncbi_client.get_transcript_info(accession)

        # Only try Ensembl if NCBI failed OR if we need exon data
        ensembl_info = None
        if not ncbi_info:
            self.log_debug(f"NCBI info not available for {accession}, trying Ensembl...")
            ensembl_info = self.ensembl_client.get_transcript_info(accession)

        # Prefer NCBI info, fall back to Ensembl
        info = ncbi_info if ncbi_info else ensembl_info

        if not info:
            self.log_warning(f"Could not fetch transcript info from NCBI or Ensembl for {accession}")
            return None

        # Get exons from Ensembl for manual CDS-to-genomic mapping fallback
        exons = []
        try:
            if not ensembl_info:
                ensembl_info = self.ensembl_client.get_transcript_info(accession)
            if ensembl_info and ensembl_info.get('ensembl_id'):
                exons = self.ensembl_client.get_transcript_exons(
                    ensembl_id=ensembl_info.get('ensembl_id')
                )
                if exons:
                    self.log_debug(f"Fetched {len(exons)} exons for {accession}")
        except Exception as e:
            self.log_debug(f"Could not fetch exons for {accession}: {e}")

        # Get sequence - try NCBI first, then Ensembl if needed
        sequence = ncbi_info.get('sequence', '') if ncbi_info else ''
        if not sequence and ensembl_info and ensembl_info.get('ensembl_id'):
            sequence = self.ensembl_client.get_transcript_sequence(
                ensembl_id=ensembl_info.get('ensembl_id'),
                accession=accession
            ) or ''

        # Build transcript object
        version = None
        base_accession = accession
        if '.' in accession:
            parts = accession.split('.')
            base_accession = parts[0]
            if len(parts) > 1:
                try:
                    version = int(parts[1])
                except ValueError:
                    version = None

        # Prefer Ensembl strand and chromosome over NCBI — the NCBI
        # GenBank XML parser doesn't reliably extract strand orientation
        # or chromosome, while Ensembl's /lookup/id endpoint always does.
        strand_str = info.get('strand', '+')
        chromosome = info.get('chromosome', '')
        if ensembl_info:
            strand_str = ensembl_info.get('strand', strand_str)
            chromosome = ensembl_info.get('chromosome', chromosome) or chromosome

        transcript = Transcript(
            accession=base_accession,
            version=version,
            gene_symbol=info.get('gene_symbol', ''),
            chromosome=chromosome,
            strand=Strand.POSITIVE if strand_str == '+' else Strand.NEGATIVE,
            exons=exons,
            sequence=sequence or ''
        )

        # Calculate CDS positions for exons
        self._calculate_exon_cds_positions(transcript)

        return transcript

    def _calculate_exon_cds_positions(self, transcript: Transcript):
        """
        Calculate CDS positions for each exon, accounting for UTR regions.

        The CDS (coding sequence) does not necessarily start at the first
        nucleotide of the first exon -- it begins at the ATG start codon,
        which may be located within the first exon (after 5'UTR) or even
        in a later exon.  Similarly, the CDS ends at the stop codon, which
        may be before the last exon ends (before 3'UTR).

        If the transcript already has cds_start/cds_end set (genomic CDS
        boundaries from NCBI/Ensembl), we use them to clip UTR regions.
        Otherwise we fall back to the old behaviour (treat all exonic
        sequence as coding).

        Args:
            transcript: Transcript with exons
        """
        if not transcript.exons:
            return

        # Sort exons by genomic position
        sorted_exons = sorted(transcript.exons, key=lambda e: e.genomic_start)

        # If on negative strand, reverse order for CDS assignment
        if transcript.strand == Strand.NEGATIVE:
            sorted_exons = list(reversed(sorted_exons))

        # Try to get genomic CDS boundaries from transcript info that was
        # fetched from NCBI/Ensembl.  These are genomic coordinates of the
        # actual coding region (ATG .. stop codon).
        genomic_cds_start = transcript.cds_start if transcript.cds_start > 0 else None
        genomic_cds_end = transcript.cds_end if transcript.cds_end > 0 else None

        # If cds_start/cds_end look like CDS positions (small numbers, < 100000)
        # rather than genomic coordinates (large numbers, > 100000), they were
        # probably already set as CDS-relative -- skip genomic clipping.
        has_genomic_cds_bounds = (
            genomic_cds_start is not None
            and genomic_cds_end is not None
            and genomic_cds_start > 100000
        )

        # Assign CDS positions
        current_cds_pos = 1
        for exon in sorted_exons:
            exon_genomic_start = exon.genomic_start
            exon_genomic_end = exon.genomic_end

            if has_genomic_cds_bounds:
                # Clip exon to CDS boundaries (remove UTR portions)
                if transcript.strand == Strand.POSITIVE:
                    coding_start = max(exon_genomic_start, genomic_cds_start)
                    coding_end = min(exon_genomic_end, genomic_cds_end)
                else:
                    # Negative strand: genomic_cds_start < genomic_cds_end
                    # but CDS direction is reversed
                    coding_start = max(exon_genomic_start, genomic_cds_start)
                    coding_end = min(exon_genomic_end, genomic_cds_end)

                if coding_start > coding_end:
                    # This exon is entirely UTR
                    exon.cds_start = None
                    exon.cds_end = None
                    continue

                exon_coding_length = coding_end - coding_start + 1
            else:
                # No genomic CDS bounds -- treat entire exon as coding
                exon_coding_length = abs(exon_genomic_end - exon_genomic_start) + 1

            exon.cds_start = current_cds_pos
            exon.cds_end = current_cds_pos + exon_coding_length - 1
            current_cds_pos += exon_coding_length

        transcript.cds_start = 1
        transcript.cds_end = current_cds_pos - 1

    def _find_exon_for_position(
        self,
        cds_position: int,
        transcript: Transcript
    ) -> Tuple[Optional[Exon], int]:
        """
        Find the exon containing a CDS position.

        Args:
            cds_position: CDS position (1-based)
            transcript: Transcript with exon data

        Returns:
            Tuple of (Exon, exon_number) or (None, 0)
        """
        if not transcript.exons:
            return None, 0

        for exon in transcript.exons:
            if exon.cds_start and exon.cds_end:
                if exon.cds_start <= cds_position <= exon.cds_end:
                    return exon, exon.number

        return None, 0

    def _find_exon_for_genomic_position(
        self,
        genomic_position: int,
        transcript: Transcript
    ) -> Tuple[Optional[Exon], int]:
        """
        Find the exon containing a genomic position.

        More reliable than _find_exon_for_position for genes with UTRs,
        because the CDS position calculation doesn't account for UTR
        regions in exons.

        Args:
            genomic_position: Genomic coordinate (1-based)
            transcript: Transcript with exon data

        Returns:
            Tuple of (Exon, exon_number) or (None, 0)
        """
        if not transcript.exons:
            return None, 0

        for exon in transcript.exons:
            if exon.genomic_start <= genomic_position <= exon.genomic_end:
                return exon, exon.number

        return None, 0

    def _manual_cds_to_genomic(
        self,
        cds_start: int,
        cds_end: int,
        transcript: Transcript
    ) -> Optional[GenomicPosition]:
        """
        Manually map CDS to genomic using transcript exon structure.

        Args:
            cds_start: CDS start position
            cds_end: CDS end position
            transcript: Transcript with exon data

        Returns:
            GenomicPosition or None
        """
        if not transcript.exons:
            return None

        genomic_start = None
        genomic_end = None

        # Find genomic positions for CDS start and end
        for exon in transcript.exons:
            if exon.cds_start is None or exon.cds_end is None:
                continue

            # Check if CDS start is in this exon
            if exon.cds_start <= cds_start <= exon.cds_end:
                offset = cds_start - exon.cds_start
                if transcript.strand == Strand.POSITIVE:
                    genomic_start = exon.genomic_start + offset
                else:
                    genomic_start = exon.genomic_end - offset

            # Check if CDS end is in this exon
            if exon.cds_start <= cds_end <= exon.cds_end:
                offset = cds_end - exon.cds_start
                if transcript.strand == Strand.POSITIVE:
                    genomic_end = exon.genomic_start + offset
                else:
                    genomic_end = exon.genomic_end - offset

        if genomic_start is not None and genomic_end is not None:
            # Ensure start < end
            if genomic_start > genomic_end:
                genomic_start, genomic_end = genomic_end, genomic_start

            return GenomicPosition(
                chromosome=transcript.chromosome,
                start=genomic_start,
                end=genomic_end,
                strand=transcript.strand,
                assembly=self.assembly
            )

        return None

    def get_flanking_region(
        self,
        genomic_pos: GenomicPosition,
        upstream_bp: int = 500,
        downstream_bp: int = 500
    ) -> GenomicPosition:
        """
        Get a region flanking a genomic position.

        Args:
            genomic_pos: Central genomic position
            upstream_bp: Base pairs upstream
            downstream_bp: Base pairs downstream

        Returns:
            Extended GenomicPosition
        """
        return GenomicPosition(
            chromosome=genomic_pos.chromosome,
            start=max(1, genomic_pos.start - upstream_bp),
            end=genomic_pos.end + downstream_bp,
            strand=genomic_pos.strand,
            assembly=genomic_pos.assembly
        )

    def get_amplicon_region(
        self,
        variant: Variant,
        min_amplicon_size: int = 200,
        max_amplicon_size: int = 500,
        min_distance_from_variant: int = 60
    ) -> Optional[GenomicPosition]:
        """
        Calculate the region where an amplicon could be placed.

        Args:
            variant: Variant with genomic position
            min_amplicon_size: Minimum amplicon size
            max_amplicon_size: Maximum amplicon size
            min_distance_from_variant: Minimum distance of primers from variant

        Returns:
            GenomicPosition for potential amplicon region
        """
        if not variant.genomic_position:
            return None

        gpos = variant.genomic_position

        # Calculate flanking region that allows for primer placement
        # Primers need to be at least min_distance_from_variant away
        # and amplicon should be within size limits

        # Maximum upstream flank: allow for reverse primer + variant + min distance
        max_upstream = (max_amplicon_size // 2) + min_distance_from_variant

        # Maximum downstream flank: same logic
        max_downstream = (max_amplicon_size // 2) + min_distance_from_variant

        return GenomicPosition(
            chromosome=gpos.chromosome,
            start=max(1, gpos.start - max_upstream),
            end=gpos.end + max_downstream,
            strand=gpos.strand,
            assembly=self.assembly
        )

    def get_exon_intron_boundaries(
        self,
        transcript: Transcript,
        exon_number: int
    ) -> Dict[str, int]:
        """
        Get the intron boundaries flanking an exon.

        Args:
            transcript: Transcript with exon data
            exon_number: Target exon number

        Returns:
            Dictionary with boundary positions
        """
        boundaries = {
            'exon_start': 0,
            'exon_end': 0,
            'upstream_intron_start': 0,
            'upstream_intron_end': 0,
            'downstream_intron_start': 0,
            'downstream_intron_end': 0
        }

        if not transcript.exons:
            return boundaries

        # Find the target exon
        target_exon = None
        prev_exon = None
        next_exon = None

        sorted_exons = sorted(transcript.exons, key=lambda e: e.genomic_start)

        for i, exon in enumerate(sorted_exons):
            if exon.number == exon_number:
                target_exon = exon
                if i > 0:
                    prev_exon = sorted_exons[i - 1]
                if i < len(sorted_exons) - 1:
                    next_exon = sorted_exons[i + 1]
                break

        if not target_exon:
            return boundaries

        boundaries['exon_start'] = target_exon.genomic_start
        boundaries['exon_end'] = target_exon.genomic_end

        # Upstream intron (5' of exon considering strand)
        if prev_exon:
            if transcript.strand == Strand.POSITIVE:
                boundaries['upstream_intron_start'] = prev_exon.genomic_end + 1
                boundaries['upstream_intron_end'] = target_exon.genomic_start - 1
            else:
                boundaries['upstream_intron_start'] = target_exon.genomic_end + 1
                boundaries['upstream_intron_end'] = prev_exon.genomic_start - 1

        # Downstream intron (3' of exon considering strand)
        if next_exon:
            if transcript.strand == Strand.POSITIVE:
                boundaries['downstream_intron_start'] = target_exon.genomic_end + 1
                boundaries['downstream_intron_end'] = next_exon.genomic_start - 1
            else:
                boundaries['downstream_intron_start'] = next_exon.genomic_end + 1
                boundaries['downstream_intron_end'] = target_exon.genomic_start - 1

        return boundaries

    def clear_cache(self):
        """Clear transcript cache."""
        self._transcript_cache.clear()
        self.log_info("Coordinate translator cache cleared")

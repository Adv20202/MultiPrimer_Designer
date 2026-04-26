"""
Primer designer using Primer3 library.
Designs primers for PCR and qPCR applications.
"""

from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass

from ..core.models import (
    Variant, Primer, PrimerPair, Amplicon, Probe, GenomicPosition,
    DesignParameters, DesignMode, PopulationVariant, DesignResult
)
from ..core.grouper import VariantGroup
from ..api.ensembl_client import EnsemblClient
from ..api.variant_db_client import VariantDBClient
from ..utils.logger import LoggerMixin
from ..utils.config import get_config

try:
    import primer3
    PRIMER3_AVAILABLE = True
    # primer3-py API changed between versions:
    #   older: primer3.bindings.design_primers()
    #   newer: primer3.design_primers()
    # Create a single reference that works with both
    if hasattr(primer3, 'bindings') and hasattr(primer3.bindings, 'design_primers'):
        _primer3_design = primer3.bindings.design_primers
    elif hasattr(primer3, 'design_primers'):
        _primer3_design = primer3.design_primers
    else:
        _primer3_design = primer3.designPrimers  # camelCase fallback
except ImportError:
    PRIMER3_AVAILABLE = False
    _primer3_design = None


@dataclass
class DesignRegion:
    """Region for primer design with sequence and constraints."""
    sequence: str
    target_start: int  # Position of target in sequence (0-based)
    target_length: int
    genomic_start: int = 0  # Genomic start coordinate for position conversion
    excluded_regions: List[Tuple[int, int]] = None  # [(start, length), ...] - deprecated, use masked_positions
    included_region: Tuple[int, int] = None  # (start, length)
    masked_positions: set = None  # Set of 0-based positions masked with 'N'
    exon_boundaries: List[Tuple[int, int, int]] = None  # [(seq_start, seq_end, exon_number), ...] for exons overlapping region

    def __post_init__(self):
        if self.excluded_regions is None:
            self.excluded_regions = []
        if self.masked_positions is None:
            self.masked_positions = set()
        if self.exon_boundaries is None:
            self.exon_boundaries = []


class PrimerDesigner(LoggerMixin):
    """
    Primer designer using Primer3 for PCR/qPCR primer design.
    Handles population variant filtering and constraint management.
    """

    def __init__(
        self,
        ensembl_client: EnsemblClient = None,
        variant_db_client: VariantDBClient = None,
        parameters: DesignParameters = None
    ):
        """
        Initialize primer designer.

        Args:
            ensembl_client: For fetching sequences
            variant_db_client: For population variant filtering
            parameters: Design parameters
        """
        if not PRIMER3_AVAILABLE:
            raise ImportError(
                "primer3-py library is required for primer design. "
                "Install it with: pip install primer3-py"
            )

        self.ensembl_client = ensembl_client or EnsemblClient()
        self.variant_db_client = variant_db_client
        self.parameters = parameters or DesignParameters()

        # Cancel flag for long operations
        self._cancel_requested = False

    def request_cancel(self):
        """Request cancellation of current operation."""
        self._cancel_requested = True

    def reset_cancel(self):
        """Reset cancellation flag."""
        self._cancel_requested = False

    def design_primers(
        self,
        variant_group: VariantGroup,
        mode: DesignMode = DesignMode.PCR,
        progress_callback=None,
        homology_result=None,
    ) -> DesignResult:
        """
        Design primers for a variant group.

        Args:
            variant_group: Group of variants to design primers for
            mode: PCR or qPCR mode
            progress_callback: Optional callback for progress updates
            homology_result: Optional HomologyResult for pseudogene discrimination

        Returns:
            DesignResult with designed amplicons
        """
        self.reset_cancel()

        if not variant_group.variants:
            return DesignResult(
                variants=[],
                amplicons=[],
                success=False,
                message="No variants in group"
            )

        try:
            self.log_info(f"[PrimerDesigner] Starting primer design for group: {variant_group.variants[0].gene_symbol if variant_group.variants else 'unknown'}")

            # Get the target region
            self.log_debug("[PrimerDesigner] Fetching design region sequence...")
            region = self._get_design_region(variant_group, mode)
            if not region:
                self.log_warning("[PrimerDesigner] Could not determine design region")
                return DesignResult(
                    variants=variant_group.variants,
                    amplicons=[],
                    success=False,
                    message="Could not determine design region"
                )

            self.log_info(f"[PrimerDesigner] Fetched {len(region.sequence)} bp sequence")

            if self._cancel_requested:
                return DesignResult(
                    variants=variant_group.variants,
                    amplicons=[],
                    success=False,
                    message="Operation cancelled"
                )

            # Get population variants if filtering enabled
            if self.parameters.filter_population_variants and self.variant_db_client:
                # Always fetch population variants for the full design region.
                # Pre-fetched _pop_variants from GUI use a smaller flank (250bp)
                # than the designer needs (max_amplicon + min_distance ≈ 560bp),
                # so we cannot safely reuse them — variants in positions 251-560bp
                # from the target would be missed, and primers could land on known SNPs.
                self.log_debug("[PrimerDesigner] Fetching population variants for design region...")
                pop_variants = self._get_population_variants(variant_group)
                self.log_info(f"[PrimerDesigner] Found {len(pop_variants)} population variants")

                region = self._apply_population_filter(region, pop_variants, region.genomic_start)

            if self._cancel_requested:
                return DesignResult(
                    variants=variant_group.variants,
                    amplicons=[],
                    success=False,
                    message="Operation cancelled"
                )

            # Extract discriminating positions from homology analysis
            discriminating_positions = {}
            homology_active = False
            num_homologous_regions = 0

            if (homology_result is not None
                    and self.parameters.use_homology_discrimination):
                from ..primer.homology_analyzer import HomologyAnalyzer
                discriminating_positions = HomologyAnalyzer.extract_discriminating_positions(
                    homology_result
                )
                secondary_hits = [
                    h for h in homology_result.hits
                    if not h.is_primary and not h.is_supplementary
                ]
                num_homologous_regions = len(secondary_hits)
                if discriminating_positions:
                    homology_active = True
                    pos_min = min(discriminating_positions.keys())
                    pos_max = max(discriminating_positions.keys())
                    self.log_info(
                        f"[PrimerDesigner] Homology discrimination active: "
                        f"{len(discriminating_positions)} discriminating positions "
                        f"(genomic range {pos_min}-{pos_max}) "
                        f"from {num_homologous_regions} homologous region(s)"
                    )
                else:
                    self.log_info(
                        "[PrimerDesigner] Homology analysis provided but "
                        "no discriminating positions found"
                    )
            elif homology_result is None:
                self.log_info(
                    "[PrimerDesigner] No homology analysis available — "
                    "primers designed without pseudogene discrimination"
                )

            # Run Primer3 — tiered approach when homology data is available
            self.log_debug("[PrimerDesigner] Running Primer3...")

            if homology_active:
                # Tiered Primer3 with homology constraints
                amplicons, tier, tier_message = self._run_primer3_tiered(
                    region, mode,
                    discriminating_positions,
                    num_homologous_regions,
                    variant_group=variant_group,
                )

                if not amplicons:
                    suggestions = self._generate_suggestions(region, mode)
                    return DesignResult(
                        variants=variant_group.variants,
                        amplicons=[],
                        success=False,
                        message="No primers found with current parameters",
                        suggestions=suggestions,
                        homology_discriminated=True,
                        num_homologous_regions=num_homologous_regions,
                        homology_tier=tier,
                        homology_tier_message=tier_message,
                        homology_result=homology_result,
                    )

                amplicons = amplicons[:self.parameters.num_primer_pairs]
                return DesignResult(
                    variants=variant_group.variants,
                    amplicons=amplicons,
                    success=True,
                    message=f"Found {len(amplicons)} primer pair(s) (homology Tier {tier})",
                    homology_discriminated=True,
                    num_homologous_regions=num_homologous_regions,
                    homology_tier=tier,
                    homology_tier_message=tier_message,
                    homology_result=homology_result,
                )

            else:
                # Standard Primer3 without homology constraints
                primer_results = self._run_primer3(region, mode)
                self.log_info(
                    f"[PrimerDesigner] Primer3 returned "
                    f"{len(primer_results) if primer_results else 0} results")

                if not primer_results:
                    suggestions = self._generate_suggestions(region, mode)
                    return DesignResult(
                        variants=variant_group.variants,
                        amplicons=[],
                        success=False,
                        message="No primers found with current parameters",
                        suggestions=suggestions,
                    )

                amplicons = self._create_amplicons(
                    primer_results, variant_group, region, mode)

                return DesignResult(
                    variants=variant_group.variants,
                    amplicons=amplicons,
                    success=len(amplicons) > 0,
                    message=f"Found {len(amplicons)} primer pair(s)",
                )

        except Exception as e:
            self.log_exception(f"Error designing primers: {e}")
            return DesignResult(
                variants=variant_group.variants,
                amplicons=[],
                success=False,
                message=f"Design error: {str(e)}"
            )

    def _get_design_region(
        self,
        variant_group: VariantGroup,
        mode: DesignMode
    ) -> Optional[DesignRegion]:
        """
        Get the sequence region for primer design.

        Args:
            variant_group: Variant group
            mode: Design mode

        Returns:
            DesignRegion or None
        """
        if not variant_group.chromosome or variant_group.genomic_start == 0:
            # Try to get position from first variant
            for variant in variant_group.variants:
                if variant.genomic_position:
                    variant_group.chromosome = variant.genomic_position.chromosome
                    variant_group.genomic_start = variant.genomic_position.start
                    variant_group.genomic_end = variant.genomic_position.end
                    break

        if not variant_group.chromosome:
            return None

        # Determine flank size based on mode
        if mode == DesignMode.QPCR:
            max_size = self.parameters.qpcr_max_amplicon_size
        else:
            max_size = self.parameters.max_amplicon_size

        # Add flanking sequence for primer placement
        flank = max_size + self.parameters.min_distance_from_variant

        region_start = max(1, variant_group.genomic_start - flank)
        region_end = variant_group.genomic_end + flank

        region_size = region_end - region_start + 1

        # Check if any variant has pre-masked sequence from GUI
        # NOTE: We intentionally do NOT use pre-masked sequence here because:
        # 1. GUI uses a smaller flank size (250bp) than PrimerDesigner needs
        # 2. Target position calculations would be wrong with mismatched sequence lengths
        # 3. We will apply our own masking based on the design parameters
        # The pre-masked sequence is used for FASTA export, not for primer design
        pre_masked_positions = set()
        for variant in variant_group.variants:
            if hasattr(variant, '_masked_positions') and variant._masked_positions:
                # We keep track of masked positions but fetch fresh sequence
                pre_masked_positions = variant._masked_positions
                self.log_info(f"[PrimerDesigner] Variant has {len(pre_masked_positions)} pre-masked positions (will apply during design)")
                break

        # Always fetch fresh sequence from Ensembl with proper flanking for primer design
        sequence = self.ensembl_client.get_genomic_sequence(
            chromosome=variant_group.chromosome,
            start=region_start,
            end=region_end
        )

        if not sequence:
            self.log_warning(f"Could not fetch sequence for {variant_group.chromosome}:{region_start}-{region_end}")
            return None

        self.log_info(f"[PrimerDesigner] Fetched {len(sequence)} bp sequence ({variant_group.chromosome}:{region_start}-{region_end})")

        # Calculate target position within sequence
        target_start = variant_group.genomic_start - region_start
        target_length = variant_group.genomic_end - variant_group.genomic_start + 1

        # Collect exon boundaries (genomic -> sequence-relative) for intron depth constraint
        exon_boundaries = []
        for variant in variant_group.variants:
            if variant.transcript and variant.transcript.exons:
                for exon in variant.transcript.exons:
                    # Convert genomic exon coords to sequence-relative coords
                    exon_seq_start = exon.genomic_start - region_start
                    exon_seq_end = exon.genomic_end - region_start
                    # Only include exons that overlap with our design region
                    if exon_seq_end >= 0 and exon_seq_start < len(sequence):
                        exon_seq_start = max(0, exon_seq_start)
                        exon_seq_end = min(len(sequence) - 1, exon_seq_end)
                        exon_boundaries.append((exon_seq_start, exon_seq_end, exon.number))
                break  # Use first variant's transcript exons
            elif variant.exon:
                # Fallback: only the variant's own exon is available
                exon_seq_start = max(0, variant.exon.genomic_start - region_start)
                exon_seq_end = min(len(sequence) - 1, variant.exon.genomic_end - region_start)
                exon_boundaries.append((exon_seq_start, exon_seq_end, variant.exon.number))
                break

        if exon_boundaries:
            self.log_debug(f"[PrimerDesigner] Found {len(exon_boundaries)} exon(s) in design region")

        return DesignRegion(
            sequence=sequence,
            target_start=target_start,
            target_length=target_length,
            genomic_start=region_start,
            masked_positions=pre_masked_positions,
            exon_boundaries=exon_boundaries
        )

    def _get_population_variants(
        self,
        variant_group: VariantGroup
    ) -> List[PopulationVariant]:
        """Get population variants in the design region."""
        if not self.variant_db_client:
            return []

        # Expand region for primer placement area
        flank = self.parameters.max_amplicon_size
        start = max(1, variant_group.genomic_start - flank)
        end = variant_group.genomic_end + flank

        return self.variant_db_client.get_variants_in_region(
            chromosome=variant_group.chromosome,
            start=start,
            end=end,
            populations=self.parameters.selected_populations,
            maf_threshold=self.parameters.maf_threshold
        )

    def _apply_population_filter(
        self,
        region: DesignRegion,
        pop_variants: List[PopulationVariant],
        region_start: int = 0
    ) -> DesignRegion:
        """
        Apply population variant filter by masking positions with 'N' in sequence.

        Args:
            region: Design region
            pop_variants: Population variants to exclude
            region_start: Genomic start position of the design region

        Returns:
            Updated DesignRegion with masked sequence
        """
        # Convert sequence to list for easier modification
        sequence = list(region.sequence)
        masked_positions = set()
        masked_count = 0

        for var in pop_variants:
            # Check if variant is above MAF threshold
            maf = var.maf_global
            for pop in self.parameters.selected_populations:
                if pop != 'global':
                    maf = max(maf, var.maf_by_population.get(pop, 0.0))

            if maf >= self.parameters.maf_threshold:
                # Convert genomic position to sequence-relative position (0-based)
                var_pos_in_seq = var.position - region_start

                # Only mask if position is within the sequence bounds
                if 0 <= var_pos_in_seq < len(sequence):
                    # Calculate variant length (for indels)
                    var_length = max(1, len(var.ref) if var.ref else 1)

                    # Mask the variant position(s) with 'N'
                    for i in range(var_pos_in_seq, min(var_pos_in_seq + var_length, len(sequence))):
                        sequence[i] = 'N'
                        masked_positions.add(i)
                    masked_count += 1

        # Update sequence with masked version and save masked positions
        region.sequence = ''.join(sequence)
        region.masked_positions = masked_positions

        # Calculate percentage of sequence masked
        pct_masked = (len(masked_positions) / len(sequence)) * 100 if sequence else 0
        self.log_info(f"[PrimerDesigner] Masked {masked_count} variant positions ({len(masked_positions)} nucleotides, {pct_masked:.1f}% of sequence) with 'N'")

        if pct_masked > 20:
            self.log_warning(f"[PrimerDesigner] High masking rate ({pct_masked:.1f}%) may prevent primer design. Consider increasing MAF threshold.")

        return region

    def _run_primer3(
        self,
        region: DesignRegion,
        mode: DesignMode,
        extra_candidates: bool = False,
        junction_list: List[int] = None,
        ok_region_list: List[List[int]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Run Primer3 for primer design.

        Args:
            region: Design region with sequence and constraints
            mode: PCR or qPCR mode
            extra_candidates: If True, request more candidates for re-ranking
            junction_list: SEQUENCE_OVERLAP_JUNCTION_LIST positions (0-based)
            ok_region_list: SEQUENCE_PRIMER_PAIR_OK_REGION_LIST entries

        Returns:
            List of primer results
        """
        # Prepare Primer3 input
        seq_args = {
            'SEQUENCE_ID': 'target',
            'SEQUENCE_TEMPLATE': region.sequence,
            'SEQUENCE_TARGET': [(region.target_start, region.target_length)]
        }

        # Build list of excluded regions
        excluded_regions = []

        # 1. Enforce minimum distance from variant by excluding a zone around the target
        if self.parameters.use_variant_distance_constraint and self.parameters.min_distance_from_variant > 0:
            min_dist = self.parameters.min_distance_from_variant
            excl_start = max(0, region.target_start - min_dist)
            excl_end = min(len(region.sequence), region.target_start + region.target_length + min_dist)
            excl_length = excl_end - excl_start
            excluded_regions.append((excl_start, excl_length))
            self.log_debug(f"[Primer3] Excluded region around target: ({excl_start}, {excl_length}) = {min_dist}bp buffer")

        # 2. Enforce min distance from exon boundary — exclude the intronic
        #    zone within N bp of each exon/intron junction (splice-site region).
        #    Primers may sit in exons freely or deep in introns, but NOT near
        #    the exon boundary on the intronic side.
        if (self.parameters.use_splice_site_constraint and
                self.parameters.min_distance_from_exon_junction > 0 and
                region.exon_boundaries):
            min_dist = self.parameters.min_distance_from_exon_junction
            seq_len = len(region.sequence)

            # Sort exon boundaries for correct intron identification
            sorted_exons = sorted(region.exon_boundaries, key=lambda x: x[0])

            # Build list of excluded intronic intervals: for each exon
            # boundary, the N bp on the INTRONIC side are forbidden.
            # Exonic positions are never excluded by this constraint.
            splice_exclusions = []
            for idx, (exon_start, exon_end, _exon_num) in enumerate(sorted_exons):
                # Left (upstream) side of exon — intron is to the LEFT
                # Exclude intronic positions [exon_start - min_dist, exon_start - 1]
                left_excl_start = max(0, exon_start - min_dist)
                left_excl_end = max(0, exon_start)  # up to (not including) exon start

                # But do NOT exclude positions that belong to a preceding exon
                if idx > 0:
                    prev_exon_end = sorted_exons[idx - 1][1]
                    left_excl_start = max(left_excl_start, prev_exon_end + 1)

                if left_excl_end > left_excl_start:
                    splice_exclusions.append((left_excl_start, left_excl_end))

                # Right (downstream) side of exon — intron is to the RIGHT
                # Exclude intronic positions [exon_end + 1, exon_end + min_dist]
                right_excl_start = min(seq_len, exon_end + 1)
                right_excl_end = min(seq_len, exon_end + 1 + min_dist)

                # But do NOT exclude positions that belong to the next exon
                if idx < len(sorted_exons) - 1:
                    next_exon_start = sorted_exons[idx + 1][0]
                    right_excl_end = min(right_excl_end, next_exon_start)

                if right_excl_end > right_excl_start:
                    splice_exclusions.append((right_excl_start, right_excl_end))

            # Check viability: if exclusions would block too much sequence,
            # auto-reduce min_dist to leave room for amplicon design
            if splice_exclusions:
                total_excluded_bp = sum(e - s for s, e in splice_exclusions)
                total_available = seq_len - total_excluded_bp
                min_required = self.parameters.min_amplicon_size + 2 * self.parameters.min_distance_from_variant

                if total_available < min_required:
                    # Reduce min_dist until enough room is available
                    reduced_dist = max(10, min_dist // 2)
                    self.log_info(
                        f"[Primer3] Splice-site exclusion zone auto-reduced: "
                        f"{min_dist}bp → {reduced_dist}bp "
                        f"(available {total_available}bp < required {min_required}bp)"
                    )
                    # Rebuild with reduced distance — recursive-safe single retry
                    splice_exclusions = []
                    for idx, (exon_start, exon_end, _exon_num) in enumerate(sorted_exons):
                        left_excl_start = max(0, exon_start - reduced_dist)
                        left_excl_end = max(0, exon_start)
                        if idx > 0:
                            left_excl_start = max(left_excl_start, sorted_exons[idx - 1][1] + 1)
                        if left_excl_end > left_excl_start:
                            splice_exclusions.append((left_excl_start, left_excl_end))

                        right_excl_start = min(seq_len, exon_end + 1)
                        right_excl_end = min(seq_len, exon_end + 1 + reduced_dist)
                        if idx < len(sorted_exons) - 1:
                            right_excl_end = min(right_excl_end, sorted_exons[idx + 1][0])
                        if right_excl_end > right_excl_start:
                            splice_exclusions.append((right_excl_start, right_excl_end))
                    min_dist = reduced_dist

            # Convert intervals to Primer3 SEQUENCE_EXCLUDED_REGION format (start, length)
            for excl_start, excl_end in splice_exclusions:
                excl_len = excl_end - excl_start
                if excl_len > 0:
                    excluded_regions.append((excl_start, excl_len))

            self.log_debug(
                f"[Primer3] Splice-site exclusion: {min_dist}bp from each exon boundary, "
                f"{len(splice_exclusions)} zone(s), {len(excluded_regions)} total excluded regions"
            )
        elif (self.parameters.use_splice_site_constraint and
              self.parameters.min_distance_from_exon_junction > 0 and
              not region.exon_boundaries):
            self.log_warning(
                "[Primer3] Splice-site exclusion constraint enabled but no exon "
                "boundaries available — constraint skipped. Exon data may not "
                "have been fetched from Ensembl for this transcript."
            )

        if excluded_regions:
            seq_args['SEQUENCE_EXCLUDED_REGION'] = excluded_regions

        if region.included_region:
            seq_args['SEQUENCE_INCLUDED_REGION'] = region.included_region

        # Homology-aware constraints
        if junction_list:
            seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = junction_list
            self.log_info(
                f"[Primer3] OVERLAP_JUNCTION_LIST: {len(junction_list)} positions "
                f"(e.g. {junction_list[:5]})")

        if ok_region_list:
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = ok_region_list
            self.log_info(
                f"[Primer3] PRIMER_PAIR_OK_REGION_LIST: {len(ok_region_list)} region pairs")

        # Global parameters
        global_args = self._get_primer3_global_args(mode)

        # Junction overlap parameters
        if junction_list:
            global_args['PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'] = 3
            global_args['PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION'] = 4

        # Request more candidates when homology re-ranking will be applied
        if extra_candidates:
            original_num = global_args.get('PRIMER_NUM_RETURN', self.parameters.num_primer_pairs)
            expanded_num = max(20, original_num * 4)
            global_args['PRIMER_NUM_RETURN'] = expanded_num
            self.log_info(
                f"[Primer3] Expanded candidate pool: {original_num} → {expanded_num} "
                f"(for homology re-ranking)"
            )

        # Log Primer3 input for debugging
        seq_len = len(region.sequence)
        n_count = region.sequence.count('N') + region.sequence.count('n')
        self.log_info(f"[Primer3] Input: seq_len={seq_len}, N_count={n_count} ({100*n_count/seq_len:.1f}%), target=({region.target_start}, {region.target_length})")
        self.log_debug(f"[Primer3] Product size range: {global_args.get('PRIMER_PRODUCT_SIZE_RANGE')}")

        # Run Primer3
        try:
            results = _primer3_design(seq_args, global_args)

            # Log Primer3 output for debugging
            num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            explain_left = results.get('PRIMER_LEFT_EXPLAIN', '')
            explain_right = results.get('PRIMER_RIGHT_EXPLAIN', '')
            explain_pair = results.get('PRIMER_PAIR_EXPLAIN', '')

            self.log_info(f"[Primer3] Output: {num_returned} primer pairs returned")
            if explain_left:
                self.log_debug(f"[Primer3] LEFT explain: {explain_left}")
            if explain_right:
                self.log_debug(f"[Primer3] RIGHT explain: {explain_right}")
            if explain_pair:
                self.log_debug(f"[Primer3] PAIR explain: {explain_pair}")

            if num_returned == 0:
                # Log more details when no primers found
                self.log_warning(f"[Primer3] No primers found. LEFT: {explain_left}")
                self.log_warning(f"[Primer3] No primers found. RIGHT: {explain_right}")
                self.log_warning(f"[Primer3] No primers found. PAIR: {explain_pair}")

            return self._parse_primer3_results(results)
        except Exception as e:
            self.log_error(f"Primer3 error: {e}")
            return []

    def _get_primer3_global_args(self, mode: DesignMode) -> Dict[str, Any]:
        """Get Primer3 global arguments based on design mode."""
        p = self.parameters

        if mode == DesignMode.QPCR:
            product_min = p.qpcr_min_amplicon_size
            product_max = p.qpcr_max_amplicon_size
            product_opt = (product_min + product_max) // 2
        else:
            product_min = p.min_amplicon_size
            product_max = p.max_amplicon_size
            product_opt = (product_min + product_max) // 2

        args = {
            # Primer size
            'PRIMER_MIN_SIZE': p.primer_min_size,
            'PRIMER_OPT_SIZE': p.primer_opt_size,
            'PRIMER_MAX_SIZE': p.primer_max_size,

            # Primer Tm
            'PRIMER_MIN_TM': p.primer_min_tm,
            'PRIMER_OPT_TM': p.primer_opt_tm,
            'PRIMER_MAX_TM': p.primer_max_tm,

            # Primer GC
            'PRIMER_MIN_GC': p.primer_min_gc,
            'PRIMER_OPT_GC_PERCENT': p.primer_opt_gc,
            'PRIMER_MAX_GC': p.primer_max_gc,

            # Product size
            'PRIMER_PRODUCT_SIZE_RANGE': [[product_min, product_max]],
            'PRIMER_PRODUCT_OPT_SIZE': product_opt,

            # Chemistry for Tm calculation
            'PRIMER_SALT_MONOVALENT': p.mv_conc,
            'PRIMER_SALT_DIVALENT': p.dv_conc,
            'PRIMER_DNTP_CONC': p.dntp_conc,
            'PRIMER_DNA_CONC': p.dna_conc,

            # Complementarity
            'PRIMER_MAX_SELF_ANY': p.max_self_complementarity,
            'PRIMER_MAX_SELF_END': p.max_end_complementarity,
            'PRIMER_PAIR_MAX_COMPL_ANY': p.max_pair_complementarity,

            # Tm difference
            'PRIMER_PAIR_MAX_DIFF_TM': p.max_tm_difference,

            # Number of primers to return
            'PRIMER_NUM_RETURN': p.num_primer_pairs,

            # Avoid poly-X runs
            'PRIMER_MAX_POLY_X': 4,

            # GC clamp
            'PRIMER_GC_CLAMP': 1,
        }

        # Add internal oligo (probe) settings for qPCR
        if mode == DesignMode.QPCR:
            args.update({
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_INTERNAL_MIN_SIZE': p.probe_min_size,
                'PRIMER_INTERNAL_MAX_SIZE': p.probe_max_size,
                'PRIMER_INTERNAL_MIN_TM': p.probe_min_tm,
                'PRIMER_INTERNAL_MAX_TM': p.probe_max_tm,
                'PRIMER_INTERNAL_OPT_TM': (p.probe_min_tm + p.probe_max_tm) / 2,
            })

        return args

    def _parse_primer3_results(self, results: Dict) -> List[Dict[str, Any]]:
        """Parse Primer3 output into structured format."""
        parsed = []

        num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

        for i in range(num_returned):
            primer_data = {
                'left_sequence': results.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
                'left_position': results.get(f'PRIMER_LEFT_{i}', (0, 0)),
                'left_tm': results.get(f'PRIMER_LEFT_{i}_TM', 0.0),
                'left_gc': results.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0.0),
                'left_self_any': results.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH', 0.0),
                'left_self_end': results.get(f'PRIMER_LEFT_{i}_SELF_END_TH', 0.0),

                'right_sequence': results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                'right_position': results.get(f'PRIMER_RIGHT_{i}', (0, 0)),
                'right_tm': results.get(f'PRIMER_RIGHT_{i}_TM', 0.0),
                'right_gc': results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0.0),
                'right_self_any': results.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH', 0.0),
                'right_self_end': results.get(f'PRIMER_RIGHT_{i}_SELF_END_TH', 0.0),

                'product_size': results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
                'pair_compl_any': results.get(f'PRIMER_PAIR_{i}_COMPL_ANY_TH', 0.0),
                'pair_compl_end': results.get(f'PRIMER_PAIR_{i}_COMPL_END_TH', 0.0),
            }

            # Internal oligo (probe) for qPCR
            if f'PRIMER_INTERNAL_{i}_SEQUENCE' in results:
                primer_data['probe_sequence'] = results.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')
                primer_data['probe_position'] = results.get(f'PRIMER_INTERNAL_{i}', (0, 0))
                primer_data['probe_tm'] = results.get(f'PRIMER_INTERNAL_{i}_TM', 0.0)
                primer_data['probe_gc'] = results.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT', 0.0)

            if primer_data['left_sequence'] and primer_data['right_sequence']:
                parsed.append(primer_data)

        return parsed

    def _create_amplicons(
        self,
        primer_results: List[Dict],
        variant_group: VariantGroup,
        region: DesignRegion,
        mode: DesignMode
    ) -> List[Amplicon]:
        """Create Amplicon objects from Primer3 results."""
        amplicons = []

        for result in primer_results:
            # Create forward primer
            left_pos = result['left_position']
            forward = Primer(
                sequence=result['left_sequence'],
                start=left_pos[0],
                end=left_pos[0] + left_pos[1] - 1,
                tm=result['left_tm'],
                gc_content=result['left_gc'],
                is_forward=True,
                self_complementarity=result['left_self_any'],
                end_complementarity=result['left_self_end']
            )

            # Create reverse primer
            right_pos = result['right_position']
            reverse = Primer(
                sequence=result['right_sequence'],
                start=right_pos[0] - right_pos[1] + 1,
                end=right_pos[0],
                tm=result['right_tm'],
                gc_content=result['right_gc'],
                is_forward=False,
                self_complementarity=result['right_self_any'],
                end_complementarity=result['right_self_end']
            )

            # Create primer pair
            primer_pair = PrimerPair(
                forward=forward,
                reverse=reverse,
                product_size=result['product_size'],
                tm_difference=abs(forward.tm - reverse.tm),
                pair_complementarity=result['pair_compl_any']
            )

            # Create probe if present (qPCR mode)
            probe = None
            if 'probe_sequence' in result:
                probe_pos = result['probe_position']
                probe = Probe(
                    sequence=result['probe_sequence'],
                    start=probe_pos[0],
                    end=probe_pos[0] + probe_pos[1] - 1,
                    tm=result['probe_tm'],
                    gc_content=result['probe_gc']
                )

            # Map exon boundaries into amplicon-relative coordinates
            exon_regions_in_amp = []
            amp_start = forward.start  # amplicon starts at forward primer start
            amp_end = reverse.end      # amplicon ends at reverse primer end
            if region.exon_boundaries:
                for exon_start, exon_end, exon_num in sorted(region.exon_boundaries, key=lambda x: x[0]):
                    # Check if exon overlaps with the amplicon
                    overlap_start = max(exon_start, amp_start)
                    overlap_end = min(exon_end, amp_end)
                    if overlap_start <= overlap_end:
                        # Convert to amplicon-relative (0-based within product)
                        rel_start = overlap_start - amp_start
                        rel_end = overlap_end - amp_start
                        exon_regions_in_amp.append((rel_start, rel_end, exon_num))

            # Create amplicon
            amplicon = Amplicon(
                primer_pair=primer_pair,
                variants=variant_group.variants,
                probe=probe,
                masked_positions=region.masked_positions,
                target_start=region.target_start,
                target_length=region.target_length,
                exon_regions_in_amplicon=exon_regions_in_amp,
            )

            amplicons.append(amplicon)

        return amplicons

    # ------------------------------------------------------------------
    #  Homology-aware tiered Primer3 design
    # ------------------------------------------------------------------

    def _genomic_to_local(
        self,
        discriminating_positions: Dict[int, int],
        region: DesignRegion,
    ) -> Dict[int, int]:
        """Convert discriminating positions from genomic to design-region coords."""
        offset = region.genomic_start
        seq_len = len(region.sequence)
        local = {}
        for gpos, count in discriminating_positions.items():
            lpos = gpos - offset
            if 0 <= lpos < seq_len:
                local[lpos] = count
        return local

    def _select_junction_positions(
        self,
        disc_local: Dict[int, int],
        max_count: int = 20,
    ) -> List[int]:
        """
        Select the strongest discriminating positions for
        SEQUENCE_OVERLAP_JUNCTION_LIST.

        Returns 0-based positions in design region, sorted by strength.
        """
        if not disc_local:
            return []
        # Sort by count descending, take top N
        sorted_pos = sorted(disc_local.items(), key=lambda x: -x[1])
        return [pos for pos, _ in sorted_pos[:max_count]]

    def _compute_discriminating_windows(
        self,
        disc_local: Dict[int, int],
        seq_len: int,
        min_window: int = 18,
        gap_tolerance: int = 10,
    ) -> List[Tuple[int, int]]:
        """
        Compute contiguous discriminating windows from local positions.

        Groups nearby discriminating positions into windows, expanding
        each by a margin suitable for primer placement.

        Args:
            disc_local: Dict of local_pos → count
            seq_len: Length of the design region sequence
            min_window: Minimum window size (= min primer length)
            gap_tolerance: Max gap between positions to merge into one window

        Returns:
            List of (start, length) tuples in design region coords.
        """
        if not disc_local:
            return []

        positions = sorted(disc_local.keys())

        # Cluster positions separated by <= gap_tolerance
        clusters: List[List[int]] = []
        current_cluster = [positions[0]]
        for pos in positions[1:]:
            if pos - current_cluster[-1] <= gap_tolerance:
                current_cluster.append(pos)
            else:
                clusters.append(current_cluster)
                current_cluster = [pos]
        clusters.append(current_cluster)

        # Convert clusters to windows with margin
        margin = 15  # extend each window by margin bp on each side
        windows = []
        for cluster in clusters:
            if len(cluster) < 2:
                continue  # need at least 2 positions for a useful window
            start = max(0, cluster[0] - margin)
            end = min(seq_len, cluster[-1] + margin + 1)
            length = end - start
            if length >= min_window:
                windows.append((start, length))

        # Merge overlapping windows
        if len(windows) > 1:
            merged = [windows[0]]
            for start, length in windows[1:]:
                prev_start, prev_len = merged[-1]
                prev_end = prev_start + prev_len
                if start <= prev_end:
                    new_end = max(prev_end, start + length)
                    merged[-1] = (prev_start, new_end - prev_start)
                else:
                    merged.append((start, length))
            windows = merged

        return windows

    def _build_ok_region_pairs(
        self,
        windows: List[Tuple[int, int]],
        target_start: int,
        target_length: int,
        seq_len: int,
    ) -> List[List[int]]:
        """
        Build SEQUENCE_PRIMER_PAIR_OK_REGION_LIST entries.

        For each window, creates a pair where:
        - One side (left or right of target) is unconstrained (-1, -1)
        - The other side is constrained to the discriminating window
        This ensures at least one primer lands in a discriminating region.

        Returns:
            List of [fwd_start, fwd_len, rev_start, rev_len] quads.
            -1 means "no constraint" for that side.
        """
        if not windows:
            return []

        target_mid = target_start + target_length // 2
        pairs = []

        for win_start, win_len in windows:
            win_mid = win_start + win_len // 2

            if win_mid < target_mid:
                # Window is left of target → constrain forward primer there
                pairs.append([win_start, win_len, -1, -1])
            else:
                # Window is right of target → constrain reverse primer there
                pairs.append([-1, -1, win_start, win_len])

        return pairs

    def _run_primer3_tiered(
        self,
        region: DesignRegion,
        mode: DesignMode,
        disc_positions_genomic: Dict[int, int],
        num_homologous_regions: int,
        variant_group=None,
    ) -> Tuple[List[Amplicon], int, str]:
        """
        Run Primer3 with tiered homology constraints.

        Tier 1: SEQUENCE_OVERLAP_JUNCTION_LIST (strict — primer must
                 overlap a discriminating position)
        Tier 2: SEQUENCE_PRIMER_PAIR_OK_REGION_LIST (moderate — primer
                 must be in a discriminating window)
        Tier 3: Expanded pool + post-hoc re-ranking (fallback)

        Returns:
            (amplicons, tier, message)
        """
        from ..primer.homology_analyzer import HomologyAnalyzer

        # Convert to local coordinates
        disc_local = self._genomic_to_local(disc_positions_genomic, region)
        if not disc_local:
            self.log_info("[PrimerDesigner] No local discriminating positions — Tier 0")
            results = self._run_primer3(region, mode)
            amplicons = self._create_amplicons(results, variant_group, region, mode) if results else []
            return amplicons, 0, ""

        # Precompute windows for Tier 2
        windows = self._compute_discriminating_windows(
            disc_local, len(region.sequence)
        )
        junctions = self._select_junction_positions(disc_local, max_count=20)

        self.log_info(
            f"[PrimerDesigner] Homology constraints: "
            f"{len(disc_local)} local positions, "
            f"{len(windows)} windows, "
            f"{len(junctions)} junction candidates"
        )

        # ---- Tier 1: OVERLAP_JUNCTION_LIST ----
        if junctions:
            self.log_info(f"[PrimerDesigner] Tier 1: trying OVERLAP_JUNCTION_LIST "
                          f"with {len(junctions)} positions")
            results = self._run_primer3(
                region, mode,
                junction_list=junctions,
            )
            if results:
                amplicons = self._create_amplicons(results, variant_group, region, mode)
                if amplicons:
                    # Score the amplicons for display purposes
                    self._score_amplicons_homology(
                        amplicons, disc_positions_genomic, region)
                    self.log_info(
                        f"[PrimerDesigner] Tier 1 SUCCESS: {len(amplicons)} pairs "
                        f"with junction overlap")
                    return (
                        amplicons,
                        1,
                        f"Tier 1 (junction overlap): {len(amplicons)} primers "
                        f"directly overlap discriminating positions from "
                        f"{num_homologous_regions} homologous region(s)",
                    )

        # ---- Tier 2: PRIMER_PAIR_OK_REGION_LIST ----
        if windows:
            ok_pairs = self._build_ok_region_pairs(
                windows, region.target_start, region.target_length,
                len(region.sequence))
            if ok_pairs:
                self.log_info(
                    f"[PrimerDesigner] Tier 2: trying OK_REGION_LIST "
                    f"with {len(ok_pairs)} region pairs")
                results = self._run_primer3(
                    region, mode,
                    ok_region_list=ok_pairs,
                )
                if results:
                    amplicons = self._create_amplicons(results, variant_group, region, mode)
                    if amplicons:
                        self._score_amplicons_homology(
                            amplicons, disc_positions_genomic, region)
                        self.log_info(
                            f"[PrimerDesigner] Tier 2 SUCCESS: {len(amplicons)} pairs "
                            f"in discriminating windows")
                        return (
                            amplicons,
                            2,
                            f"Tier 2 (discriminating windows): {len(amplicons)} primers "
                            f"placed within discriminating regions from "
                            f"{num_homologous_regions} homologous region(s)",
                        )

        # ---- Tier 3: expanded pool + re-ranking ----
        self.log_info("[PrimerDesigner] Tier 3: fallback with expanded pool + re-ranking")
        results = self._run_primer3(region, mode, extra_candidates=True)
        if results:
            amplicons = self._create_amplicons(results, variant_group, region, mode)
            if amplicons:
                amplicons = self._rerank_by_homology(
                    amplicons, disc_positions_genomic, region)
                amplicons = amplicons[:self.parameters.num_primer_pairs]
                return (
                    amplicons,
                    3,
                    f"Tier 3 (fallback): primers re-ranked by homology discrimination "
                    f"but may not fully discriminate from "
                    f"{num_homologous_regions} homologous region(s). "
                    f"Other targets may also amplify but with different product sizes.",
                )

        return [], 3, "No primers found even with relaxed homology constraints"

    def _score_amplicons_homology(
        self,
        amplicons: List[Amplicon],
        disc_positions_genomic: Dict[int, int],
        region: DesignRegion,
    ) -> None:
        """Score amplicons by homology discrimination (for display, in-place)."""
        from ..primer.homology_analyzer import HomologyAnalyzer

        genomic_offset = region.genomic_start

        for amplicon in amplicons:
            pp = amplicon.primer_pair
            fwd_score, fwd_count = HomologyAnalyzer.score_primer_discrimination(
                primer_start=genomic_offset + pp.forward.start,
                primer_length=pp.forward.length,
                is_forward=True,
                discriminating_positions=disc_positions_genomic,
            )
            rev_score, rev_count = HomologyAnalyzer.score_primer_discrimination(
                primer_start=genomic_offset + pp.reverse.start,
                primer_length=pp.reverse.length,
                is_forward=False,
                discriminating_positions=disc_positions_genomic,
            )
            pp.homology_discrimination_score = fwd_score + rev_score
            pp.fwd_discriminating_positions = fwd_count
            pp.rev_discriminating_positions = rev_count

            if fwd_count == 0 and rev_count == 0:
                pp.homology_warning = (
                    "Primers do not cover discriminating positions")
            elif fwd_count == 0:
                pp.homology_warning = (
                    "Forward primer does not cover discriminating positions")
            elif rev_count == 0:
                pp.homology_warning = (
                    "Reverse primer does not cover discriminating positions")

    def _rerank_by_homology(
        self,
        amplicons: List[Amplicon],
        discriminating_positions: Dict[int, int],
        region: DesignRegion,
    ) -> List[Amplicon]:
        """
        Re-rank amplicons by how well their primers discriminate against
        homologous regions (pseudogenes).

        Primers that cover positions where our sequence differs from
        pseudogenes receive higher scores, especially when mismatches
        fall near the 3' end of the primer.

        Args:
            amplicons: Amplicon candidates from Primer3
            discriminating_positions: Map of query_pos → homolog_mismatch_count
            region: Design region (for coordinate context)

        Returns:
            Re-ranked list of amplicons (best discrimination first)
        """
        from ..primer.homology_analyzer import HomologyAnalyzer

        scored_amplicons = []

        # discriminating_positions keys are now in GENOMIC coordinates
        # (from extract_discriminating_positions with query_genomic_start offset).
        # Primer positions from Primer3 are 0-based in the design region,
        # so we must add region.genomic_start to convert to genomic coords.
        genomic_offset = region.genomic_start

        self.log_debug(
            f"[PrimerDesigner] Homology re-rank: genomic_offset={genomic_offset}, "
            f"disc_positions range={min(discriminating_positions.keys()) if discriminating_positions else 'N/A'}"
            f"-{max(discriminating_positions.keys()) if discriminating_positions else 'N/A'}"
        )

        for amplicon in amplicons:
            pp = amplicon.primer_pair

            # Convert primer positions from design-region-relative to genomic
            fwd_genomic_start = genomic_offset + pp.forward.start
            rev_genomic_start = genomic_offset + pp.reverse.start

            # Score forward primer
            fwd_score, fwd_count = HomologyAnalyzer.score_primer_discrimination(
                primer_start=fwd_genomic_start,
                primer_length=pp.forward.length,
                is_forward=True,
                discriminating_positions=discriminating_positions,
            )

            # Score reverse primer
            rev_score, rev_count = HomologyAnalyzer.score_primer_discrimination(
                primer_start=rev_genomic_start,
                primer_length=pp.reverse.length,
                is_forward=False,
                discriminating_positions=discriminating_positions,
            )

            total_score = fwd_score + rev_score

            # Store scores in primer pair
            pp.homology_discrimination_score = total_score
            pp.fwd_discriminating_positions = fwd_count
            pp.rev_discriminating_positions = rev_count

            # Generate warning if no discrimination
            if total_score == 0:
                pp.homology_warning = (
                    "Primers do not cover any positions that differentiate "
                    "from homologous regions — may amplify pseudogene(s)"
                )
            elif fwd_count == 0:
                pp.homology_warning = (
                    "Forward primer does not cover discriminating positions"
                )
            elif rev_count == 0:
                pp.homology_warning = (
                    "Reverse primer does not cover discriminating positions"
                )

            scored_amplicons.append((total_score, amplicon))

        # Sort by discrimination score descending, preserving Primer3 rank
        # as a tiebreaker (earlier = better Primer3 quality)
        scored_amplicons.sort(key=lambda x: -x[0])

        reranked = [amp for _, amp in scored_amplicons]

        # Log re-ranking summary
        if reranked:
            best = reranked[0].primer_pair
            self.log_info(
                f"[PrimerDesigner] Homology re-ranking: best pair has "
                f"score={best.homology_discrimination_score:.1f} "
                f"(fwd={best.fwd_discriminating_positions}, "
                f"rev={best.rev_discriminating_positions})"
            )

        return reranked

    def _generate_suggestions(
        self,
        region: DesignRegion,
        mode: DesignMode
    ) -> List[str]:
        """Generate suggestions when no primers are found."""
        suggestions = []

        if mode == DesignMode.QPCR:
            suggestions.append(
                "Try increasing the maximum amplicon size (current: "
                f"{self.parameters.qpcr_max_amplicon_size} bp)"
            )
        else:
            suggestions.append(
                "Try increasing the maximum amplicon size (current: "
                f"{self.parameters.max_amplicon_size} bp)"
            )

        suggestions.append(
            f"Try relaxing Tm constraints (current: "
            f"{self.parameters.primer_min_tm}-{self.parameters.primer_max_tm}°C)"
        )

        if self.parameters.filter_population_variants:
            suggestions.append(
                "Try disabling population variant filtering or increasing MAF threshold"
            )

        if self.parameters.use_variant_distance_constraint:
            suggestions.append(
                f"Try reducing minimum distance from variant (current: "
                f"{self.parameters.min_distance_from_variant} bp)"
            )

        return suggestions

    def calculate_tm(self, sequence: str) -> float:
        """
        Calculate melting temperature for a sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Melting temperature in Celsius
        """
        return primer3.calc_tm(
            sequence,
            mv_conc=self.parameters.mv_conc,
            dv_conc=self.parameters.dv_conc,
            dntp_conc=self.parameters.dntp_conc,
            dna_conc=self.parameters.dna_conc
        )

    def check_hairpin(self, sequence: str) -> Dict[str, float]:
        """
        Check hairpin potential for a sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Dictionary with hairpin metrics
        """
        result = primer3.calc_hairpin(
            sequence,
            mv_conc=self.parameters.mv_conc,
            dv_conc=self.parameters.dv_conc,
            dntp_conc=self.parameters.dntp_conc
        )
        return {
            'tm': result.tm,
            'dg': result.dg,
            'dh': result.dh,
            'ds': result.ds
        }

    def check_dimer(self, seq1: str, seq2: str) -> Dict[str, float]:
        """
        Check dimer formation potential between two sequences.

        Args:
            seq1: First sequence
            seq2: Second sequence

        Returns:
            Dictionary with dimer metrics
        """
        result = primer3.calc_heterodimer(
            seq1, seq2,
            mv_conc=self.parameters.mv_conc,
            dv_conc=self.parameters.dv_conc,
            dntp_conc=self.parameters.dntp_conc
        )
        return {
            'tm': result.tm,
            'dg': result.dg,
            'dh': result.dh,
            'ds': result.ds
        }

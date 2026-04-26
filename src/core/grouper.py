"""
Variant grouper for identifying variants that can share a common amplicon.
Groups variants based on genomic distance and amplicon size constraints.
"""

from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

from .models import Variant, GenomicPosition, DesignParameters
from ..utils.logger import LoggerMixin


@dataclass
class VariantGroup:
    """A group of variants that can share a common amplicon."""
    group_id: int
    variants: List[Variant]
    genomic_start: int
    genomic_end: int
    chromosome: str
    total_span: int
    can_use_single_amplicon: bool

    @property
    def variant_count(self) -> int:
        return len(self.variants)

    def add_variant(self, variant: Variant):
        """Add a variant to the group and update boundaries."""
        self.variants.append(variant)
        variant.group_id = self.group_id

        if variant.genomic_position:
            self.genomic_start = min(self.genomic_start, variant.genomic_position.start)
            self.genomic_end = max(self.genomic_end, variant.genomic_position.end)
            self.total_span = self.genomic_end - self.genomic_start + 1


class VariantGrouper(LoggerMixin):
    """
    Groups variants that can share a common amplicon based on genomic proximity.
    Considers actual genomic distance (including introns), not just CDS position.
    """

    def __init__(self, parameters: DesignParameters = None):
        """
        Initialize variant grouper.

        Args:
            parameters: Design parameters including amplicon size limits
        """
        self.parameters = parameters or DesignParameters()

    def group_variants(self, variants: List[Variant]) -> List[VariantGroup]:
        """
        Group variants that can share a common amplicon.

        Args:
            variants: List of validated variants with genomic positions

        Returns:
            List of VariantGroup objects
        """
        if not variants:
            return []

        # Filter variants with valid genomic positions
        valid_variants = [v for v in variants if v.genomic_position]

        if not valid_variants:
            # Return each variant as its own group
            return self._create_individual_groups(variants)

        # Sort by chromosome and position
        sorted_variants = sorted(
            valid_variants,
            key=lambda v: (v.genomic_position.chromosome, v.genomic_position.start)
        )

        groups = []
        current_group = None

        for variant in sorted_variants:
            if not current_group:
                # Start new group
                current_group = self._create_group(len(groups), variant)
                groups.append(current_group)
            elif self._can_add_to_group(current_group, variant):
                # Add to existing group
                current_group.add_variant(variant)
                self._update_group_feasibility(current_group)
            else:
                # Start new group
                current_group = self._create_group(len(groups), variant)
                groups.append(current_group)

        # Add variants without positions as individual groups
        for variant in variants:
            if not variant.genomic_position:
                group = VariantGroup(
                    group_id=len(groups),
                    variants=[variant],
                    genomic_start=0,
                    genomic_end=0,
                    chromosome="",
                    total_span=0,
                    can_use_single_amplicon=True
                )
                variant.group_id = group.group_id
                groups.append(group)

        self.log_info(f"Grouped {len(variants)} variants into {len(groups)} groups")
        return groups

    def _create_group(self, group_id: int, variant: Variant) -> VariantGroup:
        """Create a new group with initial variant."""
        group = VariantGroup(
            group_id=group_id,
            variants=[variant],
            genomic_start=variant.genomic_position.start if variant.genomic_position else 0,
            genomic_end=variant.genomic_position.end if variant.genomic_position else 0,
            chromosome=variant.genomic_position.chromosome if variant.genomic_position else "",
            total_span=0,
            can_use_single_amplicon=True
        )
        variant.group_id = group_id
        return group

    def _can_add_to_group(self, group: VariantGroup, variant: Variant) -> bool:
        """
        Check if a variant can be added to an existing group.

        Criteria:
        1. Same chromosome
        2. Combined span fits within max amplicon size (with room for primers)
        """
        if not variant.genomic_position:
            return False

        # Must be same chromosome
        if variant.genomic_position.chromosome != group.chromosome:
            return False

        # Calculate potential new span (inclusive range)
        new_start = min(group.genomic_start, variant.genomic_position.start)
        new_end = max(group.genomic_end, variant.genomic_position.end)
        new_span = new_end - new_start + 1

        # Account for primer placement (need room for primers outside variant region)
        # Minimum distance from variant to primer
        primer_margin = self.parameters.min_distance_from_variant * 2

        # Total required amplicon size
        required_size = new_span + primer_margin

        # Check if fits within max amplicon size
        return required_size <= self.parameters.max_amplicon_size

    def _update_group_feasibility(self, group: VariantGroup):
        """Update the group's feasibility for single amplicon design."""
        if group.total_span <= 0:
            group.can_use_single_amplicon = True
            return

        # Calculate required amplicon size
        primer_margin = self.parameters.min_distance_from_variant * 2
        required_size = group.total_span + primer_margin

        group.can_use_single_amplicon = (
            required_size >= self.parameters.min_amplicon_size and
            required_size <= self.parameters.max_amplicon_size
        )

    def _create_individual_groups(self, variants: List[Variant]) -> List[VariantGroup]:
        """Create individual groups for each variant."""
        groups = []
        for i, variant in enumerate(variants):
            group = VariantGroup(
                group_id=i,
                variants=[variant],
                genomic_start=variant.genomic_position.start if variant.genomic_position else 0,
                genomic_end=variant.genomic_position.end if variant.genomic_position else 0,
                chromosome=variant.genomic_position.chromosome if variant.genomic_position else "",
                total_span=0,
                can_use_single_amplicon=True
            )
            variant.group_id = i
            groups.append(group)
        return groups

    def suggest_groupings(
        self,
        variants: List[Variant],
        max_distance: int = None
    ) -> Dict[int, List[int]]:
        """
        Suggest potential variant groupings without modifying variants.

        Args:
            variants: List of variants
            max_distance: Maximum genomic distance for grouping (optional)

        Returns:
            Dictionary mapping group ID to list of variant row numbers
        """
        if max_distance is None:
            max_distance = self.parameters.max_amplicon_size - (
                self.parameters.min_distance_from_variant * 2
            )

        # Filter variants with valid positions
        positioned = [(i, v) for i, v in enumerate(variants) if v.genomic_position]

        if not positioned:
            return {i: [v.row_number] for i, v in enumerate(variants)}

        # Sort by chromosome and position
        positioned.sort(key=lambda x: (
            x[1].genomic_position.chromosome,
            x[1].genomic_position.start
        ))

        groupings = {}
        current_group_id = 0
        current_members = []
        current_chr = None
        current_end = 0

        for idx, variant in positioned:
            pos = variant.genomic_position

            if (current_chr != pos.chromosome or
                pos.start - current_end > max_distance or
                not current_members):
                # Start new group
                if current_members:
                    groupings[current_group_id] = current_members
                    current_group_id += 1

                current_members = [variant.row_number]
                current_chr = pos.chromosome
                current_end = pos.end
            else:
                # Add to current group
                current_members.append(variant.row_number)
                current_end = max(current_end, pos.end)

        # Add last group
        if current_members:
            groupings[current_group_id] = current_members

        return groupings

    def calculate_group_statistics(self, groups: List[VariantGroup]) -> Dict:
        """
        Calculate statistics about variant groups.

        Args:
            groups: List of variant groups

        Returns:
            Dictionary with statistics
        """
        total_variants = sum(g.variant_count for g in groups)
        single_amplicon_groups = sum(1 for g in groups if g.can_use_single_amplicon)
        multi_variant_groups = sum(1 for g in groups if g.variant_count > 1)

        spans = [g.total_span for g in groups if g.total_span > 0]
        avg_span = sum(spans) / len(spans) if spans else 0

        return {
            'total_groups': len(groups),
            'total_variants': total_variants,
            'single_amplicon_feasible': single_amplicon_groups,
            'multi_variant_groups': multi_variant_groups,
            'average_group_span': avg_span,
            'max_group_span': max(spans) if spans else 0,
            'min_group_span': min(spans) if spans else 0
        }

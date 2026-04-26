"""
Visualization components for gene structure and sequence display.
"""

import tkinter as tk
from typing import List, Optional, Tuple
from ..core.models import Transcript, Exon, GenomicPosition, PopulationVariant


class GeneVisualizer:
    """
    Visualizes gene structure with exons, introns, and variant positions.
    """

    def __init__(self, canvas: tk.Canvas):
        """
        Initialize gene visualizer.

        Args:
            canvas: tkinter Canvas widget to draw on
        """
        self.canvas = canvas

        # Colors (dark-theme friendly)
        self.exon_color = "#89b4fa"
        self.intron_color = "#585b70"
        self.variant_color = "#f38ba8"
        self.primer_color = "#a6e3a1"
        self.text_color = "#cdd6f4"

        # Dimensions
        self.margin = 40
        self.exon_height = 30
        self.intron_height = 2

    def draw_gene_structure(
        self,
        transcript: Transcript,
        highlight_exon: Exon = None,
        variant_position: GenomicPosition = None,
        primer_positions: List[Tuple[int, int, bool]] = None
    ):
        """
        Draw gene structure visualization.

        Args:
            transcript: Transcript with exon information
            highlight_exon: Exon containing the variant (to highlight)
            variant_position: Position of the variant
            primer_positions: List of (start, end, is_forward) for primers
        """
        self.canvas.delete("all")

        # If no data at all, show placeholder
        if not variant_position and (not transcript or not transcript.exons):
            self._draw_placeholder()
            return

        # If no exons but we have variant position, show simplified variant-centered view
        if (not transcript or not transcript.exons) and variant_position:
            self._draw_variant_only(transcript, variant_position)
            return

        # Get canvas dimensions
        width = self.canvas.winfo_width()
        height = self.canvas.winfo_height()

        if width < 100:
            width = 800  # Default width
        if height < 50:
            height = 150  # Default height

        # Calculate scale
        exons = sorted(transcript.exons, key=lambda e: e.genomic_start)
        gene_start = exons[0].genomic_start
        gene_end = exons[-1].genomic_end
        gene_length = gene_end - gene_start

        drawable_width = width - 2 * self.margin
        scale = drawable_width / gene_length if gene_length > 0 else 1

        # Draw baseline (intron line)
        y_center = height // 2
        self.canvas.create_line(
            self.margin, y_center,
            width - self.margin, y_center,
            fill=self.intron_color, width=self.intron_height
        )

        # Draw exons
        for exon in exons:
            x1 = self.margin + (exon.genomic_start - gene_start) * scale
            x2 = self.margin + (exon.genomic_end - gene_start) * scale
            y1 = y_center - self.exon_height // 2
            y2 = y_center + self.exon_height // 2

            # Highlight if this is the variant exon
            color = self.exon_color
            if highlight_exon and exon.number == highlight_exon.number:
                color = "#9b59b6"  # Purple for highlighted exon

            self.canvas.create_rectangle(
                x1, y1, x2, y2,
                fill=color, outline="#2980b9", width=2
            )

            # Exon number
            self.canvas.create_text(
                (x1 + x2) / 2, y1 - 10,
                text=f"E{exon.number}",
                fill=self.text_color, font=("Helvetica", 8)
            )

        # Draw variant position
        if variant_position:
            var_x = self.margin + (variant_position.start - gene_start) * scale
            self.canvas.create_line(
                var_x, y_center - self.exon_height,
                var_x, y_center + self.exon_height,
                fill=self.variant_color, width=2
            )
            self.canvas.create_polygon(
                var_x - 5, y_center - self.exon_height - 10,
                var_x + 5, y_center - self.exon_height - 10,
                var_x, y_center - self.exon_height,
                fill=self.variant_color
            )
            self.canvas.create_text(
                var_x, y_center - self.exon_height - 20,
                text="Variant",
                fill=self.variant_color, font=("Helvetica", 8, "bold")
            )

        # Draw primer positions
        if primer_positions:
            for start, end, is_forward in primer_positions:
                x1 = self.margin + (start - gene_start) * scale
                x2 = self.margin + (end - gene_start) * scale
                y = y_center + self.exon_height + 10 if is_forward else y_center + self.exon_height + 25

                self.canvas.create_line(
                    x1, y, x2, y,
                    fill=self.primer_color, width=3, arrow=tk.LAST
                )
                label = "F" if is_forward else "R"
                self.canvas.create_text(
                    (x1 + x2) / 2, y + 10,
                    text=label, fill=self.primer_color, font=("Helvetica", 8, "bold")
                )

        # Draw scale bar
        self._draw_scale_bar(width, height, gene_length)

        # Draw legend
        self._draw_legend(width, height)

    def _draw_placeholder(self):
        """Draw placeholder when no data is available."""
        width = self.canvas.winfo_width() or 800
        height = self.canvas.winfo_height() or 150

        self.canvas.create_text(
            width // 2, height // 2,
            text="Select a variant to view gene structure",
            fill="#7f8c8d", font=("Helvetica", 12)
        )

    def _draw_variant_only(self, transcript: Transcript, variant_position: GenomicPosition):
        """Draw simplified view when only variant position is known (no exon data)."""
        width = self.canvas.winfo_width() or 800
        height = self.canvas.winfo_height() or 150

        if width < 100:
            width = 800
        if height < 50:
            height = 150

        y_center = height // 2

        # Draw chromosome line
        self.canvas.create_line(
            self.margin, y_center,
            width - self.margin, y_center,
            fill=self.intron_color, width=3
        )

        # Draw variant marker in center
        var_x = width // 2
        self.canvas.create_line(
            var_x, y_center - 25,
            var_x, y_center + 25,
            fill=self.variant_color, width=3
        )
        self.canvas.create_polygon(
            var_x - 8, y_center - 35,
            var_x + 8, y_center - 35,
            var_x, y_center - 25,
            fill=self.variant_color
        )

        # Variant info text
        self.canvas.create_text(
            var_x, y_center - 50,
            text="Target Variant",
            fill=self.variant_color, font=("Helvetica", 10, "bold")
        )

        # Show genomic position
        pos_text = f"Chr{variant_position.chromosome}:{variant_position.start:,}"
        self.canvas.create_text(
            var_x, y_center + 40,
            text=pos_text,
            fill=self.text_color, font=("Helvetica", 9)
        )

        # Transcript info
        if transcript:
            info_text = f"{transcript.gene_symbol} - {transcript.full_accession}"
            self.canvas.create_text(
                width // 2, 20,
                text=info_text,
                fill=self.text_color, font=("Helvetica", 11, "bold")
            )

    def _draw_scale_bar(self, width: int, height: int, gene_length: int):
        """Draw a scale bar."""
        # Determine appropriate scale
        if gene_length > 10000:
            scale_bp = 5000
            label = "5 kb"
        elif gene_length > 1000:
            scale_bp = 1000
            label = "1 kb"
        else:
            scale_bp = 100
            label = "100 bp"

        drawable_width = width - 2 * self.margin
        scale_width = (scale_bp / gene_length) * drawable_width if gene_length > 0 else 50

        x1 = width - self.margin - scale_width
        x2 = width - self.margin
        y = height - 20

        self.canvas.create_line(x1, y, x2, y, fill=self.text_color, width=2)
        self.canvas.create_line(x1, y - 3, x1, y + 3, fill=self.text_color, width=2)
        self.canvas.create_line(x2, y - 3, x2, y + 3, fill=self.text_color, width=2)
        self.canvas.create_text(
            (x1 + x2) / 2, y + 12,
            text=label, fill=self.text_color, font=("Helvetica", 8)
        )

    def _draw_legend(self, width: int, height: int):
        """Draw legend."""
        x = 10
        y = height - 20

        # Exon
        self.canvas.create_rectangle(x, y - 5, x + 15, y + 5, fill=self.exon_color)
        self.canvas.create_text(x + 25, y, text="Exon", anchor="w", font=("Helvetica", 8))

        # Variant
        x += 70
        self.canvas.create_line(x, y - 5, x, y + 5, fill=self.variant_color, width=2)
        self.canvas.create_text(x + 10, y, text="Variant", anchor="w", font=("Helvetica", 8))

        # Primer
        x += 70
        self.canvas.create_line(x, y, x + 15, y, fill=self.primer_color, width=3, arrow=tk.LAST)
        self.canvas.create_text(x + 25, y, text="Primer", anchor="w", font=("Helvetica", 8))


class SequenceVisualizer:
    """
    Visualizes DNA sequences with variant and polymorphism highlighting.
    Uses IUPAC ambiguity codes like Ensembl for displaying SNPs.
    """

    # IUPAC ambiguity codes for combining alleles
    IUPAC_CODES = {
        frozenset(['A', 'G']): 'R',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S',
        frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N',
    }

    def __init__(self, text_widget: tk.Text):
        """
        Initialize sequence visualizer.

        Args:
            text_widget: tkinter Text widget to display sequence
        """
        self.text = text_widget

        # State for collapsible sections
        self._variant_details_expanded = False
        self._variant_details_content = ""
        self._variant_details_start_mark = None

        # Configure tags - High visibility colors for dark mode compatibility
        # Font size 10 to match the rest of the interface
        # Flanking region (outside gene) - dim gray
        self.text.tag_configure("normal", foreground="#AAAAAA", font=("Courier", 10))
        # Exon region - bright white, bold (highest emphasis)
        self.text.tag_configure("exon", foreground="#FFFFFF", font=("Courier", 10, "bold"))
        # Intron region - medium brightness, distinguishable from exon and flanking
        self.text.tag_configure("intron", foreground="#CCCCCC", font=("Courier", 10))
        # Legacy alias: "transcript" maps to exon style for backward compatibility
        self.text.tag_configure("transcript", foreground="#FFFFFF", font=("Courier", 10, "bold"))
        # Target variant - red background, white text
        self.text.tag_configure("variant", foreground="white", background="#e74c3c", font=("Courier", 10, "bold"))
        # IUPAC ambiguity code colors (like Ensembl - cyan/teal background)
        self.text.tag_configure("iupac", foreground="#000000", background="#7FDBFF", font=("Courier", 10, "bold"))
        self.text.tag_configure("iupac_flank", foreground="#000000", background="#ADD8E6", font=("Courier", 10))
        # Primer/probe regions
        self.text.tag_configure("primer_fwd", foreground="white", background="#27ae60", font=("Courier", 10, "bold"))
        self.text.tag_configure("primer_rev", foreground="white", background="#3498db", font=("Courier", 10, "bold"))
        self.text.tag_configure("probe", foreground="white", background="#9b59b6", font=("Courier", 10, "bold"))
        # Position markers
        self.text.tag_configure("position", foreground="#BBBBBB", font=("Courier", 9))
        self.text.tag_configure("header", foreground="#DDDDDD", font=("Helvetica", 10, "bold"))
        # Clickable expand/collapse header
        self.text.tag_configure("expandable", foreground="#00BFFF", font=("Helvetica", 10, "bold underline"))

        # Bind click event for expandable sections
        self.text.tag_bind("expandable", "<Button-1>", self._toggle_variant_details)
        # Change cursor on hover over expandable text
        self.text.tag_bind("expandable", "<Enter>", lambda e: self.text.config(cursor="hand2"))
        self.text.tag_bind("expandable", "<Leave>", lambda e: self.text.config(cursor=""))

    def _get_iupac_code(self, ref_base: str, alt_alleles: List[str]) -> str:
        """
        Get IUPAC ambiguity code for a position with variants.

        Args:
            ref_base: Reference base
            alt_alleles: List of alternative alleles

        Returns:
            IUPAC code representing all alleles
        """
        # Collect all alleles (ref + alts), only single nucleotides
        alleles = set()
        if len(ref_base) == 1 and ref_base.upper() in 'ACGT':
            alleles.add(ref_base.upper())

        for alt in alt_alleles:
            if len(alt) == 1 and alt.upper() in 'ACGT':
                alleles.add(alt.upper())

        if len(alleles) <= 1:
            return ref_base  # No ambiguity

        # Look up IUPAC code
        allele_set = frozenset(alleles)
        return self.IUPAC_CODES.get(allele_set, 'N')

    def display_sequence(
        self,
        sequence: str,
        variant_position: int = None,
        variant_length: int = 1,
        population_variants: List[PopulationVariant] = None,
        primer_regions: List[Tuple[int, int, str]] = None,
        line_length: int = 60,
        masked_positions: set = None,
        seq_start_genomic: int = None,
        selected_populations: List[str] = None,
        transcript_info: dict = None,
        maf_display_threshold: float = 0.0,
        exon_regions: List[Tuple[int, int, int]] = None
    ):
        """
        Display a DNA sequence with IUPAC ambiguity codes for SNPs (like Ensembl).

        Args:
            sequence: DNA sequence string
            variant_position: Position of the target variant (0-based in sequence)
            variant_length: Length of the variant
            population_variants: List of population variants in the region
            primer_regions: List of (start, end, type) for primers/probes
            line_length: Characters per line
            masked_positions: Set of 0-based positions that are masked with 'N'
            seq_start_genomic: Genomic start position of the sequence (for mapping population variants)
            selected_populations: List of population codes to display MAF for
            transcript_info: Dictionary with transcript info for marking CDS region
            maf_display_threshold: MAF threshold for displaying IUPAC codes (variants below threshold show ref base)
        """
        self.text.configure(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)

        if not sequence:
            self.text.insert(tk.END, "No sequence available", "header")
            self.text.configure(state=tk.DISABLED)
            return

        # Build map of position -> list of variants for IUPAC code generation
        poly_details = {}  # Map seq_position -> list of PopulationVariant

        if masked_positions:
            for pos in masked_positions:
                if pos not in poly_details:
                    poly_details[pos] = []

        # Convert population variants to sequence positions
        if population_variants and seq_start_genomic is not None:
            for var in population_variants:
                # Calculate position in sequence (0-based)
                seq_pos = var.position - seq_start_genomic
                if 0 <= seq_pos < len(sequence):
                    if seq_pos not in poly_details:
                        poly_details[seq_pos] = []
                    poly_details[seq_pos].append(var)

        # Build primer/probe region set
        primer_positions = {}
        if primer_regions:
            for start, end, region_type in primer_regions:
                for i in range(start, end + 1):
                    primer_positions[i] = region_type

        # Build exon position set for O(1) lookup
        # exon_positions[i] = exon_number if position i is exonic, else 0
        exon_positions = {}  # pos -> exon_number
        if exon_regions:
            for ex_start, ex_end, ex_num in exon_regions:
                for pos in range(max(0, ex_start), min(len(sequence), ex_end + 1)):
                    exon_positions[pos] = ex_num

        # Determine if we have any exon data to distinguish exon/intron
        has_exon_data = bool(exon_positions)

        # Header - single line with semicolons, wraps automatically
        header_parts = [f"Length: {len(sequence)} bp"]
        if transcript_info:
            header_parts.append(f"Gene: {transcript_info.get('gene_symbol', '?')}")
            header_parts.append(f"Transcript: {transcript_info.get('accession', '?')}")
        if variant_position is not None:
            header_parts.append(f"Target pos: {variant_position + 1}")
        if seq_start_genomic:
            header_parts.append(f"Genomic: {seq_start_genomic:,}-{seq_start_genomic + len(sequence):,}")
        if population_variants:
            header_parts.append(f"SNPs: {len(population_variants)}")
        if exon_regions:
            exon_nums = sorted(set(n for _, _, n in exon_regions if n > 0))
            if exon_nums:
                exon_list = ", ".join(f"E{n}" for n in exon_nums)
                header_parts.append(f"Exons in view: {exon_list}")
        self.text.insert(tk.END, "; ".join(header_parts) + "\n", "header")

        # Legend with IUPAC codes
        self.text.insert(tk.END, "Legend: ", "header")
        self.text.insert(tk.END, "TARGET ", "variant")
        self.text.insert(tk.END, "SNP ", "iupac")
        if has_exon_data:
            self.text.insert(tk.END, "EXON ", "exon")
            self.text.insert(tk.END, "intron ", "intron")
        else:
            self.text.insert(tk.END, "TRANSCRIPT ", "transcript")
        self.text.insert(tk.END, "flanking ", "normal")
        self.text.insert(tk.END, "\n")
        self.text.insert(tk.END, "IUPAC: R=A/G  Y=C/T  S=G/C  W=A/T  K=G/T  M=A/C\n\n", "position")

        # Display sequence in lines
        for line_start in range(0, len(sequence), line_length):
            line_end = min(line_start + line_length, len(sequence))

            # Position indicator (genomic position if available)
            if seq_start_genomic:
                genomic_pos = seq_start_genomic + line_start
                pos_str = f"{genomic_pos:>10}  "
            else:
                pos_str = f"{line_start + 1:6d}  "
            self.text.insert(tk.END, pos_str, "position")

            # Sequence with highlighting and IUPAC codes
            for i in range(line_start, line_end):
                ref_base = sequence[i]
                display_base = ref_base

                # Determine base tag: exon (bold white) / intron (medium) / flanking (dim gray)
                if has_exon_data:
                    if i in exon_positions:
                        base_tag = "exon"
                        is_exonic = True
                    else:
                        base_tag = "intron"
                        is_exonic = False
                else:
                    # No exon data — fall back to legacy "transcript" tag
                    base_tag = "transcript"
                    is_exonic = True

                # Check for target variant (highest priority)
                if variant_position is not None and variant_position <= i < variant_position + variant_length:
                    tag = "variant"
                # Check for primer/probe regions (second priority)
                elif i in primer_positions:
                    region_type = primer_positions[i]
                    if region_type == "forward":
                        tag = "primer_fwd"
                    elif region_type == "reverse":
                        tag = "primer_rev"
                    elif region_type == "probe":
                        tag = "probe"
                    else:
                        tag = base_tag
                # Check for polymorphisms - display as IUPAC codes (if MAF >= threshold)
                elif i in poly_details and poly_details[i]:
                    variants_at_pos = poly_details[i]
                    # Collect all alternative alleles that pass the MAF threshold
                    alt_alleles = []
                    for var in variants_at_pos:
                        if var.alt and len(var.alt) == 1:  # Only SNPs
                            # Check if variant passes MAF display threshold
                            var_maf = var.maf_global
                            # Also check population-specific MAF if selected
                            if selected_populations and var.maf_by_population:
                                for pop in selected_populations:
                                    if pop in var.maf_by_population:
                                        var_maf = max(var_maf, var.maf_by_population.get(pop, 0))

                            # Only show IUPAC code if MAF >= threshold
                            if var_maf >= maf_display_threshold:
                                alt_alleles.append(var.alt)

                    if alt_alleles:
                        display_base = self._get_iupac_code(ref_base, alt_alleles)
                        tag = "iupac" if is_exonic else "iupac_flank"
                    else:
                        tag = base_tag
                else:
                    tag = base_tag

                self.text.insert(tk.END, display_base, tag)

            self.text.insert(tk.END, "\n")

        # Summary of population variants - collapsible section
        if population_variants:
            self.text.insert(tk.END, f"\n{'='*60}\n", "header")

            # Store variant details content for toggling
            self._build_variant_details_content(population_variants, seq_start_genomic, selected_populations)

            # Clickable header to expand/collapse
            arrow = "▶" if not self._variant_details_expanded else "▼"
            self.text.insert(tk.END, f"{arrow} ", "expandable")
            self.text.insert(tk.END, f"Population Variants Detail ({len(population_variants)} total) - Click to expand\n", "expandable")

            # Mark where details should be inserted
            self.text.mark_set("variant_details_start", tk.END)
            self.text.mark_gravity("variant_details_start", tk.LEFT)

            # If expanded, show the details
            if self._variant_details_expanded:
                self.text.insert(tk.END, self._variant_details_content)

        self.text.configure(state=tk.DISABLED)

    def _build_variant_details_content(self, population_variants: List[PopulationVariant],
                                        seq_start_genomic: int, selected_populations: List[str]):
        """Build the variant details content for the collapsible section."""
        lines = []

        # Show which populations are selected
        if selected_populations:
            pop_names = {'global': 'Global', 'nfe': 'European (NFE)', 'afr': 'African',
                       'eas': 'East Asian', 'sas': 'South Asian', 'amr': 'Latino'}
            pop_display = ', '.join(pop_names.get(p, p) for p in selected_populations)
            lines.append(f"Selected populations: {pop_display}\n\n")

        # Sort variants by position
        sorted_variants = sorted(population_variants, key=lambda v: v.position)

        def format_maf(maf):
            if maf == 0:
                return "0"
            elif maf < 0.0001:
                return f"{maf:.2e}"
            else:
                return f"{maf:.6f}"

        for var in sorted_variants[:30]:  # Show up to 30 variants
            # Calculate position in sequence
            if seq_start_genomic:
                seq_pos = var.position - seq_start_genomic + 1
                pos_display = f"seq pos {seq_pos:>4}"
            else:
                pos_display = f"pos {var.position}"

            # Get IUPAC code for this variant
            ref_base = var.ref if var.ref else '?'
            iupac = self._get_iupac_code(ref_base, [var.alt] if var.alt else [])

            maf_parts = [f"Global:{format_maf(var.maf_global)}"]
            if selected_populations and var.maf_by_population:
                for pop in selected_populations:
                    if pop != 'global' and pop in var.maf_by_population:
                        maf_parts.append(f"{pop.upper()}:{format_maf(var.maf_by_population[pop])}")

            maf_display = " ".join(maf_parts)

            lines.append(f"  {var.rsid or 'unknown':15} ({pos_display}) {ref_base}>{var.alt} [{iupac}] MAF: {maf_display}\n")

        if len(population_variants) > 30:
            lines.append(f"\n  ... and {len(population_variants) - 30} more variants\n")

        self._variant_details_content = "".join(lines)

    def _toggle_variant_details(self, event=None):
        """Toggle the visibility of variant details section."""
        self._variant_details_expanded = not self._variant_details_expanded

        self.text.configure(state=tk.NORMAL)

        try:
            # Find and update the arrow and header text
            # Search for the expandable tag
            ranges = self.text.tag_ranges("expandable")
            if ranges:
                start = ranges[0]
                # Update the arrow
                arrow = "▼" if self._variant_details_expanded else "▶"
                # Delete old arrow and insert new one
                self.text.delete(start, f"{start}+2c")
                self.text.insert(start, f"{arrow} ", "expandable")

            # Handle the content toggle
            if "variant_details_start" in self.text.mark_names():
                mark_pos = self.text.index("variant_details_start")

                if self._variant_details_expanded:
                    # Insert the details content
                    self.text.insert(mark_pos, self._variant_details_content)
                else:
                    # Remove the details content - find end and delete
                    content_len = len(self._variant_details_content)
                    if content_len > 0:
                        end_pos = f"{mark_pos}+{content_len}c"
                        self.text.delete(mark_pos, end_pos)

        finally:
            self.text.configure(state=tk.DISABLED)

    def display_amplicon(
        self,
        sequence: str,
        forward_primer: str,
        reverse_primer: str,
        probe: str = None,
        variant_position: int = None
    ):
        """
        Display an amplicon with primers highlighted.

        Args:
            sequence: Amplicon sequence
            forward_primer: Forward primer sequence
            reverse_primer: Reverse primer sequence
            probe: Probe sequence (optional)
            variant_position: Position of variant in amplicon
        """
        self.text.configure(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)

        # Find primer positions
        fwd_start = sequence.find(forward_primer)
        fwd_end = fwd_start + len(forward_primer) if fwd_start >= 0 else -1

        # Reverse primer is reverse complement
        rev_comp = self._reverse_complement(reverse_primer)
        rev_start = sequence.find(rev_comp)
        rev_end = rev_start + len(rev_comp) if rev_start >= 0 else -1

        # Find probe position
        probe_start = -1
        probe_end = -1
        if probe:
            probe_start = sequence.find(probe)
            if probe_start < 0:
                probe_start = sequence.find(self._reverse_complement(probe))
            probe_end = probe_start + len(probe) if probe_start >= 0 else -1

        # Build primer regions
        primer_regions = []
        if fwd_start >= 0:
            primer_regions.append((fwd_start, fwd_end - 1, "forward"))
        if rev_start >= 0:
            primer_regions.append((rev_start, rev_end - 1, "reverse"))
        if probe_start >= 0:
            primer_regions.append((probe_start, probe_end - 1, "probe"))

        # Display with highlighting
        self.display_sequence(
            sequence,
            variant_position=variant_position,
            primer_regions=primer_regions
        )

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of a sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        return ''.join(complement.get(b, b) for b in reversed(sequence))


# ======================================================================
#  HomologyVisualizer — BLAST-like alignment display
# ======================================================================

class HomologyVisualizer:
    """
    Renders homology analysis results in a BLAST-like format inside a
    tk.Text widget with coloured tags for matches, mismatches and gaps.
    """

    def __init__(self, text_widget: tk.Text):
        self.text = text_widget
        self._result = None          # cached HomologyResult
        self._variant_query_pos = -1  # variant position in query (0-based)

        # --- Tag configuration ---
        self.text.tag_configure(
            "hom_header", foreground="#00BFFF",
            font=("Helvetica", 10, "bold"))
        self.text.tag_configure(
            "hom_primary", foreground="#FFD700",
            font=("Helvetica", 10, "bold"))
        self.text.tag_configure(
            "hom_warning", foreground="#FF6600",
            font=("Helvetica", 10, "bold"))
        self.text.tag_configure(
            "hom_stats", foreground="#DDDDDD",
            font=("Helvetica", 10))
        self.text.tag_configure(
            "hom_match", foreground="#00FF00",
            font=("Courier", 10))
        self.text.tag_configure(
            "hom_mismatch", foreground="#FF4444", background="#330000",
            font=("Courier", 10, "bold"))
        self.text.tag_configure(
            "hom_gap", foreground="#888888",
            font=("Courier", 10))
        self.text.tag_configure(
            "hom_label", foreground="#AAAAAA",
            font=("Courier", 9))
        self.text.tag_configure(
            "hom_separator", foreground="#555555",
            font=("Courier", 10))
        # Target variant position marker
        self.text.tag_configure(
            "hom_target", foreground="#FF00FF", background="#3a003a",
            font=("Courier", 10, "bold"))
        self.text.tag_configure(
            "hom_target_label", foreground="#FF00FF",
            font=("Courier", 9, "bold"))

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def display_results(self, result, variant_query_pos: int = -1) -> None:
        """
        Display the full homology analysis result (summary view).

        Shows a short summary; individual alignments are displayed on
        demand via :meth:`display_hit`.

        Args:
            result: HomologyResult from HomologyAnalyzer
            variant_query_pos: 0-based position of the variant inside
                               the query sequence (-1 if unknown)
        """
        from ..primer.homology_analyzer import HomologyAnalyzer

        self._result = result
        self._variant_query_pos = variant_query_pos

        self.text.configure(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)

        if result.error:
            self.text.insert(tk.END, f"ERROR: {result.error}\n", "hom_warning")
            self.text.configure(state=tk.DISABLED)
            return

        # --- Summary header ---
        self.text.insert(
            tk.END,
            f"Homology Analysis: {result.query_name}\n", "hom_header")
        self.text.insert(
            tk.END,
            f"Query length: {result.query_length} bp\n", "hom_stats")
        self.text.insert(
            tk.END,
            f"Source: {HomologyAnalyzer.format_chr(result.query_chromosome)}:"
            f"{result.query_position:,}\n",
            "hom_stats")
        self.text.insert(
            tk.END,
            f"Hits found: {len(result.hits)}\n", "hom_stats")
        self.text.insert(tk.END, "=" * 72 + "\n\n", "hom_separator")

        if not result.hits:
            self.text.insert(
                tk.END,
                "No hits matching the filter criteria.\n",
                "hom_stats")
        else:
            self.text.insert(
                tk.END,
                "Select a hit in the table above to see its alignment.\n",
                "hom_stats")

        self.text.configure(state=tk.DISABLED)

    def display_hit(self, hit_index: int) -> None:
        """
        Display the alignment details for a single hit selected in the
        hits summary table.

        Args:
            hit_index: 0-based index into ``self._result.hits``
        """
        if self._result is None or hit_index >= len(self._result.hits):
            return

        hit = self._result.hits[hit_index]

        self.text.configure(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)

        self._display_single_hit(
            hit, hit_index + 1,
            self._result.query_chromosome,
            self._result.query_position)

        self.text.configure(state=tk.DISABLED)

    # ------------------------------------------------------------------
    #  Single hit
    # ------------------------------------------------------------------

    def _display_single_hit(
        self, hit, index: int, source_chr: str, source_pos: int
    ) -> None:
        """Display one hit with header + BLAST-like alignment."""
        from ..primer.homology_analyzer import HomologyAnalyzer

        # --- Hit header ---
        if hit.is_primary:
            label = "[PRIMARY]"
            tag = "hom_primary"
        else:
            label = f"[Hit #{index}]"
            tag = "hom_header"

        self.text.insert(tk.END, f"{label} ", tag)

        end_pos = hit.position + hit.aligned_length
        chr_display = HomologyAnalyzer.format_chr(hit.chromosome)
        info = (
            f"{chr_display}:{hit.position:,}-{end_pos:,} "
            f"({hit.strand}) | Identity: {hit.percent_identity}% | "
            f"Mismatches: {hit.mismatches} | "
            f"Score: {hit.alignment_score}\n"
        )
        self.text.insert(tk.END, info, "hom_stats")

        # --- Pseudogene / duplication warning ---
        src_norm = HomologyAnalyzer.normalize_chr(source_chr)
        hit_norm = HomologyAnalyzer.normalize_chr(hit.chromosome)
        is_different_pos = abs(hit.position - source_pos) > 1000

        if hit.is_primary:
            self.text.insert(
                tk.END,
                "  [TARGET] This is the expected alignment location.\n",
                "hom_stats")
        elif hit_norm == src_norm and is_different_pos:
            self.text.insert(
                tk.END,
                "  !! Potential pseudogene / duplication "
                "on SAME chromosome !!\n",
                "hom_warning")
        elif hit_norm != src_norm:
            self.text.insert(
                tk.END,
                f"  !! Homologous region on DIFFERENT chromosome "
                f"({chr_display}) !!\n",
                "hom_warning")

        # --- Variant position info ---
        vqp = self._variant_query_pos
        if vqp >= 0:
            if hit.alignment_start <= vqp < hit.alignment_end:
                self.text.insert(
                    tk.END,
                    f"  ▶ Target variant (query pos {vqp}) is INSIDE "
                    f"this alignment — marked in magenta below\n",
                    "hom_target_label")
            else:
                self.text.insert(
                    tk.END,
                    f"  ○ Target variant (query pos {vqp}) is OUTSIDE "
                    f"this alignment (aligned range: "
                    f"{hit.alignment_start}-{hit.alignment_end})\n",
                    "hom_stats")

        self.text.insert(tk.END, "\n")

        # --- BLAST-like alignment ---
        query_line, mid_line, subject_line = \
            HomologyAnalyzer.build_alignment_lines(hit)

        if not query_line:
            self.text.insert(
                tk.END,
                "  (alignment details unavailable — missing MD tag)\n",
                "hom_stats")
        else:
            # Compute variant offset inside the alignment display
            variant_display_pos = -1
            if vqp >= 0 and hit.alignment_start <= vqp < hit.alignment_end:
                variant_display_pos = vqp - hit.alignment_start

            self._render_alignment_block(
                query_line, mid_line, subject_line,
                hit.alignment_start, hit.position,
                line_length=60,
                variant_display_pos=variant_display_pos)

        self.text.insert(tk.END, "-" * 72 + "\n\n", "hom_separator")

    # ------------------------------------------------------------------
    #  Alignment block renderer
    # ------------------------------------------------------------------

    def _render_alignment_block(
        self,
        query_line: str, mid_line: str, subject_line: str,
        query_start: int, subject_start: int,
        line_length: int = 60,
        variant_display_pos: int = -1,
    ) -> None:
        """
        Render the triple (Query / Midline / Subject) in coloured blocks
        of *line_length* characters, similar to BLAST output.

        Label width is computed dynamically so that Query, Midline and
        Subject lines always have identical padding regardless of how
        many digits the position numbers need.

        Args:
            variant_display_pos: 0-based position inside the alignment
                (query_line) where the target variant falls. -1 = none.
        """
        total = len(query_line)

        # Pre-calculate the maximum position value that will appear in
        # any label so that we can set a fixed width for position numbers.
        max_q_pos = query_start + total + 1
        max_s_pos = subject_start + self._count_ref_bases(subject_line)
        pos_width = max(7, len(str(max(max_q_pos, max_s_pos))))

        # Label format:  "  Sbjct {pos}  "  →  2 + 5 + 1 + pos_width + 2
        label_len = 2 + 5 + 1 + pos_width + 2
        pad = " " * label_len

        for block_start in range(0, total, line_length):
            block_end = min(block_start + line_length, total)
            chunk_q = query_line[block_start:block_end]
            chunk_m = mid_line[block_start:block_end]
            chunk_s = subject_line[block_start:block_end]

            # Compute position labels
            q_pos = query_start + block_start + 1  # 1-based
            s_pos = subject_start + self._count_ref_bases(
                subject_line[:block_start])

            # Query line
            self.text.insert(
                tk.END, f"  Query {q_pos:>{pos_width}}  ", "hom_label")
            for ci, (qc, mc) in enumerate(zip(chunk_q, chunk_m)):
                abs_pos = block_start + ci
                if abs_pos == variant_display_pos:
                    tag = "hom_target"
                else:
                    tag = self._char_tag(mc)
                self.text.insert(tk.END, qc, tag)
            self.text.insert(tk.END, "\n")

            # Midline — same padding as Query/Sbjct labels
            self.text.insert(tk.END, pad, "hom_label")
            for ci, mc in enumerate(chunk_m):
                abs_pos = block_start + ci
                if abs_pos == variant_display_pos:
                    self.text.insert(tk.END, "▼", "hom_target_label")
                else:
                    tag = self._char_tag(mc)
                    self.text.insert(tk.END, mc, tag)
            self.text.insert(tk.END, "\n")

            # Subject line
            self.text.insert(
                tk.END, f"  Sbjct {s_pos:>{pos_width}}  ", "hom_label")
            for ci, (sc, mc) in enumerate(zip(chunk_s, chunk_m)):
                abs_pos = block_start + ci
                if abs_pos == variant_display_pos:
                    tag = "hom_target"
                else:
                    tag = self._char_tag(mc)
                self.text.insert(tk.END, sc, tag)
            self.text.insert(tk.END, "\n\n")

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _char_tag(midline_char: str) -> str:
        """Map midline character to the appropriate text tag."""
        if midline_char == "|":
            return "hom_match"
        elif midline_char == ".":
            return "hom_mismatch"
        else:
            return "hom_gap"

    @staticmethod
    def _count_ref_bases(subject_slice: str) -> int:
        """Count non-gap characters in a subject line slice."""
        return sum(1 for c in subject_slice if c != "-")

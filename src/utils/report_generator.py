"""
Report generator for primer design results.
Generates HTML reports with visualization.
"""

import html as html_module
import os
from datetime import datetime
from typing import List, Optional, Dict, Any
from pathlib import Path

from ..core.models import (
    DesignResult, Variant, Amplicon, PrimerPair, Probe, DesignParameters
)
from ..core.grouper import VariantGroup
from .logger import LoggerMixin


class ReportGenerator(LoggerMixin):
    """
    Generates HTML reports for primer design results.
    """

    def __init__(self):
        """Initialize report generator."""
        self.template_dir = Path(__file__).parent.parent / "templates"

    def generate_report(
        self,
        design_results: List[DesignResult],
        parameters: DesignParameters,
        output_path: str,
        project_name: str = "Primer Design Report",
        methodology_context: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Generate an HTML report from design results.

        Args:
            design_results: List of design results
            parameters: Design parameters used
            output_path: Path for output HTML file
            project_name: Name for the report
            methodology_context: Optional dict with extra context for methodology section

        Returns:
            True if successful
        """
        try:
            html_content = self._build_html(
                design_results, parameters, project_name,
                methodology_context=methodology_context
            )

            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_content)

            self.log_info(f"Report generated: {output_path}")
            return True

        except Exception as e:
            self.log_exception(f"Report generation error: {e}")
            return False

    def _build_html(
        self,
        design_results: List[DesignResult],
        parameters: DesignParameters,
        project_name: str,
        methodology_context: Optional[Dict[str, Any]] = None
    ) -> str:
        """Build the complete HTML report."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        # Count statistics
        total_variants = sum(len(r.variants) for r in design_results)
        successful = sum(1 for r in design_results if r.success)

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{html_module.escape(project_name)}</title>
    <style>
        {self._get_css()}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>{html_module.escape(project_name)}</h1>
            <p class="timestamp">Generated: {timestamp}</p>
        </header>

        <section class="summary">
            <h2>Summary</h2>
            <div class="stats-grid">
                <div class="stat-box">
                    <span class="stat-value">{total_variants}</span>
                    <span class="stat-label">Total Variants</span>
                </div>
                <div class="stat-box">
                    <span class="stat-value">{len(design_results)}</span>
                    <span class="stat-label">Design Groups</span>
                </div>
                <div class="stat-box success">
                    <span class="stat-value">{successful}</span>
                    <span class="stat-label">Successful Designs</span>
                </div>
                <div class="stat-box {'warning' if successful < len(design_results) else 'success'}">
                    <span class="stat-value">{successful}/{len(design_results)}</span>
                    <span class="stat-label">Success Rate</span>
                </div>
            </div>
        </section>

        <section class="parameters">
            <h2>Design Parameters</h2>
            {self._build_parameters_section(parameters)}
        </section>

        <section class="results">
            <h2>Design Results</h2>
            {self._build_results_section(design_results)}
        </section>

        <section class="methodology">
            <h2>Methodology</h2>
            {self._build_methodology_section(design_results, parameters, methodology_context)}
        </section>

        <footer>
            <p>Primer Designer [Early Beta] Report | Reference Assembly: {parameters.reference_assembly}</p>
        </footer>
    </div>

    <script>
        {self._get_javascript()}
    </script>
</body>
</html>'''

        return html

    def _get_css(self) -> str:
        """Return CSS styles for the report."""
        return '''
        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }

        header {
            background: linear-gradient(135deg, #2c3e50, #3498db);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }

        header h1 {
            margin-bottom: 10px;
        }

        .timestamp {
            opacity: 0.8;
        }

        section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }

        h2 {
            color: #2c3e50;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #3498db;
        }

        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 20px;
        }

        .stat-box {
            background: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }

        .stat-box.success {
            background: #d5f5e3;
        }

        .stat-box.warning {
            background: #fdebd0;
        }

        .stat-value {
            display: block;
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
        }

        .stat-label {
            font-size: 0.9em;
            color: #7f8c8d;
        }

        .param-table {
            width: 100%;
            border-collapse: collapse;
        }

        .param-table th, .param-table td {
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }

        .param-table th {
            background: #f8f9fa;
            font-weight: 600;
        }

        .result-card {
            border: 1px solid #ddd;
            border-radius: 8px;
            margin-bottom: 20px;
            overflow: hidden;
        }

        .result-header {
            background: #f8f9fa;
            padding: 15px;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }

        .result-header:hover {
            background: #ecf0f1;
        }

        .result-header.success {
            border-left: 4px solid #27ae60;
        }

        .result-header.failed {
            border-left: 4px solid #e74c3c;
        }

        .result-body {
            padding: 20px;
            display: none;
        }

        .result-body.active {
            display: block;
        }

        .variant-info {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 15px;
        }

        .primer-sequence {
            font-family: 'Courier New', monospace;
            background: #2c3e50;
            color: #2ecc71;
            padding: 10px 15px;
            border-radius: 5px;
            margin: 5px 0;
            word-break: break-all;
        }

        .primer-pair {
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
        }

        .primer-metrics {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 10px;
            margin-top: 10px;
        }

        .metric {
            text-align: center;
            padding: 8px;
            background: #f8f9fa;
            border-radius: 5px;
        }

        .metric-value {
            font-weight: bold;
            color: #2c3e50;
        }

        .metric-label {
            font-size: 0.8em;
            color: #7f8c8d;
        }

        .warning-box {
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 10px 15px;
            margin: 10px 0;
            border-radius: 0 5px 5px 0;
        }

        .suggestion-box {
            background: #d1ecf1;
            border-left: 4px solid #17a2b8;
            padding: 10px 15px;
            margin: 10px 0;
            border-radius: 0 5px 5px 0;
        }

        .expand-icon {
            transition: transform 0.3s;
        }

        .expanded .expand-icon {
            transform: rotate(180deg);
        }

        footer {
            text-align: center;
            padding: 20px;
            color: #7f8c8d;
        }

        .specificity-box {
            margin-top: 15px;
            padding: 12px 15px;
            border-radius: 5px;
            border-left: 4px solid #bdc3c7;
        }

        .specificity-box.specific {
            background: #d5f5e3;
            border-left-color: #27ae60;
        }

        .specificity-box.likely-specific {
            background: #d5f5e3;
            border-left-color: #2ecc71;
        }

        .specificity-box.off-target {
            background: #fadbd8;
            border-left-color: #e74c3c;
        }

        .specificity-box.pseudogene {
            background: #f5b7b1;
            border-left-color: #922b21;
        }

        .specificity-box.inconclusive {
            background: #fdebd0;
            border-left-color: #f39c12;
        }

        .specificity-box.not-checked {
            background: #f8f9fa;
            border-left-color: #bdc3c7;
        }

        .verdict-label {
            font-weight: bold;
            font-size: 1.05em;
        }

        .tool-summary {
            font-size: 0.9em;
            color: #555;
            margin-top: 4px;
        }

        .off-target-detail {
            font-size: 0.85em;
            color: #922b21;
            margin-top: 3px;
            padding-left: 15px;
        }

        .pseudogene-warning {
            font-weight: bold;
            color: #922b21;
            margin-top: 6px;
        }

        .pair-rank {
            display: inline-block;
            background: #3498db;
            color: white;
            border-radius: 50%;
            width: 24px;
            height: 24px;
            text-align: center;
            line-height: 24px;
            font-size: 0.85em;
            font-weight: bold;
            margin-right: 6px;
        }

        .pair-rank.best {
            background: #27ae60;
        }

        .homology-box {
            margin-top: 10px;
            padding: 10px 15px;
            border-radius: 5px;
            border-left: 4px solid #bdc3c7;
        }

        .homology-box.good {
            background: #d5f5e3;
            border-left-color: #27ae60;
        }

        .homology-box.moderate {
            background: #fdebd0;
            border-left-color: #f39c12;
        }

        .homology-box.weak {
            background: #fadbd8;
            border-left-color: #e74c3c;
        }

        .homology-box.none {
            background: #f5b7b1;
            border-left-color: #922b21;
        }

        .homology-label {
            font-weight: bold;
        }

        .homology-detail {
            font-size: 0.9em;
            color: #555;
            margin-top: 3px;
        }

        .tier-box {
            margin: 10px 0;
            padding: 12px 15px;
            border-radius: 5px;
            border-left: 4px solid #bdc3c7;
        }

        .tier-box.tier-1 {
            background: #d5f5e3;
            border-left-color: #27ae60;
        }

        .tier-box.tier-2 {
            background: #d6eaf8;
            border-left-color: #2980b9;
        }

        .tier-box.tier-3 {
            background: #fdebd0;
            border-left-color: #f39c12;
        }

        .tier-box.tier-0 {
            background: #f8f9fa;
            border-left-color: #bdc3c7;
        }

        .tier-label {
            font-weight: bold;
            font-size: 1.0em;
        }

        .tier-detail {
            font-size: 0.9em;
            color: #555;
            margin-top: 4px;
        }

        .homology-regions-table tbody tr:nth-child(even) {
            background: #fafafa;
        }

        .homology-regions-table tbody td {
            padding: 3px 8px;
            border-bottom: 1px solid #eee;
        }

        .methodology-text {
            line-height: 1.8;
            text-align: justify;
        }

        .methodology-text h3 {
            color: #2c3e50;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 1.1em;
        }

        .methodology-text p {
            margin-bottom: 12px;
        }

        .methodology-text ul {
            margin: 10px 0;
            padding-left: 25px;
        }

        .methodology-text li {
            margin-bottom: 4px;
        }

        .methodology-note {
            background: #fef9e7;
            border-left: 4px solid #f1c40f;
            padding: 12px 15px;
            margin-top: 20px;
            border-radius: 0 5px 5px 0;
            font-style: italic;
            color: #7d6608;
        }

        .variant-detail-table {
            width: 100%;
            border-collapse: collapse;
            margin: 10px 0;
            font-size: 0.9em;
        }

        .variant-detail-table th, .variant-detail-table td {
            padding: 6px 10px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }

        .variant-detail-table th {
            background: #f0f3f4;
            font-weight: 600;
        }

        @media print {
            .result-body {
                display: block !important;
            }
            .expand-icon {
                display: none;
            }
        }
        '''

    def _get_javascript(self) -> str:
        """Return JavaScript for interactive elements."""
        return '''
        document.querySelectorAll('.result-header').forEach(header => {
            header.addEventListener('click', () => {
                const body = header.nextElementSibling;
                body.classList.toggle('active');
                header.classList.toggle('expanded');
            });
        });
        '''

    def _build_parameters_section(self, params: DesignParameters) -> str:
        """Build HTML for parameters section."""
        return f'''
        <table class="param-table">
            <tr><th colspan="2">Amplicon Settings</th></tr>
            <tr><td>Amplicon Size Range</td><td>{params.min_amplicon_size} - {params.max_amplicon_size} bp</td></tr>
            <tr><td>qPCR Amplicon Size</td><td>{params.qpcr_min_amplicon_size} - {params.qpcr_max_amplicon_size} bp</td></tr>
            <tr><th colspan="2">Primer Settings</th></tr>
            <tr><td>Primer Size</td><td>{params.primer_min_size} - {params.primer_max_size} nt</td></tr>
            <tr><td>Primer Tm</td><td>{params.primer_min_tm}°C - {params.primer_max_tm}°C</td></tr>
            <tr><td>Primer GC%</td><td>{params.primer_min_gc}% - {params.primer_max_gc}%</td></tr>
            <tr><th colspan="2">Chemistry</th></tr>
            <tr><td>Mg2+ Concentration</td><td>{params.dv_conc} mM</td></tr>
            <tr><td>dNTP Concentration</td><td>{params.dntp_conc} mM</td></tr>
            <tr><th colspan="2">Constraints</th></tr>
            <tr><td>Min Distance from Variant</td><td>{params.min_distance_from_variant} bp</td></tr>
            <tr><td>Min Distance from Exon Junction</td><td>{params.min_distance_from_exon_junction} bp</td></tr>
            <tr><th colspan="2">Homology Discrimination</th></tr>
            <tr><td>Pseudogene Discrimination</td><td>{"Enabled" if params.use_homology_discrimination else "Disabled"}</td></tr>
        </table>
        '''

    def _build_results_section(self, results: List[DesignResult]) -> str:
        """Build HTML for results section."""
        html_parts = []

        for i, result in enumerate(results, 1):
            status_class = 'success' if result.success else 'failed'
            status_text = 'Success' if result.success else 'Failed'

            # Variant info
            variants_html = ''
            for v in result.variants:
                variants_html += f'''
                <div class="variant-info">
                    <strong>Row {v.row_number}:</strong> {v.gene_symbol} - {v.transcript_accession}<br>
                    <strong>Variant:</strong> {v.hgvs_c}
                    {f'<br><strong>Genomic:</strong> {v.genomic_position}' if v.genomic_position else ''}
                </div>
                '''

            # Sort amplicons: best specificity first
            sorted_amps = self._sort_amplicons_by_specificity(result.amplicons)
            amplicons_html = ''
            for j, amp in enumerate(sorted_amps, 1):
                amplicons_html += self._build_amplicon_html(amp, j, len(sorted_amps))

            # Homology tier info
            tier_html = self._build_tier_html(result)

            # Warnings
            warnings_html = ''
            for warning in result.warnings:
                warnings_html += f'<div class="warning-box">{html_module.escape(warning)}</div>'

            # Suggestions
            suggestions_html = ''
            for suggestion in result.suggestions:
                suggestions_html += f'<div class="suggestion-box">{html_module.escape(suggestion)}</div>'

            html_parts.append(f'''
            <div class="result-card">
                <div class="result-header {status_class}">
                    <div>
                        <strong>Group {i}</strong> - {len(result.variants)} variant(s) - {status_text}
                        <span style="color: #7f8c8d; margin-left: 10px;">{result.message}</span>
                    </div>
                    <span class="expand-icon">▼</span>
                </div>
                <div class="result-body">
                    <h4>Variants</h4>
                    {variants_html}

                    {tier_html}

                    {'<h4>Primer Pairs</h4>' + amplicons_html if amplicons_html else ''}
                    {warnings_html}
                    {suggestions_html}
                </div>
            </div>
            ''')

        return ''.join(html_parts)

    @staticmethod
    def _sort_amplicons_by_specificity(amplicons: list) -> list:
        """Sort amplicons: SPECIFIC first, then by off-target count, PSEUDOGENE last."""
        def _key(amp):
            pp = amp.primer_pair
            v = pp.specificity_verdict
            if v == "SPECIFIC":
                return (0, 0)
            elif v == "LIKELY SPECIFIC":
                return (1, 0)
            elif not v:
                return (2, 0)
            elif v == "INCONCLUSIVE":
                return (3, 0)
            elif v == "OFF-TARGET DETECTED":
                return (4, len(pp.off_target_details))
            elif v == "PSEUDOGENE RISK":
                return (1000, len(pp.off_target_details))
            return (5, 0)
        return sorted(amplicons, key=_key)

    def _build_specificity_html(self, pp: PrimerPair) -> str:
        """Build HTML for specificity verdict of a primer pair."""
        verdict = pp.specificity_verdict
        if not verdict:
            return '''
            <div class="specificity-box not-checked">
                <span class="verdict-label">Specificity: not checked</span>
                <div class="tool-summary">Enable specificity tools before designing to check for off-targets.</div>
            </div>'''

        # Map verdict to CSS class and display
        verdict_map = {
            "SPECIFIC": ("specific", "✓ SPECIFIC", "Amplifies ONLY the intended target"),
            "LIKELY SPECIFIC": ("likely-specific", "✓ LIKELY SPECIFIC", "Checked tools say specific; not all tools were available"),
            "OFF-TARGET DETECTED": ("off-target", "✗ OFF-TARGET DETECTED", "Primers may amplify additional genomic sites"),
            "PSEUDOGENE RISK": ("pseudogene", "⚠ PSEUDOGENE RISK", "Off-target on same chromosome — likely pseudogene or segmental duplication"),
            "INCONCLUSIVE": ("inconclusive", "? INCONCLUSIVE", "Could not verify specificity"),
        }
        css_class, label, description = verdict_map.get(verdict, ("inconclusive", verdict, ""))

        # Tool summary
        tool_line = ''
        if pp.specificity_tool_summary:
            tool_line = f'<div class="tool-summary">Tools: {html_module.escape(pp.specificity_tool_summary)}</div>'

        # Off-target details (show up to 5)
        off_target_lines = ''
        if pp.off_target_details:
            for detail in pp.off_target_details[:5]:
                off_target_lines += f'<div class="off-target-detail">• {html_module.escape(detail)}</div>'
            if len(pp.off_target_details) > 5:
                off_target_lines += f'<div class="off-target-detail">... and {len(pp.off_target_details) - 5} more off-target(s)</div>'

        # Pseudogene bold warning
        pseudo_line = ''
        if pp.pseudogene_risk:
            pseudo_line = '<div class="pseudogene-warning">⚠ Same-chromosome off-target — possible pseudogene. Consider redesigning this primer pair.</div>'

        return f'''
            <div class="specificity-box {css_class}">
                <span class="verdict-label">{label}</span>
                — {html_module.escape(description)}
                {tool_line}
                {off_target_lines}
                {pseudo_line}
            </div>'''

    def _build_tier_html(self, result: DesignResult) -> str:
        """Build HTML for homology tier info of a design result."""
        tier = result.homology_tier
        msg = result.homology_tier_message

        if not result.homology_discriminated and tier == 0:
            return ''

        tier_labels = {
            0: ("tier-0", "No homology constraints", "No homologous regions detected; standard primer design used."),
            1: ("tier-1", "Tier 1 — Junction Overlap", "Primers directly overlap discriminating positions (best discrimination)."),
            2: ("tier-2", "Tier 2 — Discriminating Windows", "Primers placed within discriminating regions (good discrimination)."),
            3: ("tier-3", "Tier 3 — Fallback Re-ranking", "Primers re-ranked by homology score; full discrimination not guaranteed."),
        }

        css_class, label, description = tier_labels.get(tier, ("tier-0", f"Tier {tier}", ""))

        detail_line = ''
        if msg:
            detail_line = f'<div class="tier-detail">{html_module.escape(msg)}</div>'

        # Build detailed homologous regions table
        homology_details = self._build_homologous_regions_html(result)

        return f'''
            <div class="tier-box {css_class}">
                <span class="tier-label">{label}</span>
                — {html_module.escape(description)}
                {detail_line}
                {homology_details}
            </div>'''

    def _build_homologous_regions_html(self, result: DesignResult) -> str:
        """Build HTML table showing details of homologous regions (pseudogenes/duplications).

        Extracts secondary hits from the HomologyResult to show the user exactly
        which genomic locations are homologous to the target and could potentially
        be co-amplified.
        """
        homology_result = result.homology_result
        if homology_result is None:
            return ''

        secondary_hits = [
            h for h in homology_result.hits
            if not h.is_primary and not h.is_supplementary
        ]
        if not secondary_hits:
            return ''

        # Get primary hit info for context
        primary = homology_result.primary_hit
        primary_chr = primary.chromosome if primary else homology_result.query_chromosome
        primary_pos = primary.position if primary else homology_result.query_position

        rows = ''
        for hit in secondary_hits:
            same_chr = ''
            norm_hit_chr = hit.chromosome.replace('chr', '')
            norm_primary_chr = str(primary_chr).replace('chr', '')
            if norm_hit_chr == norm_primary_chr:
                same_chr = ' <span style="color: #c0392b; font-weight: bold;">⚠ same chr</span>'

            rows += f'''
                <tr>
                    <td>{html_module.escape(hit.chromosome)}:{hit.position:,}{same_chr}</td>
                    <td>{hit.strand}</td>
                    <td>{hit.percent_identity:.1f}%</td>
                    <td>{hit.aligned_length} bp</td>
                    <td>{hit.mismatches}</td>
                    <td>{hit.bit_score:.0f}</td>
                </tr>'''

        return f'''
            <div style="margin-top: 10px;">
                <strong>Homologous regions detected ({len(secondary_hits)}):</strong>
                <div style="font-size: 0.85em; color: #666; margin-bottom: 5px;">
                    Target: {html_module.escape(str(primary_chr))}:{primary_pos:,}
                </div>
                <table class="homology-regions-table" style="width: 100%; border-collapse: collapse; font-size: 0.85em; margin-top: 5px;">
                    <thead>
                        <tr style="background: #f0f0f0; text-align: left;">
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Location</th>
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Strand</th>
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Identity</th>
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Aligned Length</th>
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Mismatches</th>
                            <th style="padding: 4px 8px; border-bottom: 2px solid #ccc;">Bit Score</th>
                        </tr>
                    </thead>
                    <tbody>
                        {rows}
                    </tbody>
                </table>
                <div style="font-size: 0.8em; color: #888; margin-top: 4px;">
                    Primers that do not discriminate from these regions may co-amplify products
                    of similar size from these locations.
                </div>
            </div>'''

    def _build_homology_html(self, pp: PrimerPair) -> str:
        """Build HTML for homology discrimination info of a primer pair."""
        score = pp.homology_discrimination_score
        fwd_n = pp.fwd_discriminating_positions
        rev_n = pp.rev_discriminating_positions
        warning = pp.homology_warning

        # No homology data at all — don't show section
        if score == 0 and not warning:
            return ''

        if score >= 6:
            css_class = "good"
            label = "GOOD"
        elif score >= 2:
            css_class = "moderate"
            label = "MODERATE"
        elif score > 0:
            css_class = "weak"
            label = "WEAK"
        else:
            css_class = "none"
            label = "NONE"

        detail = f"Score: {score:.1f} | Forward: {fwd_n} discriminating pos. | Reverse: {rev_n} discriminating pos."

        warning_line = ''
        if warning:
            warning_line = f'<div class="homology-detail" style="color: #922b21;">{html_module.escape(warning)}</div>'

        return f'''
            <div class="homology-box {css_class}">
                <span class="homology-label">Homology Discrimination: {label}</span>
                <div class="homology-detail">{detail}</div>
                {warning_line}
            </div>'''

    def _build_methodology_section(
        self,
        results: List[DesignResult],
        params: DesignParameters,
        ctx: Optional[Dict[str, Any]] = None
    ) -> str:
        """Build HTML methodology section from design results and context."""
        if ctx is None:
            ctx = {}

        assembly = params.reference_assembly
        assembly_full = "GRCh38/hg38" if assembly == "GRCh38" else "GRCh37/hg19"

        # Collect variants from results
        all_variants = []
        for r in results:
            all_variants.extend(r.variants)

        # Unique gene-transcript pairs
        gene_transcript_pairs = []
        seen_gt = set()
        for v in all_variants:
            key = (v.gene_symbol, v.transcript_accession)
            if key not in seen_gt:
                seen_gt.add(key)
                gene_transcript_pairs.append(key)

        # Context flags (from GUI or defaults)
        used_ucsc = ctx.get('used_ucsc', False)
        used_primer_blast = ctx.get('used_primer_blast', False)
        used_homology = params.use_homology_discrimination
        continuous = ctx.get('continuous_mode', False)
        max_attempts = ctx.get('max_attempts', 20)
        filter_pop = params.filter_population_variants
        pop_code = ctx.get('population', 'global')
        flanking = ctx.get('flanking_region', '250')
        blast_min_identity = ctx.get('blast_min_identity', 85.0)
        blast_min_length = ctx.get('blast_min_aligned_length', 50)

        pop_labels = {
            'global': 'global (worldwide)',
            'nfe': 'Non-Finnish European (NFE)',
            'afr': 'African/African American',
            'eas': 'East Asian',
            'sas': 'South Asian',
            'amr': 'Latino/Admixed American',
            'asj': 'Ashkenazi Jewish',
            'fin': 'Finnish',
        }
        pop_name = pop_labels.get(pop_code, pop_code)

        # Per-variant MAF thresholds (percentage)
        variant_maf_map = ctx.get('variant_maf_thresholds', {})
        default_maf_pct = ctx.get('default_maf_threshold_pct', 0.5)
        per_variant_maf = {}
        for v in all_variants:
            per_variant_maf[v.row_number] = variant_maf_map.get(
                v.row_number, default_maf_pct
            )
        unique_mafs = sorted(set(per_variant_maf.values())) if per_variant_maf else [default_maf_pct]

        # Count outcomes
        total_designed = len(results)
        successful = sum(1 for r in results if r.success)
        failed = total_designed - successful
        total_pairs = sum(len(r.amplicons) for r in results if r.success)
        specific_pairs = sum(
            sum(1 for a in r.amplicons if a.primer_pair.specificity_verdict == "SPECIFIC")
            for r in results if r.success
        )
        homology_applied = sum(1 for r in results if r.homology_discriminated)

        parts = []
        parts.append('<div class="methodology-text">')

        # -- Primer Design --
        parts.append('<h3>Primer Design</h3>')

        # Target description
        if gene_transcript_pairs:
            gene_list = ', '.join(
                f'{html_module.escape(g)} ({html_module.escape(t)})'
                for g, t in gene_transcript_pairs
            )
            parts.append(
                f'<p>PCR primers were designed to amplify regions containing '
                f'variant(s) in the following gene(s) and transcript(s): '
                f'{gene_list}. '
            )
        else:
            parts.append(
                '<p>PCR primers were designed to amplify regions containing '
                'the target variant(s). '
            )

        # Variant enumeration
        if all_variants:
            var_descs = '; '.join(
                f'{html_module.escape(v.gene_symbol)} '
                f'{html_module.escape(v.transcript_accession)}:'
                f'{html_module.escape(v.hgvs_c)}'
                for v in all_variants
            )
            parts.append(
                f'A total of {len(all_variants)} variant(s) were analysed: '
                f'{var_descs}.</p>'
            )
        else:
            parts.append('</p>')

        # Reference genome and coordinate mapping
        parts.append(
            f'<p>All genomic coordinates were referenced to the {assembly_full} '
            f'human genome assembly. Coding DNA sequence (CDS) positions '
            f'were mapped to genomic coordinates using the Ensembl REST API '
            f'(https://rest.ensembl.org) and cross-validated against NCBI '
            f'Entrez (https://eutils.ncbi.nlm.nih.gov/entrez/). Transcript '
            f'annotations were verified against the MANE (Matched Annotation '
            f'from NCBI and EBI) database to ensure use of clinically '
            f'recommended reference sequences.</p>'
        )

        # Sequence retrieval
        parts.append(
            f'<p>Genomic sequences were retrieved from the Ensembl REST API, '
            f'extending {html_module.escape(str(flanking))} bp on each side of each variant '
            f'position to provide a search space for primer placement.</p>'
        )

        # -- Population Variant Masking --
        parts.append('<h3>Population Variant Masking</h3>')
        if filter_pop:
            if len(unique_mafs) == 1:
                maf_desc = f'&ge; {unique_mafs[0]:.3f}%'
            else:
                maf_desc = (
                    f'&ge; {min(unique_mafs):.3f}%&ndash;{max(unique_mafs):.3f}% '
                    f'(per-variant thresholds applied; see Per-Variant Details below)'
                )
            parts.append(
                f'<p>Known population polymorphisms were obtained from the '
                f'Ensembl Variation API (sourced from gnomAD v4). Positions '
                f'with a minor allele frequency (MAF) {maf_desc} in '
                f'the {html_module.escape(pop_name)} population were masked with IUPAC '
                f'ambiguity code "N" in the template sequence provided to the '
                f'primer design engine. This masking ensures that primers do '
                f'not bind at polymorphic positions, which could otherwise '
                f'lead to allele-specific amplification bias or failure.</p>'
            )
        else:
            parts.append(
                '<p>Population variant masking was disabled for this design '
                'session. Primers may bind at polymorphic positions.</p>'
            )

        # -- Primer Design Engine --
        parts.append('<h3>Primer Design Engine</h3>')
        parts.append(
            f'<p>Primers were designed using the Primer3 engine (via the '
            f'primer3-py Python interface) with the following key parameters: '
            f'amplicon size {params.min_amplicon_size}&ndash;{params.max_amplicon_size}&nbsp;bp, '
            f'primer melting temperature (Tm) '
            f'{params.primer_min_tm}&ndash;{params.primer_max_tm}&nbsp;&deg;C '
            f'(optimal {params.primer_opt_tm}&nbsp;&deg;C), '
            f'primer length {params.primer_min_size}&ndash;{params.primer_max_size}&nbsp;nt '
            f'(optimal {params.primer_opt_size}&nbsp;nt), '
            f'GC content {params.primer_min_gc}&ndash;{params.primer_max_gc}% '
            f'(optimal {params.primer_opt_gc}%). '
            f'Thermodynamic calculations used the SantaLucia nearest-neighbour '
            f'model with salt corrections for '
            f'{params.mv_conc}&nbsp;mM monovalent cations, '
            f'{params.dv_conc}&nbsp;mM Mg<sup>2+</sup>, '
            f'{params.dntp_conc}&nbsp;mM dNTPs, and '
            f'{params.dna_conc}&nbsp;nM template DNA. '
            f'Maximum self-complementarity was set to {params.max_self_complementarity}, '
            f"3' end complementarity to {params.max_end_complementarity}, "
            f'and pair complementarity to {params.max_pair_complementarity} '
            f'(Primer3 thermodynamic alignment scores).</p>'
        )

        # Distance constraints
        constraint_parts = []
        if params.use_variant_distance_constraint:
            constraint_parts.append(
                f'a minimum distance of {params.min_distance_from_variant}&nbsp;bp '
                f'between the nearest primer edge and the variant position'
            )
        if params.use_splice_site_constraint:
            constraint_parts.append(
                f'a minimum distance of {params.min_distance_from_exon_junction}&nbsp;bp '
                f'from exon/intron boundaries on the intronic side to avoid '
                f'splice-site interference'
            )
        if constraint_parts:
            parts.append(
                f'<p>Primer placement was constrained to enforce '
                f'{"; and ".join(constraint_parts)}.</p>'
            )

        # -- Pseudogene Discrimination --
        parts.append('<h3>Pseudogene Discrimination</h3>')
        if used_homology:
            parts.append(
                f'<p>To minimise the risk of pseudogene co-amplification, a '
                f'BLAST+-based homology analysis was performed for each target '
                f'region. Flanking sequences were queried against the human '
                f'genome using BLASTn (word size 11, evalue &le; 1e-10, minimum '
                f'identity &ge; {blast_min_identity}%, minimum aligned length '
                f'&ge; {blast_min_length}&nbsp;nt) to identify homologous regions '
                f'(pseudogenes, segmental duplications). '
            )
            if homology_applied > 0:
                parts.append(
                    f'Homologous regions were detected for {homology_applied} of '
                    f'{total_designed} variant group(s). Primer pairs were '
                    f're-ranked to maximise the number of discriminating '
                    f"positions &mdash; nucleotides where the primer sequence "
                    f"differs from the pseudogene &mdash; preferring mismatches "
                    f"at the 3' end of the primer for maximum allele-specific "
                    f'discrimination.</p>'
                )
            else:
                parts.append(
                    f'No significant homologous regions were detected for any of '
                    f'the target loci, so standard Primer3 ranking was used.</p>'
                )
        else:
            parts.append(
                '<p>Pseudogene discrimination (homology analysis) was not '
                'enabled for this design session.</p>'
            )

        # -- Specificity Verification --
        parts.append('<h3>Specificity Verification</h3>')
        spec_tools = []
        if used_ucsc:
            spec_tools.append(
                f'UCSC In-Silico PCR (https://genome.ucsc.edu/cgi-bin/hgPcr), '
                f'which simulates PCR amplification against the full genome '
                f'assembly ({assembly_full}) to detect off-target binding sites'
            )
        if used_primer_blast:
            spec_tools.append(
                'NCBI Primer-BLAST (https://www.ncbi.nlm.nih.gov/tools/primer-blast/), '
                'which checks primer specificity against the RefSeq representative '
                'genome database'
            )
        if spec_tools:
            parts.append(
                f'<p>Primer specificity was verified using: {"; and ".join(spec_tools)}. '
            )
            if continuous:
                parts.append(
                    f'Continuous design mode was used: for each variant, up to '
                    f'{max_attempts} candidate primer pairs were generated by '
                    f'Primer3 and evaluated sequentially for specificity until '
                    f'at least {params.min_specific_pairs} verified-specific '
                    f'pair(s) were found per variant, or all candidates were '
                    f'exhausted.</p>'
                )
            else:
                parts.append(
                    f'{params.num_primer_pairs} primer pair(s) per variant were '
                    f'designed and subsequently checked for specificity.</p>'
                )
        else:
            parts.append(
                '<p>No post-design specificity verification tools were enabled. '
                'Primer specificity was not validated in silico.</p>'
            )

        # -- Results Summary --
        parts.append('<h3>Results Summary</h3>')
        parts.append(
            f'<p>Primer design was completed for {successful} of '
            f'{total_designed} variant group(s). '
        )
        if total_pairs > 0:
            parts.append(
                f'A total of {total_pairs} primer pair(s) were generated. '
            )
        if spec_tools and total_pairs > 0:
            parts.append(
                f'Of these, {specific_pairs} pair(s) were confirmed as specific '
                f'(amplifying only the intended target). '
            )
        if failed > 0:
            parts.append(
                f'Primer design failed for {failed} variant group(s), likely due '
                f'to sequence constraints (high polymorphism density, extreme GC '
                f'content, or insufficient flanking region). '
            )
        parts.append('</p>')

        # -- Per-Variant Detail Table --
        if results:
            parts.append('<h3>Per-Variant Details</h3>')
            parts.append('<table class="variant-detail-table">')
            parts.append(
                '<tr><th>Gene</th><th>Transcript</th><th>Variant</th>'
                '<th>MAF&nbsp;Threshold</th><th>Pairs</th><th>Amplicon Size</th>'
            )
            if spec_tools:
                parts.append('<th>Specific</th>')
            parts.append('<th>Homology</th></tr>')

            for r in results:
                if not r.variants:
                    continue
                v = r.variants[0]
                gene = html_module.escape(v.gene_symbol)
                transcript = html_module.escape(v.transcript_accession)
                hgvs = html_module.escape(v.hgvs_c)
                v_maf = per_variant_maf.get(v.row_number, default_maf_pct)

                if r.success and r.amplicons:
                    n_pairs = len(r.amplicons)
                    n_spec = sum(
                        1 for a in r.amplicons
                        if a.primer_pair.specificity_verdict == "SPECIFIC"
                    )
                    amp_sizes = [a.primer_pair.product_size for a in r.amplicons]
                    size_range = (f"{min(amp_sizes)}&ndash;{max(amp_sizes)}&nbsp;bp"
                                  if len(amp_sizes) > 1 else f"{amp_sizes[0]}&nbsp;bp")
                    homo = ""
                    if r.homology_discriminated:
                        homo = f"Tier {r.homology_tier}"
                    else:
                        homo = "&mdash;"

                    parts.append(f'<tr><td>{gene}</td><td>{transcript}</td>'
                                 f'<td>{hgvs}</td>'
                                 f'<td>{v_maf:.3f}%</td>'
                                 f'<td>{n_pairs}</td>'
                                 f'<td>{size_range}</td>')
                    if spec_tools:
                        parts.append(f'<td>{n_spec}</td>')
                    parts.append(f'<td>{homo}</td></tr>')
                else:
                    colspan = "5" if spec_tools else "4"
                    parts.append(f'<tr><td>{gene}</td><td>{transcript}</td>'
                                 f'<td>{hgvs}</td>'
                                 f'<td>{v_maf:.3f}%</td>'
                                 f'<td colspan="{colspan}">'
                                 f'Design failed</td></tr>')

            parts.append('</table>')

        # -- Software and Databases --
        parts.append('<h3>Software and Databases</h3>')
        parts.append(
            '<p>Primer design was performed using Primer Designer [Early Beta], '
            'a Python application integrating the following components:</p>'
            '<ul>'
            '<li>Primer3 (via primer3-py) &mdash; thermodynamic primer design engine</li>'
            '<li>Ensembl REST API &mdash; genomic sequence retrieval, variant annotation, '
            'coordinate mapping</li>'
            '<li>NCBI Entrez &mdash; transcript validation and cross-referencing</li>'
            '<li>MANE database &mdash; Matched Annotation from NCBI and EBI for '
            'clinically recommended transcripts</li>'
        )
        if filter_pop:
            parts.append(
                '<li>gnomAD v4 (via Ensembl Variation API) &mdash; population variant '
                'frequencies for primer-site masking</li>'
            )
        if used_homology:
            parts.append(
                '<li>BLAST+ (BLASTn) &mdash; local homology analysis for pseudogene '
                'detection and primer discrimination scoring</li>'
            )
        if used_ucsc:
            parts.append(
                '<li>UCSC In-Silico PCR &mdash; genome-wide specificity simulation</li>'
            )
        if used_primer_blast:
            parts.append(
                '<li>NCBI Primer-BLAST &mdash; RefSeq-based specificity verification</li>'
            )
        parts.append('</ul>')

        # Note
        parts.append(
            '<div class="methodology-note">'
            'Note: This methodology section was auto-generated based on the '
            'parameters and tools used in this design session. Verify all '
            'details before including in a publication or report.'
            '</div>'
        )

        parts.append('</div>')
        return '\n'.join(parts)

    def _build_amplicon_html(self, amp: Amplicon, index: int, total: int = 0) -> str:
        """Build HTML for a single amplicon."""
        pp = amp.primer_pair

        # Rank badge: green for #1 (best)
        rank_class = "pair-rank best" if index == 1 else "pair-rank"

        probe_html = ''
        if amp.probe:
            probe_html = f'''
            <div style="margin-top: 15px;">
                <strong>Probe:</strong>
                <div class="primer-sequence">{amp.probe.sequence}</div>
                <div class="primer-metrics">
                    <div class="metric">
                        <div class="metric-value">{amp.probe.tm:.1f}°C</div>
                        <div class="metric-label">Tm</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{amp.probe.gc_content:.1f}%</div>
                        <div class="metric-label">GC%</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{amp.probe.length}</div>
                        <div class="metric-label">Length</div>
                    </div>
                </div>
            </div>
            '''

        # Specificity section
        specificity_html = self._build_specificity_html(pp)

        return f'''
        <div class="primer-pair">
            <span class="{rank_class}">{index}</span>
            <strong>Primer Pair {index}</strong>
            <p style="color: #7f8c8d;">Product Size: {pp.product_size} bp</p>

            <div style="margin-top: 10px;">
                <strong>Forward Primer:</strong>
                <div class="primer-sequence">5'-{pp.forward.sequence}-3'</div>
                <div class="primer-metrics">
                    <div class="metric">
                        <div class="metric-value">{pp.forward.tm:.1f}°C</div>
                        <div class="metric-label">Tm</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{pp.forward.gc_content:.1f}%</div>
                        <div class="metric-label">GC%</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{pp.forward.length}</div>
                        <div class="metric-label">Length</div>
                    </div>
                </div>
            </div>

            <div style="margin-top: 15px;">
                <strong>Reverse Primer:</strong>
                <div class="primer-sequence">5'-{pp.reverse.sequence}-3'</div>
                <div class="primer-metrics">
                    <div class="metric">
                        <div class="metric-value">{pp.reverse.tm:.1f}°C</div>
                        <div class="metric-label">Tm</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{pp.reverse.gc_content:.1f}%</div>
                        <div class="metric-label">GC%</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{pp.reverse.length}</div>
                        <div class="metric-label">Length</div>
                    </div>
                </div>
            </div>

            {probe_html}

            <div style="margin-top: 15px;">
                <strong>Pair Metrics:</strong>
                <div class="primer-metrics">
                    <div class="metric">
                        <div class="metric-value">{pp.tm_difference:.1f}°C</div>
                        <div class="metric-label">ΔTm</div>
                    </div>
                    <div class="metric">
                        <div class="metric-value">{pp.pair_complementarity:.1f}</div>
                        <div class="metric-label">Pair Compl.</div>
                    </div>
                </div>
            </div>

            {specificity_html}
            {self._build_homology_html(pp)}
        </div>
        '''

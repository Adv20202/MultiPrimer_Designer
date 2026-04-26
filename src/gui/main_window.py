"""
Main GUI window for the Primer Designer application.
Built with tkinter for cross-platform compatibility.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Callable
from pathlib import Path

from ..core.models import (
    Variant, DesignParameters, DesignMode, DesignResult, ProjectState
)
from ..core.validator import VariantValidator
from ..core.coordinator import CoordinateTranslator
from ..core.grouper import VariantGrouper, VariantGroup
from ..api.ncbi_client import NCBIClient
from ..api.ensembl_client import EnsemblClient
from ..api.mane_manager import MANEManager
from ..api.variant_db_client import VariantDBClient
from ..primer.designer import PrimerDesigner
from ..utils.file_parser import FileParser
from ..utils.config import Config, get_config
from ..utils.logger import setup_logger, LoggerMixin
from ..utils.report_generator import ReportGenerator

from .visualization import SequenceVisualizer
from . import theme


class MainWindow(LoggerMixin):
    """
    Main application window with full GUI for primer design.
    """

    def __init__(self):
        """Initialize the main window."""
        self.root = tk.Tk()
        self.root.title("Primer Designer - PCR Primer Design Tool [Early Beta]")
        # Full screen width, keep current height
        screen_width = self.root.winfo_screenwidth()
        self.root.geometry(f"{screen_width}x920+0+0")
        self.root.minsize(1200, 750)

        # Apply dark theme (must be before building UI)
        theme.apply_dark_theme(self.root)

        # Initialize logger
        setup_logger()

        # Application state
        self.project_state = ProjectState()
        self.config = get_config()
        self.design_parameters = DesignParameters()

        # Processing components
        self.file_parser = FileParser()
        self.ncbi_client = NCBIClient()
        self.ensembl_client = EnsemblClient()
        # Initialize MANEManager without auto-download to prevent blocking at startup
        self.mane_manager = MANEManager(auto_download=False)
        self.variant_db_client = VariantDBClient()
        self.variant_db_client.clear_cache()  # Ensure fresh MAF data on startup
        self.validator = VariantValidator(
            self.ncbi_client, self.ensembl_client, self.mane_manager
        )
        self.coordinator = CoordinateTranslator(self.ensembl_client, self.ncbi_client)
        self.grouper = VariantGrouper(self.design_parameters)
        self.primer_designer = None  # Initialized when needed
        self.specificity_checker = None
        self.specificity_orchestrator = None  # Lazy-init on first design
        self.report_generator = ReportGenerator()

        # Threading
        self._cancel_lock = threading.Lock()
        self._cancel_requested = False
        self._current_thread: Optional[threading.Thread] = None

        # Skip-variant mechanism for continuous mode
        self._skip_variant_lock = threading.Lock()
        self._skip_variant_requested = False

        # Currently selected variant for refresh on population change
        self._selected_variant: Optional[Variant] = None

        # Pre-cached sequence data for fast display
        # Dict: row_number -> {'sequence': str, 'pop_variants': list, 'seq_start': int, 'seq_end': int}
        self._sequence_cache = {}
        self._cache_thread: Optional[threading.Thread] = None

        # Track which assembly was used when variants were loaded
        self._loaded_with_assembly: Optional[str] = None

        # Homology analysis
        self.homology_analyzer = None   # Lazy init
        self._homology_result = None    # Last result for Treeview selection
        self._homology_results_by_variant = {}  # row_number -> HomologyResult

        # Build UI
        self._build_ui()
        self._bind_events()

        self.log_info("Application initialized")

    def _build_ui(self):
        """Build the main user interface."""
        # Main container with paned window (classic tk for minsize support)
        self.main_paned = tk.PanedWindow(self.root, orient=tk.HORIZONTAL,
                                          sashwidth=6, sashrelief=tk.RAISED,
                                          bg=theme.SASH)
        self.main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel - Controls (minsize prevents label truncation)
        self.left_frame = ttk.Frame(self.main_paned)
        self.main_paned.add(self.left_frame, minsize=530, width=530, stretch='never')

        # Center panel - Results and Visualization
        self.right_frame = ttk.Frame(self.main_paned)
        self.main_paned.add(self.right_frame, stretch='always')

        # Right panel - Contextual Help
        self.help_frame = ttk.Frame(self.main_paned)
        self.main_paned.add(self.help_frame, minsize=250, width=320, stretch='middle')

        # Build left panel sections
        self._build_file_section()
        self._build_assembly_section()
        self._build_parameters_section()
        self._build_population_section()
        self._build_homology_discrimination_section()
        self._build_specificity_section()
        self._build_action_buttons()

        # Build center panel sections
        self._build_variants_table()
        self._build_visualization_panel()
        self._build_results_panel()

        # Build help panel
        self._build_help_panel()

        # Status bar
        self._build_status_bar()

        # Show initial help context
        self._update_help_context('startup')

    def _build_file_section(self):
        """Build file input section."""
        frame = ttk.LabelFrame(self.left_frame, text="Input File", padding=(10, 5))
        frame.pack(fill=tk.X, padx=5, pady=3)

        # File path entry
        self.file_path_var = tk.StringVar()
        entry_frame = ttk.Frame(frame)
        entry_frame.pack(fill=tk.X)

        self.file_entry = ttk.Entry(entry_frame, textvariable=self.file_path_var)
        self.file_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self.browse_btn = ttk.Button(entry_frame, text="Browse", command=self._browse_file)
        self.browse_btn.pack(side=tk.RIGHT, padx=(5, 0))

        # Validation options
        self.skip_api_validation = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            frame, text="Skip online validation (faster, offline mode)",
            variable=self.skip_api_validation
        ).pack(anchor=tk.W, pady=(3, 0))

        # Flanking region size
        flank_frame = ttk.Frame(frame)
        flank_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Label(flank_frame, text="Flanking region:").pack(side=tk.LEFT)
        self.flank_size_var = tk.StringVar(value="250")
        flank_entry = ttk.Entry(flank_frame, textvariable=self.flank_size_var, width=6)
        flank_entry.pack(side=tk.LEFT, padx=5)
        ttk.Label(flank_frame, text="bp each side").pack(side=tk.LEFT)

        # Load button
        self.load_btn = ttk.Button(frame, text="Load and Validate", command=self._load_file)
        self.load_btn.pack(fill=tk.X, pady=(5, 0))

    def _build_assembly_section(self):
        """Build reference assembly selection."""
        frame = ttk.LabelFrame(self.left_frame, text="Reference Genome", padding=(10, 3))
        frame.pack(fill=tk.X, padx=5, pady=2)

        self.assembly_var = tk.StringVar(value="GRCh38")
        self.assembly_var.trace_add('write', self._on_assembly_changed)

        row = ttk.Frame(frame)
        row.pack(fill=tk.X)

        ttk.Radiobutton(
            row, text="GRCh38 (HG38) - Recommended",
            variable=self.assembly_var, value="GRCh38"
        ).pack(side=tk.LEFT, padx=(0, 10))

        ttk.Radiobutton(
            row, text="GRCh37 (HG19) - Legacy",
            variable=self.assembly_var, value="GRCh37"
        ).pack(side=tk.LEFT)

    def _on_assembly_changed(self, *_args):
        """Handle assembly change — invalidate loaded variants if assembly differs."""
        new_assembly = self.assembly_var.get()

        # If variants were loaded with a different assembly, they are stale
        if (self._loaded_with_assembly is not None
                and self._loaded_with_assembly != new_assembly
                and self.project_state.variants):
            # Clear stale data
            for variant in self.project_state.variants:
                variant.genomic_position = None
                variant.exon = None
            self._sequence_cache.clear()
            self.project_state.is_validated = False
            self.design_btn.configure(state=tk.DISABLED)
            self._loaded_with_assembly = None

            messagebox.showwarning(
                "Assembly Changed",
                f"Reference genome changed to {new_assembly}.\n\n"
                "Genomic coordinates are no longer valid.\n"
                "Please click 'Load and Validate' to re-map variants."
            )

    def _build_parameters_section(self):
        """Build primer design parameters section."""
        frame = ttk.LabelFrame(self.left_frame, text="Design Parameters", padding=(10, 5))
        frame.pack(fill=tk.X, padx=5, pady=2)

        # Design mode (qPCR hidden for now — kept in code but not shown in GUI)
        self.mode_var = tk.StringVar(value="PCR")
        # mode_frame = ttk.Frame(frame)
        # mode_frame.pack(fill=tk.X, pady=(0, 5))
        # ttk.Label(mode_frame, text="Design Mode:").pack(side=tk.LEFT)
        # ttk.Radiobutton(mode_frame, text="PCR", variable=self.mode_var,
        #                 value="PCR").pack(side=tk.LEFT, padx=5)
        # ttk.Radiobutton(mode_frame, text="qPCR", variable=self.mode_var,
        #                 value="qPCR").pack(side=tk.LEFT)

        # Amplicon size
        size_frame = ttk.Frame(frame)
        size_frame.pack(fill=tk.X, pady=2)

        ttk.Label(size_frame, text="Amplicon Size:").pack(side=tk.LEFT)
        self.amp_min_var = tk.StringVar(value="100")
        self.amp_max_var = tk.StringVar(value="500")

        ttk.Entry(size_frame, textvariable=self.amp_min_var, width=6).pack(side=tk.LEFT, padx=5)
        ttk.Label(size_frame, text="-").pack(side=tk.LEFT)
        ttk.Entry(size_frame, textvariable=self.amp_max_var, width=6).pack(side=tk.LEFT, padx=5)
        ttk.Label(size_frame, text="bp").pack(side=tk.LEFT)

        # Primer Tm
        tm_frame = ttk.Frame(frame)
        tm_frame.pack(fill=tk.X, pady=2)

        ttk.Label(tm_frame, text="Primer Tm:").pack(side=tk.LEFT)
        self.tm_min_var = tk.StringVar(value="57")
        self.tm_max_var = tk.StringVar(value="64")

        ttk.Entry(tm_frame, textvariable=self.tm_min_var, width=6).pack(side=tk.LEFT, padx=5)
        ttk.Label(tm_frame, text="-").pack(side=tk.LEFT)
        ttk.Entry(tm_frame, textvariable=self.tm_max_var, width=6).pack(side=tk.LEFT, padx=5)
        ttk.Label(tm_frame, text="°C").pack(side=tk.LEFT)

        # Distance constraints
        dist_frame = ttk.Frame(frame)
        dist_frame.pack(fill=tk.X, pady=2)

        self.use_dist_constraint = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            dist_frame, text="Min distance from variant:",
            variable=self.use_dist_constraint
        ).pack(side=tk.LEFT)

        self.dist_var = tk.StringVar(value="40")
        ttk.Entry(dist_frame, textvariable=self.dist_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Label(dist_frame, text="bp").pack(side=tk.LEFT)

        # Intron constraint
        intron_frame = ttk.Frame(frame)
        intron_frame.pack(fill=tk.X, pady=2)

        self.use_intron_constraint = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            intron_frame, text="Min dist. from exon junction:",
            variable=self.use_intron_constraint
        ).pack(side=tk.LEFT)

        self.intron_dist_var = tk.StringVar(value="60")
        ttk.Entry(intron_frame, textvariable=self.intron_dist_var, width=5).pack(side=tk.LEFT, padx=5)
        ttk.Label(intron_frame, text="bp (intronic side)").pack(side=tk.LEFT)

        # Number of primer pairs
        num_pairs_frame = ttk.Frame(frame)
        num_pairs_frame.pack(fill=tk.X, pady=2)

        ttk.Label(num_pairs_frame, text="Number of primer pairs:").pack(side=tk.LEFT)
        self.num_primer_pairs_var = tk.StringVar(value="5")
        ttk.Entry(num_pairs_frame, textvariable=self.num_primer_pairs_var, width=5).pack(side=tk.LEFT, padx=5)

        # Advanced settings button
        ttk.Button(
            frame, text="Advanced Settings...",
            command=self._show_advanced_settings
        ).pack(fill=tk.X, pady=(5, 0))

    def _build_population_section(self):
        """Build population variant filtering section."""
        frame = ttk.LabelFrame(self.left_frame, text="Population Variants", padding=(5, 3))
        frame.pack(fill=tk.X, padx=5, pady=2)

        # Toggle to enable/disable population variant filtering
        self.filter_pop_variants = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            frame, text="Filter population variants",
            variable=self.filter_pop_variants
        ).pack(anchor='w', pady=(0, 2))

        # Population selection - radio buttons (single selection only)
        pop_frame = ttk.Frame(frame)
        pop_frame.pack(fill=tk.X)

        # Single variable for radio button selection
        self.selected_population = tk.StringVar(value='global')

        populations = [
            ('global', 'Global'),
            ('nfe', 'European (NFE)'),
            ('afr', 'African'),
            ('eas', 'East Asian'),
            ('sas', 'South Asian'),
            ('amr', 'Latino')
        ]

        # Display in 3 columns using radio buttons (2 rows × 3 columns)
        for i, (code, name) in enumerate(populations):
            rb = ttk.Radiobutton(pop_frame, text=name, variable=self.selected_population,
                                 value=code, command=self._on_population_change)
            row, col = divmod(i, 3)
            rb.grid(row=row, column=col, sticky='w', padx=3)

    def _build_homology_discrimination_section(self):
        """Homology discrimination checkbox is now inside _build_specificity_section."""
        pass

    def _build_specificity_section(self):
        """Build specificity checking and homology discrimination options."""
        frame = ttk.LabelFrame(self.left_frame, text="Specificity & Homology", padding=(5, 3))
        frame.pack(fill=tk.X, padx=5, pady=2)

        # --- Homology discrimination (primer design phase) ---
        self.use_homology_discrimination = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            frame,
            text="Pseudogene discrimination (auto BLAST)",
            variable=self.use_homology_discrimination
        ).pack(anchor='w')

        ttk.Separator(frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(3, 2))

        # --- Post-design specificity tools ---
        self.use_primer_blast_check = tk.BooleanVar(value=False)
        self.use_ucsc_ispcr_check = tk.BooleanVar(value=True)

        ttk.Checkbutton(
            frame, text="UCSC In-Silico PCR (online, fast)",
            variable=self.use_ucsc_ispcr_check
        ).pack(anchor='w')

        # Continuous mode controls — single line
        ttk.Separator(frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(3, 2))

        cont_row = ttk.Frame(frame)
        cont_row.pack(fill=tk.X)

        self.continuous_mode_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            cont_row, text="Continuous mode",
            variable=self.continuous_mode_var,
            command=self._on_continuous_mode_toggle
        ).pack(side=tk.LEFT)

        ttk.Label(cont_row, text="  max:").pack(side=tk.LEFT)
        self.max_continuous_var = tk.StringVar(value="20")
        self.max_continuous_entry = ttk.Entry(
            cont_row, textvariable=self.max_continuous_var,
            width=4, state=tk.NORMAL
        )
        self.max_continuous_entry.pack(side=tk.LEFT, padx=2)
        ttk.Label(cont_row, text="attempts").pack(side=tk.LEFT)

        # Minimum specific pairs row
        min_spec_row = ttk.Frame(frame)
        min_spec_row.pack(fill=tk.X)

        ttk.Label(min_spec_row, text="  Min specific pairs:").pack(side=tk.LEFT, padx=(20, 0))
        self.min_specific_pairs_var = tk.StringVar(value="3")
        self.min_specific_pairs_entry = ttk.Entry(
            min_spec_row, textvariable=self.min_specific_pairs_var,
            width=4, state=tk.NORMAL
        )
        self.min_specific_pairs_entry.pack(side=tk.LEFT, padx=2)
        ttk.Label(min_spec_row, text="(per variant)").pack(side=tk.LEFT)

    def _build_action_buttons(self):
        """Build main action buttons."""
        frame = ttk.Frame(self.left_frame)
        frame.pack(fill=tk.X, padx=5, pady=5)

        self.design_btn = ttk.Button(
            frame, text="Design Primers",
            command=self._design_primers, state=tk.DISABLED
        )
        self.design_btn.pack(fill=tk.X, pady=2)

        self.cancel_btn = ttk.Button(
            frame, text="Cancel",
            command=self._cancel_operation, state=tk.DISABLED
        )
        self.cancel_btn.pack(fill=tk.X, pady=2)

        self.skip_variant_btn = ttk.Button(
            frame, text="Skip to Next Variant",
            command=self._skip_to_next_variant, state=tk.DISABLED
        )
        self.skip_variant_btn.pack(fill=tk.X, pady=2)

        self.export_btn = ttk.Button(
            frame, text="Export Report",
            command=self._export_report, state=tk.DISABLED
        )
        self.export_btn.pack(fill=tk.X, pady=2)

    def _build_variants_table(self):
        """Build variants table in right panel."""
        frame = ttk.LabelFrame(self.right_frame, text="Variants", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=3)

        # Treeview with MAF% column for per-variant threshold
        columns = ('number', 'gene', 'transcript', 'variant', 'maf_pct', 'status')

        # Height 5 rows (increased from 3)
        self.variants_tree = ttk.Treeview(frame, columns=columns, show='headings', height=5)

        self.variants_tree.heading('number', text='#')
        self.variants_tree.heading('gene', text='Gene')
        self.variants_tree.heading('transcript', text='Transcript')
        self.variants_tree.heading('variant', text='Variant')
        self.variants_tree.heading('maf_pct', text='MAF %')
        self.variants_tree.heading('status', text='Status')

        self.variants_tree.column('number', width=40)
        self.variants_tree.column('gene', width=80)
        self.variants_tree.column('transcript', width=120)
        self.variants_tree.column('variant', width=150)
        self.variants_tree.column('maf_pct', width=60)
        self.variants_tree.column('status', width=70)

        # Scrollbars
        v_scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.variants_tree.yview)
        h_scroll = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=self.variants_tree.xview)

        self.variants_tree.configure(yscrollcommand=v_scroll.set, xscrollcommand=h_scroll.set)

        self.variants_tree.grid(row=0, column=0, sticky='nsew')
        v_scroll.grid(row=0, column=1, sticky='ns')
        h_scroll.grid(row=1, column=0, sticky='ew')

        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        # Bind selection event
        self.variants_tree.bind('<<TreeviewSelect>>', self._on_variant_select)

        # Dictionary to store per-variant MAF thresholds
        self.variant_maf_thresholds = {}  # row_number -> maf_percentage (0-100)
        self.default_maf_threshold_pct = 0.5  # Default 0.5%

    def _build_visualization_panel(self):
        """Build visualization panel with Notebook tabs: Sequence + Homology."""
        from .visualization import HomologyVisualizer

        # --- Notebook container ---
        notebook_container = ttk.Frame(self.right_frame)
        notebook_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.seq_notebook = ttk.Notebook(notebook_container)
        self.seq_notebook.pack(fill=tk.BOTH, expand=True)

        # ==============================================================
        #  Tab 1 — Sequence (existing content, unchanged logic)
        # ==============================================================
        seq_tab = ttk.Frame(self.seq_notebook, padding=5)
        self.seq_notebook.add(seq_tab, text="Sequence")

        # MAF Display Filter slider - PERCENTAGE based (0-100%)
        maf_filter_frame = ttk.Frame(seq_tab)
        maf_filter_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(maf_filter_frame, text="MAF Filter (%):").pack(side=tk.LEFT)

        # Min range entry (left side of slider) - percentage
        self.maf_range_min_var = tk.StringVar(value="0")
        self.maf_range_min_entry = ttk.Entry(maf_filter_frame, textvariable=self.maf_range_min_var, width=8)
        self.maf_range_min_entry.pack(side=tk.LEFT, padx=2)
        self.maf_range_min_entry.bind('<Return>', self._on_maf_range_change)
        self.maf_range_min_entry.bind('<FocusOut>', self._on_maf_range_change)

        # Slider for MAF display threshold - PERCENTAGE (0-1.0% default range)
        self.maf_slider_var = tk.DoubleVar(value=0.5)  # Default 0.5%
        self.maf_slider = tk.Scale(
            maf_filter_frame,
            from_=0.0,
            to=1.0,  # Percentage (0-1.0%)
            orient=tk.HORIZONTAL,
            variable=self.maf_slider_var,
            command=self._on_maf_slider_change,
            resolution=0.001,
            showvalue=False,
            length=250,
            sliderlength=12
        )
        self.maf_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)

        # Max range entry (right side of slider) - percentage (default 1.0%)
        self.maf_range_max_var = tk.StringVar(value="1.0")
        self.maf_range_max_entry = ttk.Entry(maf_filter_frame, textvariable=self.maf_range_max_var, width=8)
        self.maf_range_max_entry.pack(side=tk.LEFT, padx=2)
        self.maf_range_max_entry.bind('<Return>', self._on_maf_range_change)
        self.maf_range_max_entry.bind('<FocusOut>', self._on_maf_range_change)

        # Current value entry - percentage
        ttk.Label(maf_filter_frame, text="%").pack(side=tk.LEFT, padx=(2, 5))
        self.maf_slider_entry = ttk.Entry(maf_filter_frame, width=10)
        self.maf_slider_entry.insert(0, "0.500")
        self.maf_slider_entry.pack(side=tk.LEFT, padx=2)
        self.maf_slider_entry.bind('<Return>', self._on_maf_entry_submit)
        self.maf_slider_entry.bind('<FocusOut>', self._on_maf_entry_submit)

        # Sequence text widget
        text_frame = ttk.Frame(seq_tab)
        text_frame.pack(fill=tk.BOTH, expand=True)

        self.seq_text = tk.Text(text_frame, height=15, font=('Courier', 10), wrap=tk.NONE,
                                **theme.get_text_config('sequence'))
        seq_scroll_y = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=self.seq_text.yview)
        seq_scroll_x = ttk.Scrollbar(text_frame, orient=tk.HORIZONTAL, command=self.seq_text.xview)

        self.seq_text.configure(yscrollcommand=seq_scroll_y.set, xscrollcommand=seq_scroll_x.set)

        self.seq_text.grid(row=0, column=0, sticky='nsew')
        seq_scroll_y.grid(row=0, column=1, sticky='ns')
        seq_scroll_x.grid(row=1, column=0, sticky='ew')

        text_frame.grid_rowconfigure(0, weight=1)
        text_frame.grid_columnconfigure(0, weight=1)

        self.seq_visualizer = SequenceVisualizer(self.seq_text)

        # Keep gene_visualizer as None (not used anymore but may be referenced)
        self.gene_visualizer = None

        # ==============================================================
        #  Tab 2 — Homology Analysis
        # ==============================================================
        homology_tab = ttk.Frame(self.seq_notebook, padding=5)
        self.seq_notebook.add(homology_tab, text="Homology Analysis")
        self._build_homology_tab(homology_tab, HomologyVisualizer)

    # ------------------------------------------------------------------
    #  Homology Analysis tab
    # ------------------------------------------------------------------

    def _build_homology_tab(self, parent, HomologyVisualizerClass):
        """Build contents of the Homology Analysis tab."""
        # --- Control bar ---
        control_frame = ttk.Frame(parent)
        control_frame.pack(fill=tk.X, pady=(0, 5))

        self.homology_run_btn = ttk.Button(
            control_frame, text="Run Homology Analysis",
            command=self._run_homology_analysis)
        self.homology_run_btn.pack(side=tk.LEFT, padx=(0, 8))

        self.homology_status_var = tk.StringVar(
            value="Select a variant and click 'Run Homology Analysis'")
        ttk.Label(
            control_frame, textvariable=self.homology_status_var
        ).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        self.homology_progress = ttk.Progressbar(
            control_frame, mode='indeterminate', length=120)
        self.homology_progress.pack(side=tk.RIGHT, padx=5)

        # --- Hits summary table ---
        summary_frame = ttk.LabelFrame(parent, text="Hits summary", padding=3)
        summary_frame.pack(fill=tk.X, padx=0, pady=(0, 5))

        columns = ('chr', 'position', 'strand', 'identity', 'mismatches',
                   'score', 'length')
        self.homology_tree = ttk.Treeview(
            summary_frame, columns=columns, show='headings', height=5)

        col_cfg = [
            ('chr',         'Chr',       65),
            ('position',    'Position',  170),
            ('strand',      'Strand',    50),
            ('identity',    'Identity',  75),
            ('mismatches',  'MM',        50),
            ('score',       'Score',     65),
            ('length',      'Length',     65),
        ]
        for col_id, heading, width in col_cfg:
            self.homology_tree.heading(col_id, text=heading)
            self.homology_tree.column(col_id, width=width, minwidth=40)

        tree_scroll = ttk.Scrollbar(
            summary_frame, orient=tk.VERTICAL,
            command=self.homology_tree.yview)
        self.homology_tree.configure(yscrollcommand=tree_scroll.set)

        self.homology_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.homology_tree.bind('<<TreeviewSelect>>', self._on_homology_hit_select)

        # Row colours
        self.homology_tree.tag_configure('primary', background=theme.HOMOLOGY_PRIMARY)
        self.homology_tree.tag_configure('same_chr', background=theme.HOMOLOGY_SAME_CHR)
        self.homology_tree.tag_configure('other_chr', background=theme.HOMOLOGY_OTHER_CHR)

        # --- BLAST-like alignment text ---
        align_frame = ttk.LabelFrame(parent, text="Alignment details", padding=3)
        align_frame.pack(fill=tk.BOTH, expand=True, padx=0, pady=0)

        align_inner = ttk.Frame(align_frame)
        align_inner.pack(fill=tk.BOTH, expand=True)

        self.homology_text = tk.Text(
            align_inner, height=12, font=('Courier', 10), wrap=tk.NONE,
            **theme.get_text_config('homology'))
        hom_scroll_y = ttk.Scrollbar(
            align_inner, orient=tk.VERTICAL,
            command=self.homology_text.yview)
        hom_scroll_x = ttk.Scrollbar(
            align_inner, orient=tk.HORIZONTAL,
            command=self.homology_text.xview)
        self.homology_text.configure(
            yscrollcommand=hom_scroll_y.set,
            xscrollcommand=hom_scroll_x.set)

        self.homology_text.grid(row=0, column=0, sticky='nsew')
        hom_scroll_y.grid(row=0, column=1, sticky='ns')
        hom_scroll_x.grid(row=1, column=0, sticky='ew')

        align_inner.grid_rowconfigure(0, weight=1)
        align_inner.grid_columnconfigure(0, weight=1)

        self.homology_visualizer = HomologyVisualizerClass(self.homology_text)

    # ------------------------------------------------------------------
    #  Homology Analysis — run / display / select
    # ------------------------------------------------------------------

    def _run_homology_analysis(self):
        """Run BLAST homology analysis for the currently selected variant."""
        if not self._selected_variant:
            from tkinter import messagebox
            messagebox.showinfo(
                "Information",
                "Select a variant from the table first.")
            return

        variant = self._selected_variant
        row_num = variant.row_number
        cached = self._sequence_cache.get(row_num)

        if not cached or not cached.get('sequence'):
            from tkinter import messagebox
            messagebox.showwarning(
                "Warning",
                "Sequence not yet available. "
                "Wait for caching to complete.")
            return

        sequence = cached['sequence']

        # Disable button, start progress
        self.homology_run_btn.configure(state=tk.DISABLED)
        self.homology_progress.start(10)
        self.homology_status_var.set("Preparing analysis...")

        def _analysis_thread():
            try:
                # Lazy-init HomologyAnalyzer
                if self.homology_analyzer is None:
                    from ..primer.homology_analyzer import HomologyAnalyzer
                    self.homology_analyzer = HomologyAnalyzer()

                # Ensure BLAST ready (will download/index if needed)
                def _status(msg):
                    self.root.after(
                        0, lambda m=msg: self.homology_status_var.set(m))

                ready, error = self.homology_analyzer.ensure_blast_ready(_status)
                if not ready:
                    self.root.after(0, lambda: self._homology_error(error))
                    return

                # Build query name
                gene = variant.gene_symbol or "?"
                hgvs = variant.hgvs_c or "?"
                query_name = f"{gene}_{hgvs}"

                # Get source coordinates
                src_chr = ""
                src_pos = 0
                if variant.genomic_position:
                    src_chr = variant.genomic_position.chromosome
                    src_pos = variant.genomic_position.start

                result = self.homology_analyzer.analyze(
                    sequence=sequence,
                    query_name=query_name,
                    source_chromosome=src_chr,
                    source_position=src_pos,
                    progress_callback=_status,
                )

                # Store genomic start of query sequence for coordinate translation
                # (needed by primer designer to align discriminating positions)
                result.query_genomic_start = cached.get('seq_start', 0)

                self.root.after(0, lambda r=result: self._display_homology_results(r))

            except Exception as e:
                self.root.after(0, lambda: self._homology_error(str(e)))
            finally:
                self.root.after(
                    0, lambda: self.homology_run_btn.configure(state=tk.NORMAL))
                self.root.after(0, lambda: self.homology_progress.stop())

        thread = threading.Thread(target=_analysis_thread, daemon=True)
        thread.start()

    def _homology_error(self, message: str):
        """Display a homology analysis error."""
        from tkinter import messagebox
        self.homology_status_var.set(f"Error: {message[:80]}")
        self.homology_progress.stop()
        self.homology_run_btn.configure(state=tk.NORMAL)
        messagebox.showerror("Homology Analysis Error", message)

    def _display_homology_results(self, result):
        """Populate the hits table and alignment text from a HomologyResult."""
        from ..primer.homology_analyzer import HomologyAnalyzer

        self._homology_result = result

        # Store per-variant for primer design integration
        if self._selected_variant:
            self._homology_results_by_variant[self._selected_variant.row_number] = result

        # --- Clear & populate Treeview ---
        for item in self.homology_tree.get_children():
            self.homology_tree.delete(item)

        for i, hit in enumerate(result.hits):
            src_norm = HomologyAnalyzer.normalize_chr(result.query_chromosome)
            hit_norm = HomologyAnalyzer.normalize_chr(hit.chromosome)

            if hit.is_primary:
                row_tag = 'primary'
            elif hit_norm == src_norm:
                row_tag = 'same_chr'
            else:
                row_tag = 'other_chr'

            end_pos = hit.position + hit.aligned_length
            self.homology_tree.insert('', 'end', iid=str(i), values=(
                HomologyAnalyzer.format_chr(hit.chromosome),
                f"{hit.position:,}-{end_pos:,}",
                hit.strand,
                f"{hit.percent_identity}%",
                hit.mismatches,
                hit.alignment_score,
                f"{hit.aligned_length} bp",
            ), tags=(row_tag,))

        # --- Alignment visualisation ---
        # Calculate variant position inside query sequence (0-based).
        # The query spans [source_position - flank, source_position + flank],
        # so the variant is at the middle (query_length // 2).
        variant_query_pos = result.query_length // 2 if result.query_length > 0 else -1
        self.homology_visualizer.display_results(result, variant_query_pos)

        # --- Status ---
        n_hits = len(result.hits)
        n_other = sum(
            1 for h in result.hits
            if HomologyAnalyzer.normalize_chr(h.chromosome)
            != HomologyAnalyzer.normalize_chr(result.query_chromosome))

        if result.error:
            self.homology_status_var.set(f"Error: {result.error[:80]}")
        else:
            self.homology_status_var.set(
                f"Found {n_hits} hit(s) "
                f"({n_other} on other chromosomes)")

    def _on_homology_hit_select(self, event):
        """Display alignment details for the selected hit."""
        selection = self.homology_tree.selection()
        if not selection:
            return

        idx = int(selection[0])
        self.homology_visualizer.display_hit(idx)

    def _auto_run_homology_for_groups(self, groups):
        """
        Automatically run homology analysis for variant groups that don't
        already have a pre-computed result.  Called during primer design
        when "Pseudogene discrimination" checkbox is enabled.

        Runs sequentially (BLAST is CPU-heavy; parallel adds little benefit
        and may overload the system).  Skips groups where the result was
        already computed via the manual "Run Homology Analysis" button.
        """
        from ..primer.homology_analyzer import HomologyAnalyzer

        # Lazy-init analyzer
        if self.homology_analyzer is None:
            self.homology_analyzer = HomologyAnalyzer()

        # Ensure BLAST ready (download genome + build DB if first time)
        ready, error = self.homology_analyzer.ensure_blast_ready(
            lambda msg: self._set_status(f"Homology setup: {msg}")
        )
        if not ready:
            self.log_warning(
                f"[AutoHomology] BLAST not ready, skipping: {error}"
            )
            self._set_status(
                f"Homology analysis skipped — BLAST not available: "
                f"{error[:80]}"
            )
            return

        # Collect variants that need homology analysis
        variants_needing_analysis = []
        for group in groups:
            if not group.variants:
                continue
            first = group.variants[0]
            if first.row_number not in self._homology_results_by_variant:
                variants_needing_analysis.append(first)

        if not variants_needing_analysis:
            self.log_info("[AutoHomology] All groups already have homology results")
            return

        total = len(variants_needing_analysis)
        self.log_info(
            f"[AutoHomology] Running homology analysis for {total} variant(s)..."
        )

        for i, variant in enumerate(variants_needing_analysis):
            if self._is_cancel_requested():
                break

            row_num = variant.row_number
            cached = self._sequence_cache.get(row_num)
            if not cached or not cached.get('sequence'):
                self.log_warning(
                    f"[AutoHomology] No cached sequence for variant "
                    f"{row_num}, skipping"
                )
                continue

            gene = variant.gene_symbol or "?"
            hgvs = variant.hgvs_c or "?"
            self._set_status(
                f"Homology analysis: {i + 1}/{total} — {gene} {hgvs}..."
            )
            self._set_progress(10 + (10 * i / total))

            src_chr = ""
            src_pos = 0
            if variant.genomic_position:
                src_chr = variant.genomic_position.chromosome
                src_pos = variant.genomic_position.start

            try:
                result = self.homology_analyzer.analyze(
                    sequence=cached['sequence'],
                    query_name=f"{gene}_{hgvs}",
                    source_chromosome=src_chr,
                    source_position=src_pos,
                )
                result.query_genomic_start = cached.get('seq_start', 0)
                self._homology_results_by_variant[row_num] = result

                n_hits = len(result.hits)
                n_secondary = sum(
                    1 for h in result.hits
                    if not h.is_primary and not h.is_supplementary
                )
                self.log_info(
                    f"[AutoHomology] {gene}_{hgvs}: "
                    f"{n_hits} hit(s), {n_secondary} secondary"
                )
            except Exception as e:
                self.log_warning(
                    f"[AutoHomology] Failed for {gene}_{hgvs}: {e}"
                )

        self._set_progress(20)
        self._set_status("Homology analysis complete, designing primers...")

    def _build_results_panel(self):
        """Build results panel."""
        frame = ttk.LabelFrame(self.right_frame, text="Design Results", padding=5)
        frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=3)

        self.results_text = tk.Text(frame, height=4, wrap=tk.WORD,
                                     **theme.get_text_config('results'))
        results_scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.results_text.yview)

        self.results_text.configure(yscrollcommand=results_scroll.set)

        self.results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        results_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Configure text tags for formatting (from theme)
        for tag_name, tag_opts in theme.get_results_tags().items():
            self.results_text.tag_configure(tag_name, **tag_opts)

    def _build_status_bar(self):
        """Build status bar at bottom."""
        # Create a separator above status bar
        separator = ttk.Separator(self.root, orient=tk.HORIZONTAL)
        separator.pack(fill=tk.X, side=tk.BOTTOM, pady=(5, 0))

        self.status_frame = ttk.Frame(self.root, padding=(10, 8))
        self.status_frame.pack(fill=tk.X, side=tk.BOTTOM)

        # Animated spinner using Canvas
        self._spinner_canvas = tk.Canvas(self.status_frame, width=20, height=20,
                                          highlightthickness=0, bg=theme.BG)
        self._spinner_canvas.pack(side=tk.LEFT, padx=(10, 5))
        self._spinner_angle = 0
        self._spinner_running = False
        self._spinner_arc = None
        self._draw_spinner_idle()

        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(self.status_frame, textvariable=self.status_var)
        self.status_label.pack(side=tk.LEFT, padx=5)

        # Progress bar - make it more visible
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(
            self.status_frame, variable=self.progress_var,
            mode='determinate', length=250
        )
        self.progress_bar.pack(side=tk.RIGHT, padx=10, pady=5)

    def _draw_spinner_idle(self):
        """Draw idle spinner (gray circle)."""
        self._spinner_canvas.delete("all")
        self._spinner_canvas.create_oval(2, 2, 18, 18, outline=theme.SPINNER_IDLE, width=2)

    def _draw_spinner_frame(self):
        """Draw one frame of the animated spinner."""
        self._spinner_canvas.delete("all")
        # Background circle
        self._spinner_canvas.create_oval(2, 2, 18, 18, outline=theme.SPINNER_BG, width=2)
        # Animated arc (like iOS/Chrome spinner)
        self._spinner_canvas.create_arc(2, 2, 18, 18,
                                         start=self._spinner_angle, extent=90,
                                         outline=theme.SPINNER_ARC, width=2, style=tk.ARC)

    def _animate_spinner(self):
        """Animate the spinner."""
        if self._spinner_running:
            self._spinner_angle = (self._spinner_angle + 15) % 360
            self._draw_spinner_frame()
            self.root.after(50, self._animate_spinner)

    def _start_spinner(self):
        """Start the spinner animation."""
        if not self._spinner_running:
            self._spinner_running = True
            self._animate_spinner()

    def _stop_spinner(self):
        """Stop the spinner animation."""
        self._spinner_running = False
        self._draw_spinner_idle()

    def _build_help_panel(self):
        """Build the contextual help panel on the right side."""
        # Title
        title_frame = ttk.Frame(self.help_frame)
        title_frame.pack(fill=tk.X, padx=5, pady=(5, 0))
        ttk.Label(title_frame, text="📖 Quick Guide", style='Title.TLabel').pack(anchor='w')
        ttk.Separator(title_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(3, 0))

        # Scrollable help text
        help_text_frame = ttk.Frame(self.help_frame)
        help_text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.help_text = tk.Text(
            help_text_frame, wrap=tk.WORD, font=('Helvetica', 10),
            padx=12, pady=10,
            spacing1=2, spacing2=1, spacing3=4,
            **theme.get_text_config('help')
        )
        help_scroll = ttk.Scrollbar(help_text_frame, orient=tk.VERTICAL, command=self.help_text.yview)
        self.help_text.configure(yscrollcommand=help_scroll.set)

        self.help_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        help_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Configure text tags for formatting (from theme)
        for tag_name, tag_opts in theme.get_help_tags().items():
            self.help_text.tag_configure(tag_name, **tag_opts)

        # Make text read-only
        self.help_text.configure(state=tk.DISABLED)

    def _update_help_context(self, stage: str):
        """
        Update contextual help based on current workflow stage.

        Stages:
          - 'startup': Initial state, overview of the workflow
          - 'file_loaded': After loading and validating variants
          - 'caching': While caching sequence data
          - 'ready': Sequences cached, ready to design primers
          - 'designing': During primer design
          - 'results': After primer design is complete
        """
        help_content = self._get_help_content(stage)

        self.help_text.configure(state=tk.NORMAL)
        self.help_text.delete('1.0', tk.END)

        for tag, text in help_content:
            self.help_text.insert(tk.END, text, tag)

        self.help_text.configure(state=tk.DISABLED)
        self.help_text.see('1.0')  # Scroll to top

    def _get_help_content(self, stage: str):
        """Return help content as list of (tag, text) tuples for a given stage."""

        # Common workflow steps shown at the top of every stage
        workflow_steps = [
            ('h1', 'Workflow Overview\n'),
            ('step', '① Load & Validate → '),
            ('step', '② Review Sequences → '),
            ('step', '③ Design Primers → '),
            ('step', '④ Export Report\n\n'),
        ]

        if stage == 'startup':
            return workflow_steps + self._help_startup()
        elif stage == 'file_loaded':
            return workflow_steps + self._help_file_loaded()
        elif stage == 'caching':
            return workflow_steps + self._help_caching()
        elif stage == 'ready':
            return workflow_steps + self._help_ready()
        elif stage == 'designing':
            return workflow_steps + self._help_designing()
        elif stage == 'results':
            return workflow_steps + self._help_results()
        else:
            return workflow_steps + self._help_startup()

    def _help_startup(self):
        """Help content for the initial startup state."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Getting Started  \n\n'),

            ('h2', 'Step 1 — Prepare Your Input File\n'),
            ('body', 'Create a CSV or Excel (.xlsx) file with the following columns:\n'),
            ('bullet', '• Gene — gene symbol (e.g. NBN, BRCA1)\n'),
            ('bullet', '• Transcript — RefSeq accession (e.g. NM_002485.5)\n'),
            ('bullet', '• Position — HGVS c. notation (e.g. c.657_661delACAAA)\n\n'),

            ('tip', 'Tip: Use MANE Select transcripts when possible for best '
                    'cross-database compatibility.\n\n'),

            ('h2', 'Step 2 — Choose Reference Genome\n'),
            ('body', 'Select the reference assembly used for your variant coordinates:\n'),
            ('bullet', '• GRCh38 (HG38) — current standard, recommended for new data\n'),
            ('bullet', '• GRCh37 (HG19) — legacy; use only if your variants were annotated against HG19\n\n'),
            ('warning', 'Mismatched assemblies will produce incorrect genomic coordinates '
                       'and wrong primer sequences!\n\n'),

            ('h2', 'Step 3 — Load the File\n'),
            ('body', 'Click "Browse" to select your file, then click "Load and Validate".\n\n'),
            ('body', 'The application will:\n'),
            ('bullet', '• Parse your variant list\n'),
            ('bullet', '• Validate HGVS nomenclature\n'),
            ('bullet', '• Check MANE transcript status\n'),
            ('bullet', '• Map CDS positions to genomic coordinates\n\n'),

            ('tip', 'Tip: Check "Skip online validation" for a quick offline preview, '
                    'but online validation is strongly recommended before designing primers.\n\n'),

            ('h2', 'Flanking Region\n'),
            ('body', 'The flanking region (default: 250 bp) determines how much '
                    'genomic sequence is fetched on each side of the variant. '
                    'This provides the search space for Primer3 to find primers.\n'),
            ('bullet', '• Larger flanking → more primer candidates, but slower\n'),
            ('bullet', '• Smaller flanking → faster, but may miss optimal primers\n\n'),

            ('tip', 'Tip: 250 bp is suitable for most PCR designs.\n'),
        ]

    def _help_file_loaded(self):
        """Help content after file loading and validation."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Reviewing Variants  \n\n'),

            ('h2', 'Variants Table\n'),
            ('body', 'Your variants are now loaded in the table above. Check:\n'),
            ('bullet', '• Status column — "Valid" means the variant passed all checks\n'),
            ('bullet', '• "Warning" — minor issues; primers can still be designed\n'),
            ('bullet', '• "Error" — design may fail for this variant\n\n'),

            ('body', 'Click on any variant to view its genomic sequence in the '
                    'Sequence panel. The target variant position is highlighted.\n\n'),

            ('h2', 'Background Sequence Caching\n'),
            ('body', 'Sequences are now being downloaded and cached in the background. '
                    'The spinner in the bottom-left corner indicates this process. '
                    'You can still browse the interface while this is running.\n\n'),

            ('h2', 'What to Check Before Designing\n'),
            ('step', '1. Verify the variant table has the expected entries\n'),
            ('step', '2. Set your design parameters (left panel)\n'),
            ('step', '3. Configure MAF filtering if needed\n'),
            ('step', '4. Wait for sequence caching to complete\n'),
            ('step', '5. Click "Design Primers"\n\n'),

            ('tip', 'Tip: Use this time to adjust parameters while sequences are caching.\n'),
        ]

    def _help_caching(self):
        """Help content while sequences are being cached."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Caching Sequences  \n\n'),

            ('h2', 'What Is Happening\n'),
            ('body', 'The application is downloading genomic sequences and population '
                    'variant data for each of your variants. This data is cached '
                    'locally for fast access during primer design.\n\n'),
            ('body', 'Progress is shown in the bottom status bar (90–100%).\n\n'),

            ('h2', 'While You Wait...\n'),
            ('body', 'Review and configure your design parameters:\n\n'),

            ('h2', 'Design Parameters Explained\n'),
            ('param', 'Amplicon Size\n'),
            ('body', 'Target PCR product size range. Primer3 will try to place primers '
                    'so the product falls within this window.\n\n'),

            ('param', 'Primer Tm\n'),
            ('body', 'Melting temperature range for primers. Default 57–64°C with '
                    'optimal at 60°C works well for most standard protocols.\n\n'),

            ('param', 'Min Distance from Variant\n'),
            ('body', 'Minimum number of nucleotides between the variant and the '
                    'nearest primer edge. This ensures the variant is well inside '
                    'the amplicon, not at the primer binding site.\n'),
            ('tip', 'Tip: 40 nt is the default. Increase to 60 nt for extra safety, '
                    'or decrease to 30 nt if Primer3 cannot find solutions in a '
                    'constrained region.\n\n'),

            ('param', 'Min Distance from Exon Junction\n'),
            ('body', 'Minimum distance (in bp) that a primer must maintain from the '
                    'nearest exon/intron boundary on the intronic side. The zone '
                    'within this distance from each exon junction is excluded from '
                    'primer placement. This prevents primers from landing near '
                    'splice sites, where secondary structure and regulatory elements '
                    'may interfere with amplification.\n'),
            ('tip', 'Tip: Disabled by default. Enable and set to 60 bp to exclude '
                    'primers near splice sites (covers the branch point and '
                    'polypyrimidine tract). Increase to 100+ bp if you see non-specific '
                    'products near exon junctions.\n\n'),

            ('param', 'Number of Primer Pairs\n'),
            ('body', 'How many candidate primer pairs Primer3 should return per variant '
                    'in standard (non-continuous) mode. In continuous mode, the '
                    '"min specific pairs" setting determines how many verified-specific '
                    'pairs will be found (default: 3).\n'),
        ]

    def _help_ready(self):
        """Help content when ready to design primers."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Ready to Design  \n\n'),

            ('h2', 'Sequence Caching Complete\n'),
            ('body', 'All genomic sequences and population variants are cached.\n\n'),

            ('h2', 'Population Variant (MAF) Filtering\n'),
            ('body', 'Known polymorphic sites (SNPs from gnomAD) can be masked '
                    'with "N" in the sequence used for primer design. This prevents '
                    'primers from binding at polymorphic positions, which could cause '
                    'allele-specific amplification failure.\n\n'),

            ('param', 'MAF Filter Slider\n'),
            ('body', 'Controls the Minor Allele Frequency threshold (%) for masking. '
                    'Positions with MAF ≥ threshold will be masked with N.\n'),
            ('bullet', '• Default value: 0.5% (masks variants with MAF ≥ 0.5%)\n'),
            ('bullet', '• Lower threshold → more positions masked → stricter\n'),
            ('bullet', '• Higher threshold → fewer positions masked → more lenient\n\n'),

            ('param', 'Adjustable Slider Range (Min / Max fields)\n'),
            ('body', 'The two input fields on either side of the slider control its '
                    'minimum and maximum range. The default range is 0–1.0%, but you can '
                    'change these values to any range you need:\n'),
            ('bullet', '• For very rare variants: set range to 0–0.01%\n'),
            ('bullet', '• For common variants: set range to 0–5% or even 0–50%\n'),
            ('body', 'Simply type a new number into the Min or Max field and press Enter. '
                    'The slider will automatically adjust its scale and resolution.\n\n'),

            ('param', 'Per-Variant MAF Threshold\n'),
            ('body', 'Each variant starts with the default MAF threshold (0.5%). '
                    'You can set a different threshold for individual variants:\n'),
            ('step', '1. Select a variant in the table\n'),
            ('step', '2. Adjust the MAF slider to the desired value\n'),
            ('step', '3. Click "Apply to Selected"\n'),
            ('body', 'The MAF % column in the table shows each variant\'s current threshold. '
                    'During primer design, each variant uses its own MAF threshold '
                    'for FASTA masking.\n\n'),

            ('param', 'Population Selection\n'),
            ('body', 'Choose which population\'s allele frequencies to use for filtering:\n'),
            ('bullet', '• Global — worldwide average (recommended for general use)\n'),
            ('bullet', '• European (NFE) — non-Finnish European\n'),
            ('bullet', '• African, East Asian, South Asian, Latino\n\n'),
            ('tip', 'Tip: Choose the population matching your patient cohort for '
                    'the most relevant filtering.\n\n'),

            ('param', 'Filter Population Variants Checkbox\n'),
            ('body', 'When checked, population variants above the MAF threshold will '
                    'be masked in the FASTA sequence provided to Primer3. '
                    'Uncheck to disable masking entirely.\n\n'),

            ('h2', 'Specificity Checking (Optional)\n'),
            ('body', 'Enable specificity tools to verify that your primers amplify '
                    'ONLY the intended target and not pseudogenes or other genomic loci.\n\n'),
            ('param', 'Available Tools\n'),
            ('bullet', '• UCSC In-Silico PCR (online) — simulates PCR against the genome '
                      'assembly. Fast (a few seconds per pair). Enabled by default.\n\n'),
            ('tip', 'Tip: For diagnostic genetics, always verify specificity to catch '
                    'pseudogene co-amplification.\n\n'),

            ('h2', 'Continuous Design Mode\n'),
            ('body', 'When enabled (requires at least one specificity tool), the '
                    'application designs primers one at a time, checks each for '
                    'specificity immediately, and keeps trying until the required '
                    'number of SPECIFIC pairs is found or limits are reached.\n\n'),
            ('param', 'How It Works\n'),
            ('bullet', '• Primers are designed in batches and checked individually\n'),
            ('bullet', '• Each pair is displayed live with its specificity verdict\n'),
            ('bullet', '• Stops when the minimum number of SPECIFIC pairs is found '
                      '(default: 3 per variant)\n'),
            ('bullet', '• Stops if Primer3 keeps returning the same pairs (exhausted)\n'),
            ('bullet', '• Max attempts: limits how many unique pairs to try\n\n'),
            ('param', 'Min Specific Pairs\n'),
            ('body', 'Sets how many specific primer pairs must be found before '
                    'moving to the next variant (default: 3). This ensures you '
                    'always have backup pairs in case the top-ranked pair fails '
                    'in the lab.\n\n'),
            ('param', 'Controls During Continuous Mode\n'),
            ('bullet', '• "Skip to Next Variant" — stop current variant, move on\n'),
            ('bullet', '• "Cancel" — stop entirely, keep all results found so far\n\n'),
            ('tip', 'Tip: Use Continuous Mode for genes with pseudogenes '
                   '(e.g. CYP2D6, PMS2, SMN1) where finding specific primers is '
                   'critical. The mode automatically finds multiple specific pairs '
                   'for each variant.\n\n'),

            ('h2', 'Ready to Go\n'),
            ('highlight', 'Click "Design Primers" to start the primer design process.\n'),
        ]

    def _help_designing(self):
        """Help content during primer design."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Designing Primers  \n\n'),

            ('h2', 'Design In Progress\n'),
            ('body', 'The application is now:\n\n'),
            ('step', '1. Generating masked FASTA sequences\n'),
            ('body', '   Polymorphic positions above your MAF threshold are replaced '
                    'with "N" so Primer3 avoids them.\n\n'),
            ('step', '2. Grouping nearby variants\n'),
            ('body', '   Variants close enough on the genome are grouped to share '
                    'a single amplicon when possible.\n\n'),
            ('step', '3. Running Primer3\n'),
            ('body', '   For each variant/group, Primer3 searches the flanking sequence '
                    'for optimal primer pairs respecting all your constraints '
                    '(Tm, GC%, amplicon size, distance from variant, intron depth).\n\n'),

            ('h2', 'What Primer3 Considers\n'),
            ('bullet', '• Tm match between forward and reverse primers\n'),
            ('bullet', '• GC content within specified range\n'),
            ('bullet', '• Low self-complementarity (avoid hairpins)\n'),
            ('bullet', '• Low 3\' end complementarity (avoid primer dimers)\n'),
            ('bullet', '• Product size within your specified range\n'),
            ('bullet', '• Excluded regions (variant zone, deep intronic)\n\n'),

            ('step', '4. Checking specificity (if enabled)\n'),
            ('body', '   If you enabled specificity tools (UCSC isPCR), '
                    'each primer pair is checked for off-target '
                    'amplification after Primer3 finishes.\n\n'),

            ('step', '5. Continuous mode (if enabled)\n'),
            ('body', '   In continuous mode, design and specificity alternate:\n'),
            ('bullet', '   • Design batch → check each pair → display live → repeat\n'),
            ('bullet', '   • Use "Skip to Next Variant" to advance to the next gene\n'),
            ('bullet', '   • Results appear in real time in the panel above\n\n'),

            ('body', 'Progress is shown in the status bar at the bottom.\n'),
            ('tip', 'Tip: Click "Cancel" at any time to stop the design process.\n'),
        ]

    def _help_results(self):
        """Help content after primer design is complete."""
        return [
            ('current_step', '  ▸ YOU ARE HERE: Reviewing Results  \n\n'),
        ] + self._generate_methodology() + [
            ('h1', '\n─── Understanding the Results ───\n\n'),

            ('h2', 'Result Format\n'),
            ('body', 'Each variant shows its designed primer pairs in the '
                    'Design Results panel.\n\n'),

            ('param', 'Primer Pair Line\n'),
            ('body', 'Shows forward (F) and reverse (R) primer sequences in '
                    '5\'→3\' orientation.\n\n'),

            ('param', 'Structure Diagram\n'),
            ('body', 'A proportional ASCII visualization:\n'),
            ('bullet', '• ▶▶▶ — forward primer (length proportional to nt)\n'),
            ('bullet', '• ─── — gap between primer and target\n'),
            ('bullet', '• ◆ — target variant position\n'),
            ('bullet', '• ◀◀◀ — reverse primer\n'),
            ('body', 'Bar widths are scaled proportionally across all amplicons, '
                    'so you can visually compare product sizes.\n\n'),

            ('param', 'Distance Annotation\n'),
            ('body', 'Below the diagram: F:22nt ─87nt─ [VAR] ─95nt─ R:22nt\n'),
            ('body', 'Shows primer lengths and distances to the variant.\n\n'),

            ('h2', 'FASTA Files\n'),
            ('body', 'For each variant, a FASTA file was generated containing '
                    'the masked sequence used for primer design. Click the blue '
                    'file link to open it.\n\n'),
            ('tip', 'Tip: FASTA files include "N" at masked polymorphic positions. '
                    'This is the exact sequence Primer3 used, so you can verify '
                    'primer binding sites.\n\n'),

            ('h2', 'Design Failures\n'),
            ('body', 'If Primer3 could not find primers for a variant, suggestions '
                    'are shown in orange. Common fixes:\n'),
            ('bullet', '• Increase amplicon size range\n'),
            ('bullet', '• Reduce min distance from variant\n'),
            ('bullet', '• Widen Tm range\n'),
            ('bullet', '• Increase flanking region size\n'),
            ('bullet', '• Lower MAF threshold (less masking)\n\n'),

            ('h2', 'Export\n'),
            ('body', 'Click "Export Report" to save a full HTML report with all '
                    'primer sequences, parameters, methodology, and visualizations.\n\n'),

            ('h2', 'Specificity Verdicts\n'),
            ('body', 'If you enabled specificity tools before designing, each '
                    'primer pair shows a colored verdict:\n'),
            ('bullet', '• ✓ SPECIFIC — all enabled tools confirm only the target is amplified\n'),
            ('bullet', '• ✓ LIKELY SPECIFIC — checked tools say specific, some unavailable\n'),
            ('bullet', '• ✗ OFF-TARGET DETECTED — primers may amplify non-target sites\n'),
            ('bullet', '• ⚠ PSEUDOGENE RISK — off-target on same chromosome '
                      '(likely pseudogene or segmental duplication)\n'),
            ('bullet', '• ? INCONCLUSIVE — tools were unavailable or errored\n\n'),

            ('body', 'If no specificity tools were enabled, no verdict is shown. '
                    'You can re-design with tools enabled to check specificity.\n\n'),

            ('warning', 'For diagnostic genetics: always prefer primers marked '
                       'SPECIFIC. Avoid PSEUDOGENE RISK primers — they may produce '
                       'false results in clinical sequencing.\n\n'),

            ('highlight', 'Your primers are ready for ordering!\n'),
        ]

    def _generate_methodology(self):
        """
        Generate a dynamic methodology section based on the actual design
        parameters, variants, results, and tools used. Produces text suitable
        for inclusion in a scientific report or publication methods section.
        """
        content = []
        content.append(('h1', '─── Methodology ───\n\n'))
        content.append(('body',
            'The following methodology text has been generated automatically '
            'based on the parameters, databases, and tools used during this '
            'primer design session. It can be adapted for use in a laboratory '
            'report or publication.\n\n'))

        # -----------------------------------------------------------------
        # 1. Gather all relevant state
        # -----------------------------------------------------------------
        params = self.design_parameters
        assembly = self.assembly_var.get()
        assembly_full = "GRCh38/hg38" if assembly == "GRCh38" else "GRCh37/hg19"
        variants = self.project_state.variants or []
        results = self.project_state.design_results or []
        pop_code = self.selected_population.get()
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

        # Collect unique genes and transcripts
        gene_transcript_pairs = []
        seen_gt = set()
        for v in variants:
            key = (v.gene_symbol, v.transcript_accession)
            if key not in seen_gt:
                seen_gt.add(key)
                gene_transcript_pairs.append(key)

        # Specificity tools used
        used_ucsc = self.use_ucsc_ispcr_check.get()
        used_blast_check = self.use_primer_blast_check.get()
        used_homology = params.use_homology_discrimination
        continuous = self.continuous_mode_var.get()

        # MAF thresholds — per-variant from GUI dict, fallback to default
        per_variant_maf = {}  # row_number -> percentage
        for v in variants:
            per_variant_maf[v.row_number] = self.variant_maf_thresholds.get(
                v.row_number, self.default_maf_threshold_pct
            )
        unique_mafs = sorted(set(per_variant_maf.values()))
        filter_pop = params.filter_population_variants
        flanking = self.flank_size_var.get()

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

        # -----------------------------------------------------------------
        # 2. Build methodology text
        # -----------------------------------------------------------------
        content.append(('h2', 'Primer Design\n'))

        # Target description
        if gene_transcript_pairs:
            gene_list = ', '.join(
                f"{g} ({t})" for g, t in gene_transcript_pairs
            )
            content.append(('body',
                f'PCR primers were designed to amplify regions containing '
                f'variant(s) in the following gene(s) and transcript(s): '
                f'{gene_list}. '))
        else:
            content.append(('body',
                'PCR primers were designed to amplify regions containing '
                'the target variant(s). '))

        # Variant enumeration
        if variants:
            var_descriptions = []
            for v in variants:
                var_descriptions.append(
                    f"{v.gene_symbol} {v.transcript_accession}:{v.hgvs_c}"
                )
            content.append(('body',
                f'A total of {len(variants)} variant(s) were analysed: '
                f'{"; ".join(var_descriptions)}.\n\n'))

        # Reference genome and coordinate mapping
        content.append(('body',
            f'All genomic coordinates were referenced to the {assembly_full} '
            f'human genome assembly. Coding DNA sequence (CDS) positions '
            f'were mapped to genomic coordinates using the Ensembl REST API '
            f'(https://rest.ensembl.org) and cross-validated against NCBI '
            f'Entrez (https://eutils.ncbi.nlm.nih.gov/entrez/). Transcript '
            f'annotations were verified against the MANE (Matched Annotation '
            f'from NCBI and EBI) database to ensure use of clinically '
            f'recommended reference sequences.\n\n'))

        # Sequence retrieval
        content.append(('body',
            f'Genomic sequences were retrieved from the Ensembl REST API, '
            f'extending {flanking} bp on each side of each variant position '
            f'to provide a search space for primer placement.\n\n'))

        # Population variant filtering
        content.append(('h2', 'Population Variant Masking\n'))
        if filter_pop:
            if len(unique_mafs) == 1:
                maf_desc = f'≥ {unique_mafs[0]:.3f}%'
            else:
                maf_desc = (
                    f'≥ {min(unique_mafs):.3f}%–{max(unique_mafs):.3f}% '
                    f'(per-variant thresholds applied; see Per-Variant Details below)'
                )
            content.append(('body',
                f'Known population polymorphisms were obtained from the '
                f'Ensembl Variation API (sourced from gnomAD v4). Positions '
                f'with a minor allele frequency (MAF) {maf_desc} in '
                f'the {pop_name} population were masked with IUPAC ambiguity '
                f'code "N" in the template sequence provided to the primer '
                f'design engine. This masking ensures that primers do not '
                f'bind at polymorphic positions, which could otherwise lead '
                f'to allele-specific amplification bias or failure.\n\n'))
        else:
            content.append(('body',
                'Population variant masking was disabled for this design '
                'session. Primers may bind at polymorphic positions.\n\n'))

        # Primer design engine
        content.append(('h2', 'Primer Design Engine\n'))
        content.append(('body',
            f'Primers were designed using the Primer3 engine (via the '
            f'primer3-py Python interface) with the following key parameters: '
            f'amplicon size {params.min_amplicon_size}–{params.max_amplicon_size} bp, '
            f'primer melting temperature (Tm) {params.primer_min_tm}–{params.primer_max_tm} °C '
            f'(optimal {params.primer_opt_tm} °C), '
            f'primer length {params.primer_min_size}–{params.primer_max_size} nt '
            f'(optimal {params.primer_opt_size} nt), '
            f'GC content {params.primer_min_gc}–{params.primer_max_gc}% '
            f'(optimal {params.primer_opt_gc}%). '))
        content.append(('body',
            f'Thermodynamic calculations used the SantaLucia nearest-neighbour '
            f'model with salt corrections for '
            f'{params.mv_conc} mM monovalent cations, '
            f'{params.dv_conc} mM Mg²⁺, '
            f'{params.dntp_conc} mM dNTPs, and '
            f'{params.dna_conc} nM template DNA. '))
        content.append(('body',
            f'Maximum self-complementarity was set to {params.max_self_complementarity}, '
            f'3\' end complementarity to {params.max_end_complementarity}, '
            f'and pair complementarity to {params.max_pair_complementarity} '
            f'(Primer3 thermodynamic alignment scores).\n\n'))

        # Distance constraints
        constraint_parts = []
        if params.use_variant_distance_constraint:
            constraint_parts.append(
                f'a minimum distance of {params.min_distance_from_variant} bp '
                f'between the nearest primer edge and the variant position'
            )
        if params.use_splice_site_constraint:
            constraint_parts.append(
                f'a minimum distance of {params.min_distance_from_exon_junction} bp '
                f'from exon/intron boundaries on the intronic side to avoid '
                f'splice-site interference'
            )
        if constraint_parts:
            content.append(('body',
                f'Primer placement was constrained to enforce '
                f'{"; and ".join(constraint_parts)}.\n\n'))
        else:
            content.append(('body',
                'No additional distance constraints were applied to primer '
                'placement.\n\n'))

        # Homology discrimination
        content.append(('h2', 'Pseudogene Discrimination\n'))
        if used_homology:
            content.append(('body',
                f'To minimise the risk of pseudogene co-amplification, a '
                f'BLAST+-based homology analysis was performed for each target '
                f'region. Flanking sequences were queried against the human '
                f'genome using BLASTn (word size 11, evalue ≤ 1e-10, minimum '
                f'identity ≥ {self.config.blast.min_percent_identity}%, minimum '
                f'aligned length ≥ {self.config.blast.min_aligned_length} nt) to '
                f'identify homologous regions (pseudogenes, segmental duplications). '))
            if homology_applied > 0:
                content.append(('body',
                    f'Homologous regions were detected for {homology_applied} of '
                    f'{total_designed} variant group(s). Primer pairs were re-ranked '
                    f'to maximise the number of discriminating positions — nucleotides '
                    f'where the primer sequence differs from the pseudogene — '
                    f'preferring mismatches at the 3\' end of the primer for maximum '
                    f'allele-specific discrimination.\n\n'))
            else:
                content.append(('body',
                    f'No significant homologous regions were detected for any of '
                    f'the target loci, so standard Primer3 ranking was used.\n\n'))
        else:
            content.append(('body',
                'Pseudogene discrimination (homology analysis) was not enabled '
                'for this design session.\n\n'))

        # Specificity verification
        content.append(('h2', 'Specificity Verification\n'))
        spec_tools = []
        if used_ucsc:
            spec_tools.append(
                'UCSC In-Silico PCR (https://genome.ucsc.edu/cgi-bin/hgPcr), '
                'which simulates PCR amplification against the full genome '
                f'assembly ({assembly_full}) to detect off-target binding sites'
            )
        if used_blast_check:
            spec_tools.append(
                'NCBI Primer-BLAST (https://www.ncbi.nlm.nih.gov/tools/primer-blast/), '
                'which checks primer specificity against the RefSeq representative '
                'genome database'
            )
        if spec_tools:
            content.append(('body',
                f'Primer specificity was verified using: {"; and ".join(spec_tools)}. '))
            if continuous:
                try:
                    max_att = int(self.max_continuous_var.get())
                except ValueError:
                    max_att = 20
                content.append(('body',
                    f'Continuous design mode was used: for each variant, up to '
                    f'{max_att} candidate primer pairs were generated by Primer3 '
                    f'and evaluated sequentially for specificity until at least '
                    f'{params.min_specific_pairs} verified-specific pair(s) were '
                    f'found per variant, or all candidates were exhausted.\n\n'))
            else:
                content.append(('body',
                    f'{params.num_primer_pairs} primer pair(s) per variant were '
                    f'designed and subsequently checked for specificity.\n\n'))
        else:
            content.append(('body',
                'No post-design specificity verification tools were enabled. '
                'Primer specificity was not experimentally validated in silico.\n\n'))

        # -----------------------------------------------------------------
        # 3. Results summary
        # -----------------------------------------------------------------
        content.append(('h2', 'Results Summary\n'))
        content.append(('body',
            f'Primer design was completed for {successful} of {total_designed} '
            f'variant group(s). '))
        if total_pairs > 0:
            content.append(('body',
                f'A total of {total_pairs} primer pair(s) were generated. '))
        if spec_tools and total_pairs > 0:
            content.append(('body',
                f'Of these, {specific_pairs} pair(s) were confirmed as specific '
                f'(amplifying only the intended target). '))
        if failed > 0:
            content.append(('body',
                f'Primer design failed for {failed} variant group(s), likely due '
                f'to sequence constraints (high polymorphism density, extreme GC '
                f'content, or insufficient flanking region). '))
        content.append(('body', '\n\n'))

        # -----------------------------------------------------------------
        # 4. Per-variant detail table
        # -----------------------------------------------------------------
        if results:
            content.append(('h2', 'Per-Variant Details\n'))
            for r in results:
                if not r.variants:
                    continue
                v = r.variants[0]
                gene = v.gene_symbol
                transcript = v.transcript_accession
                hgvs = v.hgvs_c
                v_maf = per_variant_maf.get(v.row_number, self.default_maf_threshold_pct)

                if r.success and r.amplicons:
                    n_pairs = len(r.amplicons)
                    n_spec = sum(
                        1 for a in r.amplicons
                        if a.primer_pair.specificity_verdict == "SPECIFIC"
                    )
                    amp_sizes = [a.primer_pair.product_size for a in r.amplicons]
                    size_range = f"{min(amp_sizes)}–{max(amp_sizes)}" if len(amp_sizes) > 1 else str(amp_sizes[0])

                    detail = f"• {gene} ({transcript}:{hgvs}): "
                    detail += f"{n_pairs} pair(s), amplicon size {size_range} bp"
                    detail += f", MAF threshold {v_maf:.3f}%"
                    if spec_tools:
                        detail += f", {n_spec} specific"
                    if r.homology_discriminated:
                        detail += f", homology Tier {r.homology_tier}"
                        if r.homology_tier_message:
                            detail += f" ({r.homology_tier_message})"
                    detail += "\n"
                    content.append(('bullet', detail))
                else:
                    content.append(('bullet',
                        f"• {gene} ({transcript}:{hgvs}): "
                        f"MAF threshold {v_maf:.3f}%, design failed — "
                        f"{r.message}\n"))

            content.append(('body', '\n'))

        # -----------------------------------------------------------------
        # 5. Software and databases
        # -----------------------------------------------------------------
        content.append(('h2', 'Software and Databases\n'))
        content.append(('body',
            'Primer design was performed using Primer Designer [Early Beta], '
            'a Python application integrating the following components:\n'))
        content.append(('bullet',
            '• Primer3 (via primer3-py) — thermodynamic primer design engine\n'))
        content.append(('bullet',
            '• Ensembl REST API — genomic sequence retrieval, variant annotation, '
            'coordinate mapping\n'))
        content.append(('bullet',
            '• NCBI Entrez — transcript validation and cross-referencing\n'))
        content.append(('bullet',
            '• MANE database — Matched Annotation from NCBI and EBI '
            'for clinically recommended transcripts\n'))
        if filter_pop:
            content.append(('bullet',
                '• gnomAD v4 (via Ensembl Variation API) — population variant '
                'frequencies for primer-site masking\n'))
        if used_homology:
            content.append(('bullet',
                '• BLAST+ (BLASTn) — local homology analysis for pseudogene '
                'detection and primer discrimination scoring\n'))
        if used_ucsc:
            content.append(('bullet',
                '• UCSC In-Silico PCR — genome-wide specificity simulation\n'))
        if used_blast_check:
            content.append(('bullet',
                '• NCBI Primer-BLAST — RefSeq-based specificity verification\n'))
        content.append(('body', '\n'))

        content.append(('tip',
            'Note: This methodology was auto-generated based on the parameters '
            'and tools used in this session. Verify all details before including '
            'in a publication or report.\n'))

        return content

    def _bind_events(self):
        """Bind window events."""
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _browse_file(self):
        """Open file browser dialog."""
        filetypes = [
            ('Excel files', '*.xlsx'),
            ('CSV files', '*.csv'),
            ('All files', '*.*')
        ]
        filename = filedialog.askopenfilename(filetypes=filetypes)
        if filename:
            self.file_path_var.set(filename)

    def _load_file(self):
        """Load and validate input file."""
        file_path = self.file_path_var.get()
        if not file_path:
            messagebox.showerror("Error", "Please select a file first")
            return

        # Disable load button during loading
        self.load_btn.configure(state=tk.DISABLED)
        self.design_btn.configure(state=tk.DISABLED)

        skip_validation = self.skip_api_validation.get()

        # Set assembly on all API clients BEFORE any API calls
        current_assembly = self.assembly_var.get()
        self.coordinator.set_assembly(current_assembly)
        self.ensembl_client.set_assembly(current_assembly)
        self.variant_db_client.set_assembly(current_assembly)
        self.ncbi_client.set_assembly(current_assembly)

        self._set_status("Loading file...")
        self._set_progress(0)

        # Run in thread
        def load_thread():
            try:
                # Parse file
                self._set_status("Parsing file...")
                result = self.file_parser.parse_file(file_path)

                self.log_info(f"Parsed {len(result.variants)} variants with {len(result.errors)} errors")

                if result.errors:
                    self.root.after(0, lambda: self._show_parse_errors(result.errors))

                if not result.variants:
                    error_msg = "No valid variants found in file."
                    if result.errors:
                        error_msg += f"\n\nErrors encountered:\n" + "\n".join(result.errors[:3])
                    self.root.after(0, lambda: messagebox.showerror("Error", error_msg))
                    self._set_status("No variants loaded")
                    self._set_progress(0)
                    return

                self.project_state.variants = result.variants
                self._set_progress(20)

                if skip_validation:
                    # Skip API validation - just basic syntax check
                    self._set_status("Performing basic validation (offline mode)...")

                    for i, variant in enumerate(result.variants):
                        if self._is_cancel_requested():
                            break

                        # Only do basic HGVS syntax validation (no API calls)
                        hgvs_result = self.validator.hgvs_validator.validate_hgvs_syntax(variant.hgvs_c)
                        variant.add_validation_result(hgvs_result)

                        transcript_result = self.validator.hgvs_validator.validate_transcript_format(
                            variant.transcript_accession
                        )
                        variant.add_validation_result(transcript_result)

                        progress = 20 + (70 * (i + 1) / len(result.variants))
                        self._set_progress(progress)

                    self._set_status(f"Loaded {len(result.variants)} variants (offline mode)")
                else:
                    # Enable auto-download for MANE data when doing full validation
                    if self.mane_manager and not self.mane_manager._mane_data:
                        self.mane_manager.auto_download = True
                    # Full validation with API calls
                    self._set_status("Validating variants online (checking MANE transcripts, gene associations)...")

                    # --- Parallel validation ---
                    validation_errors = []
                    validation_errors_lock = threading.Lock()
                    val_completed = [0]
                    val_completed_lock = threading.Lock()
                    total_v = len(result.variants)
                    max_workers = min(self.config.max_concurrent_requests, total_v, 5)

                    def _validate_one(variant):
                        if self._is_cancel_requested():
                            return
                        try:
                            self.validator.validate_variant(variant, check_sequence=False)
                        except Exception as e:
                            error_msg = f"Validation error for {variant.gene_symbol}: {str(e)}"
                            self.log_error(error_msg)
                            with validation_errors_lock:
                                validation_errors.append(error_msg)
                        with val_completed_lock:
                            val_completed[0] += 1
                            count = val_completed[0]
                        self._set_progress(20 + (40 * count / total_v))
                        self._set_status(f"Validating variants: {count}/{total_v}")

                    with ThreadPoolExecutor(max_workers=max_workers) as executor:
                        futures = {executor.submit(_validate_one, v): v for v in result.variants}
                        for future in as_completed(futures):
                            if self._is_cancel_requested():
                                executor.shutdown(wait=False, cancel_futures=True)
                                break
                            try:
                                future.result()
                            except Exception as e:
                                v = futures[future]
                                error_msg = f"Validation error for {v.gene_symbol}: {str(e)}"
                                self.log_error(error_msg)
                                with validation_errors_lock:
                                    validation_errors.append(error_msg)

                    # --- Parallel coordinate mapping ---
                    self._set_status("Mapping genomic coordinates...")

                    mapping_errors = []
                    mapping_errors_lock = threading.Lock()
                    map_completed = [0]
                    map_completed_lock = threading.Lock()

                    def _map_one(variant):
                        if self._is_cancel_requested():
                            return None, variant
                        mapping = self.coordinator.map_cds_to_genomic(variant)
                        with map_completed_lock:
                            map_completed[0] += 1
                            count = map_completed[0]
                        self._set_progress(60 + (30 * count / total_v))
                        self._set_status(f"Mapping coordinates: {count}/{total_v}")
                        return mapping, variant

                    with ThreadPoolExecutor(max_workers=max_workers) as executor:
                        futures = {executor.submit(_map_one, v): v for v in result.variants}
                        for future in as_completed(futures):
                            if self._is_cancel_requested():
                                executor.shutdown(wait=False, cancel_futures=True)
                                break
                            try:
                                mapping, variant = future.result()
                                if mapping is None:
                                    continue
                                if mapping.success:
                                    variant.genomic_position = mapping.genomic_position
                                    variant.exon = mapping.exon
                                else:
                                    with mapping_errors_lock:
                                        mapping_errors.append(f"{variant.gene_symbol}: {mapping.error}")
                            except Exception as e:
                                v = futures[future]
                                error_msg = f"Mapping error for {v.gene_symbol}: {str(e)}"
                                self.log_error(error_msg)
                                with mapping_errors_lock:
                                    mapping_errors.append(error_msg)

                    # Show any validation or mapping errors
                    if validation_errors or mapping_errors:
                        all_errors = validation_errors + mapping_errors
                        self.root.after(0, lambda: self._show_validation_warnings(all_errors))

                self.project_state.is_validated = True
                self._loaded_with_assembly = current_assembly

                # Update UI
                self.root.after(0, self._update_variants_table)
                self.root.after(0, lambda: self.design_btn.configure(state=tk.NORMAL))
                self.root.after(0, lambda: self.load_btn.configure(state=tk.NORMAL))
                self.root.after(0, lambda: self._update_help_context('file_loaded'))

                status_msg = f"Loaded {len(result.variants)} variants — caching sequences, please wait..."
                if skip_validation:
                    status_msg = f"Loaded {len(result.variants)} variants (offline) — caching sequences, please wait..."
                self._set_status(status_msg)
                # Don't set to 100% yet - caching will continue the progress
                self._set_progress(90)

                # Start background caching of sequence data for fast display
                self.root.after(100, self._start_sequence_cache)

            except Exception as e:
                self.log_exception(f"File loading error: {e}")
                error_msg = f"Error loading file: {str(e)}\n\nPlease check:\n"
                error_msg += "- File format is correct (CSV or XLSX)\n"
                error_msg += "- Required columns exist (Gene, Transcript, Position)\n"
                if not skip_validation:
                    error_msg += "- Internet connection is available\n"
                    error_msg += "\nTip: Try enabling 'Skip online validation' for offline use"
                self.root.after(0, lambda: messagebox.showerror("Error", error_msg))
                self.root.after(0, lambda: self.load_btn.configure(state=tk.NORMAL))
                self._set_status("Error loading file")
                self._set_progress(0)
            finally:
                # Always re-enable load button
                self.root.after(0, lambda: self.load_btn.configure(state=tk.NORMAL))

        self._run_in_thread(load_thread)

    def _show_parse_errors(self, errors):
        """Show parsing errors in a dialog."""
        error_text = "\n".join(errors[:10])
        if len(errors) > 10:
            error_text += f"\n... and {len(errors) - 10} more errors"
        messagebox.showwarning("Parsing Warnings", error_text)

    def _show_validation_warnings(self, errors):
        """Show validation/mapping warnings in a dialog."""
        error_text = "Some issues were encountered during validation:\n\n"
        error_text += "\n".join(errors[:10])
        if len(errors) > 10:
            error_text += f"\n... and {len(errors) - 10} more errors"
        error_text += "\n\nYou can still proceed with primer design, but results may be affected."
        messagebox.showwarning("Validation Warnings", error_text)

    def _update_variants_table(self):
        """Update the variants table with current data."""
        # Clear existing items
        for item in self.variants_tree.get_children():
            self.variants_tree.delete(item)

        # Clear MAF thresholds dictionary
        self.variant_maf_thresholds = {}

        # Add variants with sequential numbering starting from 1
        for idx, variant in enumerate(self.project_state.variants, start=1):
            status = "Valid" if variant.is_valid else "Error"
            if variant.has_warnings:
                status = "Warning"

            # Default MAF threshold (0.001%) - format with enough decimal places
            if self.default_maf_threshold_pct < 0.01:
                default_maf = f"{self.default_maf_threshold_pct:.4f}%"
            elif self.default_maf_threshold_pct < 1.0:
                default_maf = f"{self.default_maf_threshold_pct:.3f}%"
            else:
                default_maf = f"{self.default_maf_threshold_pct:.2f}%"

            # Initialize MAF threshold for this variant
            self.variant_maf_thresholds[variant.row_number] = self.default_maf_threshold_pct

            # Store the original row_number as iid for lookup, display sequential number
            self.variants_tree.insert('', 'end', iid=str(variant.row_number), values=(
                idx,  # Sequential number starting from 1
                variant.gene_symbol,
                variant.transcript_accession,
                variant.hgvs_c,
                default_maf,  # MAF% column with default value
                status
            ))

    def _on_variant_select(self, event):
        """Handle variant selection in table."""
        selection = self.variants_tree.selection()
        if not selection:
            return

        # Get selected variant using iid (which stores original row_number)
        row_num = int(selection[0])

        variant = next((v for v in self.project_state.variants if v.row_number == row_num), None)
        if variant:
            self._selected_variant = variant
            self._show_variant_details(variant)

    def _on_population_change(self):
        """Handle population checkbox change - refresh sequence display."""
        if self._selected_variant:
            self._show_variant_details(self._selected_variant)

    def _on_maf_range_change(self, event=None):
        """Handle MAF range (min/max) change - update slider range. Values are PERCENTAGES."""
        try:
            new_min = float(self.maf_range_min_var.get())
            new_max = float(self.maf_range_max_var.get())

            # Validate range (percentages: 0-100)
            if new_min < 0:
                new_min = 0.0
                self.maf_range_min_var.set("0")
            if new_max > 100.0:
                new_max = 100.0
                self.maf_range_max_var.set("100")
            if new_min >= new_max:
                new_min = 0.0
                new_max = 50.0
                self.maf_range_min_var.set("0")
                self.maf_range_max_var.set("50")

            # Calculate appropriate resolution based on range (in percentage)
            range_size = new_max - new_min
            if range_size <= 0.01:
                resolution = 0.0001   # Very fine for tiny ranges
            elif range_size <= 0.1:
                resolution = 0.001
            elif range_size <= 1.0:
                resolution = 0.001
            elif range_size <= 10.0:
                resolution = 0.01
            else:
                resolution = 0.001    # Default fine resolution

            # Update slider range and resolution
            self.maf_slider.configure(from_=new_min, to=new_max, resolution=resolution)

            # Clamp current value to new range
            current_val = self.maf_slider_var.get()
            if current_val < new_min:
                self.maf_slider_var.set(new_min)
            elif current_val > new_max:
                self.maf_slider_var.set(new_max)

            # Update display
            self._on_maf_slider_change(self.maf_slider_var.get())
        except ValueError:
            pass

    def _on_maf_slider_change(self, value):
        """Handle MAF slider change - update entry, auto-apply to selected variant, refresh sequence."""
        try:
            maf_pct = float(value)
            # Update entry field with percentage formatting
            if maf_pct == 0:
                formatted = "0"
            elif maf_pct < 0.001:
                formatted = f"{maf_pct:.4f}"
            elif maf_pct < 0.01:
                formatted = f"{maf_pct:.3f}"
            elif maf_pct < 1.0:
                formatted = f"{maf_pct:.3f}"
            else:
                formatted = f"{maf_pct:.2f}"

            # Update the entry without triggering recursion
            self.maf_slider_entry.delete(0, tk.END)
            self.maf_slider_entry.insert(0, formatted)

            # Auto-apply MAF to the currently selected variant
            if self._selected_variant:
                row_num = self._selected_variant.row_number
                self.variant_maf_thresholds[row_num] = maf_pct
                self._update_variant_maf_in_table(row_num, maf_pct)
                self._show_variant_details(self._selected_variant)
        except ValueError:
            pass

    def _update_variant_maf_in_table(self, row_num: int, maf_pct: float):
        """Update the MAF% column for a specific variant in the table."""
        try:
            # Get current values
            item = self.variants_tree.item(str(row_num))
            values = list(item['values'])

            # Format MAF percentage with enough decimal places for small values
            if maf_pct == 0:
                maf_display = "0%"
            elif maf_pct < 0.0001:
                maf_display = f"{maf_pct:.6f}%"
            elif maf_pct < 0.01:
                maf_display = f"{maf_pct:.4f}%"
            elif maf_pct < 1.0:
                maf_display = f"{maf_pct:.3f}%"
            else:
                maf_display = f"{maf_pct:.2f}%"

            # Update MAF column (index 4)
            values[4] = maf_display

            # Update the item
            self.variants_tree.item(str(row_num), values=values)
        except Exception as e:
            self.log_warning(f"Could not update MAF in table: {e}")

    def _on_maf_entry_submit(self, event=None):
        """Handle MAF entry submit - sync with slider and auto-apply to selected variant."""
        try:
            maf_value = float(self.maf_slider_entry.get())
            # Get current slider range
            slider_min = float(self.maf_range_min_var.get())
            slider_max = float(self.maf_range_max_var.get())
            # Clamp to slider range
            maf_value = max(slider_min, min(slider_max, maf_value))
            self.maf_slider_var.set(maf_value)
            # Auto-apply and refresh display
            if self._selected_variant:
                row_num = self._selected_variant.row_number
                self.variant_maf_thresholds[row_num] = maf_value
                self._update_variant_maf_in_table(row_num, maf_value)
                self._show_variant_details(self._selected_variant)
        except ValueError:
            pass

    def _start_sequence_cache(self):
        """Start background thread to pre-cache sequence data for all variants."""
        if self._cache_thread and self._cache_thread.is_alive():
            return  # Already caching

        self._sequence_cache.clear()
        self._start_spinner()  # Start spinner for caching
        self._update_help_context('caching')

        def cache_thread():
            try:
                try:
                    flank_size = int(self.flank_size_var.get())
                except ValueError:
                    flank_size = 250

                total = len(self.project_state.variants)
                if total == 0:
                    self._set_progress(100)
                    self._set_status("Ready")
                    return

                variants_to_cache = [
                    v for v in self.project_state.variants if v.genomic_position
                ]
                if not variants_to_cache:
                    self._set_progress(100)
                    self._set_status("Ready")
                    return

                cache_completed = [0]
                cache_completed_lock = threading.Lock()
                max_workers = min(self.config.max_concurrent_requests, len(variants_to_cache), 5)

                def _cache_one(variant):
                    if self._is_cancel_requested():
                        return
                    row_num = variant.row_number
                    seq_start = variant.genomic_position.start - flank_size
                    seq_end = variant.genomic_position.end + flank_size
                    try:
                        sequence = self.ensembl_client.get_genomic_sequence(
                            variant.genomic_position.chromosome,
                            seq_start,
                            seq_end
                        )
                        pop_variants = []
                        try:
                            pop_variants = self.variant_db_client.get_variants_in_region(
                                variant.genomic_position.chromosome,
                                seq_start,
                                seq_end,
                                populations=['global'],
                                maf_threshold=0.0
                            )
                        except Exception:
                            pass
                        if sequence:
                            self._sequence_cache[row_num] = {
                                'sequence': sequence,
                                'pop_variants': pop_variants,
                                'seq_start': seq_start,
                                'seq_end': seq_end
                            }
                    except Exception as e:
                        self.log_warning(f"Could not cache sequence for variant {row_num}: {e}")

                    with cache_completed_lock:
                        cache_completed[0] += 1
                        count = cache_completed[0]
                    self._set_progress(90 + (10 * count / total))
                    self._set_status(f"Caching sequences: {count}/{total}")

                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = [executor.submit(_cache_one, v) for v in variants_to_cache]
                    for future in as_completed(futures):
                        if self._is_cancel_requested():
                            executor.shutdown(wait=False, cancel_futures=True)
                            break
                        try:
                            future.result()
                        except Exception:
                            pass

                self._set_progress(100)
                self._set_status(f"Ready - {len(self._sequence_cache)} sequences cached")
                self.root.after(0, lambda: self._update_help_context('ready'))
            finally:
                self.root.after(0, self._stop_spinner)  # Always stop spinner

        self._cache_thread = threading.Thread(target=cache_thread, daemon=True)
        self._cache_thread.start()

    def _show_variant_details(self, variant: Variant):
        """Show details for selected variant."""
        if not variant.genomic_position:
            return

        # Get flanking size from user input
        try:
            flank_size = int(self.flank_size_var.get())
        except ValueError:
            flank_size = 250

        seq_start = variant.genomic_position.start - flank_size
        seq_end = variant.genomic_position.end + flank_size

        # Try to use cached data first
        row_num = variant.row_number
        cached = self._sequence_cache.get(row_num)

        if cached and cached['seq_start'] == seq_start and cached['seq_end'] == seq_end:
            # Use cached data
            sequence = cached['sequence']
            pop_variants = cached['pop_variants']
        else:
            # Fetch fresh data
            sequence = self.ensembl_client.get_genomic_sequence(
                variant.genomic_position.chromosome,
                seq_start,
                seq_end
            )
            if not sequence:
                return

            # Get population variants - always fetch ALL variants for display
            pop_variants = []
            try:
                pop_variants = self.variant_db_client.get_variants_in_region(
                    variant.genomic_position.chromosome,
                    seq_start,
                    seq_end,
                    populations=['global'],  # Cache with global, will be filtered on display
                    maf_threshold=0.0  # Get ALL variants for visualization
                )
            except Exception as e:
                self.log_warning(f"Could not fetch population variants: {e}")

            # Cache for future use
            self._sequence_cache[row_num] = {
                'sequence': sequence,
                'pop_variants': pop_variants,
                'seq_start': seq_start,
                'seq_end': seq_end
            }

        if not sequence:
            return

        # Get selected population (single selection from radio buttons)
        selected_pops = [self.selected_population.get()]

        # Calculate variant position relative to sequence start
        variant_pos_in_seq = variant.genomic_position.start - seq_start

        # Calculate the length of the variant
        variant_length = max(1, variant.genomic_position.end - variant.genomic_position.start + 1)

        # Get transcript info for marking CDS region
        transcript_info = None
        if variant.transcript:
            transcript_info = {
                'gene_symbol': variant.transcript.gene_symbol,
                'accession': variant.transcript.full_accession,
                'cds_start': variant.transcript.cds_start,
                'cds_end': variant.transcript.cds_end
            }

        # Build exon regions in sequence coordinates for visualization
        # Each tuple: (seq_start_pos, seq_end_pos, exon_number)
        exon_regions = []
        if variant.transcript and variant.transcript.exons:
            for exon in variant.transcript.exons:
                exon_seq_start = exon.genomic_start - seq_start
                exon_seq_end = exon.genomic_end - seq_start
                # Only include exons that overlap with our displayed region
                if exon_seq_end >= 0 and exon_seq_start < len(sequence):
                    exon_seq_start = max(0, exon_seq_start)
                    exon_seq_end = min(len(sequence) - 1, exon_seq_end)
                    exon_regions.append((exon_seq_start, exon_seq_end, exon.number))
        elif variant.exon:
            # Fallback: only the variant's own exon
            exon_seq_start = max(0, variant.exon.genomic_start - seq_start)
            exon_seq_end = min(len(sequence) - 1, variant.exon.genomic_end - seq_start)
            exon_regions.append((exon_seq_start, exon_seq_end, variant.exon.number))

        # Get MAF display threshold from slider (PERCENTAGE -> convert to fraction)
        maf_display_threshold_pct = self.maf_slider_var.get()
        maf_display_threshold = maf_display_threshold_pct / 100.0  # Convert % to fraction

        self.seq_visualizer.display_sequence(
            sequence,
            variant_position=variant_pos_in_seq,
            variant_length=variant_length,
            population_variants=pop_variants,
            seq_start_genomic=seq_start,
            selected_populations=selected_pops,
            transcript_info=transcript_info,
            maf_display_threshold=maf_display_threshold,
            exon_regions=exon_regions
        )

    def _generate_masked_fasta(self, variant: Variant, maf_threshold_pct: float) -> tuple:
        """
        Generate FASTA sequence with N at positions of population variants above MAF threshold.
        Uses cached sequence data when available for speed.

        Args:
            variant: The variant to generate sequence for
            maf_threshold_pct: MAF threshold in percentage (e.g., 1.0 for 1%)

        Returns:
            Tuple of (fasta_path, masked_sequence, masked_positions, pop_variants) or (None, None, None, None) on failure
        """
        if not variant.genomic_position:
            return None, None, None, None

        try:
            flank_size = int(self.flank_size_var.get())
        except ValueError:
            flank_size = 250

        seq_start = variant.genomic_position.start - flank_size
        seq_end = variant.genomic_position.end + flank_size

        # Try to use cached data first (from _start_sequence_cache)
        row_num = variant.row_number
        cached = self._sequence_cache.get(row_num)

        if cached and cached['seq_start'] == seq_start and cached['seq_end'] == seq_end:
            # Use cached data - much faster!
            sequence = cached['sequence']
            pop_variants = cached['pop_variants']
        else:
            # Fetch fresh data (fallback)
            sequence = self.ensembl_client.get_genomic_sequence(
                variant.genomic_position.chromosome,
                seq_start,
                seq_end
            )

            if not sequence:
                return None, None, None, None

            # Get selected population
            selected_pop = self.selected_population.get()

            # Get population variants in the region
            pop_variants = []
            try:
                pop_variants = self.variant_db_client.get_variants_in_region(
                    variant.genomic_position.chromosome,
                    seq_start,
                    seq_end,
                    populations=[selected_pop],
                    maf_threshold=0.0  # Get all variants, we'll filter ourselves
                )
            except Exception as e:
                self.log_warning(f"Could not fetch population variants for masking: {e}")

            # Cache for future use
            self._sequence_cache[row_num] = {
                'sequence': sequence,
                'pop_variants': pop_variants,
                'seq_start': seq_start,
                'seq_end': seq_end
            }

        # Get selected population for filtering
        selected_pop = self.selected_population.get()

        # Convert sequence to list for modification
        seq_list = list(sequence.upper())
        masked_positions = set()
        maf_threshold_fraction = maf_threshold_pct / 100.0

        # Mask positions where population variants have MAF >= threshold
        for pv in pop_variants:
            # Get the MAF to use (population-specific or global)
            var_maf = pv.maf_global
            if selected_pop != 'global' and pv.maf_by_population:
                var_maf = max(var_maf, pv.maf_by_population.get(selected_pop, 0))

            # If MAF >= threshold, mask this position with N
            if var_maf >= maf_threshold_fraction:
                seq_pos = pv.position - seq_start
                if 0 <= seq_pos < len(seq_list):
                    # Only mask SNPs (single nucleotide positions)
                    if pv.alt and len(pv.alt) == 1 and len(pv.ref) == 1:
                        seq_list[seq_pos] = 'N'
                        masked_positions.add(seq_pos)

        masked_sequence = ''.join(seq_list)

        # Create FASTA directory if needed
        # data_cache_dir is like 'primer_designer/data/cache', go up one level for 'data'
        data_dir = Path(self.config.data_cache_dir).parent
        fasta_dir = data_dir / "fasta_sequences"
        fasta_dir.mkdir(parents=True, exist_ok=True)

        # Generate FASTA file
        safe_gene = variant.gene_symbol.replace('/', '_').replace('\\', '_')
        safe_hgvs = variant.hgvs_c.replace('>', '_').replace('<', '_').replace(':', '_')
        fasta_filename = f"{safe_gene}_{safe_hgvs}_MAF{maf_threshold_pct:.2f}pct.fasta"
        fasta_path = fasta_dir / fasta_filename

        # Write FASTA
        header = f">{variant.gene_symbol}|{variant.transcript_accession}|{variant.hgvs_c}|MAF_threshold={maf_threshold_pct}%|chr{variant.genomic_position.chromosome}:{seq_start}-{seq_end}|masked_positions={len(masked_positions)}"

        with open(fasta_path, 'w') as f:
            f.write(f"{header}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(masked_sequence), 60):
                f.write(f"{masked_sequence[i:i+60]}\n")

        self.log_info(f"Generated masked FASTA: {fasta_path} ({len(masked_positions)} positions masked)")

        return str(fasta_path), masked_sequence, masked_positions, pop_variants

    def _design_primers(self):
        """Start primer design process."""
        if not self.project_state.is_validated:
            messagebox.showerror("Error", "Please load and validate variants first")
            return

        self._set_status("Designing primers...")
        self._set_progress(0)
        self._set_cancel_requested(False)
        self._update_help_context('designing')

        # Update parameters
        self._update_design_parameters()

        def design_thread():
            try:
                # Initialize primer designer if needed
                if not self.primer_designer:
                    self.primer_designer = PrimerDesigner(
                        self.ensembl_client,
                        self.variant_db_client,
                        self.design_parameters
                    )

                # STEP 1: Generate masked FASTA sequences for each variant (parallel)
                self._set_status("Generating masked FASTA sequences...")
                fasta_files = {}  # row_number -> (fasta_path, masked_sequence, masked_positions)
                fasta_files_lock = threading.Lock()

                # Capture tkinter vars BEFORE spawning threads
                current_slider_maf = self.maf_slider_var.get()
                maf_thresholds_snapshot = dict(self.variant_maf_thresholds)

                total_variants = len(self.project_state.variants)
                fasta_completed = [0]
                fasta_completed_lock = threading.Lock()
                fasta_max_workers = min(self.config.max_concurrent_requests, total_variants, 5)

                def _generate_one_fasta(variant):
                    if self._is_cancel_requested():
                        return
                    maf_threshold_pct = maf_thresholds_snapshot.get(
                        variant.row_number, current_slider_maf
                    )
                    result = self._generate_masked_fasta(variant, maf_threshold_pct)
                    if result[0]:
                        with fasta_files_lock:
                            fasta_files[variant.row_number] = result[:3]
                        variant._masked_positions = result[2]
                        variant._masked_sequence = result[1]
                        variant._pop_variants = result[3]

                    with fasta_completed_lock:
                        fasta_completed[0] += 1
                        count = fasta_completed[0]
                    self._set_progress(5 * count / total_variants)
                    self._set_status(f"Generating FASTA: {count}/{total_variants}")

                with ThreadPoolExecutor(max_workers=fasta_max_workers) as executor:
                    futures = [executor.submit(_generate_one_fasta, v)
                               for v in self.project_state.variants]
                    for future in as_completed(futures):
                        if self._is_cancel_requested():
                            executor.shutdown(wait=False, cancel_futures=True)
                            break
                        try:
                            future.result()
                        except Exception:
                            pass

                self._set_progress(5)

                # STEP 2: Group variants
                self._set_status("Grouping variants...")
                groups = self.grouper.group_variants(self.project_state.variants)

                # Restore original input order (grouper sorts by genomic position)
                groups.sort(
                    key=lambda g: min(v.row_number for v in g.variants) if g.variants else 0
                )

                self._set_progress(10)

                # STEP 2.5: Auto-run homology analysis if checkbox enabled
                homology_ran = False
                if self.design_parameters.use_homology_discrimination:
                    self._auto_run_homology_for_groups(groups)
                    homology_ran = True

                design_progress_base = 20 if homology_ran else 10
                design_progress_range = 70 if homology_ran else 80

                # STEP 3: Design primers
                results = []
                # qPCR mode hidden in GUI for now — always defaults to PCR
                mode = DesignMode.QPCR if self.mode_var.get() == "qPCR" else DesignMode.PCR

                any_specificity = (
                    self.use_primer_blast_check.get()
                    or self.use_ucsc_ispcr_check.get()
                )
                continuous_mode = self.continuous_mode_var.get() and any_specificity

                if continuous_mode:
                    # Continuous mode: design + check + display iteratively
                    self._design_continuous(
                        groups, mode, results, fasta_files, current_slider_maf
                    )
                else:
                    # Original flow: design all, then check all, then display all
                    total_groups = len(groups)
                    for i, group in enumerate(groups):
                        if self._is_cancel_requested():
                            break

                        self._set_status(f"Designing primers for group {i + 1}/{total_groups}...")

                        if group.variants:
                            first_variant = group.variants[0]
                            maf_pct = self.variant_maf_thresholds.get(
                                first_variant.row_number, current_slider_maf
                            )
                            self.design_parameters.maf_threshold = maf_pct / 100.0
                            self.log_info(f"Setting MAF threshold for group to {maf_pct}% ({self.design_parameters.maf_threshold} fraction)")

                        # Get homology result for group's first variant (if available)
                        group_homology = None
                        if group.variants:
                            group_homology = self._homology_results_by_variant.get(
                                group.variants[0].row_number
                            )

                        result = self.primer_designer.design_primers(
                            group, mode,
                            progress_callback=lambda p, m: self._set_status(m),
                            homology_result=group_homology,
                        )
                        results.append(result)

                        if total_groups > 0:
                            progress = design_progress_base + (design_progress_range * (i + 1) / total_groups)
                            self._set_progress(progress)

                    # Check specificity only if user enabled at least one tool
                    if any_specificity:
                        self._set_status("Checking primer specificity...")
                        self._check_specificity(results)

                    self.project_state.design_results = results
                    self.project_state.is_designed = True
                    self.project_state.fasta_files = fasta_files

                    self.root.after(0, lambda: self._display_results(results))
                    self.root.after(0, lambda: self.export_btn.configure(state=tk.NORMAL))
                    self.root.after(0, lambda: self._update_help_context('results'))

                    successful = sum(1 for r in results if r.success)
                    fasta_count = len(fasta_files)
                    self._set_status(f"Design complete: {successful}/{len(results)} successful, {fasta_count} FASTA files generated")
                    self._set_progress(100)

            except Exception as e:
                self.log_exception(f"Design error: {e}")
                self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
                self._set_status("Error during design")

        self._run_in_thread(design_thread)

    def _update_design_parameters(self):
        """Update design parameters from UI."""
        try:
            self.design_parameters.min_amplicon_size = int(self.amp_min_var.get())
            self.design_parameters.max_amplicon_size = int(self.amp_max_var.get())
            self.design_parameters.primer_min_tm = float(self.tm_min_var.get())
            self.design_parameters.primer_max_tm = float(self.tm_max_var.get())
            self.design_parameters.min_distance_from_variant = int(self.dist_var.get())
            self.design_parameters.use_variant_distance_constraint = self.use_dist_constraint.get()
            self.design_parameters.min_distance_from_exon_junction = int(self.intron_dist_var.get())
            self.design_parameters.use_splice_site_constraint = self.use_intron_constraint.get()
            self.design_parameters.reference_assembly = self.assembly_var.get()

            # Number of primer pairs to return
            self.design_parameters.num_primer_pairs = int(self.num_primer_pairs_var.get())

            # Minimum specific pairs for continuous mode
            try:
                self.design_parameters.min_specific_pairs = int(self.min_specific_pairs_var.get())
            except (ValueError, AttributeError):
                self.design_parameters.min_specific_pairs = 3

            # Population variant filtering
            self.design_parameters.filter_population_variants = self.filter_pop_variants.get()

            # Selected population (single selection)
            self.design_parameters.selected_populations = [self.selected_population.get()]

            # Homology discrimination
            self.design_parameters.use_homology_discrimination = self.use_homology_discrimination.get()

            # Update coordinator and all API clients assembly
            self.coordinator.set_assembly(self.assembly_var.get())
            self.ensembl_client.set_assembly(self.assembly_var.get())
            self.variant_db_client.set_assembly(self.assembly_var.get())
            self.ncbi_client.set_assembly(self.assembly_var.get())

        except ValueError as e:
            messagebox.showerror("Error", f"Invalid parameter value: {e}")

    # ------------------------------------------------------------------
    #  Continuous design mode
    # ------------------------------------------------------------------

    def _design_continuous(self, groups, mode, results, fasta_files, current_slider_maf):
        """
        Continuous design mode: for each variant group, ask Primer3 for a large
        pool of candidates up front, then iterate through them one-by-one checking
        specificity, displaying results live, until a SPECIFIC pair is found
        or max_attempts is reached.  Runs inside the design_thread (background).

        Primer3 is deterministic — given the same input it always returns the
        same ranked list.  So we request max_attempts pairs in one call and
        walk through them sequentially.
        """
        try:
            max_attempts = int(self.max_continuous_var.get())
        except ValueError:
            max_attempts = 20

        use_primer_blast = self.use_primer_blast_check.get()
        use_ucsc = self.use_ucsc_ispcr_check.get()

        # Ensure orchestrator
        from ..primer.specificity_orchestrator import SpecificityOrchestrator
        if not self.specificity_orchestrator:
            self.specificity_orchestrator = SpecificityOrchestrator()

        # Clear results panel for live display
        self.root.after(0, lambda: self.results_text.delete('1.0', tk.END))

        # Save original num_primer_pairs to restore later
        original_num_pairs = self.design_parameters.num_primer_pairs

        total_groups = len(groups)

        for group_idx, group in enumerate(groups):
            if self._is_cancel_requested():
                break

            self._reset_skip_variant()
            self.root.after(0, lambda: self.skip_variant_btn.configure(state=tk.NORMAL))

            gene = group.variants[0].gene_symbol if group.variants else '?'
            transcript = group.variants[0].transcript_accession if group.variants else '?'
            hgvs = group.variants[0].hgvs_c if group.variants else '?'

            # Set MAF for this group
            if group.variants:
                first_variant = group.variants[0]
                maf_pct = self.variant_maf_thresholds.get(
                    first_variant.row_number, current_slider_maf
                )
                self.design_parameters.maf_threshold = maf_pct / 100.0

            # Insert variant header live
            header = f"{group_idx + 1}.{gene}; {transcript}; {hgvs}; Continuous mode\n"
            self.root.after(0, lambda t=header: self._append_live_text(t))

            # Get expected coordinates
            exp_chr, exp_start, exp_end = None, None, None
            if group.variants and group.variants[0].genomic_position:
                gp = group.variants[0].genomic_position
                exp_chr = gp.chromosome
                exp_start = gp.start
                exp_end = gp.end

            specific_amplicons = []
            non_specific_amplicons = []
            found_specific = False
            pair_counter = 0
            min_specific = self.design_parameters.min_specific_pairs

            # Request a large pool of candidates from Primer3 in one call
            self.design_parameters.num_primer_pairs = max_attempts
            self._set_status(
                f"[{group_idx+1}/{total_groups}] {gene}: "
                f"Designing up to {max_attempts} candidate pairs..."
            )

            # Get homology result for group's first variant (if available)
            cont_homology = None
            if group.variants:
                cont_homology = self._homology_results_by_variant.get(
                    group.variants[0].row_number
                )

            design_result = self.primer_designer.design_primers(
                group, mode,
                progress_callback=lambda p, m: self._set_status(m),
                homology_result=cont_homology,
            )

            if not design_result.success or not design_result.amplicons:
                msg = f"  Primer3 returned no pairs — design space exhausted.\n"
                self.root.after(0, lambda t=msg: self._append_live_text(t, 'cont_exhausted'))
            else:
                total_candidates = len(design_result.amplicons)
                tier_info = ""
                if design_result.homology_tier > 0:
                    tier_info = f" [homology Tier {design_result.homology_tier}]"
                info_msg = f"  Primer3 returned {total_candidates} candidate pair(s){tier_info}. Checking specificity...\n"
                self.root.after(0, lambda t=info_msg: self._append_live_text(t))

                # Show homologous region details if applicable
                if cont_homology and design_result.homology_tier > 0:
                    homology_detail = self._format_homology_details(cont_homology)
                    if homology_detail:
                        self.root.after(0, lambda t=homology_detail: self._append_live_text(t, 'warning'))

                # Deduplicate and iterate through candidates
                seen_sequences = set()

                for amplicon in design_result.amplicons:
                    if (self._is_cancel_requested()
                            or self._is_skip_variant_requested()
                            or len(specific_amplicons) >= min_specific):
                        break

                    pp = amplicon.primer_pair
                    pair_key = (pp.forward.sequence, pp.reverse.sequence)

                    if pair_key in seen_sequences:
                        continue
                    seen_sequences.add(pair_key)

                    pair_counter += 1

                    self._set_status(
                        f"[{group_idx+1}/{total_groups}] {gene}: "
                        f"Checking specificity for pair {pair_counter}/{total_candidates} "
                        f"(unique #{len(seen_sequences)})..."
                    )

                    combined = self.specificity_orchestrator.check_all(
                        pp,
                        expected_chromosome=exp_chr,
                        expected_start=exp_start,
                        expected_end=exp_end,
                        assembly=self.assembly_var.get(),
                        progress_callback=lambda msg: self._set_status(msg),
                        use_primer_blast=use_primer_blast,
                        use_ucsc=use_ucsc,
                    )

                    # Map result to PrimerPair
                    pp.is_specific = combined.is_specific
                    pp.specificity_warnings = combined.all_warnings
                    pp.specificity_verdict = combined.verdict.value
                    pp.pseudogene_risk = bool(combined.pseudogene_warnings)
                    pp.off_target_details = combined.off_target_summary

                    parts = []
                    for tr in combined.tool_results:
                        if tr.error == "Skipped by user":
                            continue
                        if tr.ran_successfully:
                            status = "OK" if tr.is_specific else f"{tr.off_target_count} off-target"
                        elif tr.available:
                            status = "ERROR"
                        else:
                            status = "N/A"
                        parts.append(f"{tr.tool_name}: {status}")
                    pp.specificity_tool_summary = " | ".join(parts)

                    # Display live
                    self.root.after(0,
                        lambda a=amplicon, n=pair_counter:
                            self._display_single_pair_live(a, n))

                    if combined.is_specific:
                        specific_amplicons.append(amplicon)
                        if len(specific_amplicons) >= min_specific:
                            found_specific = True
                    else:
                        non_specific_amplicons.append(amplicon)

            # End of iteration for this variant
            self.root.after(0, lambda: self.skip_variant_btn.configure(state=tk.DISABLED))

            # Summary line
            if specific_amplicons:
                n_spec = len(specific_amplicons)
                summary = f"  >>> {n_spec} SPECIFIC pair(s) found after checking {pair_counter} pair(s).\n\n"
                self.root.after(0, lambda t=summary: self._append_live_text(t, 'cont_success'))
            elif self._is_skip_variant_requested():
                summary = f"  >>> Skipped by user after checking {pair_counter} pair(s).\n\n"
                self.root.after(0, lambda t=summary: self._append_live_text(t, 'cont_skip'))
            elif self._is_cancel_requested():
                summary = f"  >>> Cancelled after checking {pair_counter} pair(s).\n\n"
                self.root.after(0, lambda t=summary: self._append_live_text(t, 'cont_skip'))
            else:
                summary = f"  >>> No specific pair found after checking {pair_counter} pair(s).\n\n"
                self.root.after(0, lambda t=summary: self._append_live_text(t, 'cont_fail'))

            # FASTA link
            if group.variants:
                row_num = group.variants[0].row_number
                fasta_info = fasta_files.get(row_num)
                if fasta_info and fasta_info[0]:
                    fasta_path = fasta_info[0]
                    self.root.after(0,
                        lambda p=fasta_path, r=row_num:
                            self._append_fasta_link_live(p, r))

            # Build DesignResult (specific first, then non-specific)
            all_amplicons = specific_amplicons + non_specific_amplicons
            result = DesignResult(
                variants=group.variants,
                amplicons=all_amplicons,
                success=len(all_amplicons) > 0,
                message=(f"Continuous: {len(specific_amplicons)} specific, "
                         f"{len(non_specific_amplicons)} non-specific "
                         f"({pair_counter} checked)"),
                homology_discriminated=design_result.homology_discriminated if design_result else False,
                num_homologous_regions=design_result.num_homologous_regions if design_result else 0,
                homology_tier=design_result.homology_tier if design_result else 0,
                homology_tier_message=design_result.homology_tier_message if design_result else "",
                homology_result=cont_homology,
            )
            results.append(result)

            progress = 10 + (85 * (group_idx + 1) / total_groups)
            self._set_progress(progress)

        # Restore original setting
        self.design_parameters.num_primer_pairs = original_num_pairs

        # Store results
        self.project_state.design_results = results
        self.project_state.is_designed = True
        self.project_state.fasta_files = fasta_files

        self.root.after(0, lambda: self.export_btn.configure(state=tk.NORMAL))
        self.root.after(0, lambda: self._update_help_context('results'))

        specific_count = sum(
            1 for r in results
            if r.success and any(
                a.primer_pair.specificity_verdict == "SPECIFIC"
                for a in r.amplicons
            )
        )
        self._set_status(
            f"Continuous design complete: {specific_count}/{len(results)} variants "
            f"have specific pairs"
        )
        self._set_progress(100)

    # ------------------------------------------------------------------
    #  Standard specificity checking (non-continuous mode)
    # ------------------------------------------------------------------

    def _check_specificity(self, results: list):
        """Check primer specificity using user-selected tools."""
        from ..primer.specificity_orchestrator import SpecificityOrchestrator

        use_primer_blast = self.use_primer_blast_check.get()
        use_ucsc = self.use_ucsc_ispcr_check.get()

        # Lazy-init the orchestrator (once per session)
        if not self.specificity_orchestrator:
            self.specificity_orchestrator = SpecificityOrchestrator()

        total_pairs = sum(len(r.amplicons) for r in results if r.success)
        checked = 0

        for result in results:
            if not result.success:
                continue
            for amplicon in result.amplicons:
                if self._is_cancel_requested():
                    return

                checked += 1
                gene = result.variants[0].gene_symbol if result.variants else '?'
                self._set_status(
                    f"Checking specificity {checked}/{total_pairs}: {gene}..."
                )
                self._set_progress(60 + (35 * checked / max(total_pairs, 1)))

                # Get expected coordinates from the variant
                exp_chr, exp_start, exp_end = None, None, None
                if result.variants and result.variants[0].genomic_position:
                    gp = result.variants[0].genomic_position
                    exp_chr = gp.chromosome
                    exp_start = gp.start
                    exp_end = gp.end

                combined = self.specificity_orchestrator.check_all(
                    amplicon.primer_pair,
                    expected_chromosome=exp_chr,
                    expected_start=exp_start,
                    expected_end=exp_end,
                    assembly=self.assembly_var.get(),
                    progress_callback=lambda msg: self._set_status(msg),
                    use_primer_blast=use_primer_blast,
                    use_ucsc=use_ucsc,
                )

                # Map combined result to PrimerPair fields
                pp = amplicon.primer_pair
                pp.is_specific = combined.is_specific
                pp.specificity_warnings = combined.all_warnings
                pp.specificity_verdict = combined.verdict.value
                pp.pseudogene_risk = bool(combined.pseudogene_warnings)
                pp.off_target_details = combined.off_target_summary

                # Build tool summary string (only show enabled tools)
                parts = []
                for tr in combined.tool_results:
                    if tr.error == "Skipped by user":
                        continue  # Don't show tools the user disabled
                    if tr.ran_successfully:
                        status = "OK" if tr.is_specific else f"{tr.off_target_count} off-target"
                    elif tr.available:
                        status = "ERROR"
                    else:
                        status = "N/A"
                    parts.append(f"{tr.tool_name}: {status}")
                pp.specificity_tool_summary = " | ".join(parts)

    def _display_results(self, results: list):
        """Display design results in the results panel."""
        self.results_text.delete('1.0', tk.END)

        # Configure hyperlink tag
        self.results_text.tag_configure('hyperlink', foreground='blue', underline=True)
        self.results_text.tag_bind('hyperlink', '<Enter>', lambda e: self.results_text.config(cursor='hand2'))
        self.results_text.tag_bind('hyperlink', '<Leave>', lambda e: self.results_text.config(cursor=''))

        # Sort results by row_number of first variant to match Variants table order
        sorted_results = sorted(results, key=lambda r: r.variants[0].row_number if r.variants else 0)

        # Find the global maximum product size across ALL results for proportional scaling
        max_product_size = 0
        for result in sorted_results:
            if result.success and result.amplicons:
                for amp in result.amplicons:
                    max_product_size = max(max_product_size, amp.primer_pair.product_size)

        # Max bar width for the longest amplicon
        MAX_BAR_WIDTH = 60

        for i, result in enumerate(sorted_results, 1):
            # Get variant info for header
            if result.variants:
                v = result.variants[0]
                gene = v.gene_symbol
                transcript = v.transcript_accession if v.transcript_accession else "?"
                hgvs = v.hgvs_c if v.hgvs_c else "?"
                row_num = v.row_number
            else:
                gene, transcript, hgvs = "?", "?", "?"
                row_num = 0

            # Main header line: number.GENE; TRANSCRIPT; VARIANT; status
            if result.success and result.amplicons:
                # Sort amplicons: best specificity first
                # Order: SPECIFIC (0), LIKELY SPECIFIC (1), no verdict (2),
                #        INCONCLUSIVE (3), OFF-TARGET by count (4+N), PSEUDOGENE RISK (1000)
                def _specificity_sort_key(amp):
                    pp = amp.primer_pair
                    v = pp.specificity_verdict
                    if v == "SPECIFIC":
                        return (0, 0)
                    elif v == "LIKELY SPECIFIC":
                        return (1, 0)
                    elif not v:
                        return (2, 0)       # not checked
                    elif v == "INCONCLUSIVE":
                        return (3, 0)
                    elif v == "OFF-TARGET DETECTED":
                        return (4, len(pp.off_target_details))
                    elif v == "PSEUDOGENE RISK":
                        return (1000, len(pp.off_target_details))
                    return (5, 0)

                sorted_amplicons = sorted(result.amplicons, key=_specificity_sort_key)
                status_text = f"Found {len(sorted_amplicons)} primer pair(s)"
                homology_tag = ""
                if result.homology_discriminated:
                    homology_tag = f" [homology-optimized, {result.num_homologous_regions} region(s)]"
                self.results_text.insert(tk.END, f"{i}.{gene}; {transcript}; {hgvs}; {status_text}{homology_tag}\n")

                # Homology tier info
                if result.homology_discriminated and result.homology_tier > 0:
                    tier = result.homology_tier
                    tier_labels = {
                        1: ("✓ Tier 1", "primers overlap discriminating positions", "specific"),
                        2: ("✓ Tier 2", "primers in discriminating windows", "specific"),
                        3: ("⚠ Tier 3", "fallback re-ranking, discrimination not guaranteed", "inconclusive"),
                    }
                    label, desc, tag = tier_labels.get(tier, (f"Tier {tier}", "", ""))
                    self.results_text.insert(tk.END, f"  Homology: {label} — {desc}\n", tag)
                    if result.homology_tier_message:
                        self.results_text.insert(tk.END, f"  {result.homology_tier_message}\n")

                    # Show detailed homologous regions
                    if result.homology_result:
                        homology_detail = self._format_homology_details(result.homology_result)
                        if homology_detail:
                            self.results_text.insert(tk.END, homology_detail, 'warning')

                # Show each primer pair with structure visualization
                for j, amp in enumerate(sorted_amplicons, 1):
                    pp = amp.primer_pair
                    # Line 1: Primer sequences
                    self.results_text.insert(tk.END,
                        f"Primer Pair {j}: F: 5'-{pp.forward.sequence}-3'; "
                        f"R: 5'-{pp.reverse.sequence}-3'\n")

                    # Line 2: ASCII art structure visualization
                    fwd_len = pp.forward.length
                    rev_len = pp.reverse.length

                    # Distance from end of forward primer to start of target
                    dist_fwd_to_target = amp.target_start - pp.forward.end - 1
                    if dist_fwd_to_target < 0:
                        dist_fwd_to_target = 0

                    # Distance from end of target to start of reverse primer
                    target_end = amp.target_start + amp.target_length - 1
                    dist_target_to_rev = pp.reverse.start - target_end - 1
                    if dist_target_to_rev < 0:
                        dist_target_to_rev = 0

                    # Scale: 1 nt = how many chars? Based on global max
                    # The longest amplicon gets MAX_BAR_WIDTH chars, shorter ones get proportionally less
                    scale = MAX_BAR_WIDTH / max_product_size if max_product_size > 0 else 1

                    fwd_bar = max(1, round(fwd_len * scale))
                    gap1_bar = max(1, round(dist_fwd_to_target * scale))
                    gap2_bar = max(1, round(dist_target_to_rev * scale))
                    rev_bar = max(1, round(rev_len * scale))

                    # Build the visual bar
                    visual = (
                        f"{'▶' * fwd_bar}"
                        f"{'─' * gap1_bar}"
                        f"◆"
                        f"{'─' * gap2_bar}"
                        f"{'◀' * rev_bar}"
                    )
                    total_bar_len = fwd_bar + gap1_bar + 1 + gap2_bar + rev_bar

                    # Build exon/intron map line (same width as visual bar)
                    exon_map = self._build_exon_map(
                        amp, pp, total_bar_len, scale
                    )

                    # Structure line with numbers
                    structure_line = f"  {visual}  {pp.product_size}bp\n"
                    if exon_map:
                        structure_line += f"  {exon_map}\n"
                    structure_line += (
                        f"  F:{fwd_len}nt ─{dist_fwd_to_target}nt─ [VAR] ─{dist_target_to_rev}nt─ R:{rev_len}nt\n"
                    )
                    self.results_text.insert(tk.END, structure_line)

                    # Specificity verdict display
                    self._insert_specificity_verdict(pp)

                # Files reference with clickable link
                fasta_info = self.project_state.fasta_files.get(row_num) if hasattr(self.project_state, 'fasta_files') else None
                if fasta_info and fasta_info[0]:
                    fasta_path = fasta_info[0]
                    self.results_text.insert(tk.END, "Files: ")
                    # Insert clickable link
                    link_tag = f'fasta_link_{row_num}'
                    self.results_text.tag_configure(link_tag, foreground='blue', underline=True)
                    self.results_text.tag_bind(link_tag, '<Button-1>', lambda e, p=fasta_path: self._open_file(p))
                    self.results_text.tag_bind(link_tag, '<Enter>', lambda e: self.results_text.config(cursor='hand2'))
                    self.results_text.tag_bind(link_tag, '<Leave>', lambda e: self.results_text.config(cursor=''))
                    self.results_text.insert(tk.END, fasta_path, link_tag)
                    self.results_text.insert(tk.END, "\n")
                else:
                    self.results_text.insert(tk.END, f"Files: (no FASTA generated)\n")

            else:
                # Failed design
                self.results_text.insert(tk.END, f"{i}.{gene}; {transcript}; {hgvs}; ", '')
                self.results_text.insert(tk.END, f"FAILED\n", 'error')
                self.results_text.insert(tk.END, f"{result.message}\n")

                # Show suggestions
                for suggestion in result.suggestions:
                    self.results_text.insert(tk.END, f"  - {suggestion}\n", 'warning')

            self.results_text.insert(tk.END, "\n")

    @staticmethod
    def _build_exon_map(amp, pp, total_bar_len: int, scale: float) -> str:
        """Build an exon/intron annotation line for the ASCII art amplicon.

        Returns a string of the same width as the visual bar, with:
          ═ for exonic positions (labeled E<n> in the centre)
          ╌ for intronic positions
        Returns empty string if no exon data is available.
        """
        exon_regions = amp.exon_regions_in_amplicon
        if not exon_regions:
            return ""

        product_size = pp.product_size
        if product_size <= 0:
            return ""

        # Build a character array for the product, then scale to bar width
        # First build at nucleotide resolution, then downsample
        # Use a list: 'E' for exon, 'I' for intron, with exon number stored separately
        nt_map = ['I'] * product_size
        nt_exon_num = [0] * product_size
        for rel_start, rel_end, exon_num in exon_regions:
            for pos in range(max(0, rel_start), min(product_size, rel_end + 1)):
                nt_map[pos] = 'E'
                nt_exon_num[pos] = exon_num

        # Downsample to bar width: for each bar character, majority vote
        bar = []
        bar_exon_num = []
        for ci in range(total_bar_len):
            # Map bar position to nucleotide range
            nt_start = int(ci * product_size / total_bar_len)
            nt_end = int((ci + 1) * product_size / total_bar_len)
            nt_end = max(nt_end, nt_start + 1)

            segment = nt_map[nt_start:nt_end]
            exon_count = segment.count('E')
            if exon_count >= len(segment) / 2:
                bar.append('E')
                # Most common exon number in segment
                nums = [nt_exon_num[p] for p in range(nt_start, min(nt_end, product_size)) if nt_map[p] == 'E']
                bar_exon_num.append(max(set(nums), key=nums.count) if nums else 0)
            else:
                bar.append('I')
                bar_exon_num.append(0)

        # Check if everything is exonic (single exon spanning whole amplicon)
        all_exon = all(c == 'E' for c in bar)
        all_intron = all(c == 'I' for c in bar)
        if all_intron:
            return ""  # No exon info to show

        # Build display string with labels
        result = list('═' * total_bar_len if all_exon else ' ' * total_bar_len)
        # Fill exon/intron characters
        for ci in range(total_bar_len):
            if bar[ci] == 'E':
                result[ci] = '═'
            else:
                result[ci] = '╌'

        # Add exon labels (E<n>) centred on each exonic stretch
        # Find contiguous exon stretches
        stretches = []  # (start_idx, end_idx, exon_number)
        i = 0
        while i < total_bar_len:
            if bar[i] == 'E':
                start = i
                num = bar_exon_num[i]
                while i < total_bar_len and bar[i] == 'E' and bar_exon_num[i] == num:
                    i += 1
                stretches.append((start, i - 1, num))
            else:
                i += 1

        for s_start, s_end, s_num in stretches:
            stretch_len = s_end - s_start + 1
            if s_num > 0:
                label = f"E{s_num}"
            else:
                label = "Ex"
            if stretch_len >= len(label) + 2:
                # Centre the label in the stretch
                mid = (s_start + s_end) // 2
                lbl_start = mid - len(label) // 2
                for ci, ch in enumerate(label):
                    if 0 <= lbl_start + ci < total_bar_len:
                        result[lbl_start + ci] = ch

        return ''.join(result)

    def _insert_specificity_verdict(self, pp):
        """Insert specificity verdict block for a primer pair into results_text."""
        verdict = pp.specificity_verdict
        if not verdict:
            self.results_text.insert(tk.END, "  Specificity: ", '')
            self.results_text.insert(tk.END, "not checked\n", 'inconclusive')
            return

        # Verdict line with color
        self.results_text.insert(tk.END, "  Specificity: ", '')
        if verdict == "SPECIFIC":
            self.results_text.insert(tk.END, "✓ SPECIFIC", 'specific')
            self.results_text.insert(tk.END, " — amplifies ONLY the intended target\n")
        elif verdict == "PSEUDOGENE RISK":
            self.results_text.insert(tk.END, "⚠ PSEUDOGENE RISK", 'pseudogene')
            self.results_text.insert(tk.END, " — may amplify pseudogene/duplication!\n")
        elif verdict == "OFF-TARGET DETECTED":
            self.results_text.insert(tk.END, "✗ OFF-TARGET DETECTED", 'not_specific')
            self.results_text.insert(tk.END, " — amplifies additional genomic sites!\n")
        elif verdict == "LIKELY SPECIFIC":
            self.results_text.insert(tk.END, "✓ LIKELY SPECIFIC", 'specific')
            self.results_text.insert(tk.END, " (not all tools were available)\n")
        elif verdict == "INCONCLUSIVE":
            self.results_text.insert(tk.END, "? INCONCLUSIVE", 'inconclusive')
            self.results_text.insert(tk.END, " — could not verify specificity\n")
        else:
            self.results_text.insert(tk.END, f"{verdict}\n")

        # Tool summary line
        if pp.specificity_tool_summary:
            self.results_text.insert(tk.END, f"  Tools: {pp.specificity_tool_summary}\n")

        # Off-target details (first 3)
        if pp.off_target_details:
            for detail in pp.off_target_details[:3]:
                self.results_text.insert(tk.END, f"    {detail}\n", 'warning')
            if len(pp.off_target_details) > 3:
                self.results_text.insert(
                    tk.END,
                    f"    ... and {len(pp.off_target_details) - 3} more off-target(s)\n",
                    'warning'
                )

        # Pseudogene bold warning
        if pp.pseudogene_risk:
            self.results_text.insert(
                tk.END,
                "  *** WARNING: Same-chromosome off-target — possible pseudogene. "
                "Consider redesigning this primer pair. ***\n",
                'pseudogene'
            )

        # Homology discrimination score
        if pp.homology_discrimination_score > 0:
            self.results_text.insert(tk.END, "  Homology discrimination: ", '')
            score = pp.homology_discrimination_score
            fwd_n = pp.fwd_discriminating_positions
            rev_n = pp.rev_discriminating_positions
            if score >= 6:
                tag = 'specific'
                label = "GOOD"
            elif score >= 2:
                tag = 'inconclusive'
                label = "MODERATE"
            else:
                tag = 'warning'
                label = "WEAK"
            self.results_text.insert(
                tk.END,
                f"{label} (score={score:.1f}, F:{fwd_n} pos, R:{rev_n} pos)\n",
                tag
            )
        elif pp.homology_warning:
            self.results_text.insert(tk.END, "  Homology discrimination: ", '')
            self.results_text.insert(
                tk.END,
                f"NONE — {pp.homology_warning}\n",
                'pseudogene'
            )

    def _open_file(self, file_path: str):
        """Open a file with the default system application."""
        import subprocess
        import platform
        try:
            if platform.system() == 'Darwin':  # macOS
                subprocess.run(['open', file_path], check=True)
            elif platform.system() == 'Windows':
                subprocess.run(['start', '', file_path], shell=True, check=True)
            else:  # Linux
                subprocess.run(['xdg-open', file_path], check=True)
        except Exception as e:
            self.log_warning(f"Could not open file {file_path}: {e}")

    def _cancel_operation(self):
        """Cancel current operation."""
        self._set_cancel_requested(True)
        self.project_state.cancel_requested = True

        if self.primer_designer:
            self.primer_designer.request_cancel()

        self._set_status("Cancelling...")
        self.cancel_btn.configure(state=tk.DISABLED)

    def _set_cancel_requested(self, value: bool):
        """Thread-safe setter for cancel flag."""
        with self._cancel_lock:
            self._cancel_requested = value

    def _is_cancel_requested(self) -> bool:
        """Thread-safe getter for cancel flag."""
        with self._cancel_lock:
            return self._cancel_requested

    # ------------------------------------------------------------------
    #  Continuous mode helpers
    # ------------------------------------------------------------------

    def _on_continuous_mode_toggle(self):
        """Enable/disable continuous mode options based on checkbox state."""
        state = tk.NORMAL if self.continuous_mode_var.get() else tk.DISABLED
        self.max_continuous_entry.configure(state=state)
        self.min_specific_pairs_entry.configure(state=state)

    def _skip_to_next_variant(self):
        """Request skipping current variant in continuous mode."""
        with self._skip_variant_lock:
            self._skip_variant_requested = True
        self._set_status("Skipping to next variant...")

    def _is_skip_variant_requested(self) -> bool:
        """Thread-safe check for skip-variant flag."""
        with self._skip_variant_lock:
            return self._skip_variant_requested

    def _reset_skip_variant(self):
        """Thread-safe reset of skip-variant flag."""
        with self._skip_variant_lock:
            self._skip_variant_requested = False

    def _append_live_text(self, text: str, tag: str = ''):
        """Append text to results_text (must be called on main thread via root.after)."""
        if tag:
            self.results_text.insert(tk.END, text, tag)
        else:
            self.results_text.insert(tk.END, text)
        self.results_text.see(tk.END)

    def _display_single_pair_live(self, amplicon, pair_number):
        """Display a single primer pair result live in results panel.
        Must be called on main thread via root.after."""
        pp = amplicon.primer_pair
        fwd_len = pp.forward.length
        rev_len = pp.reverse.length

        # Primer sequences
        self.results_text.insert(tk.END,
            f"  Pair {pair_number}: F: 5'-{pp.forward.sequence}-3'; "
            f"R: 5'-{pp.reverse.sequence}-3'\n")

        # Distances
        dist_fwd = max(0, amplicon.target_start - pp.forward.end - 1)
        target_end = amplicon.target_start + amplicon.target_length - 1
        dist_rev = max(0, pp.reverse.start - target_end - 1)

        # Exon/intron map
        if amplicon.exon_regions_in_amplicon:
            # Use a simple scale: 1 char per ~5nt, or at least 30 chars
            bar_len = max(30, min(60, pp.product_size // 5))
            scale = bar_len / pp.product_size if pp.product_size > 0 else 1
            exon_map = self._build_exon_map(amplicon, pp, bar_len, scale)
            if exon_map:
                self.results_text.insert(tk.END, f"  {exon_map}\n")

        self.results_text.insert(tk.END,
            f"  F:{fwd_len}nt —{dist_fwd}nt— [VAR] —{dist_rev}nt— R:{rev_len}nt  |  "
            f"{pp.product_size}bp\n")

        # Specificity verdict (reuse existing method)
        self._insert_specificity_verdict(pp)
        self.results_text.insert(tk.END, "\n")
        self.results_text.see(tk.END)

    def _format_homology_details(self, homology_result) -> str:
        """Format homologous region details for live display in results text."""
        if homology_result is None:
            return ''

        secondary_hits = [
            h for h in homology_result.hits
            if not h.is_primary and not h.is_supplementary
        ]
        if not secondary_hits:
            return ''

        primary = homology_result.primary_hit
        primary_chr = primary.chromosome if primary else homology_result.query_chromosome
        primary_pos = primary.position if primary else homology_result.query_position

        lines = [f"  Homologous regions ({len(secondary_hits)}) — target: {primary_chr}:{primary_pos:,}\n"]
        for hit in secondary_hits:
            same_chr_flag = ""
            norm_hit = hit.chromosome.replace('chr', '')
            norm_primary = str(primary_chr).replace('chr', '')
            if norm_hit == norm_primary:
                same_chr_flag = " [SAME CHR!]"
            lines.append(
                f"    • {hit.chromosome}:{hit.position:,} ({hit.strand}) "
                f"— {hit.percent_identity:.1f}% identity, "
                f"{hit.aligned_length}bp aligned, "
                f"{hit.mismatches} mismatches, "
                f"score={hit.bit_score:.0f}"
                f"{same_chr_flag}\n"
            )
        return ''.join(lines)

    def _append_fasta_link_live(self, fasta_path, row_num):
        """Append clickable FASTA file link. Must be called on main thread."""
        self.results_text.insert(tk.END, "  Files: ")
        link_tag = f'fasta_link_{row_num}'
        self.results_text.tag_configure(link_tag, foreground='blue', underline=True)
        self.results_text.tag_bind(link_tag, '<Button-1>',
            lambda e, p=fasta_path: self._open_file(p))
        self.results_text.tag_bind(link_tag, '<Enter>',
            lambda e: self.results_text.config(cursor='hand2'))
        self.results_text.tag_bind(link_tag, '<Leave>',
            lambda e: self.results_text.config(cursor=''))
        self.results_text.insert(tk.END, fasta_path, link_tag)
        self.results_text.insert(tk.END, "\n\n")
        self.results_text.see(tk.END)

    def _export_report(self):
        """Export results to HTML report."""
        if not self.project_state.design_results:
            messagebox.showerror("Error", "No results to export")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".html",
            filetypes=[("HTML files", "*.html")],
            initialfile="primer_design_report.html"
        )

        if filename:
            # Build methodology context for the report
            try:
                max_att = int(self.max_continuous_var.get())
            except ValueError:
                max_att = 20
            methodology_context = {
                'used_ucsc': self.use_ucsc_ispcr_check.get(),
                'used_primer_blast': self.use_primer_blast_check.get(),
                'continuous_mode': self.continuous_mode_var.get(),
                'max_attempts': max_att,
                'population': self.selected_population.get(),
                'flanking_region': self.flank_size_var.get(),
                'blast_min_identity': self.config.blast.min_percent_identity,
                'blast_min_aligned_length': self.config.blast.min_aligned_length,
                'variant_maf_thresholds': dict(self.variant_maf_thresholds),
                'default_maf_threshold_pct': self.default_maf_threshold_pct,
            }
            success = self.report_generator.generate_report(
                self.project_state.design_results,
                self.design_parameters,
                filename,
                methodology_context=methodology_context
            )

            if success:
                messagebox.showinfo("Success", f"Report saved to {filename}")
            else:
                messagebox.showerror("Error", "Failed to generate report")

    def _show_advanced_settings(self):
        """Show advanced settings dialog."""
        AdvancedSettingsDialog(self.root, self.design_parameters)

    def _set_status(self, message: str):
        """Update status bar message (thread-safe)."""
        self.root.after(0, lambda: self.status_var.set(message))

    def _set_progress(self, value: float):
        """Update progress bar (thread-safe)."""
        self.root.after(0, lambda: self.progress_var.set(value))

    def _run_in_thread(self, target):
        """Run a function in a separate thread."""
        self._set_cancel_requested(False)
        self.cancel_btn.configure(state=tk.NORMAL)
        self._start_spinner()

        def wrapper():
            try:
                target()
            finally:
                self.root.after(0, lambda: self.cancel_btn.configure(state=tk.DISABLED))
                self.root.after(0, self._stop_spinner)

        self._current_thread = threading.Thread(target=wrapper, daemon=True)
        self._current_thread.start()

    def _on_close(self):
        """Handle window close event."""
        if self._current_thread and self._current_thread.is_alive():
            self._set_cancel_requested(True)
            self._current_thread.join(timeout=2)

        self.root.destroy()

    def run(self):
        """Start the application main loop."""
        self.log_info("Starting application")
        self.root.mainloop()


class AdvancedSettingsDialog:
    """Dialog for advanced primer design settings."""

    def __init__(self, parent, parameters: DesignParameters):
        self.parameters = parameters

        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Advanced Settings")
        self.dialog.geometry("500x600")
        self.dialog.transient(parent)
        self.dialog.grab_set()

        self._build_ui()

    def _build_ui(self):
        """Build dialog UI."""
        notebook = ttk.Notebook(self.dialog)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Primer parameters tab
        primer_frame = ttk.Frame(notebook, padding=10)
        notebook.add(primer_frame, text="Primer Parameters")
        self._build_primer_params(primer_frame)

        # Chemistry tab
        chem_frame = ttk.Frame(notebook, padding=10)
        notebook.add(chem_frame, text="Chemistry")
        self._build_chemistry_params(chem_frame)

        # qPCR tab (hidden for now — kept in code but not shown in GUI)
        # qpcr_frame = ttk.Frame(notebook, padding=10)
        # notebook.add(qpcr_frame, text="qPCR Settings")
        # self._build_qpcr_params(qpcr_frame)

        # Buttons
        btn_frame = ttk.Frame(self.dialog)
        btn_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(btn_frame, text="OK", command=self._save_and_close).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=self.dialog.destroy).pack(side=tk.RIGHT)
        ttk.Button(btn_frame, text="Reset Defaults", command=self._reset_defaults).pack(side=tk.LEFT)

    def _build_primer_params(self, parent):
        """Build primer parameters section."""
        p = self.parameters

        # Primer size
        ttk.Label(parent, text="Primer Size (nt):").grid(row=0, column=0, sticky='w', pady=5)
        self.primer_min_size = tk.StringVar(value=str(p.primer_min_size))
        self.primer_max_size = tk.StringVar(value=str(p.primer_max_size))
        self.primer_opt_size = tk.StringVar(value=str(p.primer_opt_size))

        size_frame = ttk.Frame(parent)
        size_frame.grid(row=0, column=1, sticky='w')
        ttk.Entry(size_frame, textvariable=self.primer_min_size, width=5).pack(side=tk.LEFT)
        ttk.Label(size_frame, text=" - ").pack(side=tk.LEFT)
        ttk.Entry(size_frame, textvariable=self.primer_max_size, width=5).pack(side=tk.LEFT)
        ttk.Label(size_frame, text="  Opt: ").pack(side=tk.LEFT)
        ttk.Entry(size_frame, textvariable=self.primer_opt_size, width=5).pack(side=tk.LEFT)

        # GC content
        ttk.Label(parent, text="GC Content (%):").grid(row=1, column=0, sticky='w', pady=5)
        self.gc_min = tk.StringVar(value=str(p.primer_min_gc))
        self.gc_max = tk.StringVar(value=str(p.primer_max_gc))
        self.gc_opt = tk.StringVar(value=str(p.primer_opt_gc))

        gc_frame = ttk.Frame(parent)
        gc_frame.grid(row=1, column=1, sticky='w')
        ttk.Entry(gc_frame, textvariable=self.gc_min, width=5).pack(side=tk.LEFT)
        ttk.Label(gc_frame, text=" - ").pack(side=tk.LEFT)
        ttk.Entry(gc_frame, textvariable=self.gc_max, width=5).pack(side=tk.LEFT)
        ttk.Label(gc_frame, text="  Opt: ").pack(side=tk.LEFT)
        ttk.Entry(gc_frame, textvariable=self.gc_opt, width=5).pack(side=tk.LEFT)

        # Complementarity
        ttk.Label(parent, text="Max Self Complementarity:").grid(row=2, column=0, sticky='w', pady=5)
        self.self_compl = tk.StringVar(value=str(p.max_self_complementarity))
        ttk.Entry(parent, textvariable=self.self_compl, width=8).grid(row=2, column=1, sticky='w')

        ttk.Label(parent, text="Max End Complementarity:").grid(row=3, column=0, sticky='w', pady=5)
        self.end_compl = tk.StringVar(value=str(p.max_end_complementarity))
        ttk.Entry(parent, textvariable=self.end_compl, width=8).grid(row=3, column=1, sticky='w')

        ttk.Label(parent, text="Max Pair Complementarity:").grid(row=4, column=0, sticky='w', pady=5)
        self.pair_compl = tk.StringVar(value=str(p.max_pair_complementarity))
        ttk.Entry(parent, textvariable=self.pair_compl, width=8).grid(row=4, column=1, sticky='w')

        ttk.Label(parent, text="Optimal Tm (°C):").grid(row=5, column=0, sticky='w', pady=5)
        self.opt_tm = tk.StringVar(value=str(p.primer_opt_tm))
        ttk.Entry(parent, textvariable=self.opt_tm, width=8).grid(row=5, column=1, sticky='w')

        ttk.Label(parent, text="Max Tm Difference (°C):").grid(row=6, column=0, sticky='w', pady=5)
        self.tm_diff = tk.StringVar(value=str(p.max_tm_difference))
        ttk.Entry(parent, textvariable=self.tm_diff, width=8).grid(row=6, column=1, sticky='w')

    def _build_chemistry_params(self, parent):
        """Build chemistry parameters section."""
        p = self.parameters

        ttk.Label(parent, text="Monovalent Cations (mM):").grid(row=0, column=0, sticky='w', pady=5)
        self.mv_conc = tk.StringVar(value=str(p.mv_conc))
        ttk.Entry(parent, textvariable=self.mv_conc, width=8).grid(row=0, column=1, sticky='w')

        ttk.Label(parent, text="Divalent Cations - Mg2+ (mM):").grid(row=1, column=0, sticky='w', pady=5)
        self.dv_conc = tk.StringVar(value=str(p.dv_conc))
        ttk.Entry(parent, textvariable=self.dv_conc, width=8).grid(row=1, column=1, sticky='w')

        ttk.Label(parent, text="dNTP Concentration (mM):").grid(row=2, column=0, sticky='w', pady=5)
        self.dntp_conc = tk.StringVar(value=str(p.dntp_conc))
        ttk.Entry(parent, textvariable=self.dntp_conc, width=8).grid(row=2, column=1, sticky='w')

        ttk.Label(parent, text="DNA Concentration (nM):").grid(row=3, column=0, sticky='w', pady=5)
        self.dna_conc = tk.StringVar(value=str(p.dna_conc))
        ttk.Entry(parent, textvariable=self.dna_conc, width=8).grid(row=3, column=1, sticky='w')

    def _build_qpcr_params(self, parent):
        """Build qPCR parameters section."""
        p = self.parameters

        ttk.Label(parent, text="qPCR Amplicon Size (bp):").grid(row=0, column=0, sticky='w', pady=5)
        self.qpcr_min = tk.StringVar(value=str(p.qpcr_min_amplicon_size))
        self.qpcr_max = tk.StringVar(value=str(p.qpcr_max_amplicon_size))

        qpcr_frame = ttk.Frame(parent)
        qpcr_frame.grid(row=0, column=1, sticky='w')
        ttk.Entry(qpcr_frame, textvariable=self.qpcr_min, width=5).pack(side=tk.LEFT)
        ttk.Label(qpcr_frame, text=" - ").pack(side=tk.LEFT)
        ttk.Entry(qpcr_frame, textvariable=self.qpcr_max, width=5).pack(side=tk.LEFT)

        ttk.Label(parent, text="Probe Size (nt):").grid(row=1, column=0, sticky='w', pady=5)
        self.probe_min = tk.StringVar(value=str(p.probe_min_size))
        self.probe_max = tk.StringVar(value=str(p.probe_max_size))

        probe_frame = ttk.Frame(parent)
        probe_frame.grid(row=1, column=1, sticky='w')
        ttk.Entry(probe_frame, textvariable=self.probe_min, width=5).pack(side=tk.LEFT)
        ttk.Label(probe_frame, text=" - ").pack(side=tk.LEFT)
        ttk.Entry(probe_frame, textvariable=self.probe_max, width=5).pack(side=tk.LEFT)

        ttk.Label(parent, text="Probe Tm (°C):").grid(row=2, column=0, sticky='w', pady=5)
        self.probe_tm_min = tk.StringVar(value=str(p.probe_min_tm))
        self.probe_tm_max = tk.StringVar(value=str(p.probe_max_tm))

        probe_tm_frame = ttk.Frame(parent)
        probe_tm_frame.grid(row=2, column=1, sticky='w')
        ttk.Entry(probe_tm_frame, textvariable=self.probe_tm_min, width=5).pack(side=tk.LEFT)
        ttk.Label(probe_tm_frame, text=" - ").pack(side=tk.LEFT)
        ttk.Entry(probe_tm_frame, textvariable=self.probe_tm_max, width=5).pack(side=tk.LEFT)

    def _save_and_close(self):
        """Save settings and close dialog."""
        try:
            p = self.parameters

            # Primer params
            p.primer_min_size = int(self.primer_min_size.get())
            p.primer_max_size = int(self.primer_max_size.get())
            p.primer_opt_size = int(self.primer_opt_size.get())
            p.primer_min_gc = float(self.gc_min.get())
            p.primer_max_gc = float(self.gc_max.get())
            p.primer_opt_gc = float(self.gc_opt.get())
            p.max_self_complementarity = float(self.self_compl.get())
            p.max_end_complementarity = float(self.end_compl.get())
            p.max_pair_complementarity = float(self.pair_compl.get())
            p.primer_opt_tm = float(self.opt_tm.get())
            p.max_tm_difference = float(self.tm_diff.get())

            # Chemistry
            p.mv_conc = float(self.mv_conc.get())
            p.dv_conc = float(self.dv_conc.get())
            p.dntp_conc = float(self.dntp_conc.get())
            p.dna_conc = float(self.dna_conc.get())

            # qPCR (hidden for now — skip saving qPCR params)
            # p.qpcr_min_amplicon_size = int(self.qpcr_min.get())
            # p.qpcr_max_amplicon_size = int(self.qpcr_max.get())
            # p.probe_min_size = int(self.probe_min.get())
            # p.probe_max_size = int(self.probe_max.get())
            # p.probe_min_tm = float(self.probe_tm_min.get())
            # p.probe_max_tm = float(self.probe_tm_max.get())

            self.dialog.destroy()

        except ValueError as e:
            messagebox.showerror("Error", f"Invalid value: {e}")

    def _reset_defaults(self):
        """Reset to default values."""
        defaults = DesignParameters()

        self.primer_min_size.set(str(defaults.primer_min_size))
        self.primer_max_size.set(str(defaults.primer_max_size))
        self.primer_opt_size.set(str(defaults.primer_opt_size))
        self.gc_min.set(str(defaults.primer_min_gc))
        self.gc_max.set(str(defaults.primer_max_gc))
        self.gc_opt.set(str(defaults.primer_opt_gc))
        self.self_compl.set(str(defaults.max_self_complementarity))
        self.end_compl.set(str(defaults.max_end_complementarity))
        self.pair_compl.set(str(defaults.max_pair_complementarity))
        self.opt_tm.set(str(defaults.primer_opt_tm))
        self.tm_diff.set(str(defaults.max_tm_difference))
        self.mv_conc.set(str(defaults.mv_conc))
        self.dv_conc.set(str(defaults.dv_conc))
        self.dntp_conc.set(str(defaults.dntp_conc))
        self.dna_conc.set(str(defaults.dna_conc))
        # qPCR defaults (hidden for now)
        # self.qpcr_min.set(str(defaults.qpcr_min_amplicon_size))
        # self.qpcr_max.set(str(defaults.qpcr_max_amplicon_size))
        # self.probe_min.set(str(defaults.probe_min_size))
        # self.probe_max.set(str(defaults.probe_max_size))
        # self.probe_tm_min.set(str(defaults.probe_min_tm))
        # self.probe_tm_max.set(str(defaults.probe_max_tm))

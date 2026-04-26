"""
Dark theme for Primer Designer GUI.

Uses sv_ttk (Sun Valley) for a modern Windows 11 / Fluent-style dark theme
with proper checkboxes, radio buttons, and rounded widgets.
Falls back to a manual clam-based dark theme if sv_ttk is not installed.
"""

import tkinter as tk
from tkinter import ttk

# =============================================================================
#  COLOR PALETTE  (Catppuccin Mocha inspired, tuned for sv_ttk dark)
# =============================================================================

# Base colors
BG = '#1e1e2e'            # Main background (dark blue-gray)
BG_LIGHT = '#2a2a3c'      # Slightly lighter (panels, frames)
BG_PANEL = '#313145'       # Panel/card background
BG_INPUT = '#3b3b50'       # Entry/text input background
BG_HOVER = '#44445a'       # Hover state
BG_SELECTED = '#4a4a6a'    # Selected row

FG = '#cdd6f4'             # Main text (light lavender)
FG_DIM = '#9399b2'         # Secondary/dimmed text
FG_BRIGHT = '#ffffff'       # Bright text (headers)

# Accent colors
ACCENT = '#89b4fa'          # Primary accent (blue)
ACCENT_HOVER = '#b4d0fb'    # Accent hover
ACCENT_GREEN = '#a6e3a1'    # Success green
ACCENT_RED = '#f38ba8'      # Error red
ACCENT_YELLOW = '#f9e2af'   # Warning yellow
ACCENT_ORANGE = '#fab387'   # Orange
ACCENT_PURPLE = '#cba6f7'   # Purple
ACCENT_TEAL = '#94e2d5'     # Teal
ACCENT_PINK = '#f5c2e7'     # Pink

# Borders and separators
BORDER = '#45475a'          # Border color
SASH = '#585b70'            # PanedWindow sash

# Treeview
TREE_BG = '#1e1e2e'
TREE_FG = '#cdd6f4'
TREE_SELECTED_BG = '#45475a'
TREE_SELECTED_FG = '#ffffff'

# Homology tree row colors
HOMOLOGY_PRIMARY = '#1a3a1a'
HOMOLOGY_SAME_CHR = '#3a2a1a'
HOMOLOGY_OTHER_CHR = '#1a1a3a'

# Scrollbar
SCROLL_BG = '#313145'
SCROLL_TROUGH = '#1e1e2e'

# Button
BTN_BG = '#45475a'
BTN_FG = '#cdd6f4'
BTN_ACTIVE_BG = '#585b70'

# Spinner
SPINNER_IDLE = '#585b70'
SPINNER_BG = '#45475a'
SPINNER_ARC = '#89b4fa'

# Sequence/homology text (darkest)
SEQ_BG = '#11111b'

# Whether sv_ttk was successfully loaded
_sv_ttk_active = False


def apply_dark_theme(root: tk.Tk):
    """Apply dark theme to the entire application.

    Uses sv_ttk for modern look. Falls back to clam-based manual theme.
    Must be called AFTER tk.Tk() but BEFORE building widgets.
    """
    global _sv_ttk_active

    # Try sv_ttk first (modern Fluent-style theme)
    try:
        import sv_ttk
        sv_ttk.set_theme("dark")
        _sv_ttk_active = True
    except ImportError:
        # Fallback: manual clam theme
        _apply_clam_fallback(root)
        _sv_ttk_active = False

    # Common: set root bg and options for classic tk widgets
    root.configure(bg=BG)
    root.option_add('*background', BG)
    root.option_add('*foreground', FG)
    root.option_add('*PanedWindow.background', SASH)

    # Customize ttk styles on top of sv_ttk / clam
    _apply_custom_styles()


def _apply_clam_fallback(root: tk.Tk):
    """Fallback dark theme using clam base (when sv_ttk not available)."""
    style = ttk.Style()
    available = style.theme_names()
    if 'clam' in available:
        style.theme_use('clam')
    elif 'alt' in available:
        style.theme_use('alt')


def _apply_custom_styles():
    """Apply custom named styles and overrides on top of base theme."""
    style = ttk.Style()

    if _sv_ttk_active:
        # sv_ttk handles most widgets beautifully — just customize specifics

        # Treeview: override sv_ttk defaults with our palette
        style.configure('Treeview',
                        background=TREE_BG, foreground=TREE_FG,
                        fieldbackground=TREE_BG,
                        font=('Helvetica', 10))
        style.configure('Treeview.Heading',
                        font=('Helvetica', 10, 'bold'))
        style.map('Treeview',
                  background=[('selected', TREE_SELECTED_BG)],
                  foreground=[('selected', TREE_SELECTED_FG)])

        # LabelFrame label color
        style.configure('TLabelframe.Label', foreground=ACCENT)

    else:
        # Full manual styling for clam fallback
        style.configure('.', background=BG, foreground=FG, borderwidth=0,
                        focuscolor=ACCENT, font=('Helvetica', 10))

        style.configure('TFrame', background=BG)
        style.configure('TLabel', background=BG, foreground=FG)

        style.configure('TButton', background=BTN_BG, foreground=BTN_FG,
                        borderwidth=1, relief='flat', padding=(10, 4))
        style.map('TButton',
                  background=[('active', BTN_ACTIVE_BG), ('disabled', BG_LIGHT)],
                  foreground=[('disabled', FG_DIM)])

        style.configure('TLabelframe', background=BG, foreground=FG,
                        bordercolor=BORDER, relief='groove')
        style.configure('TLabelframe.Label', background=BG, foreground=ACCENT)

        style.configure('TEntry', fieldbackground=BG_INPUT, foreground=FG,
                        insertcolor=FG, bordercolor=BORDER,
                        lightcolor=BORDER, darkcolor=BORDER)
        style.map('TEntry',
                  fieldbackground=[('readonly', BG_LIGHT), ('disabled', BG_LIGHT)],
                  foreground=[('disabled', FG_DIM)])

        style.configure('TCheckbutton', background=BG, foreground=FG,
                        indicatorbackground=BG_INPUT, indicatorforeground=ACCENT)
        style.map('TCheckbutton',
                  background=[('active', BG_HOVER)],
                  indicatorbackground=[('selected', ACCENT)])

        style.configure('TRadiobutton', background=BG, foreground=FG,
                        indicatorbackground=BG_INPUT, indicatorforeground=ACCENT)
        style.map('TRadiobutton',
                  background=[('active', BG_HOVER)],
                  indicatorbackground=[('selected', ACCENT)])

        style.configure('TNotebook', background=BG, bordercolor=BORDER)
        style.configure('TNotebook.Tab', background=BG_LIGHT, foreground=FG_DIM,
                        padding=(12, 4), bordercolor=BORDER)
        style.map('TNotebook.Tab',
                  background=[('selected', BG_PANEL), ('active', BG_HOVER)],
                  foreground=[('selected', FG_BRIGHT), ('active', FG)])

        style.configure('Treeview',
                        background=TREE_BG, foreground=TREE_FG,
                        fieldbackground=TREE_BG, bordercolor=BORDER,
                        font=('Helvetica', 10))
        style.configure('Treeview.Heading',
                        background=BG_PANEL, foreground=ACCENT,
                        bordercolor=BORDER, relief='flat',
                        font=('Helvetica', 10, 'bold'))
        style.map('Treeview',
                  background=[('selected', TREE_SELECTED_BG)],
                  foreground=[('selected', TREE_SELECTED_FG)])
        style.map('Treeview.Heading',
                  background=[('active', BG_HOVER)])

        style.configure('Vertical.TScrollbar',
                        background=SCROLL_BG, troughcolor=SCROLL_TROUGH,
                        bordercolor=BORDER, arrowcolor=FG_DIM)
        style.configure('Horizontal.TScrollbar',
                        background=SCROLL_BG, troughcolor=SCROLL_TROUGH,
                        bordercolor=BORDER, arrowcolor=FG_DIM)

        style.configure('TSeparator', background=BORDER)

        style.configure('Horizontal.TProgressbar',
                        troughcolor=BG_LIGHT, background=ACCENT,
                        bordercolor=BORDER, lightcolor=ACCENT,
                        darkcolor=ACCENT)

        style.configure('Horizontal.TScale',
                        background=BG, troughcolor=BG_INPUT,
                        bordercolor=BORDER, lightcolor=ACCENT,
                        darkcolor=ACCENT)

    # ---- Custom named styles (shared for both sv_ttk and clam) ----

    style.configure('Title.TLabel', font=('Helvetica', 12, 'bold'),
                    foreground=ACCENT)
    style.configure('Header.TLabel', font=('Helvetica', 10, 'bold'),
                    foreground=FG_BRIGHT)


def get_text_config(role: str = 'default') -> dict:
    """Return config dict for tk.Text widget constructor.

    Args:
        role: One of 'default', 'help', 'results', 'sequence', 'homology'.
    """
    config = {
        'bg': BG_PANEL,
        'fg': FG,
        'insertbackground': FG,
        'selectbackground': BG_SELECTED,
        'selectforeground': FG_BRIGHT,
        'relief': tk.FLAT,
        'borderwidth': 0,
        'highlightthickness': 0,
    }

    if role == 'help':
        config['bg'] = BG_LIGHT
    elif role in ('sequence', 'homology'):
        config['bg'] = SEQ_BG
        config['fg'] = FG_DIM

    return config


def get_help_tags() -> dict:
    """Return tag configuration dict for the help text widget."""
    return {
        'h1': dict(font=('Helvetica', 13, 'bold'), foreground=ACCENT,
                    spacing1=8, spacing3=4),
        'h2': dict(font=('Helvetica', 11, 'bold'), foreground=ACCENT_TEAL,
                    spacing1=6, spacing3=3),
        'body': dict(font=('Helvetica', 10), foreground=FG),
        'tip': dict(font=('Helvetica', 10, 'italic'), foreground=ACCENT_TEAL,
                    lmargin1=15, lmargin2=15),
        'warning': dict(font=('Helvetica', 10, 'bold'), foreground=ACCENT_RED),
        'step': dict(font=('Helvetica', 10, 'bold'), foreground=ACCENT,
                     lmargin1=10, lmargin2=25),
        'bullet': dict(font=('Helvetica', 10), foreground=FG,
                       lmargin1=20, lmargin2=30),
        'param': dict(font=('Courier', 10), foreground=ACCENT_PURPLE,
                      lmargin1=20, lmargin2=30),
        'highlight': dict(font=('Helvetica', 10, 'bold'), foreground=ACCENT_YELLOW),
        'current_step': dict(font=('Helvetica', 10, 'bold'),
                             foreground=FG_BRIGHT, background=ACCENT,
                             lmargin1=5, lmargin2=5, rmargin=5,
                             spacing1=3, spacing3=3),
    }


def get_results_tags() -> dict:
    """Return tag configuration dict for the results text widget."""
    return {
        'header': dict(font=('Helvetica', 10, 'bold'), foreground=FG_BRIGHT),
        'success': dict(foreground=ACCENT_GREEN),
        'error': dict(foreground=ACCENT_RED),
        'warning': dict(foreground=ACCENT_YELLOW),
        'sequence': dict(font=('Courier', 9), foreground=ACCENT),
        'specific': dict(foreground=ACCENT_GREEN, font=('Helvetica', 9, 'bold')),
        'not_specific': dict(foreground=ACCENT_RED, font=('Helvetica', 9, 'bold')),
        'pseudogene': dict(foreground=ACCENT_RED, font=('Helvetica', 10, 'bold')),
        'inconclusive': dict(foreground=ACCENT_ORANGE, font=('Helvetica', 9)),
        'cont_success': dict(foreground=ACCENT_GREEN, font=('Courier', 10, 'bold')),
        'cont_fail': dict(foreground=ACCENT_RED, font=('Courier', 10, 'bold')),
        'cont_skip': dict(foreground=ACCENT_YELLOW, font=('Courier', 10, 'bold')),
        'cont_exhausted': dict(foreground=ACCENT_ORANGE, font=('Courier', 9, 'italic')),
    }

#!/usr/bin/env python3
"""
Primer Designer - A tool for designing PCR/qPCR primers for genetic variants.

This application provides a graphical interface for:
- Loading and validating variant data from CSV/XLSX files
- Validating HGVS notation and MANE transcripts
- Designing PCR and qPCR primers using Primer3
- Checking primer specificity
- Filtering based on population variants
- Generating HTML reports

Usage:
    python main.py

Requirements:
    - Python 3.8+
    - tkinter (usually included with Python)
    - primer3-py
    - openpyxl (for XLSX support)
"""

import sys
import os
import subprocess
import platform

# Add the src directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _try_install_tkinter_linux() -> bool:
    """Try to auto-install python3-tk via system package manager."""
    pkg_managers = [
        (["apt", "--version"], ["sudo", "apt", "install", "-y", "python3-tk"]),
        (["dnf", "--version"], ["sudo", "dnf", "install", "-y", "python3-tkinter"]),
        (["yum", "--version"], ["sudo", "yum", "install", "-y", "python3-tkinter"]),
        (["pacman", "--version"], ["sudo", "pacman", "-S", "--noconfirm", "tk"]),
    ]
    for check_cmd, install_cmd in pkg_managers:
        try:
            subprocess.run(check_cmd, capture_output=True, check=True)
            print(f"  Running: {' '.join(install_cmd)}")
            result = subprocess.run(install_cmd, text=True)
            return result.returncode == 0
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
    return False


def _print_tkinter_instructions():
    """Print manual tkinter install instructions."""
    print("\nTo install tkinter on Linux:")
    print("  sudo apt-get install python3-tk    # Debian/Ubuntu")
    print("  sudo dnf install python3-tkinter   # Fedora")
    print("  sudo yum install python3-tkinter   # CentOS/RHEL")


def check_dependencies():
    """Check that all required dependencies are available.

    Must be called BEFORE importing application modules that depend on
    tkinter/primer3 so that the user gets a clear error message instead
    of a raw ImportError traceback.
    """
    missing = []

    # Check tkinter
    try:
        import tkinter  # noqa: F401
    except ImportError:
        missing.append("tkinter")

    # Check primer3
    try:
        import primer3  # noqa: F401
    except ImportError:
        missing.append("primer3-py (pip install primer3-py)")

    # Check openpyxl (optional but recommended)
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        print("Warning: openpyxl not installed. XLSX file support will be disabled.")
        print("Install with: pip install openpyxl")

    if missing:
        print("Missing required dependencies:")
        for dep in missing:
            print(f"  - {dep}")

        if "tkinter" in missing:
            os_type = platform.system().lower()
            if os_type == "linux":
                print("\n  Attempting to auto-install tkinter...")
                installed = _try_install_tkinter_linux()
                if installed:
                    # Verify
                    try:
                        import tkinter  # noqa: F401
                        missing.remove("tkinter")
                        print("  tkinter: installed successfully!")
                    except ImportError:
                        print("  Auto-install failed.")
                        _print_tkinter_instructions()
                else:
                    _print_tkinter_instructions()
            elif os_type == "darwin":
                print("\nTo install tkinter on macOS:")
                print("  brew install python-tk")

        if missing:
            print("\nPlease install missing dependencies and try again.")
            sys.exit(1)


def clear_api_cache():
    """Clear thread-safe API caches and config singleton from previous sessions.

    This ensures a fresh start each time the application is launched,
    preventing stale data from a previous run.
    """
    try:
        from src.utils.thread_safe import ThreadSafeCache
        ThreadSafeCache.clear_all_instances()
    except Exception:
        pass

    try:
        from src.utils import config as config_module
        config_module._config = None
    except Exception:
        pass


def main():
    """Main entry point for the Primer Designer application."""
    # Check dependencies BEFORE importing any application modules
    check_dependencies()

    from src.utils.logger import setup_logger

    # Initialize logging
    logger = setup_logger()
    logger.info("Starting Primer Designer application")

    # Clear caches from previous sessions
    clear_api_cache()

    try:
        # Import GUI only after dependency check passes
        from src.gui.main_window import MainWindow

        # Create and run the main window
        app = MainWindow()
        app.run()

    except ImportError as e:
        print(f"Missing required dependency: {e}")
        print("Please run: pip install -r requirements.txt")
        sys.exit(1)

    except Exception as e:
        logger.exception(f"Application error: {e}")
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

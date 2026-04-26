#!/usr/bin/env python3
"""
Installation script for Primer Designer.
Automatically detects OS and sets up the environment.

Usage:
    python install.py

This script will:
1. Detect the operating system
2. Create a virtual environment
3. Install required Python packages
4. Download necessary data files (MANE database)
5. Check for external tools (BLAST+)
"""

import atexit
import os
import signal
import sys
import subprocess
import platform
import urllib.request
import shutil
from pathlib import Path
import json
import gzip
import re
import tarfile

# Track temporary/partial files for cleanup on interrupt
_cleanup_files = []


def _cleanup_handler(signum=None, frame=None):
    """Clean up partial downloads on Ctrl+C or unexpected exit."""
    for f in _cleanup_files:
        try:
            p = Path(f)
            if p.exists():
                p.unlink()
                print(f"Cleaned up partial file: {p}")
        except Exception:
            pass
    if signum is not None:
        print("\nInstallation interrupted by user.")
        sys.exit(1)


signal.signal(signal.SIGINT, _cleanup_handler)
atexit.register(_cleanup_handler)


class Installer:
    """Handles installation of Primer Designer."""

    def __init__(self):
        self.base_dir = Path(__file__).parent.absolute()
        self.venv_dir = self.base_dir / "venv"
        self.data_dir = self.base_dir / "data"
        self.tools_dir = self.base_dir / "tools"
        self.os_type = platform.system().lower()

        print("=" * 60)
        print("Primer Designer Installer")
        print("=" * 60)
        print(f"Operating System: {platform.system()} {platform.release()}")
        print(f"Python Version: {sys.version}")
        print(f"Installation Directory: {self.base_dir}")
        print("=" * 60)

    def run(self):
        """Run the complete installation process."""
        try:
            self.create_directories()
            self.create_virtualenv()
            self.install_dependencies()
            self.download_data()
            self.check_external_tools()
            self.setup_config()
            self.create_launcher()
            self.print_success()
        except Exception as e:
            print(f"\nError during installation: {e}")
            print("Please check the error message and try again.")
            sys.exit(1)

    def create_directories(self):
        """Create necessary directories."""
        print("\n[1/7] Creating directories...")

        directories = [
            self.data_dir,
            self.data_dir / "cache",
            self.data_dir / "fasta_sequences",
            self.data_dir / "genome" / "GRCh38",
            self.tools_dir / "blast" / "bin",
            self.base_dir / "logs"
        ]

        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
            print(f"  Created: {directory}")

    def create_virtualenv(self):
        """Create a virtual environment."""
        print("\n[2/7] Creating virtual environment...")

        # Check if venv module is available before attempting to create
        result = subprocess.run(
            [sys.executable, "-c", "import venv; import ensurepip"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print("  Python venv module is not available. Attempting auto-install...")
            if self._try_install_venv():
                # Verify it worked
                result = subprocess.run(
                    [sys.executable, "-c", "import venv; import ensurepip"],
                    capture_output=True, text=True
                )
                if result.returncode != 0:
                    print("  ERROR: Auto-install failed. Please install manually:")
                    self._print_venv_install_instructions()
                    sys.exit(1)
                print("  venv module installed successfully!")
            else:
                self._print_venv_install_instructions()
                sys.exit(1)

        if self.venv_dir.exists():
            print(f"  Virtual environment already exists at {self.venv_dir}")
            response = input("  Recreate? (y/N): ").strip().lower()
            if response == 'y':
                shutil.rmtree(self.venv_dir)
            else:
                return

        subprocess.run(
            [sys.executable, "-m", "venv", str(self.venv_dir)],
            check=True
        )
        print(f"  Created virtual environment at {self.venv_dir}")

    def _try_install_venv(self) -> bool:
        """Try to automatically install python3-venv package."""
        py_ver = f"{sys.version_info.major}.{sys.version_info.minor}"

        if self.os_type == "linux":
            # Detect package manager and try to install
            pkg_managers = [
                # (check_cmd, install_cmd)
                (["apt", "--version"], ["sudo", "apt", "install", "-y", f"python{py_ver}-venv"]),
                (["dnf", "--version"], ["sudo", "dnf", "install", "-y", "python3-virtualenv"]),
                (["yum", "--version"], ["sudo", "yum", "install", "-y", "python3-virtualenv"]),
                (["pacman", "--version"], ["sudo", "pacman", "-S", "--noconfirm", "python-virtualenv"]),
            ]
            for check_cmd, install_cmd in pkg_managers:
                try:
                    subprocess.run(check_cmd, capture_output=True, check=True)
                    print(f"  Running: {' '.join(install_cmd)}")
                    install_result = subprocess.run(install_cmd, text=True)
                    return install_result.returncode == 0
                except (FileNotFoundError, subprocess.CalledProcessError):
                    continue

        elif self.os_type == "darwin":
            # On macOS, venv should come with Python. If missing, can't auto-install.
            print("  Cannot auto-install venv on macOS.")
            return False

        return False

    def _print_venv_install_instructions(self):
        """Print manual install instructions for venv."""
        py_ver = f"{sys.version_info.major}.{sys.version_info.minor}"
        print("  ERROR: Python venv module is not available!")
        if self.os_type == "linux":
            print(f"  Install manually with:")
            print(f"    sudo apt install python{py_ver}-venv    # Debian/Ubuntu")
            print(f"    sudo dnf install python3-virtualenv     # Fedora")
            print(f"    sudo yum install python3-virtualenv     # CentOS/RHEL")
            print(f"    sudo pacman -S python-virtualenv        # Arch Linux")
        elif self.os_type == "darwin":
            print("  Reinstall Python from python.org or via: brew install python")
        else:
            print("  Reinstall Python with the 'venv' component enabled.")
        print(f"\n  After installing, re-run: python3 install.py")

    def _try_install_blast(self) -> bool:
        """Try to automatically install BLAST+ package."""
        if self.os_type == "darwin":
            # macOS: use Homebrew
            if not shutil.which('brew'):
                print("  Homebrew not found. Cannot auto-install BLAST+.")
                return False
            try:
                print("  Running: brew install blast")
                result = subprocess.run(
                    ['brew', 'install', 'blast'],
                    text=True, timeout=600
                )
                return result.returncode == 0
            except (subprocess.SubprocessError, FileNotFoundError):
                return False

        elif self.os_type == "linux":
            pkg_managers = [
                (["apt-get", "--version"], ["sudo", "apt-get", "install", "-y", "ncbi-blast+"]),
                (["yum", "--version"], ["sudo", "yum", "install", "-y", "ncbi-blast+"]),
                (["dnf", "--version"], ["sudo", "dnf", "install", "-y", "ncbi-blast+"]),
                (["pacman", "--version"], ["sudo", "pacman", "-S", "--noconfirm", "blast+"]),
            ]
            for check_cmd, install_cmd in pkg_managers:
                try:
                    subprocess.run(check_cmd, capture_output=True, check=True)
                    print(f"  Running: {' '.join(install_cmd)}")
                    result = subprocess.run(install_cmd, text=True, timeout=600)
                    return result.returncode == 0
                except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
                    continue

        elif self.os_type == "windows":
            return self._try_install_blast_windows()

        return False

    def _try_install_blast_windows(self) -> bool:
        """Download and extract BLAST+ from NCBI FTP on Windows."""
        tools_dir = self.tools_dir / "blast"
        tools_dir.mkdir(parents=True, exist_ok=True)

        # Check if already extracted in tools/
        blastn_exe = tools_dir / "bin" / "blastn.exe"
        if blastn_exe.exists():
            print(f"  BLAST+ already present at {blastn_exe}")
            self._update_config_blast_path(tools_dir / "bin")
            return True

        # Discover latest version from NCBI FTP
        ftp_url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        try:
            print("  Discovering latest BLAST+ version from NCBI...")
            with urllib.request.urlopen(ftp_url, timeout=30) as response:
                html = response.read().decode('utf-8')

            # Find the win64 tar.gz filename
            pattern = r'(ncbi-blast-[\d.+]+-x64-win64\.tar\.gz)'
            match = re.search(pattern, html)
            if not match:
                print("  Could not find Windows BLAST+ archive on NCBI FTP.")
                return False

            archive_name = match.group(1)
            download_url = ftp_url + archive_name
            archive_path = tools_dir / archive_name

            print(f"  Downloading {archive_name} (~137 MB)...")
            _cleanup_files.append(str(archive_path))
            urllib.request.urlretrieve(
                download_url, archive_path, self._download_progress
            )
            print()  # newline after progress bar

            # Extract tar.gz
            print("  Extracting BLAST+...")
            with tarfile.open(archive_path, 'r:gz') as tar:
                tar.extractall(path=tools_dir)

            # Remove archive
            archive_path.unlink(missing_ok=True)
            _cleanup_files.clear()

            # The archive extracts to ncbi-blast-X.X.X+/bin/
            # Move contents so tools/blast/bin/ has the executables
            extracted_dirs = [
                d for d in tools_dir.iterdir()
                if d.is_dir() and d.name.startswith('ncbi-blast-')
            ]
            if extracted_dirs:
                extracted = extracted_dirs[0]
                extracted_bin = extracted / "bin"
                target_bin = tools_dir / "bin"
                if target_bin.exists():
                    shutil.rmtree(target_bin)
                shutil.move(str(extracted_bin), str(target_bin))
                # Clean up extracted directory
                shutil.rmtree(extracted, ignore_errors=True)

            if blastn_exe.exists():
                self._update_config_blast_path(tools_dir / "bin")
                print(f"  BLAST+ extracted to {tools_dir / 'bin'}")
                return True
            else:
                print("  Extraction succeeded but blastn.exe not found.")
                return False

        except Exception as e:
            print(f"  Failed to download BLAST+: {e}")
            # Clean up partial downloads
            for f in tools_dir.iterdir():
                if f.suffix in ('.gz', '.tar'):
                    f.unlink(missing_ok=True)
            return False

    def _update_config_blast_path(self, bin_dir: Path):
        """Update config.json with the local BLAST+ paths."""
        config_path = self.base_dir / "config.json"
        try:
            config = json.loads(config_path.read_text()) if config_path.exists() else {}
            config.setdefault("blast", {})
            blastn = str(bin_dir / "blastn.exe") if self.os_type == "windows" else str(bin_dir / "blastn")
            makeblastdb = str(bin_dir / "makeblastdb.exe") if self.os_type == "windows" else str(bin_dir / "makeblastdb")
            config["blast"]["blastn_path"] = blastn
            config["blast"]["makeblastdb_path"] = makeblastdb
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2)
                f.write('\n')
            print(f"  Updated config.json with BLAST+ paths.")
        except Exception as e:
            print(f"  Warning: Could not update config.json: {e}")
            print("  You may need to set blast paths manually in config.json.")

    def _try_install_tkinter(self) -> bool:
        """Try to automatically install tkinter system package."""
        if self.os_type == "linux":
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
                    install_result = subprocess.run(install_cmd, text=True)
                    return install_result.returncode == 0
                except (FileNotFoundError, subprocess.CalledProcessError):
                    continue
        return False

    def _print_tkinter_install_instructions(self):
        """Print manual install instructions for tkinter."""
        print("  WARNING: tkinter is required for the graphical interface!")
        if self.os_type == "linux":
            print("  Install manually with:")
            print("    sudo apt-get install python3-tk    # Debian/Ubuntu")
            print("    sudo dnf install python3-tkinter   # Fedora")
            print("    sudo yum install python3-tkinter   # CentOS/RHEL")
            print("    sudo pacman -S tk                  # Arch Linux")
        elif self.os_type == "darwin":
            print("  Install with: brew install python-tk")
        else:
            print("  Reinstall Python with tkinter support enabled.")

    def get_pip_path(self) -> str:
        """Get the path to pip in the virtual environment."""
        if self.os_type == "windows":
            return str(self.venv_dir / "Scripts" / "pip.exe")
        else:
            return str(self.venv_dir / "bin" / "pip")

    def get_python_path(self) -> str:
        """Get the path to Python in the virtual environment."""
        if self.os_type == "windows":
            return str(self.venv_dir / "Scripts" / "python.exe")
        else:
            return str(self.venv_dir / "bin" / "python")

    def install_dependencies(self):
        """Install Python dependencies."""
        print("\n[3/7] Installing Python dependencies...")

        pip_path = self.get_pip_path()
        requirements_file = self.base_dir / "requirements.txt"

        # Upgrade pip first
        print("  Upgrading pip...")
        try:
            subprocess.run(
                [pip_path, "install", "--upgrade", "pip"],
                check=True,
                capture_output=True
            )
        except subprocess.CalledProcessError as e:
            print(f"  Warning: Pip upgrade failed ({e}). Continuing with existing pip version.")

        # Install requirements
        print("  Installing packages from requirements.txt...")
        result = subprocess.run(
            [pip_path, "install", "-r", str(requirements_file)],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"  Warning: Some packages may have failed to install:")
            print(result.stderr)
        else:
            print("  All packages installed successfully")

        # Verify critical packages
        print("  Verifying critical packages...")
        python_path = self.get_python_path()

        packages_to_verify = ["primer3", "openpyxl", "sv_ttk", "requests"]
        for package in packages_to_verify:
            try:
                subprocess.run(
                    [python_path, "-c", f"import {package}"],
                    check=True,
                    capture_output=True
                )
                print(f"    {package}: OK")
            except subprocess.CalledProcessError:
                print(f"    {package}: FAILED")

    def _discover_mane_url(self) -> str:
        """Discover the latest MANE summary file URL from NCBI FTP."""
        mane_dir_url = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
        try:
            print("  Discovering latest MANE version...")
            with urllib.request.urlopen(mane_dir_url, timeout=15) as response:
                html = response.read().decode('utf-8')

            pattern = r'MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz'
            matches = re.findall(pattern, html)
            if matches:
                versions = sorted(matches, key=lambda v: [int(x) for x in v.split('.')], reverse=True)
                latest = versions[0]
                url = f"{mane_dir_url}MANE.GRCh38.v{latest}.summary.txt.gz"
                print(f"  Found latest MANE version: v{latest}")
                return url
        except Exception as e:
            print(f"  Warning: Could not auto-discover MANE version: {e}")

        # Fallback to known working URL
        fallback = mane_dir_url + "MANE.GRCh38.v1.5.summary.txt.gz"
        print(f"  Using fallback MANE URL")
        return fallback

    def download_data(self):
        """Download necessary data files."""
        print("\n[4/7] Downloading data files...")

        mane_file = self.data_dir / "cache" / "MANE_summary.txt.gz"
        mane_extracted = self.data_dir / "cache" / "MANE_summary.txt"

        if mane_extracted.exists():
            print(f"  MANE data already exists at {mane_extracted}")
        else:
            mane_url = self._discover_mane_url()
            print(f"  Downloading MANE data from NCBI...")
            try:
                _cleanup_files.append(str(mane_file))
                urllib.request.urlretrieve(mane_url, mane_file, self._download_progress)
                print("\n  Extracting...")

                with gzip.open(mane_file, 'rb') as f_in:
                    with open(mane_extracted, 'wb') as f_out:
                        f_out.write(f_in.read())

                mane_file.unlink()  # Remove compressed file
                _cleanup_files.clear()
                print(f"  MANE data saved to {mane_extracted}")

            except Exception as e:
                print(f"  Warning: Could not download MANE data: {e}")
                print("  The application will try to download it on first run.")

    def _download_progress(self, block_num, block_size, total_size):
        """Show download progress."""
        if total_size > 0:
            percent = min(100, block_num * block_size * 100 // total_size)
            bar_length = 40
            filled = int(bar_length * percent / 100)
            bar = "=" * filled + "-" * (bar_length - filled)
            print(f"\r  [{bar}] {percent}%", end="", flush=True)

    def check_external_tools(self):
        """Check for external tools like BLAST+ and tkinter."""
        print("\n[5/7] Checking external tools and dependencies...")

        # Check for tkinter (critical for GUI)
        python_path = self.get_python_path()
        try:
            result = subprocess.run(
                [python_path, "-c", "import tkinter"],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                print("  tkinter: OK")
            else:
                print("  tkinter: NOT FOUND — attempting auto-install...")
                if self._try_install_tkinter():
                    # Verify it worked
                    result = subprocess.run(
                        [python_path, "-c", "import tkinter"],
                        capture_output=True, text=True
                    )
                    if result.returncode == 0:
                        print("  tkinter: installed successfully!")
                    else:
                        print("  tkinter: auto-install failed.")
                        self._print_tkinter_install_instructions()
                else:
                    self._print_tkinter_install_instructions()
        except Exception:
            print("  tkinter: Could not check (will be verified at runtime)")

        # Check for BLAST+
        blastn_found = False
        try:
            result = subprocess.run(
                ["blastn", "-version"],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                blastn_found = True
                version = result.stdout.strip().split('\n')[0]
                print(f"  BLAST+ (blastn): Found ({version})")
        except FileNotFoundError:
            pass

        if not blastn_found:
            # Also check local tools/ directory (Windows auto-install location)
            local_blastn = self.tools_dir / "blast" / "bin" / (
                "blastn.exe" if self.os_type == "windows" else "blastn"
            )
            if local_blastn.exists():
                try:
                    result = subprocess.run(
                        [str(local_blastn), "-version"],
                        capture_output=True, text=True
                    )
                    if result.returncode == 0:
                        blastn_found = True
                        version = result.stdout.strip().split('\n')[0]
                        print(f"  BLAST+ (blastn): Found at {local_blastn} ({version})")
                except (FileNotFoundError, subprocess.SubprocessError):
                    pass

        if not blastn_found:
            print("  BLAST+ (blastn): Not found — attempting auto-install...")
            if self._try_install_blast():
                # Verify: check PATH first, then local tools/ directory
                for blastn_cmd in ["blastn", str(self.tools_dir / "blast" / "bin" / (
                    "blastn.exe" if self.os_type == "windows" else "blastn"
                ))]:
                    try:
                        result = subprocess.run(
                            [blastn_cmd, "-version"],
                            capture_output=True, text=True
                        )
                        if result.returncode == 0:
                            version = result.stdout.strip().split('\n')[0]
                            print(f"  BLAST+ (blastn): Installed successfully ({version})")
                            blastn_found = True
                            break
                    except (FileNotFoundError, subprocess.SubprocessError):
                        continue

            if not blastn_found:
                print("  BLAST+ (blastn): Auto-install failed.")
                print("  Note: BLAST+ is recommended for homology analysis.")
                print("  Install BLAST+ manually for full functionality:")
                if self.os_type == "darwin":
                    print("    brew install blast")
                elif self.os_type == "linux":
                    print("    sudo apt-get install ncbi-blast+  # Debian/Ubuntu")
                    print("    sudo yum install ncbi-blast+      # CentOS/RHEL")
                else:
                    print("    Download from: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html")

    def setup_config(self):
        """Create/update config.json with detected paths."""
        print("\n[6/7] Setting up configuration...")

        config_path = self.base_dir / "config.json"
        blast_bin = self.tools_dir / "blast" / "bin"

        config = {
            "blast": {
                "blastn_path": str(blast_bin / ("blastn.exe" if self.os_type == "windows" else "blastn")),
                "makeblastdb_path": str(blast_bin / ("makeblastdb.exe" if self.os_type == "windows" else "makeblastdb")),
                "enabled": True
            },
            "api": {
                "ncbi_email": "",
                "max_retries": 3,
                "retry_delay": 2.0
            },
            "ensembl": {
                "delay_seconds": 3.0,
                "max_retries": 2
            }
        }

        if config_path.exists():
            with open(config_path, 'r') as f:
                existing = json.load(f)
            config.update(existing)

        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
            f.write('\n')

        print(f"  Configuration saved to {config_path}")

    def create_launcher(self):
        """Create platform-specific launcher scripts."""
        print("\n[7/7] Creating launcher scripts...")

        if self.os_type == "windows":
            self._create_windows_launcher()
        else:
            self._create_unix_launcher()

    def _create_windows_launcher(self):
        """Create Windows batch launcher with relative paths."""
        launcher_path = self.base_dir / "run_primer_designer.bat"

        content = '''@echo off
REM Launcher for Primer Designer
REM Auto-generated by install.py

set "SCRIPT_DIR=%~dp0"
set "PYTHON_PATH=%SCRIPT_DIR%venv\\Scripts\\python.exe"

if not exist "%PYTHON_PATH%" (
    echo Error: Virtual environment not found at %SCRIPT_DIR%venv
    echo Please run: python install.py
    pause
    exit /b 1
)

echo Starting Primer Designer...
"%PYTHON_PATH%" "%SCRIPT_DIR%main.py"
pause
'''
        launcher_path.write_text(content, encoding='utf-8')
        print(f"  Created: {launcher_path}")

    def _create_unix_launcher(self):
        """Create Unix shell launcher with relative paths."""
        launcher_path = self.base_dir / "run_primer_designer.sh"

        content = '''#!/bin/bash
# Launcher for Primer Designer
# Auto-generated by install.py

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_PATH="$SCRIPT_DIR/venv/bin/python"

if [ ! -f "$PYTHON_PATH" ]; then
    echo "Error: Virtual environment not found at $SCRIPT_DIR/venv"
    echo "Please run: python install.py"
    exit 1
fi

echo "Starting Primer Designer..."
"$PYTHON_PATH" "$SCRIPT_DIR/main.py"
'''
        launcher_path.write_text(content, encoding='utf-8')
        launcher_path.chmod(0o755)
        print(f"  Created: {launcher_path}")

    def print_success(self):
        """Print success message with instructions."""
        print("\n" + "=" * 60)
        print("Installation completed successfully!")
        print("=" * 60)

        print("\nTo run Primer Designer:")
        if self.os_type == "windows":
            print(f"  Double-click: run_primer_designer.bat")
            print(f"  Or from command line: {self.get_python_path()} main.py")
        else:
            print(f"  Run: ./run_primer_designer.sh")
            print(f"  Or: {self.get_python_path()} main.py")

        print("\nFor more information, see README.md")
        print("=" * 60)


def main():
    """Main entry point."""
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    main()
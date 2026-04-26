"""
BWA auto-detection, auto-installation, and genome index management.
Ensures BWA is available and genome reference index is built before
specificity checking runs.

All operations are idempotent — detect before install/download.
"""

import subprocess
import platform
import shutil
import gzip
import urllib.request
from pathlib import Path
from typing import Tuple, Optional, Callable

from ..utils.logger import LoggerMixin
from ..utils.config import get_config


class BWASetupManager(LoggerMixin):
    """
    Manages BWA installation and genome index setup.
    All methods check if the resource already exists before acting.
    """

    # NCBI FTP — GRCh38 analysis set (no alt contigs, ~800 MB compressed)
    GRCH38_FASTA_URL = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
        "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
        "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    )

    BWA_INDEX_EXTENSIONS = ['.amb', '.ann', '.bwt', '.pac', '.sa']

    def __init__(self):
        self.config = get_config()
        self.data_dir = Path(self.config.data_cache_dir).parent  # data/
        self.genome_dir = self.data_dir / "genome"

    # ------------------------------------------------------------------
    #  BWA detection / installation
    # ------------------------------------------------------------------

    def detect_bwa(self) -> Tuple[bool, str]:
        """
        Check if BWA is installed.
        Returns (is_available, bwa_path).
        """
        # 1. Configured path
        configured = self.config.bwa.bwa_path
        if configured and self._test_bwa(configured):
            return True, configured

        # 2. System PATH
        found = shutil.which('bwa')
        if found and self._test_bwa(found):
            return True, found

        # 3. Common locations (platform-specific)
        os_type = platform.system().lower()
        if os_type == 'windows':
            common_paths = []
            # BWA on Windows is typically installed via conda or manual build
            for prog_dir in [r'C:\Program Files\bwa', r'C:\Program Files (x86)\bwa']:
                prog_path = Path(prog_dir)
                if prog_path.exists():
                    candidate = prog_path / 'bwa.exe'
                    if candidate.exists():
                        common_paths.append(str(candidate))
        else:
            common_paths = ['/opt/homebrew/bin/bwa', '/usr/local/bin/bwa', '/usr/bin/bwa']

        for path in common_paths:
            if self._test_bwa(path):
                return True, path

        return False, ""

    def _test_bwa(self, path: str) -> bool:
        """Test if a given BWA path works."""
        try:
            result = subprocess.run(
                [path, 'version'],
                capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return False

    def install_bwa(self, progress_callback: Optional[Callable] = None) -> Tuple[bool, str]:
        """
        Auto-install BWA if not already present.
        Returns (success, bwa_path_or_error_message).
        """
        available, path = self.detect_bwa()
        if available:
            self.log_info(f"BWA already installed at {path}")
            return True, path

        os_type = platform.system().lower()

        if os_type == 'darwin':
            return self._install_bwa_macos(progress_callback)
        elif os_type == 'linux':
            return self._install_bwa_linux(progress_callback)
        elif os_type == 'windows':
            return False, (
                "Auto-install not supported on Windows. "
                "BWA can be installed via conda (conda install -c bioconda bwa) "
                "or downloaded from https://github.com/lh3/bwa/releases. "
                "After installation, set 'bwa_path' in config.json."
            )
        else:
            return False, (
                "Auto-install not supported on this platform. "
                "Please install BWA manually."
            )

    def _install_bwa_macos(self, cb) -> Tuple[bool, str]:
        """Install BWA via Homebrew on macOS."""
        if not shutil.which('brew'):
            return False, (
                "Homebrew not found. Install Homebrew first "
                "(https://brew.sh) or install BWA manually."
            )

        if cb:
            cb("Installing BWA via Homebrew...")

        try:
            result = subprocess.run(
                ['brew', 'install', 'bwa'],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode == 0:
                ok, path = self.detect_bwa()
                if ok:
                    self.log_info(f"BWA installed at {path}")
                    return True, path
            return False, f"brew install bwa failed: {result.stderr[:300]}"
        except subprocess.TimeoutExpired:
            return False, "BWA installation timed out (>5 min)"
        except Exception as e:
            return False, str(e)

    def _install_bwa_linux(self, cb) -> Tuple[bool, str]:
        """Install BWA via apt / yum on Linux."""
        if cb:
            cb("Installing BWA via package manager...")

        # Try apt-get (Debian / Ubuntu)
        for cmd in [
            ['sudo', 'apt-get', 'install', '-y', 'bwa'],
            ['sudo', 'yum', 'install', '-y', 'bwa'],
        ]:
            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=300
                )
                if result.returncode == 0:
                    ok, path = self.detect_bwa()
                    if ok:
                        self.log_info(f"BWA installed at {path}")
                        return True, path
            except (subprocess.SubprocessError, FileNotFoundError):
                continue

        return False, "Could not install BWA via apt-get or yum. Please install manually."

    # ------------------------------------------------------------------
    #  Genome index detection / building
    # ------------------------------------------------------------------

    def detect_genome_index(self) -> Tuple[bool, str]:
        """
        Check if BWA genome index exists.
        Returns (exists, index_base_path).
        """
        # 1. Configured path
        configured = self.config.bwa.genome_index_path
        if configured and self._verify_index(configured):
            return True, configured

        # 2. Default location
        default_base = str(self.genome_dir / "GRCh38" / "genome")
        if self._verify_index(default_base):
            return True, default_base

        return False, ""

    def _verify_index(self, base_path: str) -> bool:
        """Verify all five BWA index files exist for a given base path."""
        return all(
            Path(base_path + ext).exists()
            for ext in self.BWA_INDEX_EXTENSIONS
        )

    def setup_genome_index(
        self, progress_callback: Optional[Callable] = None
    ) -> Tuple[bool, str]:
        """
        Ensure genome index is available.
        Downloads reference FASTA and builds BWA index if needed.
        Returns (success, index_base_path_or_error).
        """
        exists, path = self.detect_genome_index()
        if exists:
            self.log_info(f"Genome index already exists at {path}")
            return True, path

        index_dir = self.genome_dir / "GRCh38"
        index_dir.mkdir(parents=True, exist_ok=True)

        fasta_path = index_dir / "genome.fa"
        index_base = str(index_dir / "genome")

        # Download FASTA if not present
        if not fasta_path.exists():
            gz_path = index_dir / "genome.fa.gz"
            if not gz_path.exists():
                if progress_callback:
                    progress_callback("Downloading GRCh38 reference genome (~800 MB)...")
                ok = self._download_file(str(gz_path), progress_callback)
                if not ok:
                    return False, "Failed to download reference genome"

            # Decompress
            if progress_callback:
                progress_callback("Decompressing reference genome...")
            try:
                self._decompress_gz(str(gz_path), str(fasta_path))
                gz_path.unlink(missing_ok=True)
            except Exception as e:
                return False, f"Decompression failed: {e}"

        # Build BWA index
        bwa_ok, bwa_path = self.detect_bwa()
        if not bwa_ok:
            return False, "BWA not installed — cannot build genome index"

        if progress_callback:
            progress_callback("Building BWA index (this may take 30–60 minutes)...")

        try:
            result = subprocess.run(
                [bwa_path, 'index', '-p', index_base, str(fasta_path)],
                capture_output=True, text=True, timeout=7200  # 2 h
            )
            if result.returncode == 0 and self._verify_index(index_base):
                self.log_info(f"BWA index built at {index_base}")
                return True, index_base
            return False, f"BWA index failed: {result.stderr[:300]}"
        except subprocess.TimeoutExpired:
            return False, "BWA indexing timed out (>2 h)"
        except Exception as e:
            return False, str(e)

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    def _download_file(self, dest: str, progress_callback: Optional[Callable] = None) -> bool:
        """Download the reference genome with progress reporting."""
        try:
            def reporthook(block_num, block_size, total_size):
                if progress_callback and total_size > 0:
                    pct = min(100, block_num * block_size * 100 // total_size)
                    progress_callback(f"Downloading genome... {pct}%")

            urllib.request.urlretrieve(self.GRCH38_FASTA_URL, dest, reporthook)
            return True
        except Exception as e:
            self.log_error(f"Genome download failed: {e}")
            return False

    @staticmethod
    def _decompress_gz(gz_path: str, out_path: str):
        """Decompress a .gz file."""
        with gzip.open(gz_path, 'rb') as f_in:
            with open(out_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

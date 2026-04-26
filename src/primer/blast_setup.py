"""
BLAST+ auto-detection, auto-installation, and database management.
Ensures blastn is available and a BLAST nucleotide database is built
from the reference genome FASTA before homology analysis runs.

All operations are idempotent -- detect before install/download.
"""

import gzip
import platform
import shutil
import subprocess
import urllib.request
import urllib.error
from pathlib import Path
from typing import Tuple, Optional, Callable

from ..utils.logger import LoggerMixin
from ..utils.config import get_config


class BlastSetupManager(LoggerMixin):
    """
    Manages BLAST+ installation and genome BLAST database setup.
    All methods check if the resource already exists before acting.
    """

    # Minimum set of extensions produced by makeblastdb (nucl)
    BLAST_DB_EXTENSIONS = ['.nhr', '.nin', '.nsq']

    def __init__(self):
        self.config = get_config()
        self.data_dir = Path(self.config.data_cache_dir).parent  # data/
        self.genome_dir = self.data_dir / "genome"

    # ------------------------------------------------------------------
    #  blastn detection / installation
    # ------------------------------------------------------------------

    def detect_blastn(self) -> Tuple[bool, str]:
        """
        Check if blastn is installed.
        Returns (is_available, blastn_path).
        """
        # 1. Configured path
        configured = self.config.blast.blastn_path
        if configured and self._test_blast(configured):
            return True, configured

        # 2. System PATH
        found = shutil.which('blastn')
        if found and self._test_blast(found):
            return True, found

        # 3. Local tools/ directory (auto-installed on Windows)
        base_dir = Path(__file__).resolve().parent.parent.parent
        os_type = platform.system().lower()
        ext = ".exe" if os_type == "windows" else ""
        local_blastn = base_dir / "tools" / "blast" / "bin" / f"blastn{ext}"
        if self._test_blast(str(local_blastn)):
            return True, str(local_blastn)

        # 4. Common locations (platform-specific)
        if os_type == 'windows':
            common_paths = [
                r'C:\Program Files\NCBI\blast-LATEST+\bin\blastn.exe',
                r'C:\Program Files (x86)\NCBI\blast-LATEST+\bin\blastn.exe',
            ]
            # Also scan for versioned directories
            for prog_dir in [r'C:\Program Files\NCBI', r'C:\Program Files (x86)\NCBI']:
                prog_path = Path(prog_dir)
                if prog_path.exists():
                    for subdir in prog_path.iterdir():
                        if subdir.is_dir() and 'blast' in subdir.name.lower():
                            candidate = subdir / 'bin' / 'blastn.exe'
                            if candidate.exists():
                                common_paths.insert(0, str(candidate))
        else:
            common_paths = [
                '/opt/homebrew/bin/blastn',
                '/usr/local/bin/blastn',
                '/usr/bin/blastn',
            ]

        for path in common_paths:
            if self._test_blast(path):
                return True, path

        return False, ""

    def detect_makeblastdb(self) -> Tuple[bool, str]:
        """
        Check if makeblastdb is installed.
        Returns (is_available, makeblastdb_path).
        """
        configured = self.config.blast.makeblastdb_path
        if configured and self._test_makeblastdb(configured):
            return True, configured

        found = shutil.which('makeblastdb')
        if found and self._test_makeblastdb(found):
            return True, found

        # Local tools/ directory (auto-installed on Windows)
        base_dir = Path(__file__).resolve().parent.parent.parent
        os_type = platform.system().lower()
        ext = ".exe" if os_type == "windows" else ""
        local_mdb = base_dir / "tools" / "blast" / "bin" / f"makeblastdb{ext}"
        if self._test_makeblastdb(str(local_mdb)):
            return True, str(local_mdb)

        if os_type == 'windows':
            common_paths = [
                r'C:\Program Files\NCBI\blast-LATEST+\bin\makeblastdb.exe',
                r'C:\Program Files (x86)\NCBI\blast-LATEST+\bin\makeblastdb.exe',
            ]
            for prog_dir in [r'C:\Program Files\NCBI', r'C:\Program Files (x86)\NCBI']:
                prog_path = Path(prog_dir)
                if prog_path.exists():
                    for subdir in prog_path.iterdir():
                        if subdir.is_dir() and 'blast' in subdir.name.lower():
                            candidate = subdir / 'bin' / 'makeblastdb.exe'
                            if candidate.exists():
                                common_paths.insert(0, str(candidate))
        else:
            common_paths = [
                '/opt/homebrew/bin/makeblastdb',
                '/usr/local/bin/makeblastdb',
                '/usr/bin/makeblastdb',
            ]

        for path in common_paths:
            if self._test_makeblastdb(path):
                return True, path

        return False, ""

    def _test_blast(self, path: str) -> bool:
        """Test if a given blastn path works."""
        try:
            result = subprocess.run(
                [path, '-version'],
                capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0 and 'blastn' in result.stdout.lower()
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return False

    def _test_makeblastdb(self, path: str) -> bool:
        """Test if a given makeblastdb path works."""
        try:
            result = subprocess.run(
                [path, '-version'],
                capture_output=True, text=True, timeout=10
            )
            return result.returncode == 0 and 'makeblastdb' in result.stdout.lower()
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return False

    def install_blast(
        self, progress_callback: Optional[Callable] = None
    ) -> Tuple[bool, str]:
        """
        Auto-install BLAST+ if not already present.
        Returns (success, blastn_path_or_error_message).
        """
        available, path = self.detect_blastn()
        if available:
            self.log_info(f"BLAST+ already installed at {path}")
            return True, path

        os_type = platform.system().lower()

        if os_type == 'darwin':
            return self._install_blast_macos(progress_callback)
        elif os_type == 'linux':
            return self._install_blast_linux(progress_callback)
        elif os_type == 'windows':
            return self._install_blast_windows(progress_callback)
        else:
            return False, (
                "Auto-install not supported on this platform. "
                "Please install BLAST+ manually: "
                "https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html"
            )

    def _install_blast_macos(self, cb) -> Tuple[bool, str]:
        """Install BLAST+ via Homebrew on macOS."""
        if not shutil.which('brew'):
            return False, (
                "Homebrew not found. Install Homebrew first "
                "(https://brew.sh) or install BLAST+ manually."
            )

        if cb:
            cb("Installing BLAST+ via Homebrew...")

        try:
            result = subprocess.run(
                ['brew', 'install', 'blast'],
                capture_output=True, text=True, timeout=600
            )
            if result.returncode == 0:
                ok, path = self.detect_blastn()
                if ok:
                    self.log_info(f"BLAST+ installed at {path}")
                    return True, path
            return False, f"brew install blast failed: {result.stderr[:300]}"
        except subprocess.TimeoutExpired:
            return False, "BLAST+ installation timed out (>10 min)"
        except Exception as e:
            return False, str(e)

    def _install_blast_linux(self, cb) -> Tuple[bool, str]:
        """Install BLAST+ via apt / yum on Linux."""
        if cb:
            cb("Installing BLAST+ via package manager...")

        for cmd in [
            ['sudo', 'apt-get', 'install', '-y', 'ncbi-blast+'],
            ['sudo', 'yum', 'install', '-y', 'ncbi-blast+'],
        ]:
            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600
                )
                if result.returncode == 0:
                    ok, path = self.detect_blastn()
                    if ok:
                        self.log_info(f"BLAST+ installed at {path}")
                        return True, path
            except (subprocess.SubprocessError, FileNotFoundError):
                continue

        return False, (
            "Could not install BLAST+ via apt-get or yum. "
            "Please install manually."
        )

    def _install_blast_windows(self, cb) -> Tuple[bool, str]:
        """Download and extract BLAST+ from NCBI FTP on Windows."""
        import tarfile
        import re as _re
        import json

        # tools/ directory next to src/
        base_dir = Path(__file__).resolve().parent.parent.parent
        tools_dir = base_dir / "tools" / "blast"
        tools_dir.mkdir(parents=True, exist_ok=True)

        # Check if already extracted
        blastn_exe = tools_dir / "bin" / "blastn.exe"
        if blastn_exe.exists():
            self._save_blast_paths_to_config(tools_dir / "bin")
            return True, str(blastn_exe)

        if cb:
            cb("Downloading BLAST+ from NCBI (~137 MB)...")

        ftp_url = (
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )
        try:
            request = urllib.request.Request(
                ftp_url,
                headers={"User-Agent": "PrimerDesigner/1.0"},
            )
            with urllib.request.urlopen(request, timeout=30) as resp:
                html = resp.read().decode('utf-8')

            pattern = r'(ncbi-blast-[\d.+]+-x64-win64\.tar\.gz)'
            match = _re.search(pattern, html)
            if not match:
                return False, "Could not find Windows BLAST+ archive on NCBI FTP."

            archive_name = match.group(1)
            download_url = ftp_url + archive_name
            archive_path = tools_dir / archive_name

            if cb:
                cb(f"Downloading {archive_name}...")

            urllib.request.urlretrieve(download_url, str(archive_path))

            if cb:
                cb("Extracting BLAST+...")

            with tarfile.open(str(archive_path), 'r:gz') as tar:
                tar.extractall(path=str(tools_dir))

            archive_path.unlink(missing_ok=True)

            # Move bin/ from extracted ncbi-blast-X.X.X+/ to tools/blast/bin/
            extracted_dirs = [
                d for d in tools_dir.iterdir()
                if d.is_dir() and d.name.startswith('ncbi-blast-')
            ]
            if extracted_dirs:
                extracted = extracted_dirs[0]
                extracted_bin = extracted / "bin"
                target_bin = tools_dir / "bin"
                if target_bin.exists():
                    shutil.rmtree(str(target_bin))
                shutil.move(str(extracted_bin), str(target_bin))
                shutil.rmtree(str(extracted), ignore_errors=True)

            if blastn_exe.exists():
                self._save_blast_paths_to_config(tools_dir / "bin")
                self.log_info(f"BLAST+ installed to {tools_dir / 'bin'}")
                return True, str(blastn_exe)

            return False, "Extraction succeeded but blastn.exe not found."

        except Exception as e:
            return False, f"Failed to download BLAST+ for Windows: {e}"

    def _save_blast_paths_to_config(self, bin_dir: Path):
        """Update config.json with local BLAST+ tool paths."""
        import json

        base_dir = Path(__file__).resolve().parent.parent.parent
        config_path = base_dir / "config.json"
        try:
            with open(config_path, 'r') as f:
                cfg = json.load(f)

            os_type = platform.system().lower()
            ext = ".exe" if os_type == "windows" else ""
            cfg.setdefault("blast", {})
            cfg["blast"]["blastn_path"] = str(bin_dir / f"blastn{ext}")
            cfg["blast"]["makeblastdb_path"] = str(bin_dir / f"makeblastdb{ext}")

            with open(config_path, 'w') as f:
                json.dump(cfg, f, indent=2)
                f.write('\n')

            self.log_info(f"Updated config.json with BLAST+ paths from {bin_dir}")
        except Exception as e:
            self.log_warning(f"Could not update config.json: {e}")

    # ------------------------------------------------------------------
    #  BLAST database detection / building
    # ------------------------------------------------------------------

    def detect_blast_db(self) -> Tuple[bool, str]:
        """
        Check if BLAST nucleotide database exists.
        Returns (exists, db_base_path).
        """
        # 1. Configured path
        configured = self.config.blast.blast_db_path
        if configured and self._verify_db(configured):
            return True, configured

        # 2. Default location (genome data dir + genome.fa)
        default_base = str(self.genome_dir / "GRCh38" / "genome")
        if self._verify_db(default_base):
            return True, default_base

        return False, ""

    def _verify_db(self, base_path: str) -> bool:
        """Verify BLAST database files exist for a given base path."""
        return all(
            Path(base_path + ext).exists()
            for ext in self.BLAST_DB_EXTENSIONS
        )

    def setup_blast_db(
        self, progress_callback: Optional[Callable] = None
    ) -> Tuple[bool, str]:
        """
        Ensure BLAST database is available.
        Builds the database from existing genome.fa if needed.
        Returns (success, db_base_path_or_error).
        """
        exists, path = self.detect_blast_db()
        if exists:
            self.log_info(f"BLAST database already exists at {path}")
            return True, path

        # Check for makeblastdb
        mdb_ok, mdb_path = self.detect_makeblastdb()
        if not mdb_ok:
            return False, (
                "makeblastdb not found -- cannot build BLAST database. "
                "Ensure BLAST+ is fully installed."
            )

        # Locate genome FASTA — download automatically if missing
        fasta_path = self.genome_dir / "GRCh38" / "genome.fa"
        if not fasta_path.exists():
            if progress_callback:
                progress_callback(
                    "Downloading GRCh38 reference genome (~900 MB compressed)..."
                )
            dl_ok, dl_msg = self._download_genome_fasta(fasta_path, progress_callback)
            if not dl_ok:
                return False, f"Could not obtain genome FASTA: {dl_msg}"

        db_base = str(self.genome_dir / "GRCh38" / "genome")

        if progress_callback:
            progress_callback(
                "Building BLAST database from genome.fa "
                "(this may take 10-15 minutes)..."
            )

        try:
            cmd = [
                mdb_path,
                '-in', str(fasta_path),
                '-dbtype', 'nucl',
                '-parse_seqids',
                '-out', db_base,
            ]
            self.log_info(f"Running: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True, text=True,
                timeout=3600  # 1 hour max
            )

            if result.returncode == 0 and self._verify_db(db_base):
                self.log_info(f"BLAST database built at {db_base}")
                return True, db_base

            return False, (
                f"makeblastdb failed (exit {result.returncode}): "
                f"{result.stderr[:500]}"
            )

        except subprocess.TimeoutExpired:
            return False, "BLAST database creation timed out (>1 hour)"
        except Exception as e:
            return False, str(e)

    # ------------------------------------------------------------------
    #  Genome FASTA auto-download
    # ------------------------------------------------------------------

    # Primary assembly FASTA (no ALT/HLA contigs — smaller, sufficient for
    # primer BLAST).  GRCh38 primary assembly from NCBI (~900 MB gzipped).
    _GENOME_URLS = [
        (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
            "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
            "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
            True,   # gzipped
        ),
        (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"
            "hg38.fa.gz",
            True,   # gzipped
        ),
    ]

    def _download_genome_fasta(
        self,
        fasta_path: Path,
        progress_callback: Optional[Callable] = None,
    ) -> Tuple[bool, str]:
        """
        Download GRCh38 reference genome FASTA.

        Tries multiple mirror URLs.  If the file is gzipped it is
        decompressed on the fly to ``fasta_path``.

        Returns (success, error_message).
        """
        fasta_path.parent.mkdir(parents=True, exist_ok=True)

        for url, is_gzipped in self._GENOME_URLS:
            try:
                self.log_info(f"Downloading genome from: {url}")
                if progress_callback:
                    progress_callback(
                        f"Downloading genome FASTA (this may take 15-30 min)..."
                    )

                request = urllib.request.Request(
                    url,
                    headers={"User-Agent": "PrimerDesigner/1.0"},
                )
                # No timeout -- genome files are ~900 MB compressed / ~3 GB
                # uncompressed; a connection timeout would interrupt large
                # downloads on slow connections.
                response = urllib.request.urlopen(request)

                # Stream-download → decompress → write
                tmp_path = fasta_path.with_suffix('.fa.tmp')
                bytes_written = 0
                chunk_size = 1024 * 1024  # 1 MB

                try:
                    if is_gzipped:
                        # Download gzipped to temp, then decompress
                        gz_path = fasta_path.with_suffix('.fa.gz')
                        with open(gz_path, 'wb') as gz_out:
                            while True:
                                chunk = response.read(chunk_size)
                                if not chunk:
                                    break
                                gz_out.write(chunk)
                                bytes_written += len(chunk)
                                if progress_callback and bytes_written % (50 * chunk_size) == 0:
                                    mb = bytes_written / (1024 * 1024)
                                    progress_callback(
                                        f"Downloading genome: {mb:.0f} MB downloaded..."
                                    )

                        if progress_callback:
                            progress_callback("Decompressing genome FASTA...")

                        # Decompress .gz → .fa
                        with gzip.open(gz_path, 'rb') as gz_in, \
                             open(tmp_path, 'wb') as fa_out:
                            while True:
                                chunk = gz_in.read(chunk_size * 4)
                                if not chunk:
                                    break
                                fa_out.write(chunk)

                        # Clean up .gz
                        gz_path.unlink(missing_ok=True)
                    else:
                        # Direct download (uncompressed)
                        with open(tmp_path, 'wb') as fa_out:
                            while True:
                                chunk = response.read(chunk_size)
                                if not chunk:
                                    break
                                fa_out.write(chunk)
                                bytes_written += len(chunk)
                                if progress_callback and bytes_written % (50 * chunk_size) == 0:
                                    mb = bytes_written / (1024 * 1024)
                                    progress_callback(
                                        f"Downloading genome: {mb:.0f} MB downloaded..."
                                    )

                    # Verify minimum size (genome should be >2.5 GB uncompressed)
                    if tmp_path.stat().st_size < 100_000_000:  # < 100 MB = broken
                        tmp_path.unlink(missing_ok=True)
                        self.log_warning(f"Downloaded file too small, trying next URL")
                        continue

                    # Atomic rename
                    tmp_path.rename(fasta_path)
                    self.log_info(
                        f"Genome FASTA downloaded to {fasta_path} "
                        f"({fasta_path.stat().st_size / 1e9:.1f} GB)"
                    )
                    return True, ""

                except Exception as e:
                    # Clean up partial files
                    tmp_path.unlink(missing_ok=True)
                    fasta_path.with_suffix('.fa.gz').unlink(missing_ok=True)
                    raise

            except urllib.error.URLError as e:
                self.log_warning(f"Download failed from {url}: {e}")
                continue
            except Exception as e:
                self.log_warning(f"Download failed from {url}: {e}")
                continue

        return False, (
            "Could not download GRCh38 genome FASTA from any mirror. "
            "Check your internet connection or download manually to: "
            f"{fasta_path}"
        )

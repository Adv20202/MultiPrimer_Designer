"""
MANE (Matched Annotation from NCBI and EMBL-EBI) transcript manager.
Handles MANE Select and MANE Plus Clinical transcript lookup.
"""

import csv
import gzip
import os
import re
import urllib.request
from typing import Optional, Dict, List, Any
from pathlib import Path

from ..core.models import MANEType
from ..utils.logger import LoggerMixin
from ..utils.config import get_config


class MANEManager(LoggerMixin):
    """
    Manager for MANE transcript data.
    Downloads and caches MANE data from NCBI/Ensembl.
    """

    # MANE data directory - we'll auto-discover the latest version
    MANE_FTP_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"

    def __init__(self, cache_dir: str = None, auto_download: bool = False):
        """
        Initialize MANE manager.

        Args:
            cache_dir: Directory to cache MANE data file
            auto_download: If True, automatically download MANE data if not available
        """
        config = get_config()
        self.cache_dir = Path(cache_dir or config.data_cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.mane_file = self.cache_dir / "MANE_summary.txt"
        self.auto_download = auto_download

        # In-memory cache: gene_symbol -> {transcript_accession -> MANE info}
        self._mane_data: Dict[str, Dict[str, Dict[str, Any]]] = {}
        self._download_attempted = False
        self._download_failed = False

        # Load data if available
        if self.mane_file.exists():
            self._load_mane_data()

    def ensure_data_available(self, progress_callback=None) -> bool:
        """
        Ensure MANE data is downloaded and loaded.

        Args:
            progress_callback: Optional callback for progress updates

        Returns:
            True if data is available
        """
        if self._mane_data:
            return True

        if not self.mane_file.exists():
            if not self._download_mane_data(progress_callback):
                return False

        return self._load_mane_data()

    def _discover_mane_file(self) -> Optional[str]:
        """
        Discover the latest MANE summary file from the FTP directory.

        Returns:
            URL to the MANE summary file or None if not found
        """
        try:
            self.log_info(f"Discovering MANE files in {self.MANE_FTP_DIR}")

            # Fetch directory listing
            with urllib.request.urlopen(self.MANE_FTP_DIR, timeout=10) as response:
                html = response.read().decode('utf-8')

            # Find all MANE.GRCh38.v*.summary.txt.gz files
            pattern = r'MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz'
            matches = re.findall(pattern, html)

            if not matches:
                self.log_warning("No MANE summary files found in directory")
                return None

            # Sort versions and get the latest
            versions = sorted(matches, key=lambda v: [int(x) for x in v.split('.')], reverse=True)
            latest_version = versions[0]

            filename = f"MANE.GRCh38.v{latest_version}.summary.txt.gz"
            url = self.MANE_FTP_DIR + filename

            self.log_info(f"Found latest MANE version: {latest_version} ({url})")
            return url

        except Exception as e:
            self.log_warning(f"Failed to discover MANE file: {e}")
            return None

    def _download_mane_data(self, progress_callback=None) -> bool:
        """
        Download MANE summary file from NCBI.
        Auto-discovers the latest version available.

        Args:
            progress_callback: Optional callback(downloaded, total) for progress

        Returns:
            True if download successful
        """
        # First, try to discover the latest file
        url = self._discover_mane_file()

        if not url:
            # Fallback: try known versions in descending order
            self.log_warning("Auto-discovery failed, trying fallback URLs")
            for fallback_ver in ["1.5", "1.4", "1.3"]:
                fallback_url = self.MANE_FTP_DIR + f"MANE.GRCh38.v{fallback_ver}.summary.txt.gz"
                try:
                    req = urllib.request.Request(fallback_url, method='HEAD')
                    urllib.request.urlopen(req, timeout=10)
                    url = fallback_url
                    self.log_info(f"Fallback MANE URL found: v{fallback_ver}")
                    break
                except Exception:
                    continue
            if not url:
                url = self.MANE_FTP_DIR + "MANE.GRCh38.v1.5.summary.txt.gz"

        try:
            self.log_info(f"Downloading MANE data from {url}")
            gz_file = self.cache_dir / "MANE_summary.txt.gz"

            # Download with progress
            def report_progress(block_num, block_size, total_size):
                if progress_callback and total_size > 0:
                    downloaded = block_num * block_size
                    progress_callback(downloaded, total_size)

            urllib.request.urlretrieve(url, gz_file, report_progress)

            # Decompress
            self.log_info("Decompressing MANE data...")
            with gzip.open(gz_file, 'rt') as f_in:
                with open(self.mane_file, 'w') as f_out:
                    f_out.write(f_in.read())

            # Remove compressed file
            os.remove(gz_file)

            self.log_info(f"MANE data saved to {self.mane_file}")
            return True

        except Exception as e:
            self.log_error(f"Failed to download MANE data: {e}")
            return False

    def _load_mane_data(self) -> bool:
        """
        Load MANE data from cached file into memory.

        Returns:
            True if loading successful
        """
        if not self.mane_file.exists():
            self.log_warning("MANE data file not found")
            return False

        try:
            self.log_info(f"Loading MANE data from {self.mane_file}")
            self._mane_data.clear()

            with open(self.mane_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')

                for row in reader:
                    gene_symbol = row.get('symbol', '').upper()
                    refseq_nuc = row.get('RefSeq_nuc', '')  # NM_xxx.x
                    mane_status = row.get('MANE_status', '')

                    if not gene_symbol or not refseq_nuc:
                        continue

                    # Determine MANE type
                    if 'MANE Select' in mane_status:
                        mane_type = MANEType.MANE_SELECT
                    elif 'MANE Plus Clinical' in mane_status:
                        mane_type = MANEType.MANE_PLUS_CLINICAL
                    else:
                        continue  # Skip non-MANE entries

                    # Extract base accession without version
                    parts = refseq_nuc.split('.')
                    base_accession = parts[0]
                    version = parts[1] if len(parts) > 1 else None

                    # Store data
                    if gene_symbol not in self._mane_data:
                        self._mane_data[gene_symbol] = {}

                    self._mane_data[gene_symbol][base_accession] = {
                        'type': mane_type,
                        'full_accession': refseq_nuc,
                        'version': version,
                        'ensembl_nuc': row.get('Ensembl_nuc', ''),
                        'ensembl_prot': row.get('Ensembl_prot', ''),
                        'refseq_prot': row.get('RefSeq_prot', ''),
                        'gene_id': row.get('GeneID', ''),
                        'hgnc_id': row.get('HGNC_ID', ''),
                        'chromosome': row.get('#NCBI_GeneID', '').split(':')[0] if row.get('#NCBI_GeneID') else '',
                    }

            self.log_info(f"Loaded MANE data for {len(self._mane_data)} genes")
            return True

        except Exception as e:
            self.log_exception(f"Failed to load MANE data: {e}")
            return False

    def get_mane_info(self, transcript_accession: str,
                      gene_symbol: str = None) -> Optional[Dict[str, Any]]:
        """
        Get MANE information for a transcript.

        Args:
            transcript_accession: RefSeq transcript accession (with or without version)
            gene_symbol: Gene symbol (optional, for verification)

        Returns:
            Dictionary with MANE info or None if not a MANE transcript
        """
        if not self._mane_data:
            # Only try to download if auto_download is enabled and not already tried
            if self.auto_download and not self._download_attempted:
                self._download_attempted = True
                if not self.ensure_data_available():
                    self._download_failed = True
                    self.log_warning("MANE data not available - validation will continue without MANE checks")
                    return None
            else:
                # Don't download, just return None
                return None

        # Extract base accession
        base_accession = transcript_accession.split('.')[0].upper()

        # If gene symbol provided, check that gene first
        if gene_symbol:
            gene_data = self._mane_data.get(gene_symbol.upper())
            if gene_data and base_accession in gene_data:
                return gene_data[base_accession]

        # Search all genes for this transcript
        for gene, transcripts in self._mane_data.items():
            if base_accession in transcripts:
                info = transcripts[base_accession].copy()
                info['gene_symbol'] = gene
                return info

        return None

    def get_mane_select(self, gene_symbol: str) -> Optional[str]:
        """
        Get the MANE Select transcript for a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            RefSeq accession or None
        """
        gene_data = self._mane_data.get(gene_symbol.upper())
        if not gene_data:
            return None

        for accession, info in gene_data.items():
            if info['type'] == MANEType.MANE_SELECT:
                return info['full_accession']

        return None

    def get_mane_clinical(self, gene_symbol: str) -> Optional[str]:
        """
        Get the MANE Plus Clinical transcript for a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            RefSeq accession or None
        """
        gene_data = self._mane_data.get(gene_symbol.upper())
        if not gene_data:
            return None

        for accession, info in gene_data.items():
            if info['type'] == MANEType.MANE_PLUS_CLINICAL:
                return info['full_accession']

        return None

    def get_all_mane_transcripts(self, gene_symbol: str) -> List[Dict[str, Any]]:
        """
        Get all MANE transcripts for a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            List of MANE transcript info dictionaries
        """
        gene_data = self._mane_data.get(gene_symbol.upper())
        if not gene_data:
            return []

        return list(gene_data.values())

    def is_mane_transcript(self, transcript_accession: str,
                           gene_symbol: str = None) -> bool:
        """
        Check if a transcript is a MANE transcript.

        Args:
            transcript_accession: RefSeq transcript accession
            gene_symbol: Gene symbol (optional)

        Returns:
            True if transcript is MANE Select or MANE Plus Clinical
        """
        return self.get_mane_info(transcript_accession, gene_symbol) is not None

    def get_mane_type(self, transcript_accession: str,
                      gene_symbol: str = None) -> MANEType:
        """
        Get the MANE type for a transcript.

        Args:
            transcript_accession: RefSeq transcript accession
            gene_symbol: Gene symbol (optional)

        Returns:
            MANEType enum value
        """
        info = self.get_mane_info(transcript_accession, gene_symbol)
        if info:
            return info['type']
        return MANEType.NOT_MANE

    def search_genes(self, query: str) -> List[str]:
        """
        Search for genes by partial name match.

        Args:
            query: Search query

        Returns:
            List of matching gene symbols
        """
        query_upper = query.upper()
        matches = [gene for gene in self._mane_data.keys()
                   if query_upper in gene]
        return sorted(matches)[:20]  # Limit results

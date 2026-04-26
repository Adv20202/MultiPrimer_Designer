"""
Configuration management for the Primer Designer application.
Handles application settings, API keys, and user preferences.
"""

import json
import os
import threading
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional, Dict, Any, List

from .logger import LoggerMixin


@dataclass
class APIConfig:
    """Configuration for external API connections."""
    ncbi_api_key: str = ""
    ncbi_email: str = ""
    request_timeout: int = 30
    max_retries: int = 3
    retry_delay: float = 1.0


@dataclass
class BlastConfig:
    """Configuration for BLAST+ aligner (homology analysis)."""
    blastn_path: str = "blastn"
    makeblastdb_path: str = "makeblastdb"
    blast_db_path: str = ""
    threads: int = 4
    evalue: float = 1e-10
    min_percent_identity: float = 85.0
    min_aligned_length: int = 50
    max_target_seqs: int = 500
    max_hsps: int = 3  # Max HSPs per subject sequence


@dataclass
class BwaConfig:
    """Configuration for BWA aligner (specificity checking)."""
    bwa_path: str = "bwa"
    genome_index_path: str = ""
    threads: int = 4
    min_seed_length: int = 10
    min_alignment_score: int = 20


@dataclass
class Primer3Config:
    """Default Primer3 configuration."""
    # These are the default values as specified in requirements
    primer_min_size: int = 18
    primer_max_size: int = 27
    primer_opt_size: int = 22
    primer_min_tm: float = 57.0
    primer_max_tm: float = 64.0
    primer_opt_tm: float = 60.0
    primer_min_gc: float = 40.0
    primer_max_gc: float = 60.0
    primer_opt_gc: float = 50.0

    # Chemistry settings
    mv_conc: float = 50.0  # Monovalent cations (Na+, K+) in mM
    dv_conc: float = 1.5   # Divalent cations (Mg2+) in mM
    dntp_conc: float = 0.6  # dNTPs in mM
    dna_conc: float = 50.0  # Template DNA in nM

    # Complementarity
    max_self_complementarity: float = 8.0
    max_end_complementarity: float = 3.0
    max_pair_complementarity: float = 8.0

    # Product size
    product_min_size: int = 100
    product_max_size: int = 500
    product_opt_size: int = 350

    # qPCR specific
    qpcr_product_min_size: int = 70
    qpcr_product_max_size: int = 150
    probe_min_size: int = 18
    probe_max_size: int = 30
    probe_min_tm: float = 65.0
    probe_max_tm: float = 70.0
    probe_opt_tm: float = 67.0


@dataclass
class PopulationConfig:
    """Configuration for population variant filtering."""
    available_populations: List[str] = field(default_factory=lambda: [
        "global",
        "afr",  # African/African American
        "ami",  # Amish
        "amr",  # Latino/Admixed American
        "asj",  # Ashkenazi Jewish
        "eas",  # East Asian
        "fin",  # Finnish
        "nfe",  # Non-Finnish European
        "sas",  # South Asian
        "oth"   # Other
    ])
    population_labels: Dict[str, str] = field(default_factory=lambda: {
        "global": "Global",
        "afr": "African/African American",
        "ami": "Amish",
        "amr": "Latino/Admixed American",
        "asj": "Ashkenazi Jewish",
        "eas": "East Asian",
        "fin": "Finnish",
        "nfe": "Non-Finnish European",
        "sas": "South Asian",
        "oth": "Other"
    })
    default_maf_threshold: float = 0.005
    default_selected_populations: List[str] = field(default_factory=lambda: ["global"])


@dataclass
class Config(LoggerMixin):
    """Main application configuration."""
    # Sub-configurations
    api: APIConfig = field(default_factory=APIConfig)
    blast: BlastConfig = field(default_factory=BlastConfig)
    bwa: BwaConfig = field(default_factory=BwaConfig)
    primer3: Primer3Config = field(default_factory=Primer3Config)
    population: PopulationConfig = field(default_factory=PopulationConfig)

    # Application settings
    default_assembly: str = "GRCh38"
    log_level: str = "INFO"
    data_cache_dir: str = ""
    max_concurrent_requests: int = 5

    # Design defaults
    min_distance_from_exon_junction: int = 60
    min_distance_from_variant: int = 40

    # File paths
    config_file_path: str = ""

    def __post_init__(self):
        """Initialize paths after object creation."""
        if not self.data_cache_dir:
            base_dir = Path(__file__).parent.parent.parent
            self.data_cache_dir = str(base_dir / "data" / "cache")
            Path(self.data_cache_dir).mkdir(parents=True, exist_ok=True)

    @classmethod
    def load(cls, config_path: Optional[str] = None) -> 'Config':
        """
        Load configuration from file.

        Args:
            config_path: Path to config file. If None, uses default location.

        Returns:
            Config instance
        """
        if config_path is None:
            base_dir = Path(__file__).parent.parent.parent
            config_path = str(base_dir / "config.json")

        config = cls()
        config.config_file_path = config_path

        if Path(config_path).exists():
            try:
                with open(config_path, 'r') as f:
                    data = json.load(f)
                config._update_from_dict(data)
                config.log_info(f"Configuration loaded from {config_path}")
            except Exception as e:
                config.log_warning(f"Failed to load config from {config_path}: {e}")

        return config

    def _update_from_dict(self, data: Dict[str, Any]):
        """Update configuration from a dictionary."""
        if 'api' in data:
            for key, value in data['api'].items():
                if hasattr(self.api, key):
                    setattr(self.api, key, value)

        if 'blast' in data:
            for key, value in data['blast'].items():
                if hasattr(self.blast, key):
                    setattr(self.blast, key, value)

        if 'bwa' in data:
            for key, value in data['bwa'].items():
                if hasattr(self.bwa, key):
                    setattr(self.bwa, key, value)

        if 'primer3' in data:
            for key, value in data['primer3'].items():
                if hasattr(self.primer3, key):
                    setattr(self.primer3, key, value)

        if 'population' in data:
            for key, value in data['population'].items():
                if hasattr(self.population, key):
                    setattr(self.population, key, value)

        # Update top-level settings
        for key in ['default_assembly', 'log_level', 'data_cache_dir',
                    'max_concurrent_requests', 'min_distance_from_exon_junction',
                    'min_distance_from_variant']:
            if key in data:
                setattr(self, key, data[key])

    def save(self, config_path: Optional[str] = None):
        """
        Save configuration to file.

        Args:
            config_path: Path to config file. If None, uses current path.
        """
        if config_path is None:
            config_path = self.config_file_path

        if not config_path:
            base_dir = Path(__file__).parent.parent.parent
            config_path = str(base_dir / "config.json")

        data = {
            'api': asdict(self.api),
            'blast': asdict(self.blast),
            'bwa': asdict(self.bwa),
            'primer3': asdict(self.primer3),
            'population': asdict(self.population),
            'default_assembly': self.default_assembly,
            'log_level': self.log_level,
            'data_cache_dir': self.data_cache_dir,
            'max_concurrent_requests': self.max_concurrent_requests,
            'min_distance_from_exon_junction': self.min_distance_from_exon_junction,
            'min_distance_from_variant': self.min_distance_from_variant
        }

        try:
            with open(config_path, 'w') as f:
                json.dump(data, f, indent=2)
            self.log_info(f"Configuration saved to {config_path}")
        except Exception as e:
            self.log_error(f"Failed to save config to {config_path}: {e}")
            raise

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'api': asdict(self.api),
            'blast': asdict(self.blast),
            'bwa': asdict(self.bwa),
            'primer3': asdict(self.primer3),
            'population': asdict(self.population),
            'default_assembly': self.default_assembly,
            'log_level': self.log_level,
            'data_cache_dir': self.data_cache_dir,
            'max_concurrent_requests': self.max_concurrent_requests,
            'min_distance_from_exon_junction': self.min_distance_from_exon_junction,
            'min_distance_from_variant': self.min_distance_from_variant
        }


# Singleton configuration instance with thread safety
_config: Optional[Config] = None
_config_lock = threading.Lock()


def get_config() -> Config:
    """Get the global configuration instance (thread-safe)."""
    global _config
    if _config is None:
        with _config_lock:
            # Double-check locking pattern
            if _config is None:
                _config = Config.load()
    return _config


def set_config(config: Config):
    """Set the global configuration instance (thread-safe)."""
    global _config
    with _config_lock:
        _config = config

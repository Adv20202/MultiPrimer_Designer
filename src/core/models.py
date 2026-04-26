"""
Core data models for the Primer Designer application.
Defines all essential domain objects used throughout the system.
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
from enum import Enum, auto


class VariantType(Enum):
    """Types of genetic variants supported by the application."""
    SUBSTITUTION = auto()
    DELETION = auto()
    INSERTION = auto()
    DUPLICATION = auto()
    DELINS = auto()
    UNKNOWN = auto()


class MANEType(Enum):
    """MANE transcript classification."""
    MANE_SELECT = "MANE Select"
    MANE_PLUS_CLINICAL = "MANE Plus Clinical"
    NOT_MANE = "Not MANE"


class DesignMode(Enum):
    """Primer design mode."""
    PCR = "PCR"
    QPCR = "qPCR"


class Strand(Enum):
    """DNA strand orientation."""
    POSITIVE = "+"
    NEGATIVE = "-"


class ValidationStatus(Enum):
    """Validation result status."""
    VALID = "valid"
    WARNING = "warning"
    ERROR = "error"


@dataclass
class ValidationResult:
    """Result of a validation check."""
    status: ValidationStatus
    message: str
    details: Optional[Dict] = None


@dataclass
class GenomicPosition:
    """Represents a position on the genome."""
    chromosome: str
    start: int
    end: int
    strand: Strand = Strand.POSITIVE
    assembly: str = "GRCh38"

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end} ({self.strand.value})"

    def overlaps(self, other: 'GenomicPosition') -> bool:
        """Check if this position overlaps with another."""
        if self.chromosome != other.chromosome:
            return False
        return not (self.end < other.start or self.start > other.end)

    def distance_to(self, other: 'GenomicPosition') -> Optional[int]:
        """Calculate distance to another position on the same chromosome."""
        if self.chromosome != other.chromosome:
            return None
        if self.overlaps(other):
            return 0
        return min(abs(self.start - other.end), abs(self.end - other.start))


@dataclass
class Exon:
    """Represents an exon in a transcript."""
    number: int
    genomic_start: int
    genomic_end: int
    cds_start: Optional[int] = None
    cds_end: Optional[int] = None

    @property
    def length(self) -> int:
        return abs(self.genomic_end - self.genomic_start) + 1


@dataclass
class Transcript:
    """Represents a transcript (NM_ identifier)."""
    accession: str  # e.g., NM_002485
    version: Optional[int] = None  # e.g., 4 for NM_002485.4
    gene_symbol: str = ""
    chromosome: str = ""
    strand: Strand = Strand.POSITIVE
    cds_start: int = 0
    cds_end: int = 0
    exons: List[Exon] = field(default_factory=list)
    sequence: str = ""
    mane_type: MANEType = MANEType.NOT_MANE

    @property
    def full_accession(self) -> str:
        """Return full accession with version if available."""
        if self.version:
            return f"{self.accession}.{self.version}"
        return self.accession

    @property
    def cds_length(self) -> int:
        """Return the length of the CDS."""
        if self.cds_start == 0 and self.cds_end == 0:
            return 0
        return max(0, self.cds_end - self.cds_start + 1)

    def get_exon_for_cds_position(self, cds_position: int) -> Optional[Exon]:
        """Get the exon containing a given CDS position."""
        current_cds_pos = 1
        for exon in self.exons:
            if exon.cds_start is None or exon.cds_end is None:
                continue
            exon_cds_length = exon.cds_end - exon.cds_start + 1
            if current_cds_pos <= cds_position <= current_cds_pos + exon_cds_length - 1:
                return exon
            current_cds_pos += exon_cds_length
        return None


@dataclass
class Gene:
    """Represents a gene."""
    symbol: str
    name: str = ""
    chromosome: str = ""
    strand: Strand = Strand.POSITIVE
    transcripts: List[Transcript] = field(default_factory=list)

    def get_mane_select(self) -> Optional[Transcript]:
        """Get the MANE Select transcript for this gene."""
        for t in self.transcripts:
            if t.mane_type == MANEType.MANE_SELECT:
                return t
        return None

    def get_mane_clinical(self) -> Optional[Transcript]:
        """Get the MANE Plus Clinical transcript for this gene."""
        for t in self.transcripts:
            if t.mane_type == MANEType.MANE_PLUS_CLINICAL:
                return t
        return None


@dataclass
class Variant:
    """Represents a genetic variant from the input file."""
    row_number: int
    gene_symbol: str
    transcript_accession: str
    hgvs_c: str  # e.g., c.657_661delACAAA

    # Parsed data
    variant_type: VariantType = VariantType.UNKNOWN
    cds_start: int = 0
    cds_end: int = 0
    ref_allele: str = ""
    alt_allele: str = ""

    # Genomic coordinates (calculated)
    genomic_position: Optional[GenomicPosition] = None

    # Validation status
    validation_results: List[ValidationResult] = field(default_factory=list)
    is_valid: bool = True

    # Associated transcript data
    transcript: Optional[Transcript] = None
    exon: Optional[Exon] = None

    # Design mode
    design_mode: DesignMode = DesignMode.PCR

    # Grouping
    group_id: Optional[int] = None

    def add_validation_result(self, result: ValidationResult):
        """Add a validation result and update validity status."""
        self.validation_results.append(result)
        if result.status == ValidationStatus.ERROR:
            self.is_valid = False

    @property
    def has_warnings(self) -> bool:
        """Check if the variant has any warnings."""
        return any(r.status == ValidationStatus.WARNING for r in self.validation_results)


@dataclass
class PopulationVariant:
    """Represents a known population variant (e.g., from gnomAD)."""
    rsid: Optional[str]
    chromosome: str
    position: int
    ref: str
    alt: str
    maf_global: float = 0.0
    maf_by_population: Dict[str, float] = field(default_factory=dict)

    def get_maf(self, population: str = "global") -> float:
        """Get MAF for a specific population."""
        if population == "global":
            return self.maf_global
        return self.maf_by_population.get(population, 0.0)


@dataclass
class Primer:
    """Represents a single primer (forward or reverse)."""
    sequence: str
    start: int  # Genomic start position
    end: int    # Genomic end position
    tm: float   # Melting temperature
    gc_content: float
    is_forward: bool

    # Quality metrics
    self_complementarity: float = 0.0
    end_complementarity: float = 0.0
    hairpin_tm: float = 0.0

    # Population variants in primer region
    population_variants: List[PopulationVariant] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def has_high_maf_variants(self, maf_threshold: float = 0.01) -> bool:
        """Check if primer region contains high MAF variants."""
        return any(v.maf_global >= maf_threshold for v in self.population_variants)


@dataclass
class Probe:
    """Represents a qPCR probe (e.g., TaqMan)."""
    sequence: str
    start: int
    end: int
    tm: float
    gc_content: float

    # Probe should cover the variant
    variant_position_in_probe: int = 0

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class PrimerPair:
    """Represents a pair of primers (forward + reverse)."""
    forward: Primer
    reverse: Primer
    product_size: int

    # Quality metrics
    tm_difference: float = 0.0
    pair_complementarity: float = 0.0

    # Specificity results
    is_specific: bool = True
    off_target_sites: List[GenomicPosition] = field(default_factory=list)
    specificity_warnings: List[str] = field(default_factory=list)

    # Detailed specificity from orchestrator (Primer-BLAST + UCSC isPCR)
    specificity_verdict: str = ""           # e.g. "SPECIFIC", "PSEUDOGENE RISK"
    specificity_tool_summary: str = ""      # e.g. "Primer-BLAST: OK | isPCR: OK"
    pseudogene_risk: bool = False
    off_target_details: List[str] = field(default_factory=list)

    # Homology discrimination — how well primers distinguish from pseudogenes
    homology_discrimination_score: float = 0.0  # Weighted score (higher = better)
    fwd_discriminating_positions: int = 0        # # of discriminating positions in fwd primer
    rev_discriminating_positions: int = 0        # # of discriminating positions in rev primer
    homology_warning: str = ""                   # Warning if no discrimination

    @property
    def mean_tm(self) -> float:
        return (self.forward.tm + self.reverse.tm) / 2


@dataclass
class Amplicon:
    """Represents a PCR amplicon covering one or more variants."""
    primer_pair: PrimerPair
    variants: List[Variant]
    sequence: str = ""

    # Optional probe for qPCR
    probe: Optional[Probe] = None

    # Region information
    genomic_position: Optional[GenomicPosition] = None

    # Target position in sequence (for visualization)
    target_start: int = 0  # 0-based position of target in design sequence
    target_length: int = 1  # Length of target region
    exons_covered: List[Exon] = field(default_factory=list)
    masked_positions: set = field(default_factory=set)  # Positions masked with 'N' due to polymorphisms

    # Exon boundaries relative to the amplicon (0-based within product)
    # Each tuple: (amplicon_start, amplicon_end, exon_number)
    exon_regions_in_amplicon: List[Tuple[int, int, int]] = field(default_factory=list)

    @property
    def is_multi_variant(self) -> bool:
        return len(self.variants) > 1


@dataclass
class DesignParameters:
    """Parameters for primer/probe design."""
    # Amplicon size
    min_amplicon_size: int = 100
    max_amplicon_size: int = 500

    # Primer parameters
    primer_min_size: int = 18
    primer_max_size: int = 27
    primer_opt_size: int = 22
    primer_min_tm: float = 57.0
    primer_max_tm: float = 64.0
    primer_opt_tm: float = 60.0
    primer_min_gc: float = 40.0
    primer_max_gc: float = 60.0
    primer_opt_gc: float = 50.0

    # Chemistry
    mv_conc: float = 50.0  # Monovalent cation concentration (mM)
    dv_conc: float = 1.5   # Divalent cation concentration (mM)
    dntp_conc: float = 0.6  # dNTP concentration (mM)
    dna_conc: float = 50.0  # DNA concentration (nM)

    # Complementarity
    max_self_complementarity: float = 8.0
    max_end_complementarity: float = 3.0
    max_pair_complementarity: float = 8.0
    max_tm_difference: float = 3.0

    # Positioning
    min_distance_from_exon_junction: int = 60  # Min bp from exon/intron boundary on intronic side (splice-site exclusion)
    min_distance_from_variant: int = 40
    use_splice_site_constraint: bool = False  # Exclude intronic positions near exon junctions
    use_variant_distance_constraint: bool = True

    # Population variant filtering
    filter_population_variants: bool = True
    maf_threshold: float = 0.005
    selected_populations: List[str] = field(default_factory=lambda: ["global"])

    # qPCR specific
    qpcr_min_amplicon_size: int = 70
    qpcr_max_amplicon_size: int = 150
    probe_min_size: int = 18
    probe_max_size: int = 30
    probe_min_tm: float = 65.0
    probe_max_tm: float = 70.0

    # Homology discrimination
    use_homology_discrimination: bool = True

    # Number of primer pairs to return
    num_primer_pairs: int = 5

    # Minimum number of specific pairs to find in continuous mode before stopping
    min_specific_pairs: int = 3

    # Reference genome
    reference_assembly: str = "GRCh38"


@dataclass
class DesignResult:
    """Result of primer design for a variant or group of variants."""
    variants: List[Variant]
    amplicons: List[Amplicon]
    success: bool
    message: str = ""
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)
    homology_discriminated: bool = False     # True if homology re-ranking was applied
    num_homologous_regions: int = 0          # Number of secondary hits used for discrimination
    homology_tier: int = 0                   # 0=none, 1=junction overlap, 2=OK regions, 3=fallback rerank
    homology_tier_message: str = ""          # Human-readable description of which tier succeeded
    homology_result: object = None           # Optional HomologyResult with secondary hit details (for reporting)


@dataclass
class ProjectState:
    """Represents the current state of a primer design project."""
    input_file_path: str = ""
    variants: List[Variant] = field(default_factory=list)
    design_parameters: DesignParameters = field(default_factory=DesignParameters)
    design_results: List[DesignResult] = field(default_factory=list)
    reference_assembly: str = "GRCh38"

    # Processing state
    is_validated: bool = False
    is_designed: bool = False

    # Cancellation flag
    cancel_requested: bool = False

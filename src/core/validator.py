"""
Validators for variant data.
Includes HGVS notation validation and comprehensive variant validation.
"""

import re
from typing import List, Optional, Tuple
from dataclasses import dataclass

from .models import (
    Variant, Transcript, VariantType, ValidationResult,
    ValidationStatus, MANEType
)
from ..utils.logger import LoggerMixin


class HGVSValidator(LoggerMixin):
    """
    Validator for HGVS c. (coding DNA) and g. (genomic) notation.
    Validates syntax according to HGVS nomenclature standards.

    Supported position prefixes:
      - c.  -- coding DNA (CDS-relative, 1-based)
      - c.* -- positions 3' of the stop codon (3'UTR)
      - c.-  -- positions 5' of the ATG start codon (5'UTR)
      - g.  -- genomic coordinates (1-based)
    """

    # Position pattern component: supports normal, 5'UTR (-N), 3'UTR (*N),
    # and intronic offsets (+/-N after main position)
    _POS = r'(-?\d+|\*\d+)'
    _OFFSET = r'(?:([+-])(\d+))?'

    # HGVS c. patterns for different variant types
    PATTERNS = {
        'substitution': re.compile(
            rf'^c\.{_POS}{_OFFSET}([ACGT])>([ACGT])$',
            re.IGNORECASE
        ),
        'deletion_single': re.compile(
            rf'^c\.{_POS}{_OFFSET}del([ACGT])?$',
            re.IGNORECASE
        ),
        'deletion_range': re.compile(
            rf'^c\.{_POS}{_OFFSET}_{_POS}{_OFFSET}del([ACGT]+)?$',
            re.IGNORECASE
        ),
        'duplication_single': re.compile(
            rf'^c\.{_POS}{_OFFSET}dup([ACGT])?$',
            re.IGNORECASE
        ),
        'duplication_range': re.compile(
            rf'^c\.{_POS}{_OFFSET}_{_POS}{_OFFSET}dup([ACGT]+)?$',
            re.IGNORECASE
        ),
        'insertion': re.compile(
            rf'^c\.{_POS}{_OFFSET}_{_POS}{_OFFSET}ins([ACGT]+)$',
            re.IGNORECASE
        ),
        'delins_single': re.compile(
            rf'^c\.{_POS}{_OFFSET}delins([ACGT]+)$',
            re.IGNORECASE
        ),
        'delins_range': re.compile(
            rf'^c\.{_POS}{_OFFSET}_{_POS}{_OFFSET}del([ACGT]*)ins([ACGT]+)$',
            re.IGNORECASE
        ),
        'repeat': re.compile(
            rf'^c\.{_POS}{_OFFSET}(?:_{_POS}{_OFFSET})?\[(\d+)\]$',
            re.IGNORECASE
        ),
    }

    # HGVS g. (genomic) patterns
    GENOMIC_PATTERNS = {
        'g_substitution': re.compile(
            r'^g\.(\d+)([ACGT])>([ACGT])$',
            re.IGNORECASE
        ),
        'g_deletion_single': re.compile(
            r'^g\.(\d+)del([ACGT])?$',
            re.IGNORECASE
        ),
        'g_deletion_range': re.compile(
            r'^g\.(\d+)_(\d+)del([ACGT]+)?$',
            re.IGNORECASE
        ),
        'g_duplication_single': re.compile(
            r'^g\.(\d+)dup([ACGT])?$',
            re.IGNORECASE
        ),
        'g_duplication_range': re.compile(
            r'^g\.(\d+)_(\d+)dup([ACGT]+)?$',
            re.IGNORECASE
        ),
        'g_insertion': re.compile(
            r'^g\.(\d+)_(\d+)ins([ACGT]+)$',
            re.IGNORECASE
        ),
        'g_delins_single': re.compile(
            r'^g\.(\d+)delins([ACGT]+)$',
            re.IGNORECASE
        ),
        'g_delins_range': re.compile(
            r'^g\.(\d+)_(\d+)del([ACGT]*)ins([ACGT]+)$',
            re.IGNORECASE
        ),
    }

    # RefSeq transcript pattern
    REFSEQ_PATTERN = re.compile(r'^NM_\d{6,}(\.\d+)?$', re.IGNORECASE)

    @staticmethod
    def _parse_hgvs_pos(pos_str: str) -> int:
        """Parse a single HGVS position token into an integer.

        Handles:
          '123'  -> 123
          '-50'  -> -50
          '*45'  -> a large sentinel value so downstream range checks
                    don't reject it; 3'UTR positions are validated
                    separately.
        """
        if pos_str.startswith('*'):
            # 3'UTR position -- return a sentinel that will always pass
            # the "within CDS" range check (handled specially later)
            return int(pos_str[1:]) + 900000000  # sentinel
        return int(pos_str)

    def validate_hgvs_syntax(self, hgvs: str) -> ValidationResult:
        """
        Validate HGVS c. or g. notation syntax.

        Args:
            hgvs: HGVS notation string

        Returns:
            ValidationResult indicating validity
        """
        if not hgvs:
            return ValidationResult(
                status=ValidationStatus.ERROR,
                message="Empty HGVS notation"
            )

        hgvs_stripped = hgvs.strip()

        # Check for g. prefix (genomic notation)
        if hgvs_stripped.lower().startswith('g.'):
            for pattern_name, pattern in self.GENOMIC_PATTERNS.items():
                if pattern.match(hgvs_stripped):
                    return ValidationResult(
                        status=ValidationStatus.VALID,
                        message=f"Valid HGVS genomic notation ({pattern_name})",
                        details={'notation_type': 'genomic'}
                    )
            return ValidationResult(
                status=ValidationStatus.ERROR,
                message=f"Invalid HGVS genomic notation format: {hgvs_stripped}"
            )

        # Check for c. prefix
        if not hgvs_stripped.lower().startswith('c.'):
            # Provide informative message for unsupported prefixes
            if hgvs_stripped.lower().startswith(('n.', 'p.', 'r.', 'm.')):
                prefix = hgvs_stripped[:2]
                return ValidationResult(
                    status=ValidationStatus.ERROR,
                    message=f"HGVS '{prefix}' notation is not supported. "
                            f"Please use 'c.' (coding DNA) or 'g.' (genomic) notation."
                )
            return ValidationResult(
                status=ValidationStatus.ERROR,
                message="HGVS notation must start with 'c.' (coding DNA) or 'g.' (genomic)"
            )

        # Try to match against known c. patterns
        for pattern_name, pattern in self.PATTERNS.items():
            if pattern.match(hgvs_stripped):
                return ValidationResult(
                    status=ValidationStatus.VALID,
                    message=f"Valid HGVS notation ({pattern_name})"
                )

        return ValidationResult(
            status=ValidationStatus.ERROR,
            message=f"Invalid HGVS notation format: {hgvs_stripped}"
        )

    def validate_transcript_format(self, transcript: str) -> ValidationResult:
        """
        Validate RefSeq transcript identifier format.

        Args:
            transcript: Transcript identifier (e.g., NM_002485.4)

        Returns:
            ValidationResult indicating validity
        """
        if not transcript:
            return ValidationResult(
                status=ValidationStatus.ERROR,
                message="Empty transcript identifier"
            )

        if self.REFSEQ_PATTERN.match(transcript):
            return ValidationResult(
                status=ValidationStatus.VALID,
                message="Valid RefSeq transcript format"
            )

        return ValidationResult(
            status=ValidationStatus.ERROR,
            message=f"Invalid RefSeq transcript format: {transcript}. Expected format: NM_XXXXXX or NM_XXXXXX.X"
        )

    def parse_position(self, hgvs: str) -> Tuple[int, int, Optional[int], Optional[int]]:
        """
        Parse CDS positions from HGVS notation.

        Supports c. and g. prefixes, including 3'UTR (*N) and 5'UTR (-N)
        positions.

        Args:
            hgvs: HGVS notation string

        Returns:
            Tuple of (start_pos, end_pos, start_offset, end_offset)
        """
        start_pos = 0
        end_pos = 0
        start_offset = None
        end_offset = None

        hgvs_stripped = hgvs.strip()

        # Handle g. (genomic) notation
        if hgvs_stripped.lower().startswith('g.'):
            for pattern_name, pattern in self.GENOMIC_PATTERNS.items():
                match = pattern.match(hgvs_stripped)
                if match:
                    groups = match.groups()
                    start_pos = int(groups[0])
                    if pattern_name.endswith('_range') or pattern_name in ['g_insertion']:
                        end_pos = int(groups[1])
                    else:
                        end_pos = start_pos
                    break
            return start_pos, end_pos, start_offset, end_offset

        # Handle c. notation
        notation = hgvs_stripped[2:] if hgvs_stripped.lower().startswith('c.') else hgvs_stripped

        # Try each pattern
        for pattern_name, pattern in self.PATTERNS.items():
            match = pattern.match(f"c.{notation}")
            if match:
                groups = match.groups()

                if pattern_name in ['substitution', 'deletion_single', 'duplication_single', 'delins_single']:
                    start_pos = end_pos = self._parse_hgvs_pos(groups[0])
                    if groups[1] and groups[2]:
                        offset_val = int(groups[2])
                        start_offset = end_offset = offset_val if groups[1] == '+' else -offset_val
                elif pattern_name == 'repeat':
                    start_pos = self._parse_hgvs_pos(groups[0])
                    if groups[1] and groups[2]:
                        offset_val = int(groups[2])
                        start_offset = offset_val if groups[1] == '+' else -offset_val
                    # Range repeat: c.123_125[4]
                    if groups[3] is not None:
                        end_pos = self._parse_hgvs_pos(groups[3])
                        if groups[4] and groups[5]:
                            offset_val = int(groups[5])
                            end_offset = offset_val if groups[4] == '+' else -offset_val
                    else:
                        end_pos = start_pos
                else:
                    # Range patterns
                    start_pos = self._parse_hgvs_pos(groups[0])
                    if groups[1] and groups[2]:
                        offset_val = int(groups[2])
                        start_offset = offset_val if groups[1] == '+' else -offset_val

                    end_pos = self._parse_hgvs_pos(groups[3])
                    if groups[4] and groups[5]:
                        offset_val = int(groups[5])
                        end_offset = offset_val if groups[4] == '+' else -offset_val

                break

        return start_pos, end_pos, start_offset, end_offset


class VariantValidator(LoggerMixin):
    """
    Comprehensive validator for variant data.
    Validates HGVS syntax, transcript existence, gene-transcript association,
    and reference sequence concordance.
    """

    def __init__(self, ncbi_client=None, ensembl_client=None, mane_manager=None):
        """
        Initialize the validator with API clients.

        Args:
            ncbi_client: NCBIClient instance for sequence fetching
            ensembl_client: EnsemblClient instance for additional validation
            mane_manager: MANEManager instance for MANE transcript lookup
        """
        self.ncbi_client = ncbi_client
        self.ensembl_client = ensembl_client
        self.mane_manager = mane_manager
        self.hgvs_validator = HGVSValidator()

    def validate_variant(self, variant: Variant, check_sequence: bool = True) -> List[ValidationResult]:
        """
        Perform comprehensive validation of a variant.

        Args:
            variant: Variant to validate
            check_sequence: Whether to check reference sequence concordance

        Returns:
            List of validation results
        """
        results = []

        # 1. Validate HGVS syntax
        hgvs_result = self.hgvs_validator.validate_hgvs_syntax(variant.hgvs_c)
        results.append(hgvs_result)
        variant.add_validation_result(hgvs_result)

        # 2. Validate transcript format
        transcript_result = self.hgvs_validator.validate_transcript_format(
            variant.transcript_accession
        )
        results.append(transcript_result)
        variant.add_validation_result(transcript_result)

        # If basic validation fails, return early
        if any(r.status == ValidationStatus.ERROR for r in results):
            return results

        # 3. Check if transcript exists and is MANE
        if self.mane_manager:
            mane_result = self._validate_mane_transcript(variant)
            results.append(mane_result)
            variant.add_validation_result(mane_result)

        # 4. Validate gene-transcript association
        if self.ncbi_client or self.ensembl_client:
            assoc_result = self._validate_gene_transcript_association(variant)
            results.append(assoc_result)
            variant.add_validation_result(assoc_result)

        # 5. Validate CDS position is within transcript length
        if variant.transcript:
            pos_result = self._validate_position_in_transcript(variant)
            results.append(pos_result)
            variant.add_validation_result(pos_result)

        # 6. Validate reference sequence concordance
        if check_sequence and variant.transcript and variant.ref_allele:
            seq_result = self._validate_reference_sequence(variant)
            results.append(seq_result)
            variant.add_validation_result(seq_result)

        return results

    def _validate_mane_transcript(self, variant: Variant) -> ValidationResult:
        """
        Validate that the transcript is a MANE transcript.
        """
        if not self.mane_manager:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message="MANE validation skipped - no MANE manager available"
            )

        # Extract base accession without version
        accession = variant.transcript_accession.split('.')[0]

        mane_info = self.mane_manager.get_mane_info(accession, variant.gene_symbol)

        if mane_info is None:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Transcript {variant.transcript_accession} is not a MANE transcript for gene {variant.gene_symbol}. "
                        f"Please verify this is the correct transcript.",
                details={'mane_type': MANEType.NOT_MANE.value}
            )

        mane_type = mane_info.get('type', MANEType.NOT_MANE)
        return ValidationResult(
            status=ValidationStatus.VALID,
            message=f"Transcript is {mane_type.value}",
            details={'mane_type': mane_type.value, 'mane_info': mane_info}
        )

    def _validate_gene_transcript_association(self, variant: Variant) -> ValidationResult:
        """
        Validate that the transcript belongs to the specified gene.
        """
        # Try NCBI first
        if self.ncbi_client:
            transcript_info = self.ncbi_client.get_transcript_info(variant.transcript_accession)
            if transcript_info:
                gene_symbol = transcript_info.get('gene_symbol', '').upper()
                if gene_symbol and gene_symbol == variant.gene_symbol.upper():
                    # Safely parse accession and version
                    parts = variant.transcript_accession.split('.')
                    base_accession = parts[0]
                    version = None
                    if len(parts) > 1:
                        try:
                            version = int(parts[1])
                        except ValueError:
                            version = None

                    # Import Strand from models
                    from .models import Strand

                    variant.transcript = Transcript(
                        accession=base_accession,
                        version=version,
                        gene_symbol=gene_symbol,
                        chromosome=transcript_info.get('chromosome', ''),
                        strand=Strand.POSITIVE if transcript_info.get('strand', '+') == '+' else Strand.NEGATIVE,
                        cds_start=transcript_info.get('cds_start', 0),
                        cds_end=transcript_info.get('cds_end', 0),
                        sequence=transcript_info.get('sequence', '')
                    )
                    return ValidationResult(
                        status=ValidationStatus.VALID,
                        message=f"Transcript {variant.transcript_accession} belongs to gene {variant.gene_symbol}"
                    )
                elif gene_symbol:
                    return ValidationResult(
                        status=ValidationStatus.ERROR,
                        message=f"Transcript {variant.transcript_accession} belongs to gene {gene_symbol}, not {variant.gene_symbol}"
                    )

        # Try Ensembl as fallback
        if self.ensembl_client:
            transcript_info = self.ensembl_client.get_transcript_info(variant.transcript_accession)
            if transcript_info:
                gene_symbol = transcript_info.get('gene_symbol', '').upper()
                if gene_symbol and gene_symbol == variant.gene_symbol.upper():
                    return ValidationResult(
                        status=ValidationStatus.VALID,
                        message=f"Transcript {variant.transcript_accession} belongs to gene {variant.gene_symbol}"
                    )
                elif gene_symbol:
                    return ValidationResult(
                        status=ValidationStatus.ERROR,
                        message=f"Transcript {variant.transcript_accession} belongs to gene {gene_symbol}, not {variant.gene_symbol}"
                    )

        return ValidationResult(
            status=ValidationStatus.WARNING,
            message=f"Could not verify gene-transcript association for {variant.transcript_accession}"
        )

    def _validate_position_in_transcript(self, variant: Variant) -> ValidationResult:
        """
        Validate that the variant position is within the transcript CDS length.

        5'UTR (negative positions) and 3'UTR (*N / sentinel positions) are
        allowed with a WARNING -- they are valid HGVS but outside CDS.
        """
        if not variant.transcript:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message="Cannot validate position - transcript data not loaded"
            )

        cds_length = variant.transcript.cds_length
        if cds_length <= 0:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message="Cannot validate position - CDS length unknown"
            )

        # Detect 3'UTR sentinel positions (set by _parse_hgvs_pos for *N)
        is_3utr = variant.cds_start > 900000000 or variant.cds_end > 900000000

        # Detect 5'UTR positions (negative CDS positions)
        is_5utr = variant.cds_start < 0 or variant.cds_end < 0

        if is_3utr:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Variant is in the 3'UTR (downstream of stop codon). "
                        f"Primer design will use genomic coordinates.",
                details={'region': '3_prime_utr'}
            )

        if is_5utr:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Variant is in the 5'UTR (upstream of start codon, position {variant.cds_start}). "
                        f"Primer design will use genomic coordinates.",
                details={'region': '5_prime_utr'}
            )

        # Check if position is within CDS
        if variant.cds_start < 1 or variant.cds_end > cds_length:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Position {variant.cds_start}-{variant.cds_end} is outside CDS range (1-{cds_length}). "
                        f"The variant may be in a UTR region."
            )

        return ValidationResult(
            status=ValidationStatus.VALID,
            message=f"Position {variant.cds_start}-{variant.cds_end} is within CDS (1-{cds_length})"
        )

    def _validate_reference_sequence(self, variant: Variant) -> ValidationResult:
        """
        Validate that the reference allele matches the transcript sequence.
        """
        if not variant.transcript or not variant.transcript.sequence:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message="Cannot validate reference sequence - sequence not available"
            )

        if not variant.ref_allele:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message="No reference allele specified - cannot validate sequence"
            )

        sequence = variant.transcript.sequence

        # Extract the reference sequence at the variant position
        # Note: CDS positions are 1-based, Python strings are 0-based
        start_idx = variant.cds_start - 1
        end_idx = variant.cds_end

        if start_idx < 0 or end_idx > len(sequence):
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Position out of range for sequence extraction"
            )

        actual_ref = sequence[start_idx:end_idx].upper()
        expected_ref = variant.ref_allele.upper()

        if actual_ref == expected_ref:
            return ValidationResult(
                status=ValidationStatus.VALID,
                message=f"Reference allele '{expected_ref}' matches transcript sequence"
            )
        else:
            return ValidationResult(
                status=ValidationStatus.WARNING,
                message=f"Reference allele mismatch: expected '{expected_ref}' but found '{actual_ref}' at position {variant.cds_start}-{variant.cds_end}. "
                        f"The variant may still be processed, but please verify the notation.",
                details={
                    'expected': expected_ref,
                    'found': actual_ref,
                    'position': f"{variant.cds_start}-{variant.cds_end}"
                }
            )

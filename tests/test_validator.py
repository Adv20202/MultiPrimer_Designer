"""
Unit tests for the variant validator module.
"""

import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.core.validator import HGVSValidator, VariantValidator
from src.core.models import Variant, ValidationStatus


class TestHGVSValidator(unittest.TestCase):
    """Tests for HGVS notation validation."""

    def setUp(self):
        self.validator = HGVSValidator()

    def test_valid_substitution(self):
        """Test valid substitution notation."""
        result = self.validator.validate_hgvs_syntax("c.657A>G")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_deletion_single(self):
        """Test valid single nucleotide deletion."""
        result = self.validator.validate_hgvs_syntax("c.657del")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_deletion_range(self):
        """Test valid range deletion."""
        result = self.validator.validate_hgvs_syntax("c.657_661delACAAA")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_duplication(self):
        """Test valid duplication."""
        result = self.validator.validate_hgvs_syntax("c.657dupA")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_insertion(self):
        """Test valid insertion."""
        result = self.validator.validate_hgvs_syntax("c.657_658insACG")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_delins(self):
        """Test valid deletion-insertion."""
        result = self.validator.validate_hgvs_syntax("c.657_661delACAAinsGT")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_missing_c_prefix(self):
        """Test notation without c. prefix."""
        result = self.validator.validate_hgvs_syntax("657A>G")
        self.assertEqual(result.status, ValidationStatus.ERROR)

    def test_invalid_format(self):
        """Test invalid notation format."""
        result = self.validator.validate_hgvs_syntax("c.invalid")
        self.assertEqual(result.status, ValidationStatus.ERROR)

    def test_empty_notation(self):
        """Test empty notation."""
        result = self.validator.validate_hgvs_syntax("")
        self.assertEqual(result.status, ValidationStatus.ERROR)

    def test_valid_transcript_format(self):
        """Test valid RefSeq transcript format."""
        result = self.validator.validate_transcript_format("NM_002485.4")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_valid_transcript_without_version(self):
        """Test valid RefSeq transcript without version."""
        result = self.validator.validate_transcript_format("NM_002485")
        self.assertEqual(result.status, ValidationStatus.VALID)

    def test_invalid_transcript_format(self):
        """Test invalid transcript format."""
        result = self.validator.validate_transcript_format("INVALID")
        self.assertEqual(result.status, ValidationStatus.ERROR)

    def test_parse_position_substitution(self):
        """Test parsing position from substitution."""
        start, end, start_off, end_off = self.validator.parse_position("c.657A>G")
        self.assertEqual(start, 657)
        self.assertEqual(end, 657)

    def test_parse_position_range(self):
        """Test parsing position from range deletion."""
        start, end, start_off, end_off = self.validator.parse_position("c.657_661delACAAA")
        self.assertEqual(start, 657)
        self.assertEqual(end, 661)


class TestVariantValidator(unittest.TestCase):
    """Tests for comprehensive variant validation."""

    def setUp(self):
        self.validator = VariantValidator()

    def test_validate_basic_variant(self):
        """Test basic variant validation (without API calls)."""
        variant = Variant(
            row_number=1,
            gene_symbol="BRCA1",
            transcript_accession="NM_007294.4",
            hgvs_c="c.116G>A"
        )

        results = self.validator.validate_variant(variant, check_sequence=False)

        # Should have at least HGVS and transcript format validation
        self.assertTrue(len(results) >= 2)


if __name__ == '__main__':
    unittest.main()

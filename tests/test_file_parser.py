"""
Unit tests for the file parser module.
"""

import unittest
import tempfile
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.utils.file_parser import FileParser
from src.core.models import VariantType


class TestFileParser(unittest.TestCase):
    """Tests for file parsing functionality."""

    def setUp(self):
        self.parser = FileParser()

    def test_parse_csv_comma_delimiter(self):
        """Test parsing CSV with comma delimiter."""
        content = """Gene,Transcript,Position
BRCA1,NM_007294.4,c.116G>A
CHEK2,NM_007194.4,c.470T>C"""

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(len(result.variants), 2)
                self.assertEqual(result.variants[0].gene_symbol, "BRCA1")
                self.assertEqual(result.variants[0].transcript_accession, "NM_007294.4")
                self.assertEqual(result.variants[0].hgvs_c, "c.116G>A")

            finally:
                os.unlink(f.name)

    def test_parse_csv_tab_delimiter(self):
        """Test parsing CSV with tab delimiter."""
        content = "Gene\tTranscript\tPosition\nBRCA1\tNM_007294.4\tc.116G>A"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(len(result.variants), 1)
                self.assertEqual(result.variants[0].gene_symbol, "BRCA1")

            finally:
                os.unlink(f.name)

    def test_parse_csv_semicolon_delimiter(self):
        """Test parsing CSV with semicolon delimiter."""
        content = "Gene;Transcript;Position\nBRCA1;NM_007294.4;c.116G>A"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(len(result.variants), 1)
                self.assertEqual(result.variants[0].gene_symbol, "BRCA1")

            finally:
                os.unlink(f.name)

    def test_parse_variant_type_substitution(self):
        """Test parsing substitution variant."""
        content = "Gene,Transcript,Position\nBRCA1,NM_007294.4,c.116G>A"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(result.variants[0].variant_type, VariantType.SUBSTITUTION)
                self.assertEqual(result.variants[0].cds_start, 116)
                self.assertEqual(result.variants[0].cds_end, 116)
                self.assertEqual(result.variants[0].ref_allele, "G")
                self.assertEqual(result.variants[0].alt_allele, "A")

            finally:
                os.unlink(f.name)

    def test_parse_variant_type_deletion(self):
        """Test parsing deletion variant."""
        content = "Gene,Transcript,Position\nNBN,NM_002485.4,c.657_661delACAAA"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(result.variants[0].variant_type, VariantType.DELETION)
                self.assertEqual(result.variants[0].cds_start, 657)
                self.assertEqual(result.variants[0].cds_end, 661)
                self.assertEqual(result.variants[0].ref_allele, "ACAAA")

            finally:
                os.unlink(f.name)

    def test_parse_html_entities(self):
        """Test parsing with HTML entities (e.g., &gt;)."""
        content = "Gene,Transcript,Position\nBRCA1,NM_007294.4,c.116G&gt;A"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            f.flush()

            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(result.variants[0].hgvs_c, "c.116G>A")

            finally:
                os.unlink(f.name)

    def test_missing_file(self):
        """Test handling of missing file."""
        result = self.parser.parse_file("/nonexistent/file.csv")

        self.assertEqual(len(result.variants), 0)
        self.assertTrue(len(result.errors) > 0)

    def test_unsupported_format(self):
        """Test handling of unsupported file format."""
        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            try:
                result = self.parser.parse_file(f.name)

                self.assertEqual(len(result.variants), 0)
                self.assertTrue(len(result.errors) > 0)

            finally:
                os.unlink(f.name)


if __name__ == '__main__':
    unittest.main()

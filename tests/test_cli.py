"""Run with: python3 tests/test_cli.py /absolute/path/to/IEContig"""

import pathlib
import subprocess
import sys
import tempfile
import unittest


BINARY = str(pathlib.Path(sys.argv.pop(1)).resolve())
NT = "ATGGCTTGTGATGAATTTGGTCATATTAAACTG"
AA = "MACDEFGHIKL"


class CliTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="iecontig-test-")
        self.addCleanup(self.directory.cleanup)
        self.root = pathlib.Path(self.directory.name)
        self.input = self.root / "input.fasta"
        self.reference = self.root / "reference.fasta"
        self.output = self.root / "output.fasta"

    def run_cli(self, sequence=NT, reference=AA, input_type="NT", output_type="AA",
                threshold="0", input_path=None):
        self.input.write_text("" if sequence is None else f">contig\n{sequence}\n")
        self.reference.write_text("" if reference is None else f">reference\n{reference}\n")
        return subprocess.run(
            [BINARY, "-i", str(self.input if input_path is None else input_path), "-t", input_type,
             "-I", str(self.reference), "-o", str(self.output),
             "-y", output_type, "-s", threshold],
            capture_output=True, text=True, timeout=30,
        )

    def assert_output(self, result, sequence):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), f">contig\n{sequence}\n")

    def assert_input_error(self, **kwargs):
        # Invalid input must fail normally before touching an existing result.
        self.output.write_text("existing result\n")
        result = self.run_cli(**kwargs)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stderr)
        self.assertNotIn("terminate called", result.stderr)
        self.assertEqual(self.output.read_text(), "existing result\n")
        return result

    def test_nt_to_aa(self):
        self.assert_output(self.run_cli(), AA)

    def test_lowercase_input(self):
        self.assert_output(self.run_cli(sequence=NT.lower()), AA)

    def test_lowercase_reference(self):
        self.assert_output(self.run_cli(reference=AA.lower()), AA)

    def test_empty_input(self):
        result = self.assert_input_error(sequence=None)
        self.assertIn("Input file contains no sequences", result.stderr)

    def test_empty_reference(self):
        result = self.assert_input_error(reference=None)
        self.assertIn("Reference file contains no sequences", result.stderr)

    def test_empty_sequence_records(self):
        for kwargs in ({"sequence": ""}, {"reference": ""},
                       {"sequence": "---"}, {"reference": "-.+"}):
            with self.subTest(kwargs=kwargs):
                self.assert_input_error(**kwargs)

    def test_invalid_type(self):
        self.assert_input_error(input_type="BAD")

    def test_invalid_threshold(self):
        self.assert_input_error(threshold="100")

    def test_missing_required_option(self):
        result = subprocess.run([BINARY, "-t", "NT"], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("required", result.stderr)
        self.assertNotIn("terminate called", result.stderr)

    def test_missing_input_file(self):
        result = self.assert_input_error(input_path=self.root / "missing.fasta")
        self.assertIn("missing.fasta", result.stderr)

    def test_aa_input(self):
        # Flanks distinguish direct protein alignment from six-frame translation.
        self.assert_output(self.run_cli(sequence="PPP" + AA + "PPP", input_type="AA"), AA)

    def test_lowercase_aa_input(self):
        self.assert_output(self.run_cli(sequence=AA.lower(), reference=AA.lower(), input_type="AA"), AA)

    def test_aa_to_nt_rejected(self):
        result = self.assert_input_error(sequence=AA, input_type="AA", output_type="NT")
        self.assertIn("AA input supports only AA output", result.stderr)

    def test_six_reading_frames(self):
        complement = str.maketrans("ACGT", "TGCA")
        for reverse in (False, True):
            for offset in range(3):
                sequence = "A" * offset + NT
                if reverse:
                    sequence = sequence.translate(complement)[::-1]
                for output_type, expected in (("NT", NT), ("AA", AA)):
                    with self.subTest(reverse=reverse, offset=offset, output_type=output_type):
                        self.assert_output(self.run_cli(sequence=sequence, output_type=output_type), expected)

    def test_multiple_references_and_inputs(self):
        sequence = NT + "\n>second\n" + NT.lower()
        reference = "PPPPPPPPPPP\n>matching_reference\n" + AA
        result = self.run_cli(sequence=sequence, reference=reference)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), f">contig\n{AA}\n>second\n{AA}\n")

    def test_short_nt_input(self):
        result = self.run_cli(sequence="AT")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), "")

    def test_malformed_record(self):
        self.input.write_text("@broken\nATGC\n+\nA\n")
        self.reference.write_text(f">reference\n{AA}\n")
        result = subprocess.run(
            [BINARY, "-i", str(self.input), "-t", "NT", "-I", str(self.reference),
             "-o", str(self.output), "-y", "AA", "-s", "0"],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("Failed to read sequence", result.stderr)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()

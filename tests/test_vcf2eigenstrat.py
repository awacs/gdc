"""End-to-end regression tests; run with python -m unittest discover -s tests."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


SCRIPT = Path(__file__).resolve().parents[1] / "vcf2eigenstrat.py"
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
"""


def row(pos, ref="A", alt="C", info="AA=A", gt="0|1", chrom="chr1", name="."):
    return f"{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{gt}\n"


class ConverterTests(unittest.TestCase):
    def convert(self, records, *args, filetype="vcf", success=True):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "input.vcf"
            source.write_text(HEADER + "".join(records))
            if filetype != "vcf":
                destination = directory / ("input." + filetype)
                with pysam.VariantFile(str(source)) as reader:
                    mode = "wb" if filetype == "bcf" else "wz"
                    with pysam.VariantFile(str(destination), mode, header=reader.header) as writer:
                        for record in reader:
                            writer.write(record)
                source = destination
            output = directory / "out"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "-v", str(source), "-o", str(output), *args],
                capture_output=True, text=True,
            )
            if success:
                self.assertEqual(result.returncode, 0, result.stderr)
            else:
                self.assertNotEqual(result.returncode, 0)
            contents = {
                suffix: output.with_suffix("." + suffix).read_text().splitlines()
                for suffix in ("snp", "geno", "ind")
            }
            return result, contents

    def test_split_and_unsplit_multiallelic_positions_in_all_formats(self):
        records = [
            row(1),
            row(10, info="AC=0;AA=A", gt="0|0"),
            row(10, alt="G", info="AC=0;AA=A", gt="0|0"),
            row(20, alt="C,G", gt="1|2"),
            row(30, alt="T", gt="1|0"),
            row(40),
            row(40, alt="G"),  # Also exercise a split position at EOF.
        ]
        for filetype in ("vcf", "vcf.gz", "bcf"):
            with self.subTest(filetype=filetype):
                result, files = self.convert(records, "--phased", filetype=filetype)
                self.assertEqual(files["snp"], [
                    "chr1:1\tchr1\t0.0\t1\tA\tC",
                    "chr1:30\tchr1\t0.0\t30\tA\tT",
                ])
                self.assertEqual(files["geno"], ["01", "10"])
                self.assertEqual(files["ind"], ["sample\tU\tPOP"])
                self.assertIn("Excluded 5 multiallelic", result.stdout)

    def test_polarization_and_existing_ids_preserved(self):
        _, files = self.convert([
            row(1, info="AA=c|||", name="rsExample"),
            row(2, alt="G", gt="0|0"),
        ], "--aa", "--phased")
        self.assertEqual(files["snp"], [
            "rsExample\tchr1\t0.0\t1\tC\tA",
            "chr1:2\tchr1\t0.0\t2\tA\tG",
        ])
        self.assertEqual(files["geno"], ["10", "00"])

    def test_standard_genotypes_and_same_position_on_different_chromosomes(self):
        _, files = self.convert([
            row(1, gt="0/0"), row(2, gt="0/1"), row(3, gt="1/1"),
            row(4, gt="./."), row(1, alt="G", chrom="chr2"),
        ])
        self.assertEqual(files["geno"], ["2", "1", "0", "9", "1"])
        self.assertEqual(files["snp"][-1], "chr2:1\tchr2\t0.0\t1\tA\tG")

    def test_multiallelic_filter_precedes_aa_and_indel_filters(self):
        result, files = self.convert([
            row(1), row(10), row(10, alt="G", info="."),
            row(20), row(20, alt="AT"), row(30, alt="AT"),
        ], "--aa")
        self.assertEqual(len(files["snp"]), 1)
        self.assertIn("Excluded 4 multiallelic", result.stdout)
        self.assertIn("Excluded 1 indel", result.stdout)
        self.assertNotIn("unknownAA", result.stdout)

    def test_identical_duplicates_are_not_multiallelic(self):
        result, files = self.convert([row(1), row(1)])
        self.assertEqual(len(files["snp"]), 2)
        self.assertIn("Excluded 0 records total", result.stdout)

    def test_unsorted_positions_rejected(self):
        result, _ = self.convert([row(2), row(1)], success=False)
        self.assertIn("Input must be coordinate-sorted", result.stderr)

    def test_repeated_contig_rejected(self):
        result, _ = self.convert([
            row(1), row(1, chrom="chr2"), row(2),
        ], success=False)
        self.assertIn("Input must be coordinate-sorted", result.stderr)


if __name__ == "__main__":
    unittest.main()

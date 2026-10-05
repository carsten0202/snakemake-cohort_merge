import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
BCFTOOLS = shutil.which("bcftools")


def write_vcf(path, chromosome, records, samples):
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1>",
        "##contig=<ID=chrX>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                *samples,
            ]
        ),
    ]

    for position, identifier, ref, alt, genotypes in records:
        lines.append(
            "\t".join(
                [
                    chromosome,
                    str(position),
                    identifier,
                    ref,
                    alt,
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    *genotypes,
                ]
            )
        )

    path.write_text("\n".join(lines) + "\n")


def compress_vcf(source, destination):
    subprocess.run(
        [
            BCFTOOLS,
            "view",
            "--output-type",
            "z",
            "--output",
            str(destination),
            str(source),
        ],
        check=True,
        text=True,
        capture_output=True,
    )


def build_test_project(root, mismatch=False):
    cohort_a = root / "cohort_a"
    cohort_b = root / "cohort_b"
    output = root / "output"

    cohort_a.mkdir()
    cohort_b.mkdir()

    chr1_a = cohort_a / "chr1.vcf"
    chr_x_a = cohort_a / "chrX.vcf"
    chr1_b = cohort_b / "chr1.vcf"
    chr_x_b = cohort_b / "chrX.vcf"

    write_vcf(
        chr1_a,
        "chr1",
        [
            (100, "rs1", "A", "G", ["0/1", "0/0"]),
            (200, "rs2", "C", "T", ["1/1", "0/1"]),
        ],
        ["duplicate", "cohort_a_sample"],
    )
    write_vcf(
        chr_x_a,
        "chrX",
        [(100, "rsX1", "A", "C", ["0/1", "0/0"])],
        ["duplicate", "cohort_a_sample"],
    )

    cohort_b_alt = "T" if mismatch else "G"
    write_vcf(
        chr1_b,
        "chr1",
        [
            (100, "rs1", "A", cohort_b_alt, ["0/0", "0/1"]),
            (200, "rs2", "C", "T", ["0/1", "1/1"]),
        ],
        ["duplicate", "cohort_b_sample"],
    )
    write_vcf(
        chr_x_b,
        "chrX",
        [(100, "rsX1", "A", "C", ["0/0", "0/1"])],
        ["duplicate", "cohort_b_sample"],
    )

    source_vcfs = [
        cohort_a / "chr1.vcf.gz",
        cohort_a / "chrX.vcf.gz",
        cohort_b / "chr1.filtered.vcf.gz",
        cohort_b / "chr23.filtered.vcf.gz",
    ]

    compress_vcf(chr1_a, source_vcfs[0])
    compress_vcf(chr_x_a, source_vcfs[1])
    compress_vcf(chr1_b, source_vcfs[2])
    compress_vcf(chr_x_b, source_vcfs[3])

    config = {
        "super_cohort": "Glostrup",
        "output_directory": str(output),
        "chromosomes": ["chr1", "chrX"],
        "cohorts": {
            "cohort_a": {
                "directory": str(cohort_a),
            },
            "cohort_b": {
                "directory": str(cohort_b),
                "pattern": "{chromosome}.filtered.vcf.gz",
                "chromosome_aliases": {
                    "chrX": "chr23",
                },
            },
        },
    }

    config_path = root / "config.json"
    config_path.write_text(json.dumps(config))

    return config_path, output, source_vcfs


def run_workflow(config):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--configfile",
            str(config),
            "--cores",
            "4",
        ],
        cwd=REPOSITORY,
        text=True,
        capture_output=True,
    )


@unittest.skipUnless(BCFTOOLS, "bcftools is required for integration tests")
class WorkflowIntegrationTest(unittest.TestCase):
    def test_complete_workflow(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, source_vcfs = build_test_project(
                Path(temporary_directory)
            )

            result = run_workflow(config)
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            for vcf in source_vcfs:
                self.assertTrue(Path(f"{vcf}.tbi").is_file())

            chromosome_outputs = [
                output / "by_chromosome" / "Glostrup.chr1.vcf.gz",
                output / "by_chromosome" / "Glostrup.chrX.vcf.gz",
            ]
            for vcf in chromosome_outputs:
                self.assertTrue(vcf.is_file())
                self.assertTrue(Path(f"{vcf}.tbi").is_file())

            final_vcf = output / "Glostrup.vcf.gz"
            self.assertTrue(final_vcf.is_file())
            self.assertTrue(Path(f"{final_vcf}.tbi").is_file())

            samples = subprocess.run(
                [BCFTOOLS, "query", "--list-samples", str(final_vcf)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout.splitlines()
            self.assertEqual(
                samples,
                [
                    "duplicate",
                    "cohort_a_sample",
                    "2:duplicate",
                    "cohort_b_sample",
                ],
            )

            records = subprocess.run(
                [
                    BCFTOOLS,
                    "query",
                    "--format",
                    "%CHROM\t%POS\t%REF\t%ALT\n",
                    str(final_vcf),
                ],
                check=True,
                text=True,
                capture_output=True,
            ).stdout.splitlines()
            self.assertEqual(
                records,
                [
                    "chr1\t100\tA\tG",
                    "chr1\t200\tC\tT",
                    "chrX\t100\tA\tC",
                ],
            )

            index_stats = subprocess.run(
                [BCFTOOLS, "index", "--stats", str(final_vcf)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout.splitlines()
            self.assertEqual(index_stats, ["chr1\t.\t2", "chrX\t.\t1"])

            self.assertTrue(
                (output / "qc" / "variant_sets" / "chr1.sha256.tsv").is_file()
            )
            self.assertTrue(
                (output / "qc" / "variant_sets" / "chrX.sha256.tsv").is_file()
            )

    def test_mismatched_variants_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                mismatch=True,
            )

            result = run_workflow(config)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Variant records differ between cohorts",
                result.stdout + result.stderr,
            )
            self.assertFalse((output / "Glostrup.vcf.gz").exists())


if __name__ == "__main__":
    unittest.main()

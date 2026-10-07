import json
import os
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


def build_test_project(root, mismatch=False, empty_intersection=False):
    cohort_a = root / "cohort_a"
    cohort_b = root / "cohort_b"
    output = root / "output"

    cohort_a.mkdir()
    cohort_b.mkdir()

    chr1_a = cohort_a / "chr1.vcf"
    chr_x_a = cohort_a / "chrX.vcf"
    chr1_b = cohort_b / "chr1.vcf"
    chr_x_b = cohort_b / "chrX.vcf"

    chr1_a_records = [
        (100, "rs1", "A", "G", ["0/1", "0/0"]),
        (200, "rs2", "C", "T", ["1/1", "0/1"]),
    ]
    if mismatch:
        chr1_a_records.append((300, "rs3", "G", "A", ["0/1", "0/0"]))
    write_vcf(
        chr1_a,
        "chr1",
        chr1_a_records,
        ["duplicate", "cohort_a_sample"],
    )
    write_vcf(
        chr_x_a,
        "chrX",
        [(100, "rsX1", "A", "C", ["0/1", "0/0"])],
        ["duplicate", "cohort_a_sample"],
    )

    cohort_b_first_alt = "T" if mismatch or empty_intersection else "G"
    cohort_b_second_alt = "G" if empty_intersection else "T"
    write_vcf(
        chr1_b,
        "chr1",
        [
            (100, "rs1", "A", cohort_b_first_alt, ["0/0", "0/1"]),
            (200, "rs2", "C", cohort_b_second_alt, ["0/1", "1/1"]),
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


def write_fake_analysis_tools(bin_directory, call_log):
    plink2 = bin_directory / "plink2"
    plink2.write_text(
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
prefix = Path(args[args.index("--out") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
for suffix in (".bed", ".bim", ".fam", ".log"):
    Path(f"{prefix}{suffix}").touch()
with open(os.environ["FAKE_TOOL_CALLS"], "a") as handle:
    handle.write("plink2 " + " ".join(args) + "\\n")
"""
    )
    plink2.chmod(0o755)

    quick_pcair = bin_directory / "quick-pcair"
    quick_pcair.write_text(
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
prefix = Path(args[args.index("--output") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
suffixes = (
    ".gds",
    ".KINGkinship.png",
    ".KINGkinship.tsv",
    ".ldprune.snpids.txt",
    ".eigenvalues.tsv",
    ".eigenvectors.tsv",
    ".pcair_1v2.png",
    ".pcair_3v4.png",
    ".pcrelate.RData",
    ".pcrelate_1v2.png",
)
for suffix in suffixes:
    Path(f"{prefix}{suffix}").touch()
with open(os.environ["FAKE_TOOL_CALLS"], "a") as handle:
    handle.write("quick-pcair " + " ".join(args) + "\\n")
"""
    )
    quick_pcair.chmod(0o755)

    return {
        **os.environ,
        "FAKE_TOOL_CALLS": str(call_log),
        "PATH": f"{bin_directory}{os.pathsep}{os.environ['PATH']}",
    }


def run_workflow(config, targets=None, env=None):
    command = [
        sys.executable,
        "-m",
        "snakemake",
        "--configfile",
        str(config),
        "--cores",
        "4",
    ]
    if targets:
        command.extend(str(target) for target in targets)

    return subprocess.run(
        command,
        cwd=REPOSITORY,
        env=env,
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

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )
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

            chr1_report = (
                output / "qc" / "variant_intersection" / "chr1.tsv"
            )
            chr_x_report = (
                output / "qc" / "variant_intersection" / "chrX.tsv"
            )
            self.assertEqual(
                chr1_report.read_text().splitlines(),
                [
                    "cohort\tinput_records\tshared_records\t"
                    "excluded_records\texcluded_fraction",
                    "cohort_a\t2\t2\t0\t0.000000",
                    "cohort_b\t2\t2\t0\t0.000000",
                ],
            )
            self.assertTrue(chr_x_report.is_file())
            self.assertFalse(
                (output / "intermediate" / "intersection" / "chr1").exists()
            )
            self.assertFalse(
                (output / "intermediate" / "intersection" / "chrX").exists()
            )

    def test_nonshared_variants_are_excluded(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                mismatch=True,
            )

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            final_vcf = output / "Glostrup.vcf.gz"
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
                    "chr1\t200\tC\tT",
                    "chrX\t100\tA\tC",
                ],
            )

            report = output / "qc" / "variant_intersection" / "chr1.tsv"
            self.assertEqual(
                report.read_text().splitlines(),
                [
                    "cohort\tinput_records\tshared_records\t"
                    "excluded_records\texcluded_fraction",
                    "cohort_a\t3\t1\t2\t0.666667",
                    "cohort_b\t2\t1\t1\t0.500000",
                ],
            )

    def test_empty_intersection_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                empty_intersection=True,
            )

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "No variants are shared by all cohorts for chr1",
                result.stdout + result.stderr,
            )
            self.assertFalse((output / "Glostrup.vcf.gz").exists())

    def test_population_structure_workflow(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config, output, _ = build_test_project(root)
            bin_directory = root / "bin"
            bin_directory.mkdir()
            call_log = root / "analysis_calls.log"
            env = write_fake_analysis_tools(bin_directory, call_log)

            result = run_workflow(config, env=env)
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            plink_prefix = output / "analysis" / "plink" / "Glostrup"
            for suffix in (".bed", ".bim", ".fam", ".log"):
                self.assertTrue(Path(f"{plink_prefix}{suffix}").is_file())

            pcair_prefix = output / "analysis" / "pcair" / "Glostrup"
            pcair_suffixes = (
                ".gds",
                ".KINGkinship.png",
                ".KINGkinship.tsv",
                ".ldprune.snpids.txt",
                ".eigenvalues.tsv",
                ".eigenvectors.tsv",
                ".pcair_1v2.png",
                ".pcair_3v4.png",
                ".pcrelate.RData",
                ".pcrelate_1v2.png",
            )
            for suffix in pcair_suffixes:
                self.assertTrue(Path(f"{pcair_prefix}{suffix}").is_file())

            calls = call_log.read_text().splitlines()
            self.assertEqual(len(calls), 2)
            self.assertIn("--vcf-require-gt", calls[0])
            self.assertIn("--double-id", calls[0])
            self.assertIn("--autosome", calls[0])
            self.assertIn("--snps-only just-acgt", calls[0])
            self.assertIn("--max-alleles 2", calls[0])
            self.assertIn(f"--plink {plink_prefix}", calls[1])
            self.assertIn(f"--output {pcair_prefix}", calls[1])


if __name__ == "__main__":
    unittest.main()

import csv
import json
import runpy
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


REPOSITORY = Path(__file__).resolve().parents[2]
BCFTOOLS = shutil.which("bcftools")
RSCRIPT = shutil.which("Rscript")


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


def write_symbolic_vcf(path, samples, flexible_numbers):
    if flexible_numbers:
        af_number, ds_number, gp_number, hds_number = "A", "A", "G", "."
    else:
        af_number, ds_number, gp_number, hds_number = "1", "1", "3", "2"

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1>",
        f"##INFO=<ID=AF,Number={af_number},Type=Float,"
        'Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"##FORMAT=<ID=DS,Number={ds_number},Type=Float,"
        'Description="Dosage">',
        f"##FORMAT=<ID=GP,Number={gp_number},Type=Float,"
        'Description="Genotype probabilities">',
        f"##FORMAT=<ID=HDS,Number={hds_number},Type=Float,"
        'Description="Haplotype dosage">',
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
    sample_values = [
        "0|0:0,0:1,0,0:0" if index % 2 == 0 else "0|1:0,1:0,1,0:1"
        for index in range(len(samples))
    ]
    lines.append(
        "\t".join(
            [
                "chr1",
                "39119",
                "chr1:39119:T:<CN2>",
                "T",
                "<CN2>",
                ".",
                "PASS",
                "AF=0.03",
                "GT:HDS:GP:DS",
                *sample_values,
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


def build_test_project(
    root,
    mismatch=False,
    empty_intersection=False,
    exclusion_entries=None,
):
    cohort_a = root / "cohort_a"
    cohort_b = root / "cohort_b"
    output_root = root / "output"
    output = output_root / "Glostrup"

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
        "output_directory": str(output_root),
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

    if exclusion_entries is not None:
        exclusions = root / "cohort_b_exclusions.txt"
        exclusions.write_text("\n".join(exclusion_entries) + "\n")
        config["cohorts"]["cohort_b"]["exclude_samples"] = str(exclusions)

    config_path = root / "config.json"
    config_path.write_text(json.dumps(config))

    return config_path, output, source_vcfs


def write_fake_analysis_tools(bin_directory, call_log):
    plink2 = bin_directory / "plink2"
    plink2.write_text(
        """#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

args = sys.argv[1:]
prefix = Path(args[args.index("--out") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
for suffix in (".bed", ".bim", ".fam", ".log"):
    Path(f"{prefix}{suffix}").touch()
samples = subprocess.run(
    ["bcftools", "query", "--list-samples", args[args.index("--vcf") + 1]],
    check=True, text=True, capture_output=True,
).stdout
Path(f"{prefix}.fam").write_text(samples)
with open(os.environ["FAKE_TOOL_CALLS"], "a") as handle:
    handle.write("plink2 " + " ".join(args) + "\\n")
"""
    )
    plink2.chmod(0o755)

    quick_pcair = bin_directory / "quick-pcair"
    quick_pcair.write_text(
        """#!/usr/bin/env python3
import os
import subprocess
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
king_lines = (
    " duplicate cohort_a_sample 2:duplicate cohort_b_sample\\n"
    "duplicate 0.5 0.354 0.177 0.0884\\n"
    "cohort_a_sample 0.354 0.5 0.0442 0.02\\n"
    "2:duplicate 0.177 0.0442 0.5 0.01\\n"
    "cohort_b_sample 0.0884 0.02 0.01 0.5\\n"
).splitlines()
fam = args[args.index("--plink") + 1] + ".fam"
samples = Path(fam).read_text().splitlines()
all_samples = king_lines[0].split()
indices = [all_samples.index(sample) for sample in samples]
king_rows = [" " + " ".join(samples)]
for index in indices:
    fields = king_lines[index + 1].split()
    king_rows.append(" ".join([fields[0]] + [fields[i + 1] for i in indices]))
Path(f"{prefix}.KINGkinship.tsv").write_text("\\n".join(king_rows) + "\\n")
subprocess.run(
    [
        "Rscript",
        "-e",
        (
            "args <- commandArgs(TRUE); "
            "mypcrelate <- list(kinBtwn=data.frame("
            "ID1=c('duplicate','duplicate','duplicate','cohort_a_sample',"
            "'cohort_a_sample','2:duplicate'),"
            "ID2=c('cohort_a_sample','2:duplicate','cohort_b_sample',"
            "'2:duplicate','cohort_b_sample','cohort_b_sample'),"
            "kin=c(0.36,0.17,0.04,0.03,0.05,0.01))); "
            "if (nzchar(args[3])) mypcrelate$kinBtwn$kin <- "
            "as.numeric(strsplit(args[3], ',', fixed=TRUE)[[1]]); "
            "ids <- readLines(args[2]); "
            "mypcrelate$kinBtwn <- subset(mypcrelate$kinBtwn, ID1 %in% ids & ID2 %in% ids); "
            "save(mypcrelate, file=args[1])"
        ),
        f"{prefix}.pcrelate.RData",
        fam,
        os.environ.get("FAKE_PCRELATE_KINSHIPS", ""),
    ],
    check=True,
)
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


@unittest.skipUnless(
    BCFTOOLS and RSCRIPT,
    "bcftools and Rscript are required for integration tests",
)
class WorkflowIntegrationTest(unittest.TestCase):
    def test_manifest_handles_existing_prefixed_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config, output, source_vcfs = build_test_project(root)
            data = json.loads(config.read_text())
            data["deduplication"] = {"enabled": True}
            for priority, cohort in enumerate(data["cohorts"].values(), 1):
                cohort["priority"] = priority
            config.write_text(json.dumps(data))
            for chromosome, destination in zip(("chr1", "chrX"), source_vcfs[:2]):
                source = root / "cohort_a" / f"{chromosome}.vcf"
                source.write_text(source.read_text().replace("cohort_a_sample", "2:duplicate"))
                compress_vcf(source, destination)
            qc = output / "qc" / "deduplication"
            result = run_workflow(config, [qc / "sample_manifest.tsv"])
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (qc / "sample_manifest.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(row["merged_sample_id"], row["original_sample_id"], row["cohort"]) for row in rows],
                [("duplicate", "duplicate", "cohort_a"),
                 ("2:duplicate", "2:duplicate", "cohort_a"),
                 ("2:2:duplicate", "duplicate", "cohort_b"),
                 ("cohort_b_sample", "cohort_b_sample", "cohort_b")],
            )

    def test_deduplication(self):
        scenarios = [
            # Cohort-a twins: seeded selection, without affecting unrelated samples.
            ("within_cohort", "0.36,0.17,0.04,0.03,0.05,0.01", (1, 2), None,
             {"duplicate", "cohort_a_sample"}, 1, True),
            # Cohort b wins despite appearing second in merge order; exact boundary.
            ("cross_cohort", "0.353999,0.354,0.04,0.03,0.05,0.01", (2, 1), None,
             {"duplicate", "2:duplicate"}, 1, True),
            # A--C--B is one component; choose between the two cohort-a members.
            ("chain", "0.1,0.354,0.04,0.36,0.05,0.01", (1, 2), None,
             {"duplicate", "cohort_a_sample", "2:duplicate"}, 2, False),
            # KING's duplicate classification and non-finite PC-Relate values do not remove samples.
            ("no_duplicates", "NA,Inf,0.04,NaN,0.05,-Inf", (1, 2), None,
             set(), 0, True),
            ("existing_exclusions", "0.1,0.36,0.04,0.03,0.05,0.01", (2, 1),
             ["cohort_b_sample"], {"duplicate", "2:duplicate"}, 1, True),
        ]
        for name, kinships, priorities, exclusions, members, count, connected in scenarios:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                config, output, _ = build_test_project(root, exclusion_entries=exclusions)
                data = json.loads(config.read_text())
                data["deduplication"] = {"enabled": True, "random_seed": 42}
                # The report cutoff must not hide deduplication candidates.
                data["relatedness"] = {"minimum_kinship": 0.5}
                for cohort, priority in zip(data["cohorts"].values(), priorities):
                    cohort["priority"] = priority
                config.write_text(json.dumps(data))
                bin_directory = root / "bin"
                bin_directory.mkdir()
                env = write_fake_analysis_tools(bin_directory, root / "calls.log")
                env["FAKE_PCRELATE_KINSHIPS"] = kinships
                result = run_workflow(config, env=env)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

                qc = output / "qc" / "deduplication"
                with (qc / "sample_manifest.tsv").open() as handle:
                    manifest = list(csv.DictReader(handle, delimiter="\t"))
                mapping = {row["merged_sample_id"]: row for row in manifest}
                self.assertEqual(mapping["2:duplicate"]["original_sample_id"], "duplicate")
                self.assertEqual(mapping["2:duplicate"]["cohort"], "cohort_b")
                self.assertEqual("cohort_b_sample" in mapping, exclusions is None)
                with (qc / "sample_decisions.tsv").open() as handle:
                    decisions = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual({row["merged_sample_id"] for row in decisions}, members)
                self.assertTrue(all(
                    row["fully_pairwise_connected"] == str(connected).lower()
                    for row in decisions
                ))
                removals = (qc / "removed_samples.txt").read_text().splitlines()
                self.assertEqual(len(removals), count)
                if decisions:
                    kept = [row for row in decisions if row["decision"] == "keep"]
                    self.assertEqual(len(kept), 1)
                    winner = kept[0]
                    self.assertEqual(int(winner["priority"]), min(
                        int(row["priority"]) for row in decisions
                    ))
                    self.assertTrue(all(
                        row["retained_representative"] == winner["merged_sample_id"]
                        for row in decisions
                    ))
                    if name in ("cross_cohort", "existing_exclusions"):
                        self.assertEqual(winner["merged_sample_id"], "2:duplicate")
                    if name in ("within_cohort", "chain"):
                        self.assertEqual(winner["selection_reason"], "seeded_random_selection")

                # Repeat selection with reversed input rows: seed and canonical ordering
                # must preserve both the representative and the decision report.
                expected_decisions = (qc / "sample_decisions.tsv").read_text()
                for filename in ("sample_manifest.tsv", "duplicate_pairs.tsv"):
                    path = qc / filename
                    lines = path.read_text().splitlines()
                    path.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")
                runpy.run_path(str(REPOSITORY / "workflow/scripts/select_duplicate_representatives.py"),
                    run_name="__main__", init_globals={"snakemake": SimpleNamespace(
                        input=SimpleNamespace(manifest=qc / "sample_manifest.tsv", pairs=qc / "duplicate_pairs.tsv"),
                        output=SimpleNamespace(decisions=qc / "sample_decisions.tsv", removals=qc / "removed_samples.txt"),
                        params=SimpleNamespace(random_seed=42), log=[root / "selection.log"],
                    )})
                self.assertEqual((qc / "sample_decisions.tsv").read_text(), expected_decisions)

                original = output / "Glostrup.vcf.gz"
                filtered = output / "Glostrup.deduplicated.vcf.gz"
                samples = subprocess.run(
                    [BCFTOOLS, "query", "--list-samples", str(filtered)],
                    check=True, text=True, capture_output=True,
                ).stdout.splitlines()
                self.assertEqual(samples, [row["merged_sample_id"] for row in manifest
                                           if row["merged_sample_id"] not in removals])
                def query(vcf, selection):
                    return subprocess.run(
                        [BCFTOOLS, "query", *selection, "--format", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
                        check=True, text=True, capture_output=True,
                    ).stdout
                self.assertEqual(query(filtered, []), query(original, ["--samples", ",".join(samples)]))
                stats = subprocess.run(
                    [BCFTOOLS, "index", "--stats", str(filtered)],
                    check=True, text=True, capture_output=True,
                ).stdout.splitlines()
                self.assertEqual(stats, ["chr1\t.\t2", "chrX\t.\t1"])

    def test_deduplication_rejects_invalid_priorities(self):
        for priorities, message in [
            ((None, 2), "requires a priority for every cohort"),
            ((1, 1), "requires unique cohort priorities"),
            ((0, 2), "0 is less than the minimum"),
            ((True, 2), "not of type 'integer'"),
        ]:
            with self.subTest(priorities=priorities), tempfile.TemporaryDirectory() as directory:
                config, _, _ = build_test_project(Path(directory))
                data = json.loads(config.read_text())
                data["deduplication"] = {"enabled": True}
                for cohort, priority in zip(data["cohorts"].values(), priorities):
                    if priority is not None:
                        cohort["priority"] = priority
                config.write_text(json.dumps(data))
                result = run_workflow(config)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stdout + result.stderr)

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
                output / "by_chromosome" / "chr1.vcf.gz",
                output / "by_chromosome" / "chrX.vcf.gz",
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
            merge_inputs = output / "intermediate" / "merge_input"
            self.assertEqual(list(merge_inputs.rglob("*.bcf")), [])
            self.assertEqual(list(merge_inputs.rglob("*.bcf.csi")), [])

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

    def test_conflicting_annotations_are_removed_before_merge(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config, output, source_vcfs = build_test_project(root)

            cohort_a_vcf = root / "cohort_a.symbolic.vcf"
            cohort_b_vcf = root / "cohort_b.symbolic.vcf"
            write_symbolic_vcf(
                cohort_a_vcf,
                ["duplicate", "cohort_a_sample"],
                flexible_numbers=False,
            )
            write_symbolic_vcf(
                cohort_b_vcf,
                ["duplicate", "cohort_b_sample"],
                flexible_numbers=True,
            )
            source_vcfs[0].unlink()
            source_vcfs[2].unlink()
            compress_vcf(cohort_a_vcf, source_vcfs[0])
            compress_vcf(cohort_b_vcf, source_vcfs[2])

            merged_vcf = output / "by_chromosome" / "chr1.vcf.gz"
            result = run_workflow(config, [merged_vcf])
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            header = subprocess.run(
                [BCFTOOLS, "view", "--header-only", str(merged_vcf)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout
            self.assertIn("##FORMAT=<ID=GT,", header)
            self.assertNotIn("##FORMAT=<ID=DS,", header)
            self.assertNotIn("##FORMAT=<ID=GP,", header)
            self.assertNotIn("##FORMAT=<ID=HDS,", header)
            self.assertNotIn("##INFO=<ID=AF,", header)

            record = subprocess.run(
                [BCFTOOLS, "view", "--no-header", str(merged_vcf)],
                check=True,
                text=True,
                capture_output=True,
            ).stdout.strip().split("\t")
            self.assertEqual(
                record[0:5],
                ["chr1", "39119", "chr1:39119:T:<CN2>", "T", "<CN2>"],
            )
            self.assertEqual(record[7], ".")
            self.assertEqual(record[8], "GT")
            self.assertEqual(record[9:], ["0|0", "0|1", "0|0", "0|1"])

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

    def test_samples_are_excluded_from_configured_cohort(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                exclusion_entries=[
                    "cohort_b_sample barcode-123 PCA",
                    "unknown_sample barcode-456 PCA",
                ],
            )

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            final_vcf = output / "Glostrup.vcf.gz"
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
                ],
            )

            preparation_log = output / "logs" / "sample_exclusions.prepare.log"
            self.assertIn(
                "sample 'unknown_sample' is not present in cohort 'cohort_b'; "
                "skipping",
                preparation_log.read_text(),
            )
            self.assertFalse(
                (output / "intermediate" / "sample_exclusions.txt").exists()
            )
            self.assertFalse(
                (
                    output
                    / "intermediate"
                    / "unfiltered"
                    / "Glostrup.bcf"
                ).exists()
            )

    def test_ambiguous_excluded_sample_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                exclusion_entries=["duplicate barcode-123 PCA"],
            )

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Cannot globally exclude sample 'duplicate' from cohort "
                "'cohort_b': the identifier is also present in cohort(s) "
                "cohort_a",
                result.stdout + result.stderr,
            )
            self.assertFalse((output / "Glostrup.vcf.gz").exists())

    def test_unknown_exclusions_leave_samples_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config, output, _ = build_test_project(
                Path(temporary_directory),
                exclusion_entries=["unknown_sample barcode-456 PCA"],
            )

            result = run_workflow(
                config,
                [output / "Glostrup.vcf.gz.tbi"],
            )
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            samples = subprocess.run(
                [
                    BCFTOOLS,
                    "query",
                    "--list-samples",
                    str(output / "Glostrup.vcf.gz"),
                ],
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

            relatedness_report = Path(f"{pcair_prefix}.related_pairs.tsv")
            self.assertEqual(
                relatedness_report.read_text().splitlines(),
                [
                    "sample_1\tsample_2\tking_kinship\tking_relationship\t"
                    "pcrelate_kinship\tpcrelate_relationship",
                    "duplicate\tcohort_a_sample\t0.354\tduplicate_or_mz_twin\t"
                    "0.36\tduplicate_or_mz_twin",
                    "duplicate\t2:duplicate\t0.177\tfirst_degree\t0.17\t"
                    "second_degree",
                    "duplicate\tcohort_b_sample\t0.0884\tsecond_degree\t0.04\t"
                    "distant_or_unrelated",
                    "cohort_a_sample\tcohort_b_sample\t0.02\t"
                    "distant_or_unrelated\t0.05\tthird_degree",
                    "cohort_a_sample\t2:duplicate\t0.0442\tthird_degree\t0.03\t"
                    "distant_or_unrelated",
                ],
            )

            calls = call_log.read_text().splitlines()
            self.assertFalse((output / "Glostrup.deduplicated.vcf.gz").exists())
            self.assertFalse((output / "qc" / "deduplication").exists())
            self.assertEqual(len(calls), 2)
            self.assertIn("--vcf-require-gt", calls[0])
            self.assertIn("--double-id", calls[0])
            self.assertIn("--autosome", calls[0])
            self.assertIn("--snps-only just-acgt", calls[0])
            self.assertIn("--max-alleles 2", calls[0])
            self.assertIn(f"--plink {plink_prefix}", calls[1])
            self.assertIn(f"--output {pcair_prefix}", calls[1])

    def test_empty_relatedness_report_has_header(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config, output, _ = build_test_project(root)
            config_data = json.loads(config.read_text())
            config_data["relatedness"] = {"minimum_kinship": 0.5}
            config.write_text(json.dumps(config_data))

            bin_directory = root / "bin"
            bin_directory.mkdir()
            call_log = root / "analysis_calls.log"
            env = write_fake_analysis_tools(bin_directory, call_log)

            result = run_workflow(config, env=env)
            if result.returncode != 0:
                self.fail(result.stdout + result.stderr)

            report = (
                output
                / "analysis"
                / "pcair"
                / "Glostrup.related_pairs.tsv"
            )
            self.assertEqual(
                report.read_text(),
                "sample_1\tsample_2\tking_kinship\tking_relationship\t"
                "pcrelate_kinship\tpcrelate_relationship\n",
            )


if __name__ == "__main__":
    unittest.main()

import subprocess
from pathlib import Path


cohorts = list(snakemake.params.cohorts)
exclusion_cohorts = list(snakemake.params.exclusion_cohorts)
vcfs = list(snakemake.input.vcfs)
exclusion_files = list(snakemake.input.exclusions)

sample_sets = {}
for cohort, vcf in zip(cohorts, vcfs):
    result = subprocess.run(
        ["bcftools", "query", "--list-samples", str(vcf)],
        check=True,
        text=True,
        capture_output=True,
    )
    sample_sets[cohort] = set(result.stdout.splitlines())

output = Path(snakemake.output[0])
log = Path(snakemake.log[0])
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)

selected = []
seen = set()
warnings = []

for cohort, exclusion_file in zip(exclusion_cohorts, exclusion_files):
    lines = Path(exclusion_file).read_text().splitlines()
    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        if not fields:
            continue

        sample = fields[0]
        if sample not in sample_sets[cohort]:
            warnings.append(
                f"{exclusion_file}:{line_number}: sample {sample!r} is not present "
                f"in cohort {cohort!r}; skipping"
            )
            continue

        other_cohorts = [
            other
            for other, samples in sample_sets.items()
            if other != cohort and sample in samples
        ]
        if other_cohorts:
            conflicts = ", ".join(other_cohorts)
            raise ValueError(
                f"Cannot globally exclude sample {sample!r} from cohort {cohort!r}: "
                f"the identifier is also present in cohort(s) {conflicts}"
            )

        if sample not in seen:
            selected.append(sample)
            seen.add(sample)

output.write_text("".join(f"{sample}\n" for sample in selected))
log.write_text("".join(f"WARNING: {warning}\n" for warning in warnings))

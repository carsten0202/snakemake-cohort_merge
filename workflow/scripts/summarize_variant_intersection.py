import subprocess
from pathlib import Path


def record_count(vcf):
    result = subprocess.run(
        ["bcftools", "index", "--nrecords", str(vcf)],
        check=True,
        text=True,
        capture_output=True,
    )
    return int(result.stdout.strip())


input_counts = [record_count(vcf) for vcf in snakemake.input.vcfs]
shared_counts = [record_count(vcf) for vcf in snakemake.params.intersection_vcfs]

if len(set(shared_counts)) != 1:
    details = ", ".join(
        f"{cohort}={count}"
        for cohort, count in zip(snakemake.params.cohorts, shared_counts)
    )
    raise ValueError(f"Intersection VCF record counts differ: {details}")

shared_count = shared_counts[0]
if shared_count == 0:
    raise ValueError(
        "No variants are shared by all cohorts for "
        f"{snakemake.wildcards.chromosome}"
    )

report = Path(snakemake.output.report)
report.parent.mkdir(parents=True, exist_ok=True)
with report.open("w") as output:
    output.write(
        "cohort\tinput_records\tshared_records\texcluded_records\t"
        "excluded_fraction\n"
    )
    for cohort, input_count in zip(snakemake.params.cohorts, input_counts):
        excluded_count = input_count - shared_count
        excluded_fraction = excluded_count / input_count
        output.write(
            f"{cohort}\t{input_count}\t{shared_count}\t{excluded_count}\t"
            f"{excluded_fraction:.6f}\n"
        )

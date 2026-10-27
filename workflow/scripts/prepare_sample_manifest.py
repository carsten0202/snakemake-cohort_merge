import csv
import re
import subprocess
from pathlib import Path


def sample_ids(vcf):
    result = subprocess.run(
        ["bcftools", "query", "--list-samples", str(vcf)],
        check=True, text=True, capture_output=True,
    )
    samples = result.stdout.splitlines()
    if len(samples) != len(set(samples)):
        raise ValueError(f"Duplicate sample IDs within {vcf}")
    return samples


def main(snakemake):
    # bcftools merge preserves input sample order, including renamed samples.
    originals = []
    for index, (cohort, priority, vcf) in enumerate(zip(
        snakemake.params.cohorts, snakemake.params.priorities, snakemake.input.vcfs
    ), 1):
        originals.extend(
            (sample, cohort, priority, index) for sample in sample_ids(vcf)
        )
    merged = sample_ids(snakemake.input.merged)
    if len(originals) != len(merged):
        raise ValueError("Merged sample count does not match the cohort sample lists")

    manifest = {}
    for merged_id, (original, cohort, priority, index) in zip(merged, originals):
        # Prefixes may repeat when an input already contains a prefixed ID.
        if not re.fullmatch(rf"(?:{index}:)*{re.escape(original)}", merged_id):
            raise ValueError(
                f"Cannot map merged ID {merged_id!r} to {cohort}: {original!r}"
            )
        manifest[merged_id] = (merged_id, original, cohort, priority)

    retained = sample_ids(snakemake.input.final)
    retained_set = set(retained)
    unknown = retained_set - manifest.keys()
    if unknown:
        raise ValueError(f"Final VCF contains unmapped samples: {sorted(unknown)}")
    if [sample for sample in merged if sample in retained_set] != retained:
        raise ValueError("Final VCF sample order differs from the merged header")

    output = Path(snakemake.output.manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["merged_sample_id", "original_sample_id", "cohort", "priority"])
        writer.writerows(manifest[sample] for sample in retained)
    Path(snakemake.log[0]).write_text(
        f"Mapped {len(retained)} retained samples from {len(originals)} input samples\n"
    )


if __name__ == "__main__":
    main(snakemake)

import subprocess
from pathlib import Path


result = subprocess.run(
    ["bcftools", "query", "--list-samples", str(snakemake.input.analysis)],
    check=True,
    text=True,
    capture_output=True,
)
analysis_samples = result.stdout.splitlines()
if len(analysis_samples) != len(set(analysis_samples)):
    raise ValueError("Analysis BCF contains duplicate sample IDs")

removals = set()
for filename in snakemake.input.removals:
    removals.update(Path(filename).read_text().splitlines())
unknown = removals - set(analysis_samples)
if unknown:
    raise ValueError(f"Removal list contains unknown analysis samples: {sorted(unknown)}")
retained = [sample for sample in analysis_samples if sample not in removals]
if not retained:
    raise ValueError("No samples remain for deployment")

Path(snakemake.output.samples).write_text("".join(f"{sample}\n" for sample in retained))
# Conversion uses --double-id, so both FID and IID equal the merged sample ID.
Path(snakemake.output.plink_keep).write_text(
    "".join(f"{sample}\t{sample}\n" for sample in retained)
)
Path(snakemake.log[0]).write_text(
    f"Retained {len(retained)} of {len(analysis_samples)} analysis samples\n"
)

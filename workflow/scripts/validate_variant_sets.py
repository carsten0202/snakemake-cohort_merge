import hashlib
import subprocess
from pathlib import Path


def variant_digest(vcf):
    command = [
        "bcftools",
        "query",
        "--format",
        "%CHROM\t%POS\t%REF\t%ALT\n",
        str(vcf),
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    digest = hashlib.sha256()
    record_count = 0

    assert process.stdout is not None
    for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
        digest.update(chunk)
        record_count += chunk.count(b"\n")

    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)

    return record_count, digest.hexdigest()


results = [
    (cohort, str(vcf), *variant_digest(vcf))
    for cohort, vcf in zip(snakemake.params.cohorts, snakemake.input.vcfs)
]

expected = results[0][2:]
mismatches = [result for result in results[1:] if result[2:] != expected]
if mismatches:
    details = "\n".join(
        f"{cohort}: {records} records, sha256={digest}"
        for cohort, _, records, digest in results
    )
    raise ValueError(f"Variant records differ between cohorts:\n{details}")

report = Path(snakemake.output.report)
report.parent.mkdir(parents=True, exist_ok=True)
with report.open("w") as output:
    output.write("cohort\tvcf\trecords\tsha256\n")
    for cohort, vcf, records, digest in results:
        output.write(f"{cohort}\t{vcf}\t{records}\t{digest}\n")

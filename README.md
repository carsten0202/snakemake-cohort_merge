# Cohort Merge

Snakemake workflow for merging imputed cohort VCF files into a single
mega-cohort VCF and performing kinship-aware population structure analysis.

The workflow intersects and merges cohorts independently by chromosome, then
concatenates the chromosome-level results. Only exact variants present in every
cohort are retained. Do not commit real cohort data or other sensitive inputs
to this repository.

After merging, the workflow converts autosomal biallelic SNP hard calls to
PLINK BED format and runs quick-pcair. The analysis calculates KING-robust
kinship, performs LD pruning, runs PC-AiR, and runs PC-Relate.

## Local development

Create and activate the development environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
```

Run the synthetic integration tests:

```bash
python -m unittest discover -s tests/integration -v
```

The tests require `bcftools` to be available on `PATH`. PLINK and quick-pcair
are represented by local test doubles, so their production modules are not
required. The tests use temporary, synthetic VCFs, including cohorts with
differing variants, and do not access real cohort data.

The example configuration contains placeholder paths and therefore cannot
construct a complete workflow DAG. It can still be used for configuration and
syntax checks:

```bash
snakemake \
    --configfile config/super_cohorts/Glostrup.example.yaml \
    --list-rules
```

## Configuration

`config/defaults.yaml` defines shared defaults, including:

- Canonical chromosome order.
- The default VCF filename pattern.
- Environment modules.
- Initial thread and memory requirements.

A super-cohort configuration supplies its name, shared output root, and cohort
input directories. The workflow writes all artifacts beneath a directory named
for the super-cohort. Values in the super-cohort configuration recursively
override the defaults.

Real super-cohort configurations under `config/super_cohorts/` are ignored by
Git. Tracked `*.example.yaml` files document their expected structure without
exposing server paths.

Each cohort uses `input.pattern` unless it defines its own `pattern`. Patterns
must contain `{chromosome}` and end in `.vcf.gz`.

A cohort can optionally provide a file of samples to remove from the final
super-cohort:

```yaml
cohorts:
  cohort_a:
    directory: /path/to/cohort_a

  cohort_b:
    directory: /path/to/cohort_b
    exclude_samples: /path/to/cohort_b/excluded_samples.txt
```

The exclusion file is whitespace-delimited and has no header. The first column
must contain the VCF sample identifier; additional columns, such as a barcode
or exclusion reason, are ignored. Blank lines are also ignored.

A cohort whose X chromosome file is named with `chr23`, while the VCF itself
still uses `chrX`, can define:

```yaml
chromosome_aliases:
  chrX: chr23
```

## Input requirements

Each cohort must provide one BGZF-compressed VCF per configured chromosome.

Input VCFs must:

- Use compatible genome builds and contig names.
- Contain non-overlapping sample sets, except for known or suspected duplicates.
- Use `chrX` internally even when the filename contains `chr23`.
- Be sorted and suitable for TBI indexing.
- Contain hard genotype calls in `FORMAT/GT` for population structure analysis.

The workflow validates exclusion identifiers against the first configured
chromosome of every cohort. An identifier absent from its configured cohort is
reported in `logs/sample_exclusions.prepare.log` and skipped. An identifier
present in both its configured cohort and another cohort is rejected because a
global exclusion would be ambiguous.

Each input VCF should have a matching `.vcf.gz.tbi`. Missing or stale indexes
are generated beside the source VCF, so cohort directories must be writable.

## Variant intersection

Before merging a chromosome, the workflow uses `bcftools isec` to select
variants present in every cohort. Records must match exactly on `CHROM`, `POS`,
`REF`, and `ALT`; matching positions with different alleles are not retained.

The filtered per-cohort VCFs are temporary and are removed after a successful
merge. The workflow fails if a chromosome has no variants shared by every
cohort.

Before merging, the workflow creates temporary BCF copies with all INFO fields
and all FORMAT fields except `GT` removed. This prevents incompatible
cohort-specific annotation definitions from affecting the merge. The original
intersection files are used for the retained QC counts below.

Retained QC reports contain the input, shared, and excluded record counts for
each cohort:

```text
<output_directory>/<super_cohort>/qc/variant_intersection/<chromosome>.tsv
```

## Sample duplicates

`bcftools merge --force-samples` retains exact duplicate sample IDs. Later
occurrences receive bcftools prefixes based on input order, for example:

```text
duplicate
2:duplicate
```

Biological duplicates that use different sample IDs remain unchanged. Both
types are retained for the downstream KING-robust kinship analysis.

Cohort order in the super-cohort YAML determines bcftools input order and
duplicate-sample prefixes.

When sample exclusions are configured, the workflow combines validated IDs
from all cohort exclusion files and removes them once after chromosome
concatenation. General duplicate-ID support remains enabled, but an ID selected
for exclusion must be unique to its configured cohort.

## Retained VCF fields

The chromosome and final VCFs retain variant coordinates, IDs, alleles, quality,
filters, and hard genotype calls in `FORMAT/GT`. Cohort-specific INFO fields and
other FORMAT fields such as `DS`, `GP`, and `HDS` are intentionally discarded:
they can have incompatible definitions and do not describe the combined cohort.
Symbolic and non-SNP records remain in the merged VCF, but INFO annotations that
describe structural intervals, including `END`, are not retained. The merged
VCF is therefore intended for the documented hard-call population analysis,
not as a structural-variant output. The downstream PLINK conversion retains
only biallelic A/C/G/T SNPs.

## Population structure analysis

PLINK 2 converts the combined VCF to BED/BIM/FAM files before quick-pcair is
run. Conversion retains autosomal, biallelic A/C/G/T SNPs and requires the VCF
`GT` field. Dosage fields such as `DS` and `GP` are not imported, and PLINK does
not derive new hard calls from them. This is intentional because PLINK BED and
quick-pcair's GDS input represent hard genotype calls rather than dosages.

VCF files do not carry PLINK pedigree information. Conversion therefore uses
`--double-id`, assigning every sample its own family ID while retaining the VCF
sample name as its individual ID. This prevents quick-pcair from passing one
shared family ID for the entire cohort to SNPRelate's KING implementation.

quick-pcair performs its own LD pruning with SNPRelate. Version 1.1.0 uses a
10 Mb sliding window and a correlation threshold of `sqrt(0.1)`, corresponding
to an r-squared threshold of 0.1. A separate PLINK LD-pruning step is therefore
not performed.

The quick-pcair analysis includes:

- KING-robust kinship estimation.
- PC-AiR using the LD-pruned SNP set.
- PC-Relate using the first two PC-AiR components.

The default quick-pcair memory request is 256 GB because its kinship matrix can
grow quadratically with sample count. Override the rule resources in the
super-cohort configuration when the cohort requires a larger allocation.

## Outputs

For a super-cohort named `Glostrup`, the workflow produces:

```text
<output_directory>/
`-- Glostrup/
    |-- Glostrup.vcf.gz
    |-- Glostrup.vcf.gz.tbi
    |-- analysis/
    |   |-- pcair/
    |   |   |-- Glostrup.gds
    |   |   |-- Glostrup.KINGkinship.png
    |   |   |-- Glostrup.KINGkinship.tsv
    |   |   |-- Glostrup.ldprune.snpids.txt
    |   |   |-- Glostrup.eigenvalues.tsv
    |   |   |-- Glostrup.eigenvectors.tsv
    |   |   |-- Glostrup.pcair_1v2.png
    |   |   |-- Glostrup.pcair_3v4.png
    |   |   |-- Glostrup.pcrelate.RData
    |   |   `-- Glostrup.pcrelate_1v2.png
    |   `-- plink/
    |       |-- Glostrup.bed
    |       |-- Glostrup.bim
    |       |-- Glostrup.fam
    |       `-- Glostrup.log
    |-- by_chromosome/
    |   |-- chr1.vcf.gz
    |   |-- chr1.vcf.gz.tbi
    |   `-- ...
    |-- logs/
    |   |-- concatenate.log
    |   |-- plink_conversion.log
    |   |-- quick_pcair.log
    |   |-- sample_exclusions.log       (when exclusions are configured)
    |   |-- sample_exclusions.prepare.log
    |   |-- intersect/
    |   |   |-- chr1.log
    |   |   `-- ...
    |   |-- prepare_merge/
    |   |   |-- chr1/
    |   |   `-- ...
    |   `-- merge/
    |       |-- chr1.log
    |       `-- ...
    `-- qc/
        `-- variant_intersection/
            |-- chr1.tsv
            `-- ...
```

Chromosome-level outputs are retained after final concatenation.

## Server execution

Production rules load these environment modules in order:

```text
perl
gsl/2.5
bcftools/1.21
plink/2.0-alpha-6.2
gcc/13.2.0
openjdk/20.0.0
R/4.4.2
quick-pcair/1.1.0
```

Each rule loads only the modules it needs. The quick-pcair rule loads
`gcc/13.2.0`, `openjdk/20.0.0`, `R/4.4.2`, and `quick-pcair/1.1.0` in
prerequisite order.

Run Glostrup with:

```bash
snakemake \
    --configfile config/super_cohorts/Glostrup.yaml \
    --use-envmodules \
    --cores all
```

Snakemake purges inherited modules before loading the modules declared by each
rule. A scheduler-specific execution profile can be added once the server
scheduler is known.

## Layout

```text
config/defaults.yaml                        Shared workflow defaults
config/super_cohorts/                       Super-cohort configurations
tests/integration/test_workflow.py          Synthetic integration tests
workflow/Snakefile                          Workflow entry point
workflow/rules/bcftools.smk                 BCFtools workflow rules
workflow/rules/pcair.smk                    PLINK and quick-pcair rules
workflow/schemas/config.schema.yaml         Configuration schema
workflow/scripts/prepare_sample_exclusions.py
                                            Sample exclusion validator
workflow/scripts/summarize_variant_intersection.py
                                            Intersection QC reporter
requirements-dev.txt                        Local development dependencies
```

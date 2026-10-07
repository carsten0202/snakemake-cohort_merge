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

A super-cohort configuration supplies its name, output directory, and cohort
input directories. Values in the super-cohort configuration recursively
override the defaults.

Real super-cohort configurations under `config/super_cohorts/` are ignored by
Git. Tracked `*.example.yaml` files document their expected structure without
exposing server paths.

Each cohort uses `input.pattern` unless it defines its own `pattern`. Patterns
must contain `{chromosome}` and end in `.vcf.gz`.

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

Each input VCF should have a matching `.vcf.gz.tbi`. Missing or stale indexes
are generated beside the source VCF, so cohort directories must be writable.

## Variant intersection

Before merging a chromosome, the workflow uses `bcftools isec` to select
variants present in every cohort. Records must match exactly on `CHROM`, `POS`,
`REF`, and `ALT`; matching positions with different alleles are not retained.

The filtered per-cohort VCFs are temporary and are removed after a successful
merge. The workflow fails if a chromosome has no variants shared by every
cohort.

Retained QC reports contain the input, shared, and excluded record counts for
each cohort:

```text
<output_directory>/qc/variant_intersection/<chromosome>.tsv
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

## INFO fields

The workflow currently preserves bcftools' default INFO merge behavior.

Standard bcftools defaults aggregate selected fields such as `DP` and `DP4`;
other cohort-level annotations may be inherited from the first input rather
than recalculated for the mega-cohort. In particular, imputation-quality and
allele-frequency annotations should not be assumed to represent the combined
cohort until an explicit recalculation step is added.

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

The default quick-pcair memory request is 64 GB because its kinship matrix can
grow quadratically with sample count. Override the rule resources in the
super-cohort configuration when the cohort requires a larger allocation.

## Outputs

For a super-cohort named `Glostrup`, the workflow produces:

```text
<output_directory>/
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
|   |-- Glostrup.chr1.vcf.gz
|   |-- Glostrup.chr1.vcf.gz.tbi
|   `-- ...
|-- logs/
|   |-- concatenate.log
|   |-- plink_conversion.log
|   |-- quick_pcair.log
|   |-- intersect/
|   |   |-- chr1.log
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
R/4.4.2
quick-pcair/1.1.0
```

Each rule loads only the modules it needs. The quick-pcair rule loads
`R/4.4.2` before `quick-pcair/1.1.0`, as required by the module definition.

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
workflow/scripts/summarize_variant_intersection.py
                                            Intersection QC reporter
requirements-dev.txt                        Local development dependencies
```

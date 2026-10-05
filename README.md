# Cohort Merge

Snakemake workflow for merging imputed cohort VCF files into a single
mega-cohort VCF.

The workflow merges cohorts independently by chromosome, validates that their
variant records match, and concatenates the chromosome-level results. Do not
commit real cohort data or other sensitive inputs to this repository.

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

The tests require `bcftools` to be available on `PATH`. They use temporary,
synthetic VCFs and do not access real cohort data.

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

- Use the same ordered `CHROM`, `POS`, `REF`, and `ALT` records.
- Contain non-overlapping sample sets, except for known or suspected duplicates.
- Use `chrX` internally even when the filename contains `chr23`.
- Be sorted and suitable for TBI indexing.

Each input VCF should have a matching `.vcf.gz.tbi`. Missing or stale indexes
are generated beside the source VCF, so cohort directories must be writable.

## Variant validation

Before merging a chromosome, the workflow streams `CHROM`, `POS`, `REF`, and
`ALT` from every cohort VCF and calculates a SHA-256 digest.

The merge is rejected if cohorts differ in record count, variant identity, or
record order. Validation reports are written to:

```text
<output_directory>/qc/variant_sets/<chromosome>.sha256.tsv
```

## Sample duplicates

`bcftools merge --force-samples` retains exact duplicate sample IDs. Later
occurrences receive bcftools prefixes based on input order, for example:

```text
duplicate
2:duplicate
```

Biological duplicates that use different sample IDs remain unchanged. Both
types can therefore be handled together by a later kinship workflow.

Cohort order in the super-cohort YAML determines bcftools input order and
duplicate-sample prefixes.

## INFO fields

The workflow currently preserves bcftools' default INFO merge behavior.

Standard bcftools defaults aggregate selected fields such as `DP` and `DP4`;
other cohort-level annotations may be inherited from the first input rather
than recalculated for the mega-cohort. In particular, imputation-quality and
allele-frequency annotations should not be assumed to represent the combined
cohort until an explicit recalculation step is added.

## Outputs

For a super-cohort named `Glostrup`, the workflow produces:

```text
<output_directory>/
|-- Glostrup.vcf.gz
|-- Glostrup.vcf.gz.tbi
|-- by_chromosome/
|   |-- Glostrup.chr1.vcf.gz
|   |-- Glostrup.chr1.vcf.gz.tbi
|   `-- ...
|-- logs/
|   |-- concatenate.log
|   `-- merge/
|       |-- chr1.log
|       `-- ...
`-- qc/
    `-- variant_sets/
        |-- chr1.sha256.tsv
        `-- ...
```

Chromosome-level outputs are retained after final concatenation.

## Server execution

Production rules load these environment modules in order:

```text
perl
gsl/2.5
bcftools/1.21
```

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
workflow/schemas/config.schema.yaml         Configuration schema
workflow/scripts/validate_variant_sets.py   Variant-set validator
requirements-dev.txt                        Local development dependencies
```

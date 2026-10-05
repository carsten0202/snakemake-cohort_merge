# Cohort Merge

Snakemake workflow for merging individual cohort VCF files into a single
mega-cohort VCF.

The workflow currently provides its configuration scaffold but no data
processing rules. Do not commit real cohort data or other sensitive inputs to
this repository.

## Local development

Create and activate the development environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
```

Validate the workflow from the repository root:

```bash
snakemake --configfile config/super_cohorts/Glostrup.example.yaml --dry-run
```

## Configuration

`config/defaults.yaml` defines shared defaults, including the canonical
chromosomes and VCF filename pattern. A super-cohort config supplies its name,
output directory, and cohort input directories. Values in the super-cohort
config recursively override the defaults.

Real super-cohort configs under `config/super_cohorts/` are ignored by Git.
Tracked `*.example.yaml` files document their expected structure without
exposing server paths.

Run Glostrup with its populated config:

```bash
snakemake --configfile config/super_cohorts/Glostrup.yaml
```

Each cohort uses `input.pattern` unless it defines its own `pattern`. Patterns
must contain `{chromosome}` and end in `.vcf.gz`. A cohort whose X chromosome
file is named with `chr23` can define:

```yaml
chromosome_aliases:
  chrX: chr23
```

## Server execution

Production tools will be supplied through environment modules rather than
Conda environments. Rules requiring modules will use Snakemake's `envmodules`
directive and the workflow will be executed with `--use-envmodules`.

## Layout

```text
config/defaults.yaml                  Shared workflow defaults
config/super_cohorts/                 Super-cohort configurations
workflow/Snakefile                    Workflow entry point
workflow/rules/                       Modular Snakemake rule files
workflow/schemas/config.schema.yaml   Configuration schema
workflow/scripts/                     Supporting scripts
requirements-dev.txt                  Local development dependencies
```

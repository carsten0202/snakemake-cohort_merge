# Cohort Merge

Snakemake workflow for merging individual cohort VCF files into a single
mega-cohort VCF.

The workflow is currently an empty scaffold. Do not commit real cohort data or
other sensitive inputs to this repository.

## Local development

Create and activate the development environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
```

Validate the workflow from the repository root:

```bash
snakemake --snakefile workflow/Snakefile --dry-run
```

## Server execution

Production tools will be supplied through environment modules rather than
Conda environments. Rules requiring modules will use Snakemake's `envmodules`
directive and the workflow will be executed with `--use-envmodules`.

## Layout

```text
config/config.yaml       Workflow configuration
workflow/Snakefile       Workflow entry point
workflow/rules/          Modular Snakemake rule files
workflow/scripts/        Supporting scripts
requirements-dev.txt     Local development dependencies
```

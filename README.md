# Cohort Merge

Snakemake workflow for merging imputed cohort VCF files into a mega-cohort
and performing kinship-aware population structure analysis.

The workflow intersects and merges cohorts independently by chromosome. Only
exact variants present in every cohort are retained. Final chromosome VCFs,
PLINK files, and analysis results are bundled under `deliverables/` for copying
to the deployed site. Do not commit real cohort data or other sensitive inputs
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

The tests require `bcftools` and `Rscript` to be available on `PATH`. PLINK and
quick-pcair are represented by local test doubles, so their production modules
are not required. The tests use temporary, synthetic VCFs, including cohorts
with differing variants, and do not access real cohort data.

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
- Initial thread and memory requests.

A super-cohort configuration supplies its name, shared output root, and cohort
input directories. The workflow writes all artifacts beneath a directory named
for the super-cohort. Values in the super-cohort configuration recursively
override the defaults.

The resource defaults are scheduling starting points rather than hardware
requirements. Override them for a deployment based on its executor, available
hardware, and observed resource usage. A rule's `threads` value is the total
number of CPUs available to the job. For BCFtools commands, the workflow
derives the `--threads` argument from this allocation because BCFtools counts
compression workers separately from its main thread. `bcftools isec` creates a
compression worker pool for each selected output, so its workers are divided
across the cohort outputs.

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
Compressed VCF and BCF files produced by the workflow are indexed by their
producing BCFtools command, avoiding a separate indexing pass.

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
concatenation for analysis. Those exclusions also apply to every delivered
chromosome VCF and the final PLINK dataset. General duplicate-ID support remains
enabled, but an ID selected for exclusion must be unique to its configured cohort.

### Cohort-prioritized deduplication

Optionally retain one representative from each duplicate/monozygotic-twin group
after PC-Relate. Enable the stage and assign a unique positive integer priority
to every cohort; **lower numbers take precedence**:

```yaml
deduplication:
  enabled: true
  minimum_kinship: 0.354
  random_seed: 42

cohorts:
  cohort_a:
    directory: /path/to/cohort_a
    priority: 1
  cohort_b:
    directory: /path/to/cohort_b
    priority: 2
```

Deduplication defaults to disabled. When enabled, missing or repeated cohort
priorities are rejected. Priority does not alter YAML merge order or sample-ID
prefixes. The example super-cohort configuration enables this stage.

The workflow extracts sample lists from each cohort's first configured
chromosome and maps them, in merge order, to the actual chromosome-level merged
header. It validates original IDs and bcftools-generated prefixes, then retains
only samples still present after the existing configured exclusions. Unmapped
or inconsistent sample headers cause an error.

Qualifying pairs are exported directly from `mypcrelate$kinBtwn` in the
PC-Relate RData file, independently of `relatedness.minimum_kinship`. Only
finite PC-Relate estimates at or above `deduplication.minimum_kinship` create
connections; KING estimates do not determine removal. Missing and non-finite
PC-Relate estimates do not create connections.

Each connected component is one duplicate/twin group. For example, A--B and
B--C form one group even if A--C is below the threshold or absent. The report
flags groups that are not fully pairwise connected. This intentionally retains
exactly one representative per connected group, rather than one per qualifying
pair. Samples outside duplicate groups are retained.

Within each group, the workflow first selects the highest-priority cohort.
If multiple members belong to that cohort, it chooses one using the configured
random seed. Groups and candidates are sorted before selection, making results
reproducible for the same inputs and seed regardless of PC-Relate row order.
Changing the seed or group membership may change the representative.

Enabling the stage applies duplicate removals to the delivered chromosome VCFs
and PLINK files, and adds these default workflow targets:

```text
<output_directory>/<super_cohort>/
    qc/deduplication/
        sample_manifest.tsv
        duplicate_pairs.tsv
        sample_decisions.tsv
        removed_samples.txt
```

- `sample_manifest.tsv`: merged sample ID, original sample ID, cohort, priority;
  includes all samples present before deduplication.
- `duplicate_pairs.tsv`: qualifying sample pairs and their PC-Relate kinship.
- `sample_decisions.tsv`: one row per duplicate-group member, with group ID,
  sample IDs, cohort, priority, keep/remove decision, retained representative,
  selection reason, and `fully_pairwise_connected` flag.
- `removed_samples.txt`: headerless list of merged sample IDs removed from the VCF.

When there are no qualifying pairs, the pair and decision reports contain only
headers, the removal list is empty, and all samples are retained. Filtering
preserves variants and the retained samples' genotype calls. PCs and kinship
describe the analysis sample set after configured exclusions and before
deduplication; they are not recalculated automatically. The pre-deduplication
genotype inputs are temporary.

Stage logs are `logs/sample_manifest.log`, `logs/duplicate_pairs.log`,
`logs/duplicate_selection.log`, and `logs/final_samples.log`. Resources can be
overridden through `resources.prepare_sample_manifest`,
`resources.export_duplicate_pairs`, `resources.select_duplicate_representatives`,
and `resources.prepare_final_samples`. Final genotype export uses
`resources.export_chromosome_vcf` and `resources.export_final_plink`.

## Retained VCF fields

The delivered chromosome VCFs retain variant coordinates, IDs, alleles, quality,
filters, and hard genotype calls in `FORMAT/GT`. Cohort-specific INFO fields and
other FORMAT fields such as `DS`, `GP`, and `HDS` are intentionally discarded:
they can have incompatible definitions and do not describe the combined cohort.
Symbolic and non-SNP records remain in the merged VCF, but INFO annotations that
describe structural intervals, including `END`, are not retained. The merged
VCF is therefore intended for the documented hard-call population analysis,
not as a structural-variant output. The downstream PLINK conversion retains
only biallelic A/C/G/T SNPs.

## Population structure analysis

PLINK 2 converts a temporary combined BCF to analysis BED/BIM/FAM files before
quick-pcair is run. Conversion retains autosomal, biallelic A/C/G/T SNPs and requires the
`GT` field. Dosage fields such as `DS` and `GP` are not imported, and PLINK does
not derive new hard calls from them. This is intentional because PLINK BED and
quick-pcair's GDS input represent hard genotype calls rather than dosages.

VCF files do not carry PLINK pedigree information. Conversion therefore uses
`--double-id`, assigning every sample its own family ID while retaining the VCF
sample name as its individual ID. This prevents quick-pcair from passing one
shared family ID for the entire cohort to SNPRelate's KING implementation.

After sample selection, PLINK 2 subsets the analysis dataset using explicit
FID/IID pairs and writes the final dataset to `deliverables/plink/`. Its sample
IDs and order match every delivered chromosome VCF. Its variants remain the
autosomal, biallelic SNP subset, while VCFs include all configured chromosomes
and retained variant types. The analysis PLINK dataset and combined BCF are
temporary and are removed after their consumers finish successfully.

quick-pcair performs its own LD pruning with SNPRelate. Version 1.1.0 uses a
10 Mb sliding window and a correlation threshold of `sqrt(0.1)`, corresponding
to an r-squared threshold of 0.1. A separate PLINK LD-pruning step is therefore
not performed.

The quick-pcair analysis includes:

- KING-robust kinship estimation.
- PC-AiR using the LD-pruned SNP set.
- PC-Relate using the first two PC-AiR components.

The default quick-pcair memory request is 524,288 MB, approximately 512 GiB.
KING and PC-Relate resource requirements can grow quadratically with sample
count, so this default does not guarantee that every cohort will fit. Override
the rule resources in the super-cohort configuration when observed usage or
deployment limits require it.
The default quick-pcair runtime request is 11,520 minutes, or eight days.

## Related sample report

The workflow combines the KING-robust and ancestry-adjusted PC-Relate estimates
into:

```text
<output_directory>/<super_cohort>/deliverables/pcair/<super_cohort>.related_pairs.tsv
```

The report contains one row per unordered sample pair and the columns
`sample_1`, `sample_2`, `king_kinship`, `king_relationship`,
`pcrelate_kinship`, and `pcrelate_relationship`. A pair is included when either
kinship estimate is at least `relatedness.minimum_kinship`. The default is
`0.0442`, corresponding to third-degree relatives or closer. Override it in a
super-cohort configuration with:

```yaml
relatedness:
  minimum_kinship: 0.0442
```

Relationship classes are assigned independently to each estimate using these
kinship boundaries:

- `>= 0.354`: duplicate or monozygotic twin.
- `>= 0.177`: first degree.
- `>= 0.0884`: second degree.
- `>= 0.0442`: third degree.
- `< 0.0442`: distant or unrelated.

The classes are inferred degree ranges, not definitive pedigree assignments.
For example, kinship alone does not distinguish a parent-child pair from full
siblings. KING and PC-Relate can assign different classes because PC-Relate
adjusts its estimates for ancestry.

## Outputs

For a super-cohort named `Glostrup`, the workflow produces:

```text
<output_directory>/
`-- Glostrup/
    |-- deliverables/
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
    |   |   |-- Glostrup.pcrelate_1v2.png
    |   |   `-- Glostrup.related_pairs.tsv
    |   |-- plink/
    |   |   |-- Glostrup.bed
    |   |   |-- Glostrup.bim
    |   |   |-- Glostrup.fam
    |   |   `-- Glostrup.log
    |   `-- vcf/
    |       |-- chr1.vcf.gz
    |       |-- chr1.vcf.gz.tbi
    |       |-- chrX.vcf.gz
    |       |-- chrX.vcf.gz.tbi
    |       `-- ...
    |-- logs/
    |   |-- concatenate.log
    |   |-- plink_conversion.log
    |   |-- plink_export.log
    |   |-- final_samples.log
    |   |-- quick_pcair.log
    |   |-- sample_exclusions.log       (when exclusions are configured)
    |   |-- sample_exclusions.prepare.log
    |   |-- export_vcf/
    |   |   |-- chr1.log
    |   |   `-- ...
    |   |-- intersect/
    |   |   |-- chr1.log
    |   |   `-- ...
    |   |-- prepare_merge/
    |   |   |-- chr1/
    |   |   `-- ...
    |   `-- merge/
    |       |-- chr1.log
    |       `-- ...
    |-- qc/
    |   |-- final_samples.txt
    |   |-- deduplication/               (when enabled)
    |   |   |-- sample_manifest.tsv
    |   |   |-- duplicate_pairs.tsv
    |   |   |-- sample_decisions.tsv
    |   |   `-- removed_samples.txt
    |   `-- variant_intersection/
    |       |-- chr1.tsv
    |       `-- ...
    `-- intermediate/                   (temporary processing artifacts)
```

Copy `deliverables/` to deploy the genotype datasets and analysis outputs.
Chromosome filenames use configured canonical names, including `chrX.vcf.gz`.

`qc/final_samples.txt` is a retained, headerless list of merged sample IDs in
original merge order, after configured exclusions and optional deduplication.
This one sample list controls both chromosome VCF export and PLINK subsetting.
If no samples remain or a duplicate-removal ID is absent from the analysis
dataset, the workflow fails rather than silently exporting inconsistent samples.

| Output | Sample set | Variant set |
| --- | --- | --- |
| `deliverables/vcf/chr*.vcf.gz` | Final retained samples | Exact-intersection variants for that chromosome |
| `deliverables/plink/` | Same final samples, in the same order | Autosomal biallelic A/C/G/T SNPs |
| `deliverables/pcair/` | After configured exclusions, before deduplication | Existing analysis inputs and estimates |

Combined VCFs are not retained. The workflow uses temporary merged chromosome
VCFs, combined BCFs, and analysis PLINK files under `intermediate/`, then removes
them after successful use. The default workflow targets include all delivered
chromosome VCFs and indexes, final PLINK files, analysis outputs, the final
sample list, and applicable QC reports. Files from older runs in the previous
output layout are not automatically moved or deleted.

## Execution

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

Scheduler configuration is intentionally kept outside this repository because
executor plugins, accounts, partitions, and submission policies are
deployment-specific. When using a scheduler, supply the appropriate external
Snakemake profile:

```bash
snakemake \
    --profile /path/to/profile \
    --configfile config/super_cohorts/Glostrup.yaml \
    --use-envmodules
```

The profile should map Snakemake's `threads`, `mem_mb`, and `runtime` resources
to the corresponding scheduler requests.

For direct execution without a scheduler:

```bash
snakemake \
    --configfile config/super_cohorts/Glostrup.yaml \
    --use-envmodules \
    --cores all
```

Snakemake purges inherited modules before loading the modules declared by each
rule.

## Layout

```text
config/defaults.yaml                        Shared workflow defaults
config/super_cohorts/                       Super-cohort configurations
tests/integration/test_workflow.py          Synthetic integration tests
workflow/Snakefile                          Workflow entry point
workflow/rules/bcftools.smk                 BCFtools workflow rules
workflow/rules/pcair.smk                    PLINK and quick-pcair rules
workflow/rules/deduplication.smk            Cohort-prioritized deduplication rules
workflow/schemas/config.schema.yaml         Configuration schema
workflow/scripts/prepare_sample_exclusions.py
                                             Sample exclusion validator
workflow/scripts/prepare_sample_manifest.py
                                             Merged sample provenance mapper
workflow/scripts/prepare_final_samples.py    Final sample list and PLINK keep pairs
workflow/scripts/export_duplicate_pairs.R   PC-Relate duplicate-pair exporter
workflow/scripts/select_duplicate_representatives.py
                                             Seeded cohort-prioritized selection
workflow/scripts/summarize_variant_intersection.py
                                            Intersection QC reporter
requirements-dev.txt                        Local development dependencies
```

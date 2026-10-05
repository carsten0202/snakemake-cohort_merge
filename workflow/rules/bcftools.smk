BCFTOOLS_MODULES = tuple(config["modules"]["bcftools"])

MERGED_VCF_PATTERN = str(
    Path(OUTPUT_DIRECTORY)
    / "by_chromosome"
    / f"{SUPER_COHORT}.{{chromosome}}.vcf.gz"
)

VARIANT_REPORT_PATTERN = str(
    Path(OUTPUT_DIRECTORY) / "qc" / "variant_sets" / "{chromosome}.sha256.tsv"
)


def cohort_vcf(cohort, chromosome):
    cohort_config = config["cohorts"][cohort]
    input_chromosome = cohort_config.get("chromosome_aliases", {}).get(
        chromosome, chromosome
    )
    pattern = cohort_config.get("pattern", config["input"]["pattern"])
    filename = pattern.format(chromosome=input_chromosome)
    return str(Path(cohort_config["directory"]) / filename)


def chromosome_vcfs(wildcards):
    return [cohort_vcf(cohort, wildcards.chromosome) for cohort in COHORTS]


def chromosome_indexes(wildcards):
    return [f"{vcf}.tbi" for vcf in chromosome_vcfs(wildcards)]


rule index_vcf:
    input:
        "{vcf_path}.vcf.gz"
    output:
        "{vcf_path}.vcf.gz.tbi"
    threads:
        config["resources"]["index_vcf"]["threads"]
    resources:
        mem_mb=config["resources"]["index_vcf"]["mem_mb"]
    envmodules:
        *BCFTOOLS_MODULES
    shell:
        """
        bcftools index \
            --tbi \
            --threads {threads} \
            --output {output:q} \
            {input:q}
        """


rule validate_variant_sets:
    input:
        vcfs=chromosome_vcfs
    output:
        report=VARIANT_REPORT_PATTERN
    params:
        cohorts=COHORTS
    threads:
        config["resources"]["validate_variant_sets"]["threads"]
    resources:
        mem_mb=config["resources"]["validate_variant_sets"]["mem_mb"]
    envmodules:
        *BCFTOOLS_MODULES
    script:
        "../scripts/validate_variant_sets.py"


rule merge_chromosome:
    input:
        vcfs=chromosome_vcfs,
        indexes=chromosome_indexes,
        variant_report=VARIANT_REPORT_PATTERN,
    output:
        vcf=MERGED_VCF_PATTERN
    threads:
        config["resources"]["merge_chromosome"]["threads"]
    resources:
        mem_mb=config["resources"]["merge_chromosome"]["mem_mb"]
    log:
        str(Path(OUTPUT_DIRECTORY) / "logs" / "merge" / "{chromosome}.log")
    envmodules:
        *BCFTOOLS_MODULES
    shell:
        """
        bcftools merge \
            --force-samples \
            --output-type z \
            --threads {threads} \
            --output {output.vcf:q} \
            {input.vcfs:q} \
            2> {log:q}
        """


rule concatenate_chromosomes:
    input:
        vcfs=expand(MERGED_VCF_PATTERN, chromosome=CHROMOSOMES),
        indexes=expand(f"{MERGED_VCF_PATTERN}.tbi", chromosome=CHROMOSOMES),
    output:
        vcf=FINAL_VCF
    threads:
        config["resources"]["concatenate_chromosomes"]["threads"]
    resources:
        mem_mb=config["resources"]["concatenate_chromosomes"]["mem_mb"]
    log:
        str(Path(OUTPUT_DIRECTORY) / "logs" / "concatenate.log")
    envmodules:
        *BCFTOOLS_MODULES
    shell:
        """
        bcftools concat \
            --output-type z \
            --threads {threads} \
            --output {output.vcf:q} \
            {input.vcfs:q} \
            2> {log:q}
        """

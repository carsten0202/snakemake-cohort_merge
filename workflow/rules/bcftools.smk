BCFTOOLS_MODULES = tuple(config["modules"]["bcftools"])

MERGED_VCF_PATTERN = str(
    Path(OUTPUT_DIRECTORY)
    / "by_chromosome"
    / "{chromosome}.vcf.gz"
)

INTERSECTION_DIRECTORY_PATTERN = str(
    Path(OUTPUT_DIRECTORY) / "intermediate" / "intersection" / "{chromosome}"
)

INTERSECTION_REPORT_PATTERN = str(
    Path(OUTPUT_DIRECTORY)
    / "qc"
    / "variant_intersection"
    / "{chromosome}.tsv"
)

SAMPLE_EXCLUSION_LIST = str(
    Path(OUTPUT_DIRECTORY) / "intermediate" / "sample_exclusions.txt"
)

UNFILTERED_BCF = str(
    Path(OUTPUT_DIRECTORY)
    / "intermediate"
    / "unfiltered"
    / f"{SUPER_COHORT}.bcf"
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


def intersection_vcfs(wildcards):
    directory = Path(
        INTERSECTION_DIRECTORY_PATTERN.format(chromosome=wildcards.chromosome)
    )
    return [str(directory / f"{index:04d}.vcf.gz") for index in range(len(COHORTS))]


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


rule intersect_variants:
    input:
        vcfs=chromosome_vcfs,
        indexes=chromosome_indexes,
    output:
        intersection=temp(directory(INTERSECTION_DIRECTORY_PATTERN))
    params:
        nfiles=len(COHORTS),
        write=",".join(str(index) for index in range(1, len(COHORTS) + 1)),
    threads:
        config["resources"]["intersect_variants"]["threads"]
    resources:
        mem_mb=config["resources"]["intersect_variants"]["mem_mb"]
    log:
        str(Path(OUTPUT_DIRECTORY) / "logs" / "intersect" / "{chromosome}.log")
    envmodules:
        *BCFTOOLS_MODULES
    shell:
        """
        bcftools isec \
            --collapse none \
            --nfiles ={params.nfiles} \
            --write {params.write} \
            --output-type z \
            --write-index=tbi \
            --threads {threads} \
            --prefix {output.intersection:q} \
            {input.vcfs:q} \
            2> {log:q}
        """


rule summarize_intersection:
    input:
        vcfs=chromosome_vcfs,
        indexes=chromosome_indexes,
        intersection=INTERSECTION_DIRECTORY_PATTERN,
    output:
        report=INTERSECTION_REPORT_PATTERN
    params:
        cohorts=COHORTS,
        intersection_vcfs=intersection_vcfs,
    threads:
        config["resources"]["summarize_intersection"]["threads"]
    resources:
        mem_mb=config["resources"]["summarize_intersection"]["mem_mb"]
    envmodules:
        *BCFTOOLS_MODULES
    script:
        "../scripts/summarize_variant_intersection.py"


rule merge_chromosome:
    input:
        intersection=INTERSECTION_DIRECTORY_PATTERN,
        report=INTERSECTION_REPORT_PATTERN,
    output:
        vcf=MERGED_VCF_PATTERN
    params:
        vcfs=intersection_vcfs
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
            {params.vcfs:q} \
            2> {log:q}
        """


if HAS_SAMPLE_EXCLUSIONS:
    rule prepare_sample_exclusions:
        input:
            exclusions=[
                config["cohorts"][cohort]["exclude_samples"]
                for cohort in EXCLUSION_COHORTS
            ],
            vcfs=[cohort_vcf(cohort, CHROMOSOMES[0]) for cohort in COHORTS],
        output:
            samples=temp(SAMPLE_EXCLUSION_LIST)
        params:
            cohorts=COHORTS,
            exclusion_cohorts=EXCLUSION_COHORTS,
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "sample_exclusions.prepare.log")
        envmodules:
            *BCFTOOLS_MODULES
        script:
            "../scripts/prepare_sample_exclusions.py"


    rule concatenate_chromosomes:
        input:
            vcfs=expand(MERGED_VCF_PATTERN, chromosome=CHROMOSOMES),
            indexes=expand(f"{MERGED_VCF_PATTERN}.tbi", chromosome=CHROMOSOMES),
        output:
            bcf=temp(UNFILTERED_BCF),
            index=temp(f"{UNFILTERED_BCF}.csi"),
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
                --output-type b \
                --write-index=csi \
                --threads {threads} \
                --output {output.bcf:q} \
                {input.vcfs:q} \
                2> {log:q}
            """


    rule exclude_samples:
        input:
            bcf=UNFILTERED_BCF,
            index=f"{UNFILTERED_BCF}.csi",
            samples=SAMPLE_EXCLUSION_LIST,
        output:
            vcf=FINAL_VCF
        threads:
            config["resources"]["exclude_samples"]["threads"]
        resources:
            mem_mb=config["resources"]["exclude_samples"]["mem_mb"]
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "sample_exclusions.log")
        envmodules:
            *BCFTOOLS_MODULES
        shell:
            """
            if [[ -s {input.samples:q} ]]; then
                bcftools view \
                    --force-samples \
                    --samples-file ^{input.samples:q} \
                    --output-type z \
                    --threads {threads} \
                    --output {output.vcf:q} \
                    {input.bcf:q} \
                    2> {log:q}
            else
                bcftools view \
                    --output-type z \
                    --threads {threads} \
                    --output {output.vcf:q} \
                    {input.bcf:q} \
                    2> {log:q}
            fi
            """
else:
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

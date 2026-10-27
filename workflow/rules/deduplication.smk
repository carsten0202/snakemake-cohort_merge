if DEDUPLICATION_ENABLED:
    rule prepare_sample_manifest:
        input:
            vcfs=[cohort_vcf(cohort, CHROMOSOMES[0]) for cohort in COHORTS],
            merged=MERGED_VCF_PATTERN.format(chromosome=CHROMOSOMES[0]),
            final=FINAL_VCF,
        output:
            manifest=SAMPLE_MANIFEST,
        params:
            cohorts=COHORTS,
            priorities=[config["cohorts"][cohort]["priority"] for cohort in COHORTS],
        threads:
            config["resources"]["prepare_sample_manifest"]["threads"]
        resources:
            mem_mb=config["resources"]["prepare_sample_manifest"]["mem_mb"]
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "sample_manifest.log")
        envmodules:
            *BCFTOOLS_MODULES
        script:
            "../scripts/prepare_sample_manifest.py"


    rule export_duplicate_pairs:
        input:
            pcrelate_data=f"{PCAIR_PREFIX}.pcrelate.RData",
            manifest=SAMPLE_MANIFEST,
        output:
            pairs=DUPLICATE_PAIRS,
        params:
            minimum_kinship=config["deduplication"]["minimum_kinship"],
        threads:
            config["resources"]["export_duplicate_pairs"]["threads"]
        resources:
            mem_mb=config["resources"]["export_duplicate_pairs"]["mem_mb"]
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "duplicate_pairs.log")
        envmodules:
            *PCAIR_MODULES
        script:
            "../scripts/export_duplicate_pairs.R"


    rule select_duplicate_representatives:
        input:
            manifest=SAMPLE_MANIFEST,
            pairs=DUPLICATE_PAIRS,
        output:
            decisions=SAMPLE_DECISIONS,
            removals=DUPLICATE_REMOVALS,
        params:
            random_seed=config["deduplication"]["random_seed"],
        threads:
            config["resources"]["select_duplicate_representatives"]["threads"]
        resources:
            mem_mb=config["resources"]["select_duplicate_representatives"]["mem_mb"]
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "duplicate_selection.log")
        script:
            "../scripts/select_duplicate_representatives.py"


    rule deduplicate_vcf:
        input:
            vcf=FINAL_VCF,
            index=f"{FINAL_VCF}.tbi",
            removals=DUPLICATE_REMOVALS,
        output:
            vcf=DEDUPLICATED_VCF,
            index=f"{DEDUPLICATED_VCF}.tbi",
        params:
            workers=lambda wildcards, threads: bcftools_workers(threads),
        threads:
            config["resources"]["deduplicate_vcf"]["threads"]
        resources:
            mem_mb=config["resources"]["deduplicate_vcf"]["mem_mb"]
        log:
            str(Path(OUTPUT_DIRECTORY) / "logs" / "deduplication.log")
        envmodules:
            *BCFTOOLS_MODULES
        shell:
            """
            selection=()
            if [[ -s {input.removals:q} ]]; then
                selection=(--samples-file ^{input.removals:q})
            fi
            bcftools view \
                "${{selection[@]}}" \
                --no-update \
                --output-type z \
                --write-index=tbi \
                --threads {params.workers} \
                --output {output.vcf:q} \
                {input.vcf:q} \
                2> {log:q}
            """


    ruleorder: deduplicate_vcf > index_vcf

PLINK_MODULES = tuple(config["modules"]["plink"])
PCAIR_MODULES = tuple(config["modules"]["quick_pcair"])


rule convert_vcf_to_plink:
    input:
        vcf=FINAL_VCF,
        index=f"{FINAL_VCF}.tbi",
    output:
        bed=f"{PLINK_PREFIX}.bed",
        bim=f"{PLINK_PREFIX}.bim",
        fam=f"{PLINK_PREFIX}.fam",
        plink_log=f"{PLINK_PREFIX}.log",
    params:
        prefix=PLINK_PREFIX,
    threads:
        config["resources"]["convert_vcf_to_plink"]["threads"]
    resources:
        mem_mb=config["resources"]["convert_vcf_to_plink"]["mem_mb"]
    log:
        str(Path(OUTPUT_DIRECTORY) / "logs" / "plink_conversion.log")
    envmodules:
        *PLINK_MODULES
    shell:
        """
        plink2 \
            --vcf {input.vcf:q} \
            --vcf-require-gt \
            --double-id \
            --autosome \
            --snps-only just-acgt \
            --max-alleles 2 \
            --make-bed \
            --threads {threads} \
            --memory {resources.mem_mb} \
            --out {params.prefix:q} \
            > {log:q} 2>&1
        """


rule run_quick_pcair:
    input:
        bed=f"{PLINK_PREFIX}.bed",
        bim=f"{PLINK_PREFIX}.bim",
        fam=f"{PLINK_PREFIX}.fam",
    output:
        gds=f"{PCAIR_PREFIX}.gds",
        king_plot=f"{PCAIR_PREFIX}.KINGkinship.png",
        king_matrix=f"{PCAIR_PREFIX}.KINGkinship.tsv",
        ld_snps=f"{PCAIR_PREFIX}.ldprune.snpids.txt",
        eigenvalues=f"{PCAIR_PREFIX}.eigenvalues.tsv",
        eigenvectors=f"{PCAIR_PREFIX}.eigenvectors.tsv",
        pcair_1v2=f"{PCAIR_PREFIX}.pcair_1v2.png",
        pcair_3v4=f"{PCAIR_PREFIX}.pcair_3v4.png",
        pcrelate_data=f"{PCAIR_PREFIX}.pcrelate.RData",
        pcrelate_plot=f"{PCAIR_PREFIX}.pcrelate_1v2.png",
    params:
        plink_prefix=PLINK_PREFIX,
        output_prefix=PCAIR_PREFIX,
    threads:
        config["resources"]["run_quick_pcair"]["threads"]
    resources:
        mem_mb=config["resources"]["run_quick_pcair"]["mem_mb"]
    log:
        str(Path(OUTPUT_DIRECTORY) / "logs" / "quick_pcair.log")
    envmodules:
        *PCAIR_MODULES
    shell:
        """
        rm -f {output:q}

        quick-pcair \
            --plink {params.plink_prefix:q} \
            --output {params.output_prefix:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """

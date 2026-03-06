rule variant_filtration:
    input:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
        **REF_FILES,
    output:
        vcf="results/vcfs/filtered.vcf.gz",
        tbi="results/vcfs/filtered.vcf.gz.tbi",
    params:
        filter_args=get_gatk_hard_filter_args(),
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/variant_filtration.txt"
    log:
        "logs/variant_filtration.txt"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            --output {output.vcf} \
            {params.filter_args} \
            --create-output-variant-index \
            --invalidate-previous-filters true \
            &> {log}
        """

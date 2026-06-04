rule variant_filtration:
    input:
        vcf=RAW_VCF,
        idx=RAW_VCF_INDEX,
        **REF_FILES,
    output:
        vcf=FILTERED_VCF,
        idx=FILTERED_VCF_INDEX,
    params:
        filter_args=get_gatk_hard_filter_args(),
        index_args=BCFTOOLS_INDEX_ARGS,
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
            --create-output-variant-index false \
            --invalidate-previous-filters true \
            &> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """

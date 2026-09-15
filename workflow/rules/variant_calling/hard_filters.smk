if LONG_CONTIG_MODE:

    # GATK/htsjdk has no CSI reader for VCFs: it looks only for a .tbi and
    # refuses block-compressed input without one ("An index is required but was
    # not found ... Support for unindexed block-compressed files has been
    # temporarily disabled"). tabix cannot index a contig longer than 2^29, so
    # there is no index that would satisfy GATK here. Filter the uncompressed
    # work VCF instead -- the same path every other GATK rule takes in this mode
    # -- then compress and CSI-index the result.
    #
    # This keeps results/vcfs/work/raw.vcf (temp) on disk until filtering
    # finishes, which costs roughly the uncompressed call set for large cohorts.

    rule variant_filtration:
        input:
            vcf=RAW_VCF_WORK,
            idx=RAW_VCF_WORK_INDEX,
            **REF_FILES,
        output:
            vcf=temp(FILTERED_VCF_WORK),
            idx=temp(FILTERED_VCF_WORK_INDEX),
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
                --invalidate-previous-filters true \
                &> {log}
            """


    rule compress_filtered_vcf:
        input:
            vcf=FILTERED_VCF_WORK,
            idx=FILTERED_VCF_WORK_INDEX,
        output:
            vcf=FILTERED_VCF,
            idx=FILTERED_VCF_INDEX,
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/bcftools.yaml"
        benchmark:
            "benchmarks/compress_filtered_vcf.txt"
        log:
            "logs/compress_filtered_vcf.txt"
        shell:
            """
            bcftools view -Oz -o {output.vcf} {input.vcf} 2> {log}
            bcftools index {params.index_args} {output.vcf} 2>> {log}
            """

else:

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

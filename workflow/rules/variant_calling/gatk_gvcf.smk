# Shared by both GATK paths (gatk.smk without intervals, gatk_intervals.smk
# with them).


def external_gvcf_input(wildcards):
    if not sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(
            f"Sample {wildcards.sample} does not have input_type 'gvcf'"
        )
    return {"gvcf": get_final_gvcf(wildcards.sample)}


# Outputs of the rule that finishes a sample's GATK gVCF (gatk_haplotypecaller
# or concat_interval_gvcfs). In long-contig mode that is an uncompressed work
# gVCF, which GATK reads and archive_gatk_gvcf compresses; otherwise it is the
# archived gVCF itself.
if LONG_CONTIG_MODE:
    GATK_SAMPLE_GVCF = {
        "gvcf": temp("results/gvcfs/work/{sample}.g.vcf"),
        "idx": temp("results/gvcfs/work/{sample}.g.vcf.idx"),
    }
else:
    GATK_SAMPLE_GVCF = {
        "gvcf": get_archive_gvcf("{sample}"),
        "idx": get_archive_gvcf_index("{sample}"),
    }


if LONG_CONTIG_MODE:

    rule normalize_external_gvcf_for_gatk:
        input:
            unpack(external_gvcf_input),
        output:
            gvcf=temp("results/gvcfs/work/external/{sample}.g.vcf"),
            idx=temp("results/gvcfs/work/external/{sample}.g.vcf.idx"),
        conda:
            "../../envs/gatk.yaml"
        benchmark:
            "benchmarks/normalize_external_gvcf_for_gatk/{sample}.txt"
        log:
            "logs/normalize_external_gvcf_for_gatk/{sample}.txt"
        shell:
            """
            bcftools view -O v -o {output.gvcf} {input.gvcf} 2> {log}
            gatk IndexFeatureFile -I {output.gvcf} &>> {log}
            """


    rule archive_gatk_gvcf:
        input:
            gvcf=lambda wc: get_gatk_work_gvcf(wc.sample),
            idx=lambda wc: get_gatk_work_gvcf_index(wc.sample),
        output:
            gvcf=get_archive_gvcf("{sample}"),
            idx=get_archive_gvcf_index("{sample}"),
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/bcftools.yaml"
        benchmark:
            "benchmarks/archive_gatk_gvcf/{sample}.txt"
        log:
            "logs/archive_gatk_gvcf/{sample}.txt"
        shell:
            """
            bcftools view -Oz -o {output.gvcf} {input.gvcf} 2> {log}
            bcftools index {params.index_args} {output.gvcf} 2>> {log}
            """

def haplotype_caller_input(wildcards):
    if sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(
            f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller"
        )

    return {
        **get_indexed_final_bam_input(wildcards.sample),
        **REF_FILES,
    }


<<<<<<< HEAD
rule gatk_haplotypecaller:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
        min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
    threads: 1
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_haplotypecaller/{sample}.txt"
    log:
        "logs/gatk_haplotypecaller/{sample}.txt"
    shell:
        """
        gatk HaplotypeCaller \
            --java-options '-Xmx{resources.mem_mb_reduced}m' \
            -R {input.ref} \
            -I {input.bam} \
            --read-index {input.bam_index} \
            -O {output.gvcf} \
            -ploidy {params.ploidy} \
            --native-pair-hmm-threads {threads} \
            --emit-ref-confidence GVCF \
            --min-pruning {params.min_pruning} \
            --min-dangling-branch-length {params.min_dangling} \
            &> {log}
        """
=======
def external_gvcf_input(wildcards):
    if not sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(
            f"Sample {wildcards.sample} does not have input_type 'gvcf'"
        )
    return {"gvcf": get_final_gvcf(wildcards.sample)}


if LONG_CONTIG_MODE:

    rule gatk_haplotypecaller:
        input:
            unpack(haplotype_caller_input),
        output:
            gvcf=temp("results/gvcfs/work/{sample}.g.vcf"),
            idx=temp("results/gvcfs/work/{sample}.g.vcf.idx"),
        params:
            ploidy=config["variant_calling"]["ploidy"],
            min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
            min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
        threads: 1
        conda:
            "../../envs/gatk.yaml"
        benchmark:
            "benchmarks/gatk_haplotypecaller/{sample}.txt"
        log:
            "logs/gatk_haplotypecaller/{sample}.txt"
        shell:
            """
            gatk HaplotypeCaller \
                --java-options '-Xmx{resources.mem_mb_reduced}m' \
                -R {input.ref} \
                -I {input.bam} \
                -O {output.gvcf} \
                -ploidy {params.ploidy} \
                --native-pair-hmm-threads {threads} \
                --emit-ref-confidence GVCF \
                --min-pruning {params.min_pruning} \
                --min-dangling-branch-length {params.min_dangling} \
                &> {log}
            """


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

else:

    rule gatk_haplotypecaller:
        input:
            unpack(haplotype_caller_input),
        output:
            gvcf="results/gvcfs/{sample}.g.vcf.gz",
            idx=get_compressed_vcf_index("results/gvcfs/{sample}.g.vcf.gz"),
        params:
            ploidy=config["variant_calling"]["ploidy"],
            min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
            min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
        threads: 1
        conda:
            "../../envs/gatk.yaml"
        benchmark:
            "benchmarks/gatk_haplotypecaller/{sample}.txt"
        log:
            "logs/gatk_haplotypecaller/{sample}.txt"
        shell:
            """
            gatk HaplotypeCaller \
                --java-options '-Xmx{resources.mem_mb_reduced}m' \
                -R {input.ref} \
                -I {input.bam} \
                -O {output.gvcf} \
                -ploidy {params.ploidy} \
                --native-pair-hmm-threads {threads} \
                --emit-ref-confidence GVCF \
                --min-pruning {params.min_pruning} \
                --min-dangling-branch-length {params.min_dangling} \
                &> {log}
            """
>>>>>>> 7b09cc8 (implement long-contig option that avoids tbi indexes for variant call files)

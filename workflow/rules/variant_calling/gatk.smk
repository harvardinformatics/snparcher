def haplotype_caller_input(wildcards):
    input_type = get_sample_input_type(wildcards.sample)

    if input_type == "gvcf":
        raise ValueError(
            f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller"
        )

    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


rule gatk_haplotypecaller:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
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
            --java-options '{params.java_opts}' \
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

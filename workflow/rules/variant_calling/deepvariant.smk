def deepvariant_input(wildcards):
    if sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call deepvariant")

    return {
        **get_indexed_final_bam_input(wildcards.sample),
        **REF_FILES,
    }


if LONG_CONTIG_MODE:

    rule deepvariant_call:
        input:
            unpack(deepvariant_input),
        output:
            gvcf=temp("results/deepvariant/{sample}.g.vcf"),
            vcf=temp("results/deepvariant/{sample}.vcf"),
        params:
            model_type=config["variant_calling"]["deepvariant"]["model_type"],
        threads: config["variant_calling"]["deepvariant"]["num_shards"]
        container:
            "docker://google/deepvariant:1.10.0"
        benchmark:
            "benchmarks/deepvariant_call/{sample}.txt"
        log:
            "logs/deepvariant_call/{sample}.txt"
        shell:
            """
            mkdir -p results/deepvariant/{wildcards.sample}
            /opt/deepvariant/bin/run_deepvariant \
                --ref {input.ref} \
                --reads {input.bam} \
                --output_vcf {output.vcf} \
                --output_gvcf {output.gvcf} \
                --model_type {params.model_type} \
                --num_shards {threads} \
                --intermediate_results_dir results/deepvariant/{wildcards.sample} \
                &> {log}
            """


    rule archive_deepvariant_gvcf:
        input:
            gvcf="results/deepvariant/{sample}.g.vcf",
        output:
            gvcf="results/gvcfs/{sample}.g.vcf.gz",
            idx=get_compressed_vcf_index("results/gvcfs/{sample}.g.vcf.gz"),
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/bcftools.yaml"
        benchmark:
            "benchmarks/archive_deepvariant_gvcf/{sample}.txt"
        log:
            "logs/archive_deepvariant_gvcf/{sample}.txt"
        shell:
            """
            bcftools view -Oz -o {output.gvcf} {input.gvcf} 2> {log}
            bcftools index {params.index_args} {output.gvcf} 2>> {log}
            """

else:

    rule deepvariant_call:
        input:
            unpack(deepvariant_input),
        output:
            gvcf="results/gvcfs/{sample}.g.vcf.gz",
            tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
            vcf=temp("results/deepvariant/{sample}.vcf.gz"),
            vcf_tbi=temp("results/deepvariant/{sample}.vcf.gz.tbi"),
        params:
            model_type=config["variant_calling"]["deepvariant"]["model_type"],
        threads: config["variant_calling"]["deepvariant"]["num_shards"]
        container:
            "docker://google/deepvariant:1.10.0"
        benchmark:
            "benchmarks/deepvariant_call/{sample}.txt"
        log:
            "logs/deepvariant_call/{sample}.txt"
        shell:
            """
            mkdir -p results/deepvariant/{wildcards.sample}
            /opt/deepvariant/bin/run_deepvariant \
                --ref {input.ref} \
                --reads {input.bam} \
                --output_vcf {output.vcf} \
                --output_gvcf {output.gvcf} \
                --model_type {params.model_type} \
                --num_shards {threads} \
                --intermediate_results_dir results/deepvariant/{wildcards.sample} \
                &> {log}
            """

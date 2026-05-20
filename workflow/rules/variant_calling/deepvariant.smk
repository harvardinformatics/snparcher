def deepvariant_input(wildcards):
    if sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call deepvariant")

    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


def glnexus_joint_input(wc):
    gvcfs = [get_final_gvcf(s) for s in SAMPLES_ALL]
    return {
        "gvcfs": gvcfs,
        "tbis": [g + ".tbi" for g in gvcfs],
    }


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


rule glnexus_joint:
    input:
        unpack(glnexus_joint_input),
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        db="results/vcfs/GLnexus.DB",
    threads: 1
    conda:
        "../../envs/glnexus.yaml"
    benchmark:
        "benchmarks/glnexus_joint.txt"
    log:
        "logs/glnexus_joint.txt"
    shell:
        """
        glnexus_db="{params.db}"
        rm -rf "$glnexus_db"
        trap 'rm -rf "$glnexus_db"' EXIT

        glnexus_cli --dir "$glnexus_db" --config DeepVariant --threads {threads} {input.gvcfs} 2> {log} \
            | bcftools view -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

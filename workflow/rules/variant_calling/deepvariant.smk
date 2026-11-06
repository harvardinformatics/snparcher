from snakemake.settings.types import DeploymentMethod


def apptainer_enabled():
    settings = workflow.deployment_settings
    return settings is not None and DeploymentMethod.APPTAINER in settings.deployment_method


def deepvariant_input(wildcards):
    if sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call deepvariant")
    # run_deepvariant exists only inside the container image, and Snakemake
    # ignores the container directive unless Apptainer is enabled. Fail while
    # building the DAG rather than after mapping has finished.
    if not apptainer_enabled():
        raise WorkflowError(
            "variant_calling.tool 'deepvariant' runs DeepVariant from its container image, "
            "so Snakemake must run with Apptainer enabled. Add --use-apptainer to the "
            "command line, or 'use-apptainer: True' to your workflow profile."
        )

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
            # Intermediate tfrecords can take tens of GB per sample. Keep them in
            # a private directory under the working directory, which Apptainer
            # binds (a custom tmpdir may not be), and remove it when the job ends.
            WORK=$(mktemp -d "results/deepvariant/{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$WORK"' EXIT
            /opt/deepvariant/bin/run_deepvariant \
                --ref {input.ref} \
                --reads {input.bam} \
                --output_vcf {output.vcf} \
                --output_gvcf {output.gvcf} \
                --model_type {params.model_type} \
                --num_shards {threads} \
                --intermediate_results_dir "$WORK" \
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
            # Intermediate tfrecords can take tens of GB per sample. Keep them in
            # a private directory under the working directory, which Apptainer
            # binds (a custom tmpdir may not be), and remove it when the job ends.
            WORK=$(mktemp -d "results/deepvariant/{wildcards.sample}.XXXXXX")
            trap 'rm -rf "$WORK"' EXIT
            /opt/deepvariant/bin/run_deepvariant \
                --ref {input.ref} \
                --reads {input.bam} \
                --output_vcf {output.vcf} \
                --output_gvcf {output.gvcf} \
                --model_type {params.model_type} \
                --num_shards {threads} \
                --intermediate_results_dir "$WORK" \
                &> {log}
            """

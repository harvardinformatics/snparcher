def parabricks_haplotypecaller_input(wildcards):
    input_type = get_sample_input_type(wildcards.sample)

    if input_type == "gvcf":
        raise ValueError(
            f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotypecaller"
        )

    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


rule parabricks_haplotypecaller:
    input:
        unpack(parabricks_haplotypecaller_input),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        image=config["variant_calling"]["parabricks"]["container_image"],
        extra_args=config["variant_calling"]["parabricks"]["extra_args"],
        ploidy=config["variant_calling"]["ploidy"],
        num_cpu_threads=config["variant_calling"]["parabricks"]["num_cpu_threads"],
    threads: config["variant_calling"]["parabricks"]["num_cpu_threads"]
    resources:
        gpus=config["variant_calling"]["parabricks"]["num_gpus"]
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/parabricks_haplotypecaller/{sample}.txt"
    log:
        "logs/parabricks_haplotypecaller/{sample}.txt"
    shell:
        """
        apptainer exec --nv {params.image} \
            pbrun haplotypecaller \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf} \
            --gvcf \
            --tmp-dir {resources.tmpdir} \
            --num-gpus {resources.gpus} \
            --num-cpu-threads {params.num_cpu_threads} \
            --ploidy {params.ploidy} \
            {params.extra_args} \
            &> {log}
        tabix -p vcf {output.gvcf} 2>> {log}
        """

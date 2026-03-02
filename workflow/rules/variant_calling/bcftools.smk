def bcftools_call_input(wc):
    bams = [get_final_bam(s) for s in SAMPLES_WITH_BAM]
    return {
        "bams": bams,
        "bais": [f"{bam}.bai" for bam in bams],
        **REF_FILES,
    }


rule bcftools_call:
    input:
        unpack(bcftools_call_input),
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        min_mapq=config["variant_calling"]["bcftools"]["min_mapq"],
        min_baseq=config["variant_calling"]["bcftools"]["min_baseq"],
        max_depth=config["variant_calling"]["bcftools"]["max_depth"],
        ploidy=config["variant_calling"]["ploidy"],
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_call.txt"
    log:
        "logs/bcftools_call.txt"
    shell:
        """
        bcftools mpileup \
            -f {input.ref} \
            -q {params.min_mapq} \
            -Q {params.min_baseq} \
            -d {params.max_depth} \
            -Ou {input.bams} 2> {log} \
        | bcftools call \
            -m \
            --ploidy {params.ploidy} \
            -v \
            -Oz \
            -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

# Region checkpoint and helpers are in contig_regions.smk.


def get_bcftools_region_vcfs(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand("results/vcfs/regions/{region_id}.vcf.gz", region_id=region_ids)


def get_bcftools_region_vcf_tbis(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand(
        get_compressed_vcf_index("results/vcfs/regions/{region_id}.vcf.gz"),
        region_id=region_ids,
    )


rule bcftools_call:
    input:
        unpack(bcftools_call_input),
        regions_tsv="results/vcfs/regions/regions.tsv",
    output:
        vcf=temp("results/vcfs/regions/{region_id}.vcf.gz"),
        idx=temp(get_compressed_vcf_index("results/vcfs/regions/{region_id}.vcf.gz")),
    params:
        min_mapq=config["variant_calling"]["bcftools"]["min_mapq"],
        min_baseq=config["variant_calling"]["bcftools"]["min_baseq"],
        max_depth=config["variant_calling"]["bcftools"]["max_depth"],
        ploidy=config["variant_calling"]["ploidy"],
        contig=lambda wc: get_bcftools_region_name(wc.region_id),
        index_args=BCFTOOLS_INDEX_ARGS,
    threads: 1
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_call/{region_id}.txt"
    log:
        "logs/bcftools_call/{region_id}.txt"
    shell:
        """
        bcftools mpileup \
            -f {input.ref} \
            -q {params.min_mapq} \
            -Q {params.min_baseq} \
            -d {params.max_depth} \
            -r {params.contig} \
            --threads {threads} \
            -Ou {input.bams} 2> {log} \
        | bcftools call \
            -m \
            --ploidy {params.ploidy} \
            --threads {threads} \
            -v \
            -Oz \
            -o {output.vcf} - 2>> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """


rule bcftools_concat_regions:
    input:
        vcfs=get_bcftools_region_vcfs,
        tbis=get_bcftools_region_vcf_tbis,
    output:
        vcf=temp(RAW_VCF),
        idx=temp(RAW_VCF_INDEX),
    params:
        index_args=BCFTOOLS_INDEX_ARGS,
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_concat_regions.txt"
    log:
        "logs/bcftools_concat_regions.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} - 2>> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """

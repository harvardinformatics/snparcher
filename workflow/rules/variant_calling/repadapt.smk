# RepAdapt's calling model (https://github.com/RepAdapt/nextflow_snp_calling_linux),
# run on snpArcher's final BAMs. Region checkpoint and helpers are in
# contig_regions.smk.
#
# Differences from RepAdapt's literal command, neither of which changes calls:
# - mpileup -q 10 instead of -q 5 stands in for RepAdapt's BAM-level
#   `samtools view -q 10`; on BAMs that already went through that filter the
#   two select the same reads.
# - --ploidy is passed; bcftools 1.16 accepts only its predefined aliases, so
#   ploidy is validated as 1 or 2 in common.smk. Diploid calls are identical
#   to omitting it.
#
# bcftools is pinned to 1.16, RepAdapt's version: 1.16 writes INFO/MQ as an
# Integer and newer versions as a Float, which can move sites across the
# MQ < 30 filter. -Q, -d and --threads are deliberately left at RepAdapt's
# (bcftools 1.16's) defaults: min base quality 1, max depth 250.

localrules: repadapt_region_list


REPADAPT_REGION_VCF = "results/vcfs/regions/repadapt/{region_id}.vcf.gz"


def get_repadapt_region_vcfs(wc):
    return expand(REPADAPT_REGION_VCF, region_id=get_bcftools_region_ids(wc))


rule repadapt_call:
    input:
        unpack(bcftools_call_input),
        regions_tsv="results/vcfs/regions/regions.tsv",
    output:
        vcf=temp(REPADAPT_REGION_VCF),
    wildcard_constraints:
        region_id=r"L\d{6}",
    params:
        ploidy=config["variant_calling"]["ploidy"],
        contig=lambda wc: get_bcftools_region_name(wc.region_id),
    threads: 1
    conda:
        "../../envs/repadapt/bcftools.yaml"
    benchmark:
        "benchmarks/repadapt_call/{region_id}.txt"
    log:
        "logs/repadapt_call/{region_id}.txt"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -r {params.contig} {input.bams} \
            -q 10 -I -a FMT/AD,FMT/DP 2> {log} \
        | bcftools call -G - -f GQ -mv --ploidy {params.ploidy} \
            -Oz -o {output.vcf} - 2>> {log}
        """


rule repadapt_region_list:
    """Write the region VCFs, in .fai order, to a file for bcftools concat -f.
    A file list keeps fragmented references within command-length limits."""
    input:
        vcfs=get_repadapt_region_vcfs,
    output:
        txt=temp("results/vcfs/regions/repadapt/regions.list"),
    run:
        Path(output.txt).write_text("".join(f"{vcf}\n" for vcf in input.vcfs))


rule repadapt_concat_regions:
    input:
        vcfs=get_repadapt_region_vcfs,
        vcf_list="results/vcfs/regions/repadapt/regions.list",
    output:
        # Not temp: the raw call set is kept alongside the filtered one.
        vcf=RAW_VCF,
        idx=RAW_VCF_INDEX,
    params:
        index_args=BCFTOOLS_INDEX_ARGS,
    conda:
        "../../envs/repadapt/bcftools.yaml"
    benchmark:
        "benchmarks/repadapt_concat_regions.txt"
    log:
        "logs/repadapt_concat_regions.txt"
    shell:
        # Regions are whole contigs in .fai order, so a plain concat is sorted.
        """
        bcftools concat -f {input.vcf_list} -Oz -o {output.vcf} 2> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """


if APPLY_REPADAPT_FILTER:

    # RepAdapt drops records matching `AC=AN || MQ < 30`. Here they are
    # soft-filtered instead, so the PASS records are exactly RepAdapt's call
    # set. bcftools sets PASS on records whose FILTER was '.', and `-m +`
    # appends LowMQ to AllHomAlt (a lone PASS is replaced, not appended to).
    rule repadapt_filter:
        input:
            vcf=RAW_VCF,
            idx=RAW_VCF_INDEX,
        output:
            vcf=FILTERED_VCF,
            idx=FILTERED_VCF_INDEX,
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/repadapt/bcftools.yaml"
        benchmark:
            "benchmarks/repadapt_filter.txt"
        log:
            "logs/repadapt_filter.txt"
        shell:
            """
            bcftools filter -s AllHomAlt -e 'AC=AN' -Ou {input.vcf} 2> {log} \
            | bcftools filter -m + -s LowMQ -e 'MQ<30' -Oz -o {output.vcf} - 2>> {log}
            bcftools index {params.index_args} {output.vcf} 2>> {log}
            """

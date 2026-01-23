localrules: create_db_mapfile

def get_java_mem(wildcards, resources):
    """Generate Java options string with heap size from resources.mem_mb.

    Calculates Java heap size (90% of mem_mb) to leave room for JVM overhead.

    Args:
        wildcards: Snakemake wildcards object (unused, required for params signature)
        resources: Snakemake resources object

    Returns:
        Java options string, e.g. "-Xmx7200m"

    Example:
        rule example:
            resources:
                mem_mb=8000,
            params:
                java_mem=get_java_mem,
            shell:
                "gatk Tool --java-options {params.java_mem} ..."
    """
    import math

    mem_mb = getattr(resources, "mem_mb")
    heap_mb = math.floor(mem_mb * 0.9)
    return f"-Xmx{heap_mb}m"


def haplotype_caller_input(wildcards):
    input_type = samples_df.loc[wildcards.sample, "input_type"]
    ref_files = REF_FILES
    if input_type == "bam":
        bam = samples_df.loc[wildcards.sample, "input"]
    else:  # fastq, srr
        bam = f"results/mapping/{wildcards.sample}.bam"
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


def get_interval_gvcfs(wc):
    """Get interval list numbers from checkpoint output."""
    checkpoint_output = checkpoints.create_gvcf_intervals.get(**wc).output[0]
    with checkpoint_output.open() as f:
        lines = [l.strip() for l in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]
    intervals = [f.replace("-scattered.interval_list", "") for f in list_files]

    return {
        "gvcfs": expand(
            "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz",
            sample=wc.sample,
            interval=intervals,
        ),
        "tbis": expand(
            "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz.tbi",
            sample=wc.sample,
            interval=intervals,
        ),
    }

def get_interval_vcfs(wc):
    """Get interval VCF list numbers from checkpoint output."""
    checkpoint_output = checkpoints.create_db_intervals.get(**wc).output[0]
    with checkpoint_output.open() as f:
        lines = [l.strip() for l in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]
    intervals = [f.replace("-scattered.interval_list", "") for f in list_files]

    return {
        "vcfs": expand(
            "results/vcfs/intervals/filters_applied_{interval}.vcf.gz",
            interval=intervals,
        ),
        "tbis": expand(
            "results/vcfs/intervals/filters_applied_{interval}.vcf.gz.tbi",
            interval=intervals,
        ),
    }


def get_gvcfs_for_db(wc):
    """Get all sample gVCFs for genomics DB import."""
    gvcfs = []
    tbis = []
    for sample in samples_df.index:
        input_type = samples_df.loc[sample, "input_type"]
        if input_type == "gvcf":
            gvcf = samples_df.loc[sample, "input"]
        else:
            gvcf = f"results/gvcfs/{sample}.g.vcf.gz"
        gvcfs.append(gvcf)
        tbis.append(gvcf + ".tbi")
    return {
        "gvcfs": gvcfs,
        "tbis": tbis,
        "interval": f"results/intervals/db/{wc.interval}-scattered.interval_list",
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }

rule gatk_haplotypecaller_interval:
    input:
        unpack(haplotype_caller_input),
    output:
        vcf=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz"),
        tbi=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz.tbi"),
    params:
        java_mem=lambda wc, res: get_java_mem,
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
        min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_haplotypecaller/{sample}/{interval}.benchmark.txt"
    log:
        "logs/gatk_haplotypecaller/{sample}/{interval}.log.txt",
    shell:
        """
        gatk HaplotypeCaller \
        --java-options {params.java_mem} \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -L {input.interval} \
        -ploidy {params.ploidy} \
        --emit-ref-confidence GVCF \
        --min-pruning {params.min_pruning} \
        --min-dangling-branch-length {params.min_dangling} &> {log}
        """

rule concat_interval_gvcfs:
    input:
        unpack(get_interval_gvcfs)
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/concat_interval_gvcfs/{sample}.txt"
    log:
        "logs/concat_interval_gvcfs/{sample}.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.gvcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir} -Oz -o {output.gvcf} - 2>> {log}
        tabix -p vcf {output.gvcf} 2>> {log}
        """

rule create_db_mapfile:
    input:
        gvcfs=lambda wc: [
            samples_df.loc[s, "input"] if samples_df.loc[s, "input_type"] == "gvcf"
            else f"results/gvcfs_norm/{s}.g.vcf.gz"
            for s in samples_df.index
        ],
    output:
        mapfile="results/genomics_db/mapfile.txt",
    run:
        with open(output.mapfile, "w") as f:
            for gvcf in input.gvcfs:
                sample = os.path.basename(gvcf).replace(".g.vcf.gz", "")
                print(sample, gvcf, sep="\t", file=f)


rule gatk_genomics_db_import:
    input:
        unpack(get_gvcfs_for_db),
    output:
        db=temp(directory("results/gatk_genomics_db/L{interval}")),
        tar="results/gatk_genomics_db/L{interval}.tar",
    params:
        java_mem=get_java_mem,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genomics_db_import/{interval}.txt"
    log:
        "logs/gatk_genomics_db_import/{interval}.txt"
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '{params.java_mem}' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            -L {input.interval} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            &> {log}
        tar -cf {output.tar} {output.db} &>> {log}
        """

rule gatk_variant_filtration:
    input:
        unpack(REF_FILES),
        vcf="results/vcfs/intervals/L{interval}.vcf.gz",
        tbi="results/vcfs/intervals/L{interval}.vcf.gz.tbi",
    output:
        vcf=temp("results/vcfs/intervals/filters_applied_{interval}.vcf.gz"),
        tbi=temp("results/vcfs/intervals/filters_applied_{interval}.vcf.gz.tbi"),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_variant_filtration/{interval}.txt"
    log:
        "logs/gatk_variant_filtration/{interval}.txt"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            --output {output.vcf} \
            --filter-name "RPRS_filter" \
            --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" \
            --filter-name "FS_SOR_filter" \
            --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \
            --filter-name "MQ_filter" \
            --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
            --filter-name "QUAL_filter" \
            --filter-expression "QUAL < 30.0" \
            --create-output-variant-index \
            --invalidate-previous-filters true \
            &> {log}
        """


rule concat_interval_vcfs:
    input:
        unpack(get_interval_vcfs),
    output:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/concat_interval_vcfs/benchmark.txt"
    log:
        "logs/concat_interval_vcfs/log.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir} -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """
localrules: create_db_mapfile


def _java_opts_from_resources(resources, default_mem_mb=4096):
    mem_mb = getattr(resources, "mem_mb", default_mem_mb)
    try:
        mem_mb = int(float(mem_mb))
    except (TypeError, ValueError):
        mem_mb = default_mem_mb
    return f"-Xmx{int(mem_mb * 0.9)}m"


def get_db_intervals(wc):
    """Get DB interval IDs from checkpoint output."""
    intervals_file = "results/intervals/db/intervals.txt"

    if exists(intervals_file):
        with open(intervals_file) as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        checkpoint_output = checkpoints.create_db_intervals.get(**wc).output[0]
        with open(checkpoint_output) as f:
            lines = [l.strip() for l in f.readlines()]

    list_files = [os.path.basename(x) for x in lines]
    return [f.replace("-scattered.interval_list", "") for f in list_files]


def get_interval_vcfs(wc):
    """Get filtered interval VCF files."""
    intervals = get_db_intervals(wc)
    return expand(
        "results/vcfs/intervals/filtered_{interval}.vcf.gz",
        interval=intervals,
    )


def get_interval_vcf_tbis(wc):
    """Get filtered interval VCF index files."""
    intervals = get_db_intervals(wc)
    return expand(
        "results/vcfs/intervals/filtered_{interval}.vcf.gz.tbi",
        interval=intervals,
    )


def get_gvcfs_for_db(wc):
    gvcfs = [get_final_gvcf(s) for s in SAMPLES_ALL]
    return {
        "gvcfs": gvcfs,
        "tbis": [g + ".tbi" for g in gvcfs],
        "interval": f"results/intervals/db/{wc.interval}-scattered.interval_list",
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }


rule create_db_mapfile:
    input:
        gvcfs=[get_final_gvcf(s) for s in SAMPLES_ALL],
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
        java_mem=lambda wildcards, resources: _java_opts_from_resources(resources),
    resources:
        mem_mb=4096,
    conda:
        "../../envs/gatk.yaml"
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
        tar -cf {output.tar} {output.db} >> {log} 2>&1
        """


rule gatk_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db/L{interval}.tar",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/intervals/L{interval}.vcf.gz"),
        tbi=temp("results/vcfs/intervals/L{interval}.vcf.gz.tbi"),
    params:
        java_opts=lambda wildcards, resources: _java_opts_from_resources(resources),
        het_prior=config["variant_calling"]["gatk"]["het_prior"],
        db=subpath(input.db, strip_suffix=".tar"),
    resources:
        mem_mb=4096,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genotype_gvcfs/{interval}.txt"
    log:
        "logs/gatk_genotype_gvcfs/{interval}.txt"
    shell:
        """
        tar -xf {input.db}
        gatk GenotypeGVCFs \
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            --heterozygosity {params.het_prior} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} \
            &> {log}
        """


rule gatk_variant_filtration:
    input:
        vcf="results/vcfs/intervals/L{interval}.vcf.gz",
        tbi="results/vcfs/intervals/L{interval}.vcf.gz.tbi",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/intervals/filtered_{interval}.vcf.gz"),
        tbi=temp("results/vcfs/intervals/filtered_{interval}.vcf.gz.tbi"),
    params:
        filter_args=get_gatk_hard_filter_args(),
    conda:
        "../../envs/gatk.yaml"
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
            {params.filter_args} \
            --create-output-variant-index \
            --invalidate-previous-filters true \
            &> {log}
        """


rule concat_interval_vcfs:
    input:
        vcfs=get_interval_vcfs,
        tbis=get_interval_vcf_tbis,
    output:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
    conda:
        "../../envs/bcftools.yaml"
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

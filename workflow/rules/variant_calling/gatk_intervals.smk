localrules: create_db_mapfile

def haplotype_caller_input(wildcards):
    input_type = samples_df.loc[wildcards.sample, "input_type"]
    
    if input_type == "gvcf":
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller")
    
    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        "interval": f"results/intervals/gvcf/{wildcards.interval}-scattered.interval_list",
        **REF_FILES,
    }


def get_interval_gvcfs(wc):
    """Get interval gVCF files for a sample."""
    intervals_file = "results/intervals/gvcf/intervals.txt"
    
    if exists(intervals_file):
        # Read directly from file if it exists
        with open(intervals_file) as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        # Fall back to checkpoint mechanism to trigger creation
        checkpoint_output = checkpoints.create_gvcf_intervals.get(**wc).output[0]
        with open(checkpoint_output) as f:
            lines = [l.strip() for l in f.readlines()]
    
    list_files = [os.path.basename(x) for x in lines]
    intervals = [f.replace("-scattered.interval_list", "") for f in list_files]
    return expand(
        "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz",
        sample=wc.sample,
        interval=intervals,
    )


def get_interval_gvcf_tbis(wc):
    """Get interval gVCF index files for a sample."""
    intervals_file = "results/intervals/gvcf/intervals.txt"
    
    if exists(intervals_file):
        with open(intervals_file) as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        checkpoint_output = checkpoints.create_gvcf_intervals.get(**wc).output[0]
        with open(checkpoint_output) as f:
            lines = [l.strip() for l in f.readlines()]
    
    list_files = [os.path.basename(x) for x in lines]
    intervals = [f.replace("-scattered.interval_list", "") for f in list_files]
    return expand(
        "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz.tbi",
        sample=wc.sample,
        interval=intervals,
    )


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
    print(f"get_gvcfs_for_db: gvcfs={gvcfs}")
    return {
        "gvcfs": gvcfs,
        "tbis": [g + ".tbi" for g in gvcfs],
        "interval": f"results/intervals/db/{wc.interval}-scattered.interval_list",
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }

rule gatk_haplotypecaller_interval:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz"),
        tbi=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz.tbi"),
    params:
        java_mem=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
        min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
    resources:
        mem_mb=4096,
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
        gvcfs=get_interval_gvcfs,
        tbis=get_interval_gvcf_tbis,
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    conda:
        "../../envs/bcftools.yaml"
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
        java_mem=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
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
        tar -cf {output.tar} {output.db} &>> {log}
        """

rule gatk_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db/L{interval}.tar",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/intervals/L{interval}.vcf.gz"),
        tbi=temp("results/vcfs/intervals/L{interval}.vcf.gz.tbi"),
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
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
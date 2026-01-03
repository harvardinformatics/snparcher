"""
GATK variant calling rules for snpArcher v2 (interval-based mode).

Used for large genomes where parallelization improves performance.
Requires interval rules from call/intervals.smk.
"""


def _get_all_gvcfs_for_db(wildcards=None):
    if config.get("normalize_gvcfs", True):
        return expand("results/gvcfs_norm/{sample}.g.vcf.gz", sample=SAMPLE_IDS)
    return expand("results/gvcfs/{sample}.g.vcf.gz", sample=SAMPLE_IDS)


def _get_all_gvcf_indices_for_db(wildcards=None):
    if config.get("normalize_gvcfs", True):
        return expand("results/gvcfs_norm/{sample}.g.vcf.gz.tbi", sample=SAMPLE_IDS)
    return expand("results/gvcfs/{sample}.g.vcf.gz.tbi", sample=SAMPLE_IDS)


def _get_interval_gvcfs(wildcards):
    """Get interval GVCFs after checkpoint."""
    checkpoint_output = checkpoints.intervals_gvcf.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        lines = [line.strip() for line in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]
    list_numbers = [f.replace("-scattered.interval_list", "") for f in list_files]
    return expand(
        "results/interval_gvcfs/{sample}/{l}.raw.g.vcf.gz",
        sample=wildcards.sample,
        l=list_numbers,
    )


def _get_interval_gvcfs_idx(wildcards):
    return [f + ".tbi" for f in _get_interval_gvcfs(wildcards)]


def _get_interval_vcfs(wildcards):
    """Get interval VCFs after checkpoint."""
    checkpoint_output = checkpoints.intervals_db.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        lines = [line.strip() for line in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]
    list_numbers = [f.replace("-scattered.interval_list", "") for f in list_files]
    return expand("results/vcfs/intervals/filtered_L{l}.vcf.gz", l=list_numbers)


def _get_interval_vcfs_idx(wildcards):
    return [f + ".tbi" for f in _get_interval_vcfs(wildcards)]


localrules:
    call_create_db_mapfile_intervals,


rule call_haplotypecaller_interval:
    """Generate GVCF for a single sample and interval using GATK HaplotypeCaller."""
    input:
        unpack(get_final_bam),
        ref=f"results/reference/{REF_FILE}",
        fai=f"results/reference/{REF_FILE}.fai",
        dictf=f"results/reference/{REF_NAME}.dict",
        interval="results/intervals/gvcf_intervals/{l}-scattered.interval_list",
    output:
        gvcf="results/interval_gvcfs/{sample}/{l}.raw.g.vcf.gz",
        tbi="results/interval_gvcfs/{sample}/{l}.raw.g.vcf.gz.tbi",
    params:
        minPrun=config.get("minP", 1),
        minDang=config.get("minD", 1),
        ploidy=config.get("ploidy", 2),
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_haplotypecaller_interval/{sample}/{l}.txt"
    benchmark:
        "benchmarks/call_haplotypecaller_interval/{sample}_{l}.txt"
    threads: 4
    resources:
        mem_mb=8000,
        mem_mb_reduced=7200,
    shell:
        """
        gatk HaplotypeCaller \
            --java-options "-Xmx{resources.mem_mb_reduced}m" \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -L {input.interval} \
            -ploidy {params.ploidy} \
            --emit-ref-confidence GVCF \
            --min-pruning {params.minPrun} \
            --min-dangling-branch-length {params.minDang} \
            &> {log}
        """


rule call_concat_gvcfs:
    input:
        gvcfs=_get_interval_gvcfs,
        tbis=_get_interval_gvcfs_idx,
    output:
        gvcf=temp("results/gvcfs/{sample}.g.vcf.gz"),
        tbi=temp("results/gvcfs/{sample}.g.vcf.gz.tbi"),
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/call_concat_gvcfs/{sample}.txt"
    benchmark:
        "benchmarks/call_concat_gvcfs/{sample}.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.gvcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir} -Oz -o {output.gvcf} - 2>> {log}
        tabix -p vcf {output.gvcf} 2>> {log}
        """


rule call_bcftools_norm:
    """Normalize multiallelic sites in GVCFs."""
    input:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
    output:
        gvcf="results/gvcfs_norm/{sample}.g.vcf.gz",
        tbi="results/gvcfs_norm/{sample}.g.vcf.gz.tbi",
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/call_bcftools_norm/{sample}.txt"
    benchmark:
        "benchmarks/call_bcftools_norm/{sample}.txt"
    shell:
        """
        bcftools norm -m +any -Oz -o {output.gvcf} {input.gvcf} 2> {log}
        tabix -p vcf {output.gvcf} 2>> {log}
        """


rule call_create_db_mapfile_intervals:
    input:
        gvcfs=_get_all_gvcfs_for_db,
    output:
        db_mapfile="results/genomics_db/DB_mapfile.txt",
    run:
        with open(output.db_mapfile, "w") as f:
            for file_path in input.gvcfs:
                sample_name = os.path.basename(file_path).replace(".g.vcf.gz", "")
                print(sample_name, file_path, sep="\t", file=f)


rule call_genomicsdb_import_interval:
    input:
        gvcfs=_get_all_gvcfs_for_db(),
        tbis=_get_all_gvcf_indices_for_db(),
        interval="results/intervals/db_intervals/{l}-scattered.interval_list",
        db_mapfile="results/genomics_db/DB_mapfile.txt",
    output:
        db=temp(directory("results/genomics_db/DB_L{l}")),
        tar=temp("results/genomics_db/DB_L{l}.tar"),
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_genomicsdb_import_interval/{l}.txt"
    benchmark:
        "benchmarks/call_genomicsdb_import_interval/{l}.txt"
    threads: 4
    resources:
        mem_mb=16000,
        mem_mb_reduced=14400,
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '-Xmx{resources.mem_mb_reduced}m -Xms{resources.mem_mb_reduced}m' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            -L {input.interval} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} &> {log}
        
        tar -cf {output.tar} {output.db}
        """


rule call_genotype_gvcfs_interval:
    """Joint genotyping from GenomicsDB for a single interval."""
    input:
        db="results/genomics_db/DB_L{l}.tar",
        ref=f"results/reference/{REF_FILE}",
        fai=f"results/reference/{REF_FILE}.fai",
        dictf=f"results/reference/{REF_NAME}.dict",
    output:
        vcf=temp("results/vcfs/intervals/L{l}.vcf.gz"),
        vcfidx=temp("results/vcfs/intervals/L{l}.vcf.gz.tbi"),
    params:
        het=config.get("het_prior", 0.005),
        db=lambda wc, input: input.db[:-4],
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_genotype_gvcfs_interval/{l}.txt"
    benchmark:
        "benchmarks/call_genotype_gvcfs_interval/{l}.txt"
    threads: 4
    resources:
        mem_mb=16000,
        mem_mb_reduced=14400,
    shell:
        """
        tar -xf {input.db}
        gatk GenotypeGVCFs \
            --java-options '-Xmx{resources.mem_mb_reduced}m -Xms{resources.mem_mb_reduced}m' \
            -R {input.ref} \
            --heterozygosity {params.het} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} &> {log}
        """


rule call_filter_vcf_interval:
    """Apply GATK recommended hard filters to interval VCF."""
    input:
        vcf="results/vcfs/intervals/L{l}.vcf.gz",
        vcfidx="results/vcfs/intervals/L{l}.vcf.gz.tbi",
        ref=f"results/reference/{REF_FILE}",
        fai=f"results/reference/{REF_FILE}.fai",
        dictf=f"results/reference/{REF_NAME}.dict",
    output:
        vcf=temp("results/vcfs/intervals/filtered_L{l}.vcf.gz"),
        vcfidx=temp("results/vcfs/intervals/filtered_L{l}.vcf.gz.tbi"),
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_filter_vcf_interval/{l}.txt"
    benchmark:
        "benchmarks/call_filter_vcf_interval/{l}.txt"
    resources:
        mem_mb=8000,
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
            --invalidate-previous-filters true &> {log}
        """


rule call_gather_vcfs:
    input:
        vcfs=_get_interval_vcfs,
        tbis=_get_interval_vcfs_idx,
    output:
        vcf=f"results/{REF_NAME}_raw.vcf.gz",
        vcfidx=f"results/{REF_NAME}_raw.vcf.gz.tbi",
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/call_gather_vcfs.txt"
    benchmark:
        "benchmarks/call_gather_vcfs.txt"
    threads: 4
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir} -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

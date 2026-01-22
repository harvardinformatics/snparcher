"""
GATK variant calling rules for snpArcher v2 (non-interval mode).
"""


def _get_all_gvcfs():
    return expand("results/gvcfs/{sample}.g.vcf.gz", sample=SAMPLE_IDS)


def _get_all_gvcf_indices():
    return expand("results/gvcfs/{sample}.g.vcf.gz.tbi", sample=SAMPLE_IDS)


localrules:
    call_create_db_mapfile,


rule call_haplotypecaller:
    """Generate GVCF for a single sample using GATK HaplotypeCaller."""
    wildcard_constraints:
        sample="|".join(BAM_REQUIRED_SAMPLES) if BAM_REQUIRED_SAMPLES else "NOMATCH",
    input:
        unpack(get_final_bam),
        unpack(get_ref_bundle),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        min_prun=config["variant_calling"]["gatk"]["min_pruning"],
        min_dang=config["variant_calling"]["gatk"]["min_dangling"],
        ploidy=config["variant_calling"]["gatk"]["ploidy"],
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_haplotypecaller/{sample}.txt"
    benchmark:
        "benchmarks/call_haplotypecaller/{sample}.txt"
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
        gatk HaplotypeCaller \
            --java-options "-Xmx{params.java_mem}m" \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ploidy {params.ploidy} \
            --emit-ref-confidence GVCF \
            --min-pruning {params.min_prun} \
            --min-dangling-branch-length {params.min_dang} \
            &> {log}
        """


rule call_create_db_mapfile:
    input:
        gvcfs=_get_all_gvcfs(),
    output:
        db_mapfile="results/genomics_db/DB_mapfile.txt",
    run:
        with open(output.db_mapfile, "w") as f:
            for file_path in input.gvcfs:
                sample_name = os.path.basename(file_path).replace(".g.vcf.gz", "")
                print(sample_name, file_path, sep="\t", file=f)


rule call_prepare_db_intervals:
    """Create interval list for GenomicsDBImport from reference FAI."""
    input:
        unpack(get_ref_bundle),
    output:
        intervals="results/genomics_db/db_intervals.list",
    run:
        with open(output.intervals, "w") as out:
            with open(input.fai, "r") as f:
                for line in f:
                    fields = line.strip().split()
                    chrom, end = fields[0], fields[1]
                    print(f"{chrom}:1-{end}", file=out)


rule call_genomicsdb_import:
    input:
        gvcfs=_get_all_gvcfs(),
        tbis=_get_all_gvcf_indices(),
        db_mapfile="results/genomics_db/DB_mapfile.txt",
        intervals="results/genomics_db/db_intervals.list",
    output:
        db=temp(directory("results/genomics_db/DB")),
        tar=temp("results/genomics_db/DB.tar"),
    params:
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_genomicsdb_import.txt"
    benchmark:
        "benchmarks/call_genomicsdb_import.txt"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options "-Xmx{params.java_mem}m" \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            -L {input.intervals} \
            --merge-input-intervals \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} &> {log}
        
        tar -cf {output.tar} {output.db} &>> {log}
        """


rule call_genotype_gvcfs:
    """Joint genotyping using GenomicsDB."""
    input:
        unpack(get_ref_bundle),
        db="results/genomics_db/DB.tar",
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        vcfidx=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        het=config["variant_calling"]["gatk"]["het_prior"],
        db=lambda wc, input: input.db[:-4],
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_genotype_gvcfs.txt"
    benchmark:
        "benchmarks/call_genotype_gvcfs.txt"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        tar -xf {input.db} &> {log}
        gatk GenotypeGVCFs \
            --java-options "-Xmx{params.java_mem}m" \
            -R {input.ref} \
            --heterozygosity {params.het} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} &>> {log}
        """


rule call_filter_vcf:
    """Apply GATK recommended hard filters to VCF."""
    input:
        unpack(get_ref_bundle),
        vcf="results/vcfs/raw.vcf.gz",
        vcfidx="results/vcfs/raw.vcf.gz.tbi",
    output:
        vcf=temp("results/vcfs/filtered.vcf.gz"),
        vcfidx=temp("results/vcfs/filtered.vcf.gz.tbi"),
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_filter_vcf.txt"
    benchmark:
        "benchmarks/call_filter_vcf.txt"
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


rule call_sort_vcf:
    """Sort and index final VCF."""
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        vcfidx="results/vcfs/filtered.vcf.gz.tbi",
    output:
        vcf=f"results/{REF_NAME}_raw.vcf.gz",
        vcfidx=f"results/{REF_NAME}_raw.vcf.gz.tbi",
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/call_sort_vcf.txt"
    benchmark:
        "benchmarks/call_sort_vcf.txt"
    threads: 4
    shell:
        """
        bcftools sort -Oz -o {output.vcf} {input.vcf} 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

localrules: create_db_mapfile


def haplotype_caller_input(wildcards):
    input_type = samples_df.loc[wildcards.sample, "input_type"]
    
    if input_type == "gvcf":
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller")
    
    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


def get_gvcfs_for_db(wc):
    return {
        "gvcfs": [get_final_gvcf(s) for s in SAMPLES_ALL],
        "tbis": [get_final_gvcf(s) + ".tbi" for s in SAMPLES_ALL],
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }


rule gatk_haplotypecaller:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
        min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_haplotypecaller/{sample}.txt"
    log:
        "logs/gatk_haplotypecaller/{sample}.txt"
    shell:
        """
        gatk HaplotypeCaller \
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ploidy {params.ploidy} \
            --emit-ref-confidence GVCF \
            --min-pruning {params.min_pruning} \
            --min-dangling-branch-length {params.min_dangling} \
            &> {log}
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
        db=temp(directory("results/gatk_genomics_db")),
        tar="results/gatk_genomics_db.tar",
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genomics_db_import.txt"
    log:
        "logs/gatk_genomics_db_import.txt"
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '{params.java_opts}' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            -L {input.ref_idx} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            &> {log}
        tar -cf {output.tar} {output.db} &>> {log}
        """


rule gatk_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db.tar",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        het_prior=config["variant_calling"]["gatk"]["het_prior"],
        db=lambda wc, input: subpath(input.db, strip_suffix=".tar"),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genotype_gvcfs.txt"
    log:
        "logs/gatk_genotype_gvcfs.txt"
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
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
        **REF_FILES,
    output:
        vcf="results/vcfs/filtered.vcf.gz",
        tbi="results/vcfs/filtered.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_variant_filtration.txt"
    log:
        "logs/gatk_variant_filtration.txt"
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

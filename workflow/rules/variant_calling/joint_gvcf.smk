localrules: create_db_mapfile


def get_gvcfs_for_db(wc):
    inputs = {
        "gvcfs": get_joint_gvcf_paths(),
        "tbis": get_joint_gvcf_indexes(),
        "db_mapfile": "results/genomics_db/mapfile.txt",
        **REF_FILES,
    }
    if GATK_LONG_CONTIG_MODE:
        inputs.update(
            {
                "archive_gvcfs": get_generated_gvcf_archives(),
                "archive_indexes": get_generated_gvcf_archive_indexes(),
            }
        )
    return inputs


def create_db_mapfile_input(wc):
    inputs = {"gvcfs": get_joint_gvcf_paths()}
    if GATK_LONG_CONTIG_MODE:
        inputs.update(
            {
                "gvcf_indexes": get_joint_gvcf_indexes(),
                "archive_gvcfs": get_generated_gvcf_archives(),
                "archive_indexes": get_generated_gvcf_archive_indexes(),
            }
        )
    return inputs


rule create_db_mapfile:
    input:
        unpack(create_db_mapfile_input),
    output:
        mapfile="results/genomics_db/mapfile.txt",
    run:
        write_joint_gvcf_mapfile(output.mapfile)


rule joint_genomics_db_import:
    input:
        unpack(get_gvcfs_for_db),
    output:
        db=temp(directory("results/gatk_genomics_db")),
        tar="results/gatk_genomics_db.tar",
    params:
        interval_tools=INTERVAL_LIST_TOOLS,
        merge_contig_threshold=GENOMICSDB_MERGE_CONTIG_THRESHOLD,
    threads: 1
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/joint_genomics_db_import.txt"
    log:
        "logs/joint_genomics_db_import.txt"
    shell:
        """
        : > {log}
        export TILEDB_DISABLE_FILE_LOCKING=1
        MERGE_CONTIGS_ARG=$(python {params.interval_tools} genomicsdb-merge-contigs-arg \
            --input {input.ref_fai} \
            --threshold {params.merge_contig_threshold} 2>> {log})
        gatk GenomicsDBImport \
            --java-options '-Xmx{resources.mem_mb_reduced}m' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals $MERGE_CONTIGS_ARG \
            --reader-threads {threads} \
            -L {input.ref_fai} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            &>> {log}
        tar -cf {output.tar} {output.db} >> {log} 2>&1
        """


if LONG_CONTIG_MODE:

    rule joint_genotype_gvcfs:
        input:
            db="results/gatk_genomics_db.tar",
            **REF_FILES,
        output:
            vcf=temp(RAW_VCF_WORK),
            idx=temp(RAW_VCF_WORK_INDEX),
        params:
            het_prior=config["variant_calling"]["gatk"]["het_prior"],
            db_rel=lambda wc, input: subpath(input.db, strip_suffix=".tar"),
        conda:
            "../../envs/gatk.yaml"
        benchmark:
            "benchmarks/joint_genotype_gvcfs.txt"
        log:
            "logs/joint_genotype_gvcfs.txt"
        shell:
            """
            EXTRACT_DIR=$(mktemp -d {resources.tmpdir}/joint_genotype_gvcfs.XXXXXX)
            trap 'rm -rf "$EXTRACT_DIR"' EXIT
            tar -xf {input.db} -C "$EXTRACT_DIR"
            gatk GenotypeGVCFs \
                --java-options '-Xmx{resources.mem_mb_reduced}m' \
                -R {input.ref} \
                --heterozygosity {params.het_prior} \
                --genomicsdb-shared-posixfs-optimizations true \
                -V gendb://"$EXTRACT_DIR/{params.db_rel}" \
                -O {output.vcf} \
                --tmp-dir {resources.tmpdir} \
                &> {log}
            """


    rule compress_joint_raw_vcf:
        input:
            vcf=RAW_VCF_WORK,
            idx=RAW_VCF_WORK_INDEX,
        output:
            vcf=temp(RAW_VCF),
            idx=temp(RAW_VCF_INDEX),
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/bcftools.yaml"
        benchmark:
            "benchmarks/compress_joint_raw_vcf.txt"
        log:
            "logs/compress_joint_raw_vcf.txt"
        shell:
            """
            bcftools view -Oz -o {output.vcf} {input.vcf} 2> {log}
            bcftools index {params.index_args} {output.vcf} 2>> {log}
            """

else:

    rule joint_genotype_gvcfs:
        input:
            db="results/gatk_genomics_db.tar",
            **REF_FILES,
        output:
            vcf=temp(RAW_VCF),
            tbi=temp(RAW_VCF_INDEX),
        params:
            het_prior=config["variant_calling"]["gatk"]["het_prior"],
            db_rel=lambda wc, input: subpath(input.db, strip_suffix=".tar"),
        conda:
            "../../envs/gatk.yaml"
        benchmark:
            "benchmarks/joint_genotype_gvcfs.txt"
        log:
            "logs/joint_genotype_gvcfs.txt"
        shell:
            """
            EXTRACT_DIR=$(mktemp -d {resources.tmpdir}/joint_genotype_gvcfs.XXXXXX)
            trap 'rm -rf "$EXTRACT_DIR"' EXIT
            tar -xf {input.db} -C "$EXTRACT_DIR"
            gatk GenotypeGVCFs \
                --java-options '-Xmx{resources.mem_mb_reduced}m' \
                -R {input.ref} \
                --heterozygosity {params.het_prior} \
                --genomicsdb-shared-posixfs-optimizations true \
                -V gendb://"$EXTRACT_DIR/{params.db_rel}" \
                -O {output.vcf} \
                --tmp-dir {resources.tmpdir} \
                &> {log}
            """

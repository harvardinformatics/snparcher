localrules: create_db_mapfile


def _java_opts_from_resources(resources, default_mem_mb=4096):
    mem_mb = getattr(resources, "mem_mb", default_mem_mb)
    try:
        mem_mb = int(float(mem_mb))
    except (TypeError, ValueError):
        mem_mb = default_mem_mb
    return f"-Xmx{int(mem_mb * 0.9)}m"


def get_gvcfs_for_db(wc):
    return {
        "gvcfs": get_joint_gvcf_paths(),
        "tbis": get_joint_gvcf_tbis(),
        "db_mapfile": "results/genomics_db/mapfile.txt",
        **REF_FILES,
    }


rule create_db_mapfile:
    input:
        gvcfs=get_joint_gvcf_paths(),
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
        java_opts=lambda wildcards, resources: _java_opts_from_resources(resources),
    threads: 1
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/joint_genomics_db_import.txt"
    log:
        "logs/joint_genomics_db_import.txt"
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '{params.java_opts}' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            --reader-threads {threads} \
            -L {input.ref_fai} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            &> {log}
        tar -cf {output.tar} {output.db} >> {log} 2>&1
        """


rule joint_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db.tar",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        java_opts=lambda wildcards, resources: _java_opts_from_resources(resources),
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
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            --heterozygosity {params.het_prior} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://"$EXTRACT_DIR/{params.db_rel}" \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} \
            &> {log}
        """

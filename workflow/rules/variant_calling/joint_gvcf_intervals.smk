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
    """Get unfiltered interval VCF files."""
    intervals = get_db_intervals(wc)
    return expand(
        "results/vcfs/intervals/L{interval}.vcf.gz",
        interval=intervals,
    )


def get_interval_vcf_tbis(wc):
    """Get unfiltered interval VCF index files."""
    intervals = get_db_intervals(wc)
    return expand(
        "results/vcfs/intervals/L{interval}.vcf.gz.tbi",
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


def get_concat_batch_size():
    """Get batch size for staged bcftools concat operations."""
    size = int(config["variant_calling"]["gatk"]["concat_batch_size"])
    if size < 2:
        raise ValueError("variant_calling.gatk.concat_batch_size must be >= 2")
    return size


def get_concat_max_rounds():
    """Get max allowed rounds for staged concat operations."""
    rounds = int(config["variant_calling"]["gatk"]["concat_max_rounds"])
    if rounds < 1:
        raise ValueError("variant_calling.gatk.concat_max_rounds must be >= 1")
    return rounds


def _ceil_div(a, b):
    return (a + b - 1) // b


def get_stage_chunk_counts(num_files):
    """Compute chunk counts for each staged concat round."""
    if num_files < 1:
        raise ValueError("Staged concat requires at least one input file")

    batch_size = get_concat_batch_size()
    max_rounds = get_concat_max_rounds()

    counts = []
    current = num_files
    while current > 1:
        if len(counts) >= max_rounds:
            raise ValueError(
                f"Staged concat exceeded variant_calling.gatk.concat_max_rounds={max_rounds}. "
                "Increase concat_batch_size or concat_max_rounds."
            )
        n_chunks = _ceil_div(current, batch_size)
        counts.append(n_chunks)
        current = n_chunks
    return counts


def staged_vcf_path(wc, round_idx, chunk_idx):
    return f"results/vcfs/staged/r{round_idx}/c{chunk_idx}.vcf.gz"


def get_stage_inputs(base_files, round_idx, chunk_idx, wc, path_builder):
    """Return input files for one staged concat chunk."""
    stage_counts = get_stage_chunk_counts(len(base_files))
    if not stage_counts:
        raise ValueError("Requested staged concat inputs with only one base file")

    if round_idx < 1 or round_idx > len(stage_counts):
        raise ValueError(
            f"Invalid staged concat round {round_idx}. Valid rounds are 1..{len(stage_counts)}"
        )

    n_chunks = stage_counts[round_idx - 1]
    if chunk_idx < 0 or chunk_idx >= n_chunks:
        raise ValueError(
            f"Invalid staged concat chunk {chunk_idx} for round {round_idx}. "
            f"Valid chunks are 0..{n_chunks - 1}"
        )

    if round_idx == 1:
        round_inputs = list(base_files)
    else:
        prev_chunks = stage_counts[round_idx - 2]
        round_inputs = [path_builder(wc, round_idx - 1, i) for i in range(prev_chunks)]

    batch_size = get_concat_batch_size()
    start = chunk_idx * batch_size
    end = min(start + batch_size, len(round_inputs))
    selected = round_inputs[start:end]

    if not selected:
        raise ValueError(
            f"No files selected for staged concat round={round_idx}, chunk={chunk_idx}"
        )
    return selected


def get_final_stage_file(base_files, wc, path_builder):
    """Return final staged output path (or original file if no staging is needed)."""
    if not base_files:
        raise ValueError("No files available for concat")

    stage_counts = get_stage_chunk_counts(len(base_files))
    if not stage_counts:
        return base_files[0]

    final_round = len(stage_counts)
    return path_builder(wc, final_round, 0)


def get_interval_vcf_stage_inputs(wc):
    return get_stage_inputs(
        base_files=get_interval_vcfs(wc),
        round_idx=int(wc.round),
        chunk_idx=int(wc.chunk),
        wc=wc,
        path_builder=staged_vcf_path,
    )


def get_interval_vcf_stage_tbis(wc):
    stage_inputs = get_interval_vcf_stage_inputs(wc)
    return [f"{vcf}.tbi" for vcf in stage_inputs]


def get_final_interval_vcf_stage_file(wc):
    return get_final_stage_file(
        base_files=get_interval_vcfs(wc),
        wc=wc,
        path_builder=staged_vcf_path,
    )


def get_final_interval_vcf_stage_tbi(wc):
    return get_final_interval_vcf_stage_file(wc) + ".tbi"


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
    threads: 1
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
            --reader-threads {threads} \
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


rule concat_interval_vcfs:
    input:
        vcf=get_final_interval_vcf_stage_file,
        tbi=get_final_interval_vcf_stage_tbi,
    output:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
    benchmark:
        "benchmarks/concat_interval_vcfs/benchmark.txt"
    log:
        "logs/concat_interval_vcfs/log.txt"
    shell:
        """
        mv {input.vcf} {output.vcf} 2> {log}
        mv {input.tbi} {output.tbi} 2>> {log}
        """


rule concat_interval_vcfs_stage:
    input:
        vcfs=get_interval_vcf_stage_inputs,
        tbis=get_interval_vcf_stage_tbis,
    output:
        vcf=temp("results/vcfs/staged/r{round}/c{chunk}.vcf.gz"),
        tbi=temp("results/vcfs/staged/r{round}/c{chunk}.vcf.gz.tbi"),
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    log:
        "logs/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

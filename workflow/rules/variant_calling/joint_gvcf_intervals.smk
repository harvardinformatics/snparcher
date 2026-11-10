localrules: create_db_mapfile

import shlex

wildcard_constraints:
    round=r"\d+",
    chunk=r"\d+",


INTERVAL_VCF_PATTERN = (
    "results/vcfs/intervals/L{interval}.vcf"
    if LONG_CONTIG_MODE
    else "results/vcfs/intervals/L{interval}.vcf.gz"
)


STAGED_VCF_PATTERN = (
    "results/vcfs/work/staged/r{round}/c{chunk}.vcf"
    if LONG_CONTIG_MODE
    else "results/vcfs/staged/r{round}/c{chunk}.vcf.gz"
)


def get_db_intervals(wc):
    """Get DB interval IDs."""
    return interval_list_ids("db", wc)


def get_interval_vcfs(wc):
    """Get unfiltered interval VCF files."""
    intervals = get_db_intervals(wc)
    return expand(INTERVAL_VCF_PATTERN, interval=intervals)


def get_gvcfs_for_db(wc):
    inputs = {
        "gvcfs": get_joint_gvcf_paths(),
        "tbis": get_joint_gvcf_indexes(),
        "interval": f"results/intervals/db/{wc.interval}-scattered.interval_list",
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }
    if GATK_LONG_CONTIG_MODE:
        inputs.update(
            {
                "archive_gvcfs": get_generated_gvcf_archives(),
                "archive_indexes": get_generated_gvcf_archive_indexes(),
            }
        )
    return inputs


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
    return STAGED_VCF_PATTERN.format(round=round_idx, chunk=chunk_idx)


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


def format_picard_vcf_inputs(files):
    """Format VCF inputs as repeated Picard SortVcf I= arguments."""
    return " ".join(f"I={shlex.quote(str(path))}" for path in files)


# Long-contig mode gathers interval VCFs with picard SortVcf, not a streaming
# gather: GenotypeGVCFs writes a multi-contig DB shard in contig-name order,
# not reference order. Giving GenotypeGVCFs the shard's intervals (-L, with
# --merge-input-intervals --only-output-calls-starting-in-intervals, as GATK's
# WARP pipeline does) writes reference order, which would allow picard
# GatherVcfs instead.
#
# SortVcf holds MAX_RECORDS_IN_RAM records on the heap before spilling to
# TMP_DIR, and a record grows with the sample count: about 110 bytes per sample
# for GATK genotype fields, measured. Picard's default of 500,000 records needs
# roughly 65 MB of heap per sample, so the gather rules cap the batch at half
# the heap, using this conservative per-record estimate. The heap is only known
# when the job runs, so the division happens in the shell.
SORTVCF_RECORD_BYTES = 2048


SORTVCF_SAMPLE_BYTES = 250


def sortvcf_record_bytes(n_samples):
    """Estimate the heap one SortVcf record with n_samples occupies."""
    return SORTVCF_RECORD_BYTES + SORTVCF_SAMPLE_BYTES * n_samples


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
    return [get_vcf_index(vcf) for vcf in stage_inputs]


def get_interval_vcf_stage_picard_inputs(wc):
    return format_picard_vcf_inputs(get_interval_vcf_stage_inputs(wc))


def get_final_interval_vcf_stage_file(wc):
    return get_final_stage_file(
        base_files=get_interval_vcfs(wc),
        wc=wc,
        path_builder=staged_vcf_path,
    )


def get_final_interval_vcf_stage_tbi(wc):
    return get_vcf_index(get_final_interval_vcf_stage_file(wc))


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


rule gatk_genomics_db_import:
    input:
        unpack(get_gvcfs_for_db),
    output:
        db=temp(directory("results/gatk_genomics_db/L{interval}")),
        tar="results/gatk_genomics_db/L{interval}.tar",
    params:
        interval_tools=INTERVAL_LIST_TOOLS,
        merge_contig_threshold=GENOMICSDB_MERGE_CONTIG_THRESHOLD,
    threads: 1
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genomics_db_import/{interval}.txt"
    log:
        "logs/gatk_genomics_db_import/{interval}.txt"
    shell:
        """
        : > {log}
        export TILEDB_DISABLE_FILE_LOCKING=1
        MERGE_CONTIGS_ARG=$(python {params.interval_tools} genomicsdb-merge-contigs-arg \
            --input {input.interval} \
            --threshold {params.merge_contig_threshold} 2>> {log})
        gatk GenomicsDBImport \
            --java-options '-Xmx{resources.mem_mb_reduced}m' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals $MERGE_CONTIGS_ARG \
            --reader-threads {threads} \
            -L {input.interval} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            >> {log} 2>&1
        tar -cf {output.tar} {output.db} >> {log} 2>&1
        """


rule gatk_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db/L{interval}.tar",
        **REF_FILES,
    output:
        vcf=temp(INTERVAL_VCF_PATTERN),
        idx=temp(get_vcf_index(INTERVAL_VCF_PATTERN)),
    params:
        het_prior=config["variant_calling"]["gatk"]["het_prior"],
        db_rel=subpath(input.db, strip_suffix=".tar"),
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genotype_gvcfs/{interval}.txt"
    log:
        "logs/gatk_genotype_gvcfs/{interval}.txt"
    shell:
        """
        EXTRACT_DIR=$(mktemp -d {resources.tmpdir}/gatk_genotype_gvcfs.{wildcards.interval}.XXXXXX)
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


rule concat_interval_vcfs_stage:
    input:
        vcfs=get_interval_vcf_stage_inputs,
        tbis=get_interval_vcf_stage_tbis,
    output:
        vcf=temp(STAGED_VCF_PATTERN),
        idx=temp(get_vcf_index(STAGED_VCF_PATTERN)),
    params:
        index_args=BCFTOOLS_INDEX_ARGS,
        picard_inputs=get_interval_vcf_stage_picard_inputs,
        sortvcf_record_bytes=sortvcf_record_bytes(len(SAMPLES_ALL)),
        long_mode=LONG_CONTIG_MODE,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    log:
        "logs/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    shell:
        """
        if [ "{params.long_mode}" = "True" ]; then
            # Keep half the heap's worth of records in RAM, at most Picard's
            # default of 500,000 (see SORTVCF_RECORD_BYTES).
            max_records=$(( {resources.mem_mb_reduced} * 1048576 / 2 / {params.sortvcf_record_bytes} ))
            max_records=$(( max_records > 500000 ? 500000 : (max_records < 1 ? 1 : max_records) ))
            picard SortVcf \
                -Xmx{resources.mem_mb_reduced}m \
                {params.picard_inputs} \
                O={output.vcf} \
                TMP_DIR={resources.tmpdir} \
                MAX_RECORDS_IN_RAM=$max_records \
                CREATE_INDEX=false \
                > {log} 2>&1
            gatk IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        else
            bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
                | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} - 2>> {log}
            bcftools index {params.index_args} {output.vcf} 2>> {log}
        fi
        """


rule concat_interval_vcfs:
    input:
        vcf=get_final_interval_vcf_stage_file,
        tbi=get_final_interval_vcf_stage_tbi,
    output:
        vcf=temp(RAW_VCF_WORK) if LONG_CONTIG_MODE else raw_vcf_output(RAW_VCF),
        idx=temp(RAW_VCF_WORK_INDEX) if LONG_CONTIG_MODE else raw_vcf_output(RAW_VCF_INDEX),
    benchmark:
        "benchmarks/concat_interval_vcfs/benchmark.txt"
    log:
        "logs/concat_interval_vcfs/log.txt"
    shell:
        """
        mv {input.vcf} {output.vcf} 2> {log}
        mv {input.tbi} {output.idx} 2>> {log}
        touch {output.vcf} {output.idx}
        """


if LONG_CONTIG_MODE:

    rule compress_interval_raw_vcf:
        input:
            vcf=RAW_VCF_WORK,
            idx=RAW_VCF_WORK_INDEX,
        output:
            vcf=raw_vcf_output(RAW_VCF),
            idx=raw_vcf_output(RAW_VCF_INDEX),
        params:
            index_args=BCFTOOLS_INDEX_ARGS,
        conda:
            "../../envs/bcftools.yaml"
        benchmark:
            "benchmarks/compress_interval_raw_vcf.txt"
        log:
            "logs/compress_interval_raw_vcf.txt"
        shell:
            """
            bcftools view -Oz -o {output.vcf} {input.vcf} 2> {log}
            bcftools index {params.index_args} {output.vcf} 2>> {log}
            """

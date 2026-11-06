def get_db_interval_count(wildcards):
    """Calculate number of DB intervals based on sample count and config."""
    num_samples = len(SAMPLES_ALL)
    scatter_factor = config["intervals"]["db_scatter_factor"]
    num_gvcf_intervals = config["intervals"]["num_gvcf_intervals"]
    return max(int(scatter_factor * num_samples * num_gvcf_intervals), 1)


INTERVAL_LIST_TOOLS = str(Path(workflow.basedir, "scripts/interval_list_tools.py"))


def interval_settings(kind):
    """Settings that determine the interval lists in results/intervals/{kind}."""
    cfg = config["intervals"]
    settings = {
        "reference_name": REF_NAME,
        "reference_source": config["reference"]["source"],
        "min_nmer": cfg["min_nmer"],
        "min_contig_length": cfg["min_contig_length"],
    }
    if kind == "gvcf":
        settings["num_gvcf_intervals"] = cfg["num_gvcf_intervals"]
    else:
        settings.update(
            scatter_count=get_db_interval_count(None),
            max_intervals_per_shard=cfg["db_max_intervals_per_shard"],
            max_contigs_per_shard=cfg["db_max_contigs_per_shard"],
            merge_contigs_threshold=GENOMICSDB_MERGE_CONTIG_THRESHOLD,
        )
    return [f"{key}={value}" for key, value in settings.items()]


def reference_contigs_sha1(fai):
    """SHA-1 of the contig names and lengths, as `cut -f1,2 <fai> | sha1sum` gives."""
    with open(fai) as handle:
        text = "".join("\t".join(line.split("\t")[:2]) + "\n" for line in handle)
    return hashlib.sha1(text.encode()).hexdigest()


def changed_interval_settings(kind):
    """Return "key: old -> new" for each setting that differs from the recorded one."""
    path = Path(f"results/intervals/{kind}/settings.txt")
    if not path.exists():
        # Not built yet, or built before snpArcher 2.3, which recorded nothing.
        return []
    recorded = dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)
    current = dict(item.split("=", 1) for item in interval_settings(kind))
    fai = Path(REF_FILES["ref_fai"])
    if fai.exists():
        current["reference_contigs_sha1"] = reference_contigs_sha1(fai)
    else:
        recorded.pop("reference_contigs_sha1", None)
    return [
        f"{key}: {recorded.get(key)} -> {current.get(key)}"
        for key in sorted(recorded.keys() | current.keys())
        if recorded.get(key) != current.get(key)
    ]


@cache
def check_interval_settings():
    """Refuse to reuse interval lists built with other settings or another reference.

    The input functions read existing lists straight from disk (see
    get_interval_gvcfs), so Snakemake never reruns the checkpoint that wrote
    them when a setting changes.
    """
    problems = []
    for kind in ("gvcf", "db"):
        changed = changed_interval_settings(kind)
        if changed:
            problems.append(f"results/intervals/{kind}:\n    " + "\n    ".join(changed))
    if not problems:
        return
    message = (
        "These interval lists were built with different settings, and snpArcher "
        "would reuse them unchanged:\n  " + "\n  ".join(problems) + "\n"
    )
    if any("scatter_count" in problem for problem in problems):
        message += (
            "scatter_count is db_scatter_factor x samples x num_gvcf_intervals, so it "
            "also changes when samples are added or removed.\n"
        )
    raise WorkflowError(
        message + "Restore the previous settings, or delete the listed directories to "
        "rebuild them (every job that uses those intervals then reruns)."
    )


rule picard_intervals:
    input:
        **REF_FILES,
    output:
        intervals=temp("results/intervals/picard.interval_list"),
    conda:
        "../envs/gatk.yaml"
    params:
        min_nmer=config["intervals"]["min_nmer"],
    benchmark:
        "benchmarks/intervals_picard.txt"
    log:
        "logs/intervals_picard.txt",
    shell:
        """
        picard ScatterIntervalsByNs \
            -Xmx{resources.mem_mb_reduced}m \
            REFERENCE={input.ref} \
            OUTPUT={output.intervals} \
            MAX_TO_MERGE={params.min_nmer} \
            OUTPUT_TYPE=ACGT \
            &> {log}
        """


rule filter_picard_intervals:
    input:
        intervals="results/intervals/picard.interval_list",
    output:
        intervals="results/intervals/filtered.interval_list",
    params:
        min_contig_length=config["intervals"]["min_contig_length"],
        interval_tools=INTERVAL_LIST_TOOLS,
    benchmark:
        "benchmarks/intervals_filter.txt"
    log:
        "logs/intervals_filter.txt",
    shell:
        """
        python {params.interval_tools} filter \
            --input {input.intervals} \
            --output {output.intervals} \
            --min-contig-length {params.min_contig_length} \
            &> {log}
        """


checkpoint create_gvcf_intervals:
    input:
        **REF_FILES,
        intervals="results/intervals/filtered.interval_list",
    output:
        fof="results/intervals/gvcf/intervals.txt",
        out_dir=directory("results/intervals/gvcf"),
    params:
        scatter=config["intervals"]["num_gvcf_intervals"],
        settings=interval_settings("gvcf"),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/intervals_gvcf.txt"
    log:
        "logs/intervals_gvcf.txt",
    shell:
        """
        gatk SplitIntervals \
            --java-options '-Xmx{resources.mem_mb_reduced}m' \
            -R {input.ref} \
            -L {input.intervals} \
            -O {output.out_dir} \
            --scatter-count {params.scatter} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        ls -1 {output.out_dir}/*-scattered.interval_list > {output.fof}
        # Record what the lists were built from; see check_interval_settings.
        {{
            printf 'reference_contigs_sha1=%s\\n' "$(cut -f1,2 {input.ref_fai} | sha1sum | cut -d' ' -f1)"
            printf '%s\\n' {params.settings:q}
        }} > {output.out_dir}/settings.txt
        """


checkpoint create_db_intervals:
    input:
        **REF_FILES,
        intervals="results/intervals/filtered.interval_list",
    output:
        fof="results/intervals/db/intervals.txt",
        out_dir=directory("results/intervals/db"),
    params:
        scatter=get_db_interval_count,
        interval_tools=INTERVAL_LIST_TOOLS,
        max_intervals=config["intervals"]["db_max_intervals_per_shard"],
        max_contigs=config["intervals"]["db_max_contigs_per_shard"],
        merge_contig_threshold=GENOMICSDB_MERGE_CONTIG_THRESHOLD,
        settings=interval_settings("db"),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/intervals_db.txt"
    log:
        "logs/intervals_db.txt",
    shell:
        """
        mkdir -p {output.out_dir}
        raw_dir=$(mktemp -d {output.out_dir}/gatk_split_raw.XXXXXX)
        gatk SplitIntervals \
            --java-options '-Xmx{resources.mem_mb_reduced}m' \
            -R {input.ref} \
            -L {input.intervals} \
            -O "$raw_dir" \
            --scatter-count {params.scatter} \
            --subdivision-mode INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        python {params.interval_tools} split-db \
            --input-dir "$raw_dir" \
            --output-dir {output.out_dir} \
            --fof {output.fof} \
            --max-intervals-per-shard {params.max_intervals} \
            --max-contigs-per-shard {params.max_contigs} \
            --merge-contigs-threshold {params.merge_contig_threshold} \
            >> {log} 2>&1
        # Record what the lists were built from; see check_interval_settings.
        {{
            printf 'reference_contigs_sha1=%s\\n' "$(cut -f1,2 {input.ref_fai} | sha1sum | cut -d' ' -f1)"
            printf '%s\\n' {params.settings:q}
        }} > {output.out_dir}/settings.txt
        """

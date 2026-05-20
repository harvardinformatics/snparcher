def get_db_interval_count(wildcards):
    """Calculate number of DB intervals based on sample count and config."""
    num_samples = len(SAMPLES_ALL)
    scatter_factor = config["intervals"]["db_scatter_factor"]
    num_gvcf_intervals = config["intervals"]["num_gvcf_intervals"]
    return max(int(scatter_factor * num_samples * num_gvcf_intervals), 1)


INTERVAL_LIST_TOOLS = str(Path(workflow.basedir, "scripts/interval_list_tools.py"))


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
        """

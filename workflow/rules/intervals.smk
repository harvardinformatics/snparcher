def get_db_interval_count(wildcards):
    """Calculate number of DB intervals based on sample count and config."""
    num_samples = len(samples_df)
    scatter_factor = config["intervals"]["db_scatter_factor"]
    num_gvcf_intervals = config["intervals"]["num_gvcf_intervals"]
    return max(int(scatter_factor * num_samples * num_gvcf_intervals), 1)


rule picard_intervals:
    input:
        **REF_FILES,
    output:
        intervals=temp("results/intervals/picard.interval_list"),
    conda:
        "../envs/gatk.yaml"
    params:
        min_nmer=config["intervals"]["min_nmer"],
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
    resources:
        mem_mb=4096,
    benchmark:
        "benchmarks/intervals_picard.txt"
    log:
        "logs/intervals_picard.txt"
    shell:
        """
        picard ScatterIntervalsByNs \
            {params.java_opts} \
            REFERENCE={input.ref} \
            OUTPUT={output.intervals} \
            MAX_TO_MERGE={params.min_nmer} \
            OUTPUT_TYPE=ACGT \
            &> {log}
        """

checkpoint create_gvcf_intervals:
    input:
        intervals="results/intervals/picard.interval_list",
        **REF_FILES,
    output:
        fof="results/intervals/gvcf/intervals.txt",
        out_dir=directory("results/intervals/gvcf"),
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        scatter=config["intervals"]["num_gvcf_intervals"],
    resources:
        mem_mb=4096,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/intervals_gvcf.txt"
    log:
        "logs/intervals_gvcf.txt"
    shell:
        """
        gatk SplitIntervals \
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            -L {input.intervals} \
            -O {output.out_dir} \
            --scatter-count {params.scatter} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        ls {output.out_dir}/*-scattered.interval_list > {output.fof}
        """


checkpoint create_db_intervals:
    input:
        intervals="results/intervals/picard.interval_list",
        **REF_FILES,
    output:
        fof="results/intervals/db/intervals.txt",
        out_dir=directory("results/intervals/db"),
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        scatter=get_db_interval_count,
    resources:
        mem_mb=4096,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/intervals_db.txt"
    log:
        "logs/intervals_db.txt"
    shell:
        """
        gatk SplitIntervals \
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            -L {input.intervals} \
            -O {output.out_dir} \
            --scatter-count {params.scatter} \
            --subdivision-mode INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        ls {output.out_dir}/*-scattered.interval_list > {output.fof}
        """
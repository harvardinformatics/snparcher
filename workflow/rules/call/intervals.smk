"""
Interval generation rules for snpArcher v2.

Creates interval files for parallelizing GVCF generation and GenomicsDB import.
"""


def _get_db_interval_count(wildcards=None):
    """Calculate number of GenomicsDB intervals based on sample count."""
    num_samples = len(SAMPLE_IDS)
    num_intervals = max(
        int(config["intervals"]["db_scatter_factor"] * num_samples * config["intervals"]["num_gvcf_intervals"]),
        1,
    )
    return num_intervals


rule intervals_picard:
    """Generate interval list by splitting reference on N stretches."""
    input:
        unpack(get_ref_bundle),
    output:
        intervals=temp("results/intervals/picard_interval_list.list"),
    params:
        min_nmer=int(config["intervals"]["min_nmer"]),
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/intervals_picard.txt"
    benchmark:
        "benchmarks/intervals_picard.txt"
    resources:
        mem_mb=8000,
    shell:
        """
        picard ScatterIntervalsByNs \
            -Xmx{params.java_mem}m \
            REFERENCE={input.ref} \
            OUTPUT={output.intervals} \
            MAX_TO_MERGE={params.min_nmer} \
            OUTPUT_TYPE=ACGT \
            &> {log}
        """


rule intervals_format:
    """Convert Picard interval list to simple format."""
    input:
        intervals="results/intervals/picard_interval_list.list",
    output:
        intervals="results/intervals/master_interval_list.list",
    run:
        with open(output.intervals, "w") as out:
            with open(input.intervals, "r") as inp:
                for line in inp:
                    if not line.startswith("@"):
                        fields = line.strip().split("\t")
                        chrom, start, end = fields[0], fields[1], fields[2]
                        print(f"{chrom}:{start}-{end}", file=out)


checkpoint intervals_gvcf:
    """Split intervals for parallelizing HaplotypeCaller.

    Uses BALANCING_WITHOUT_INTERVAL_SUBDIVISION mode to keep intervals intact.
    """
    input:
        unpack(get_ref_bundle),
        intervals="results/intervals/master_interval_list.list",
    output:
        fof="results/intervals/gvcf_intervals/intervals.txt",
        out_dir=directory("results/intervals/gvcf_intervals"),
    params:
        max_intervals=config["intervals"]["num_gvcf_intervals"],
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/intervals_gvcf.txt"
    benchmark:
        "benchmarks/intervals_gvcf.txt"
    resources:
        mem_mb=8000,
    shell:
        """
        gatk SplitIntervals \
            --java-options "-Xmx{params.java_mem}m" \
            -L {input.intervals} \
            -O {output.out_dir} \
            -R {input.ref} \
            -scatter {params.max_intervals} \
            -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof} 2>> {log}
        """


checkpoint intervals_db:
    """Split intervals for parallelizing GenomicsDBImport.

    Uses INTERVAL_SUBDIVISION mode and scales with sample count.
    """
    input:
        unpack(get_ref_bundle),
        intervals="results/intervals/master_interval_list.list",
    output:
        fof="results/intervals/db_intervals/intervals.txt",
        out_dir=directory("results/intervals/db_intervals"),
    params:
        max_intervals=_get_db_interval_count,
        java_mem=get_java_mem,
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/intervals_db.txt"
    benchmark:
        "benchmarks/intervals_db.txt"
    resources:
        mem_mb=8000,
    shell:
        """
        gatk SplitIntervals \
            --java-options "-Xmx{params.java_mem}m" \
            -L {input.intervals} \
            -O {output.out_dir} \
            -R {input.ref} \
            -scatter {params.max_intervals} \
            -mode INTERVAL_SUBDIVISION \
            --interval-merging-rule OVERLAPPING_ONLY \
            &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof} 2>> {log}
        """

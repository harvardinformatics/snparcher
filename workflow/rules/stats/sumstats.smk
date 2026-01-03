"""
Summary statistics rules for snpArcher v2.

Collects BAM statistics, fastp output, and generates summary reports.

Rules:
- stats_bam_sumstats: Compute BAM coverage and alignment stats
- stats_sentieon_bam: Compute additional Sentieon-specific BAM stats
- stats_collect_fastp: Aggregate fastp output across units
- stats_collect_sumstats: Generate final summary statistics file
"""


rule stats_bam_sumstats:
    """Compute BAM coverage and alignment summary statistics."""
    input:
        unpack(get_final_bam),
        ref=f"results/reference/{REF_FILE}",
    output:
        cov="results/summary_stats/{sample}_coverage.txt",
        alnSum="results/summary_stats/{sample}_AlnSumMets.txt",
    conda:
        "../../envs/fastq2bam.yml"
    log:
        "logs/stats_bam_sumstats/{sample}.txt"
    benchmark:
        "benchmarks/stats_bam_sumstats/{sample}.txt"
    shell:
        """
        samtools coverage --output {output.cov} {input.bam} 2> {log}
        samtools flagstat -O tsv {input.bam} > {output.alnSum} 2>> {log}
        """


rule stats_sentieon_bam:
    """Compute Sentieon-specific BAM statistics."""
    input:
        unpack(get_final_bam),
        ref=f"results/reference/{REF_FILE}",
        fai=f"results/reference/{REF_FILE}.fai",
    output:
        insert_file="results/summary_stats/{sample}_insert_metrics.txt",
        qd="results/summary_stats/{sample}_qd_metrics.txt",
        gc="results/summary_stats/{sample}_gc_metrics.txt",
        gc_summary="results/summary_stats/{sample}_gc_summary.txt",
        mq="results/summary_stats/{sample}_mq_metrics.txt",
    params:
        lic=config.get("sentieon_lic", ""),
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/stats_sentieon_bam/{sample}.txt"
    benchmark:
        "benchmarks/stats_sentieon_bam/{sample}.txt"
    threads: 4
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} \
            -t {threads} -i {input.bam} \
            --algo MeanQualityByCycle {output.mq} \
            --algo QualDistribution {output.qd} \
            --algo GCBias --summary {output.gc_summary} {output.gc} \
            --algo InsertSizeMetricAlgo {output.insert_file} \
            2> {log}
        """


rule stats_collect_fastp:
    """Aggregate fastp output across all units for a sample."""
    input:
        collect_fastp_stats_input,
    output:
        "results/summary_stats/{sample}_fastp.out",
    run:
        combine_fastp_files(input, output)


rule stats_collect_sumstats:
    """Generate final summary statistics file."""
    input:
        unpack(get_sumstats_input),
    output:
        f"results/summary_stats/{REF_NAME}_bam_sumstats.txt",
    run:
        FractionReadsPassFilter, NumReadsPassFilter = collectFastpOutput(input.fastpFiles)
        aln_metrics = collectAlnSumMets(input.alnSumMetsFiles)
        SeqDepths, CoveredBases = collectCoverageMetrics(input.coverageFiles)
        
        if config.get("sentieon", False):
            median_inserts, median_insert_std = collect_inserts(input.insert_files)
            printBamSumStats(
                SeqDepths,
                CoveredBases,
                aln_metrics,
                FractionReadsPassFilter,
                NumReadsPassFilter,
                output[0],
                median_inserts,
                median_insert_std,
            )
        else:
            printBamSumStats(
                SeqDepths,
                CoveredBases,
                aln_metrics,
                FractionReadsPassFilter,
                NumReadsPassFilter,
                output[0],
            )

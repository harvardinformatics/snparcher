# The sentieon mapping pipeline, used with variant_calling.tool: sentieon:
# fastp, Sentieon bwa mem and per-library Dedup, then the shared merges in
# rules/mapping/common.smk, plus Sentieon's BAM metrics for the QC report.

include: "fastp.smk"


rule sentieon_map:
    input:
        unpack(bwa_mem_input),
    output:
        bam=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam"),
    params:
        rg=get_read_group,
        lic=config["variant_calling"]["sentieon"]["license"],
    threads: 8
    conda:
        "../../envs/sentieon.yaml"
    benchmark:
        "benchmarks/sentieon_map/{sample}/{library}/{input_unit}.txt"
    log:
        "logs/sentieon_map/{sample}/{library}/{input_unit}.txt"
    shell:
        """
        export MALLOC_CONF=lg_dirty_mult:-1
        export SENTIEON_LICENSE={params.lic}
        sentieon bwa mem -M -R {params.rg} -t {threads} -K 10000000 {input.ref} {input.r1} {input.r2} 2> {log} \
            | sentieon util sort --bam_compression 1 -r {input.ref} -o {output.bam} -t {threads} --sam2bam -i - 2>> {log}
        """


rule sentieon_dedup_library:
    input:
        unpack(dedup_library_input),
    output:
        bam=temp("results/bams/library_markdup/{sample}/{library}.bam"),
        score=temp("results/bams/library_markdup/{sample}/{library}_score.txt"),
        metrics=temp("results/bams/library_markdup/{sample}/{library}_metrics.txt"),
    params:
        lic=config["variant_calling"]["sentieon"]["license"],
    threads: 4
    conda:
        "../../envs/sentieon.yaml"
    benchmark:
        "benchmarks/sentieon_dedup/{sample}/{library}.txt"
    log:
        "logs/sentieon_dedup/{sample}/{library}.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -t {threads} -i {input.bam} \
            --algo LocusCollector --fun score_info {output.score} \
            2> {log}
        sentieon driver -t {threads} -i {input.bam} \
            --algo Dedup --score_info {output.score} --metrics {output.metrics} \
            --bam_compression 1 {output.bam} \
            2>> {log}
        rm -f {output.bam}.bai
        """


rule sentieon_bam_stats:
    input:
        bam=lambda wc: get_final_bam(wc.sample),
        **REF_FILES,
    output:
        insert_metrics="results/qc_metrics/sentieon/{sample}_insert_metrics.txt",
        qd="results/qc_metrics/sentieon/{sample}_qd_metrics.txt",
        gc="results/qc_metrics/sentieon/{sample}_gc_metrics.txt",
        gc_summary="results/qc_metrics/sentieon/{sample}_gc_summary.txt",
        mq="results/qc_metrics/sentieon/{sample}_mq_metrics.txt",
    params:
        lic=config["variant_calling"]["sentieon"]["license"],
    threads: 4
    conda:
        "../../envs/sentieon.yaml"
    benchmark:
        "benchmarks/sentieon_bam_stats/{sample}.txt"
    log:
        "logs/sentieon_bam_stats/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver \
            -r {input.ref} \
            -t {threads} \
            -i {input.bam} \
            --algo MeanQualityByCycle {output.mq} \
            --algo QualDistribution {output.qd} \
            --algo GCBias --summary {output.gc_summary} {output.gc} \
            --algo InsertSizeMetricAlgo {output.insert_metrics} \
            2> {log}
        """

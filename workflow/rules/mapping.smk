# mapping.smk

def bwa_mem_input(wildcards):
    """Get input fastqs for alignment."""
    return {
        "r1": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}_1.fastq.gz",
        "r2": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}_2.fastq.gz",
        **REF_FILES,
    }


def get_read_group(wildcards):
    """Generate read group string for BWA."""
    return f"'@RG\\tID:{wildcards.library}\\tSM:{wildcards.sample}\\tLB:{wildcards.library}\\tPL:ILLUMINA'"


def merge_bams_input(wildcards):
    """Get all library BAMs for a sample to merge."""
    libraries = get_sample_libraries(wildcards.sample)
    return {
        "bams": [f"results/bams/raw/{wildcards.sample}/{lib}.bam" for lib in libraries],
    }


def dedup_input(wildcards):
    """Get input BAM for deduplication."""
    if sample_has_multiple_libraries(wildcards.sample):
        bam = f"results/bams/merged/{wildcards.sample}.bam"
    else:
        lib = get_sample_libraries(wildcards.sample)[0]
        bam = f"results/bams/raw/{wildcards.sample}/{lib}.bam"

    return {
        "bam": bam,
        "bai": bam + ".bai",
    }


if USE_SENTIEON:

    rule sentieon_map:
        input:
            unpack(bwa_mem_input),
        output:
            bam=temp("results/bams/raw/{sample}/{library}.bam"),
            bai=temp("results/bams/raw/{sample}/{library}.bam.bai"),
        params:
            rg=get_read_group,
            lic=config["sentieon"]["license"],
        threads: 8
        conda:
            "../envs/sentieon.yaml"
        benchmark:
            "benchmarks/sentieon_map/{sample}/{library}.txt"
        log:
            "logs/sentieon_map/{sample}/{library}.txt"
        shell:
            """
            export MALLOC_CONF=lg_dirty_mult:-1
            export SENTIEON_LICENSE={params.lic}
            sentieon bwa mem -M -R {params.rg} -t {threads} -K 10000000 {input.ref} {input.r1} {input.r2} 2> {log} \
                | sentieon util sort --bam_compression 1 -r {input.ref} -o {output.bam} -t {threads} --sam2bam -i - 2>> {log}
            samtools index {output.bam} 2>> {log}
            """

    rule sentieon_dedup:
        input:
            unpack(dedup_input),
        output:
            bam="results/bams/markdup/{sample}.bam",
            bai="results/bams/markdup/{sample}.bam.bai",
            score=temp("results/bams/markdup/{sample}_score.txt"),
            metrics=temp("results/bams/markdup/{sample}_metrics.txt"),
        params:
            lic=config["sentieon"]["license"],
        threads: 4
        conda:
            "../envs/sentieon.yaml"
        benchmark:
            "benchmarks/sentieon_dedup/{sample}.txt"
        log:
            "logs/sentieon_dedup/{sample}.txt"
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
            """

    rule sentieon_bam_stats:
        input:
            bam=lambda wc: get_final_bam(wc.sample),
            **REF_FILES,
        output:
            insert="results/qc_metrics/sentieon/{sample}_insert_metrics.txt",
            qd="results/qc_metrics/sentieon/{sample}_qd_metrics.txt",
            gc="results/qc_metrics/sentieon/{sample}_gc_metrics.txt",
            gc_summary="results/qc_metrics/sentieon/{sample}_gc_summary.txt",
            mq="results/qc_metrics/sentieon/{sample}_mq_metrics.txt",
        params:
            lic=config["sentieon"]["license"],
        threads: 4
        conda:
            "../envs/sentieon.yaml"
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
                --algo InsertSizeMetricAlgo {output.insert} \
                2> {log}
            """

else:

    rule bwa_mem:
        input:
            unpack(bwa_mem_input),
        output:
            bam=temp("results/bams/raw/{sample}/{library}.bam"),
            bai=temp("results/bams/raw/{sample}/{library}.bam.bai"),
        params:
            rg=get_read_group,
        threads: 8
        conda:
            "../envs/samtools.yaml"
        benchmark:
            "benchmarks/bwa_mem/{sample}/{library}.txt"
        log:
            "logs/bwa_mem/{sample}/{library}.txt"
        shell:
            """
            bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} \
                | samtools sort -o {output.bam} - 2>> {log}
            samtools index {output.bam} 2>> {log}
            """

    rule markdup:
        input:
            unpack(dedup_input),
        output:
            bam="results/bams/markdup/{sample}.bam",
            bai="results/bams/markdup/{sample}.bam.bai",
        threads: 4
        conda:
            "../envs/sambamba.yaml"
        benchmark:
            "benchmarks/markdup/{sample}.txt"
        log:
            "logs/markdup/{sample}.txt"
        shell:
            """
            sambamba markdup -t {threads} {input.bam} {output.bam} 2> {log}
            """


rule merge_bams:
    input:
        unpack(merge_bams_input),
    output:
        bam=temp("results/bams/merged/{sample}.bam"),
        bai=temp("results/bams/merged/{sample}.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_bams/{sample}.txt"
    log:
        "logs/merge_bams/{sample}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule bam_stats:
    input:
        bam=lambda wc: get_final_bam(wc.sample),
    output:
        coverage=temp("results/qc_metrics/bam/{sample}_coverage.txt"),
        flagstat=temp("results/qc_metrics/bam/{sample}_flagstat.txt"),
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/bam_stats/{sample}.txt"
    log:
        "logs/bam_stats/{sample}.txt"
    shell:
        """
        samtools coverage {input.bam} -o {output.coverage} 2> {log}
        samtools flagstat -O tsv {input.bam} > {output.flagstat} 2>> {log}
        """
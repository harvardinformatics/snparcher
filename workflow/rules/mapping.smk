# mapping.smk

def bwa_mem_input(wildcards):
    """Get input fastqs for alignment."""
    return {
        "r1": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}/{wildcards.input_unit}_1.fastq.gz",
        "r2": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}/{wildcards.input_unit}_2.fastq.gz",
        **REF_FILES,
    }


def get_read_group(wildcards):
    """Generate read group string for BWA."""
    return (
        f"'@RG\\tID:{wildcards.library}.{wildcards.input_unit}"
        f"\\tSM:{wildcards.sample}\\tLB:{wildcards.library}\\tPL:ILLUMINA'"
    )


def get_library_rows(sample, library):
    """Get all rows for a sample/library pair."""
    sample_rows = get_sample_rows(sample)
    library_rows = sample_rows[sample_rows["library_id"] == library]
    if library_rows.empty:
        raise ValueError(f"No rows found for sample '{sample}', library '{library}'")
    return library_rows


def merge_library_bams_input(wildcards):
    """Get all per-input BAMs for a library."""
    library_rows = get_library_rows(wildcards.sample, wildcards.library)
    input_units = library_rows["input_unit"].tolist()
    return {
        "bams": [
            f"results/bams/raw/{wildcards.sample}/{wildcards.library}/{input_unit}.bam"
            for input_unit in input_units
        ],
    }


def dedup_library_input(wildcards):
    """Get input BAM for library-level duplicate marking."""
    bam = f"results/bams/library/{wildcards.sample}/{wildcards.library}.bam"
    return {
        "bam": bam,
        "bai": bam + ".bai",
    }


def merge_dedup_libraries_input(wildcards):
    """Get deduplicated library BAMs for sample-level merge."""
    if not get_sample_mark_duplicates(wildcards.sample):
        raise ValueError(
            f"Sample '{wildcards.sample}' has mark_duplicates=False, "
            "but merge_dedup_libraries was requested."
        )
    libraries = get_sample_libraries(wildcards.sample)
    return {
        "bams": [
            f"results/bams/library_markdup/{wildcards.sample}/{lib}.bam"
            for lib in libraries
        ],
    }


def merge_library_level_bams_input(wildcards):
    """Get library BAMs for sample-level merge when duplicate marking is disabled."""
    if get_sample_mark_duplicates(wildcards.sample):
        raise ValueError(
            f"Sample '{wildcards.sample}' has mark_duplicates=True, "
            "but merge_library_level_bams was requested."
        )
    libraries = get_sample_libraries(wildcards.sample)
    return {
        "bams": [f"results/bams/library/{wildcards.sample}/{lib}.bam" for lib in libraries],
    }


if USE_SENTIEON:

    rule sentieon_map:
        input:
            unpack(bwa_mem_input),
        output:
            bam=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam"),
            bai=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam.bai"),
        params:
            rg=get_read_group,
            lic=config["variant_calling"]["sentieon"]["license"],
        threads: 8
        conda:
            "../envs/sentieon.yaml"
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
            samtools index {output.bam} 2>> {log}
            """

    rule sentieon_dedup_library:
        input:
            unpack(dedup_library_input),
        output:
            bam=temp("results/bams/library_markdup/{sample}/{library}.bam"),
            bai=temp("results/bams/library_markdup/{sample}/{library}.bam.bai"),
            score=temp("results/bams/library_markdup/{sample}/{library}_score.txt"),
            metrics=temp("results/bams/library_markdup/{sample}/{library}_metrics.txt"),
        params:
            lic=config["variant_calling"]["sentieon"]["license"],
        threads: 4
        conda:
            "../envs/sentieon.yaml"
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
            lic=config["variant_calling"]["sentieon"]["license"],
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
            bam=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam"),
            bai=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam.bai"),
        params:
            rg=get_read_group,
        threads: 8
        conda:
            "../envs/samtools.yaml"
        benchmark:
            "benchmarks/bwa_mem/{sample}/{library}/{input_unit}.txt"
        log:
            "logs/bwa_mem/{sample}/{library}/{input_unit}.txt"
        shell:
            """
            bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} \
                | samtools sort -o {output.bam} - 2>> {log}
            samtools index {output.bam} 2>> {log}
            """

    rule markdup_library:
        input:
            unpack(dedup_library_input),
        output:
            bam=temp("results/bams/library_markdup/{sample}/{library}.bam"),
            bai=temp("results/bams/library_markdup/{sample}/{library}.bam.bai"),
        threads: 4
        conda:
            "../envs/sambamba.yaml"
        benchmark:
            "benchmarks/markdup/{sample}/{library}.txt"
        log:
            "logs/markdup/{sample}/{library}.txt"
        shell:
            """
            sambamba markdup -t {threads} {input.bam} {output.bam} 2> {log}
            """


rule merge_library_bams:
    input:
        unpack(merge_library_bams_input),
    output:
        bam=temp("results/bams/library/{sample}/{library}.bam"),
        bai=temp("results/bams/library/{sample}/{library}.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_library_bams/{sample}/{library}.txt"
    log:
        "logs/merge_library_bams/{sample}/{library}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule merge_dedup_libraries:
    input:
        unpack(merge_dedup_libraries_input),
    output:
        bam="results/bams/markdup/{sample}.bam",
        bai="results/bams/markdup/{sample}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_dedup_libraries/{sample}.txt"
    log:
        "logs/merge_dedup_libraries/{sample}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule merge_library_level_bams:
    input:
        unpack(merge_library_level_bams_input),
    output:
        bam="results/bams/merged/{sample}.bam",
        bai="results/bams/merged/{sample}.bam.bai",
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

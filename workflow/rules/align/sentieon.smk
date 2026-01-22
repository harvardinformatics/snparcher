"""
Sentieon alignment rules for snpArcher v2.

Requires Sentieon license configured in config['sentieon_lic'].

Rules:
- align_sentieon_map: Map reads with Sentieon BWA
- align_sentieon_merge: Merge BAMs for multi-unit samples
- align_sentieon_dedup: Mark duplicates with Sentieon LocusCollector/Dedup
"""


rule align_sentieon_map:
    """Map paired-end reads using Sentieon BWA."""
    input:
        unpack(get_ref_with_bwa),
        r1="results/filtered_fastqs/{sample}/{unit}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{unit}_2.fastq.gz",
    output:
        bam=temp("results/bams/preMerge/{sample}/{unit}.bam"),
        bai=temp("results/bams/preMerge/{sample}/{unit}.bam.bai"),
    params:
        rg=get_read_group,
        lic=config["variant_calling"]["sentieon"]["license"],
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/align_sentieon_map/{sample}/{unit}.txt",
    benchmark:
        "benchmarks/align_sentieon_map/{sample}_{unit}.txt"
    threads: 8
    shell:
        """
        export MALLOC_CONF=lg_dirty_mult:-1
        export SENTIEON_LICENSE={params.lic}
        sentieon bwa mem -M -R {params.rg} -t {threads} -K 10000000 \
            {input.ref} {input.r1} {input.r2} 2> {log} \
            | sentieon util sort --bam_compression 1 -r {input.ref} \
                -o {output.bam} -t {threads} --sam2bam -i - 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """


rule align_sentieon_merge:
    """Merge multiple BAM files for a single sample."""
    input:
        get_unit_bams,
    output:
        bam=temp("results/bams/merged/{sample}.bam"),
        bai=temp("results/bams/merged/{sample}.bam.bai"),
    conda:
        "../../envs/fastq2bam.yml"
    log:
        "logs/align_sentieon_merge/{sample}.txt",
    benchmark:
        "benchmarks/align_sentieon_merge/{sample}.txt"
    threads: 4
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule align_sentieon_dedup:
    """Mark duplicate reads with Sentieon LocusCollector and Dedup."""
    input:
        unpack(get_dedup_input),
    output:
        bam="results/bams/{sample}.bam",
        bai="results/bams/{sample}.bam.bai",
        score=temp("results/summary_stats/{sample}/sentieon_dedup_score.txt"),
        metrics=temp("results/summary_stats/{sample}/sentieon_dedup_metrics.txt"),
    params:
        lic=config["sentieon_lic"],
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/align_sentieon_dedup/{sample}.txt",
    benchmark:
        "benchmarks/align_sentieon_dedup/{sample}.txt"
    threads: 8
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -t {threads} -i {input.bam} \
            --algo LocusCollector --fun score_info {output.score} &> {log}
        sentieon driver -t {threads} -i {input.bam} \
            --algo Dedup --score_info {output.score} \
            --metrics {output.metrics} --bam_compression 1 {output.bam} &>> {log}
        """

"""
BWA alignment rules for snpArcher v2.

Rules:
- align_bwa_map: Map reads to reference with BWA-MEM
- align_merge_bams: Merge multiple BAMs for same sample
- align_dedup: Mark duplicates with sambamba
"""


rule align_bwa_map:
    """Map paired-end reads to reference genome using BWA-MEM."""
    input:
        ref=f"results/reference/{REF_FILE}",
        r1="results/filtered_fastqs/{sample}/{unit}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{unit}_2.fastq.gz",
        indexes=expand(
            f"results/reference/{REF_FILE}.{{ext}}",
            ext=["sa", "pac", "bwt", "ann", "amb", "fai"],
        ),
    output:
        bam=temp("results/bams/preMerge/{sample}/{unit}.bam"),
        bai=temp("results/bams/preMerge/{sample}/{unit}.bam.bai"),
    params:
        rg=get_read_group,
    conda:
        "../../envs/fastq2bam.yml"
    log:
        "logs/align_bwa_map/{sample}/{unit}.txt"
    benchmark:
        "benchmarks/align_bwa_map/{sample}_{unit}.txt"
    threads: 8
    shell:
        """
        bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} \
            | samtools sort -o {output.bam} - \
            && samtools index {output.bam} {output.bai}
        """


rule align_merge_bams:
    """Merge multiple BAM files for a single sample."""
    input:
        get_unit_bams,
    output:
        bam=temp("results/bams/merged/{sample}.bam"),
        bai=temp("results/bams/merged/{sample}.bam.bai"),
    conda:
        "../../envs/fastq2bam.yml"
    log:
        "logs/align_merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/align_merge_bams/{sample}.txt"
    threads: 4
    shell:
        """
        samtools merge -@ {threads} {output.bam} {input} 2> {log} \
            && samtools index {output.bam}
        """


rule align_dedup:
    """Mark duplicate reads with sambamba."""
    input:
        unpack(get_dedup_input),
    output:
        bam="results/bams/{sample}.bam",
        bai="results/bams/{sample}.bam.bai",
    conda:
        "../../envs/sambamba.yml"
    log:
        "logs/align_dedup/{sample}.txt"
    benchmark:
        "benchmarks/align_dedup/{sample}.txt"
    threads: 4
    shell:
        """
        sambamba markdup -t {threads} {input.bam} {output.bam} 2> {log}
        """

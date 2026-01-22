"""
FASTQ processing rules for snpArcher v2.

Rules:
- fastq_download_srr: Download FASTQ files from SRA
- fastq_sort_reads: Sort reads by name (optional)
- fastq_fastp: Quality filtering and adapter trimming
"""


rule fastq_download_srr:
    """Download FASTQ files from SRA using prefetch/fasterq-dump or ENA."""
    output:
        r1=temp("results/fastq/{sample}/{srr}_1.fastq.gz"),
        r2=temp("results/fastq/{sample}/{srr}_2.fastq.gz"),
    params:
        outdir="results/fastq/{sample}/",
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/fastq_download_srr/{sample}/{srr}.txt"
    benchmark:
        "benchmarks/fastq_download_srr/{sample}_{srr}.txt"
    shell:
        """
        set +e
        # Delete existing prefetch file in case of previous run failure
        rm -rf {wildcards.srr}
        
        # Attempt to get SRA file from NCBI (prefetch) or ENA (wget)
        prefetch --max-size 1T {wildcards.srr} &> {log}
        prefetchExit=$?
        
        if [[ $prefetchExit -ne 0 ]]; then
            ffq --ftp {wildcards.srr} \
                | grep -Eo '"url": "[^"]*"' \
                | grep -o '"[^"]*"$' \
                | grep "fastq" \
                | xargs curl --remote-name-all --output-dir {params.outdir} &>> {log}
        else
            fasterq-dump {wildcards.srr} -O {params.outdir} -e {threads} -t {resources.tmpdir} &>> {log}
            pigz -p {threads} {params.outdir}{wildcards.srr}*.fastq &>> {log}
        fi
        
        rm -rf {wildcards.srr}
        """


rule fastq_sort_reads:
    """Sort reads by name for improved compression and downstream processing."""
    input:
        unpack(get_reads),
    output:
        r1=temp("results/sorted_reads/{sample}/{unit}_1.fastq.gz"),
        r2=temp("results/sorted_reads/{sample}/{unit}_2.fastq.gz"),
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/fastq_sort_reads/{sample}/{unit}.txt"
    benchmark:
        "benchmarks/fastq_sort_reads/{sample}_{unit}.txt"
    shell:
        """
        sortbyname.sh in={input.r1} out={output.r1} &> {log}
        sortbyname.sh in={input.r2} out={output.r2} &>> {log}
        """


rule fastq_fastp:
    """Quality filter and trim adapters from paired-end reads."""
    input:
        unpack(get_reads_fastp),
    output:
        r1="results/filtered_fastqs/{sample}/{unit}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{unit}_2.fastq.gz",
        json="results/fastp/{sample}/{unit}.fastp.json",
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/fastq_fastp/{sample}/{unit}.txt"
    benchmark:
        "benchmarks/fastq_fastp/{sample}_{unit}.txt"
    threads: 4
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            -j {output.json} -h /dev/null \
            &> {log}
        """


def _get_fastp_jsons(wildcards):
    """Get all fastp JSON files for a sample's units."""
    units = get_sample_units(wildcards.sample)
    return expand(
        "results/fastp/{sample}/{unit}.fastp.json",
        sample=wildcards.sample,
        unit=units,
    )


rule fastq_merge_fastp_stats:
    """Merge fastp stats across all units for a sample."""
    input:
        _get_fastp_jsons,
    output:
        "results/fastp/{sample}.fastp.json",
    run:
        import json

        total_before = 0
        total_after = 0

        for fn in input:
            with open(fn) as f:
                data = json.load(f)
            total_before += data["summary"]["before_filtering"]["total_reads"]
            total_after += data["summary"]["after_filtering"]["total_reads"]

        merged = {
            "summary": {
                "before_filtering": {"total_reads": total_before},
                "after_filtering": {"total_reads": total_after},
            }
        }

        with open(output[0], "w") as f:
            json.dump(merged, f)


# =============================================================================
# External Input Symlink Rules
# For samples with pre-aligned BAMs or pre-computed GVCFs
# =============================================================================
rule link_provided_bam:
    """Symlink user-provided BAM to results directory for BAM input type samples."""
    wildcard_constraints:
        sample="|".join(BAM_SAMPLES) if BAM_SAMPLES else "NOMATCH",
    input:
        unpack(get_provided_bam),
    output:
        bam="results/bams/{sample}.bam",
        bai="results/bams/{sample}.bam.bai",
    log:
        "logs/link_provided_bam/{sample}.txt",
    shell:
        """
        ln -sf $(realpath {input.bam}) {output.bam} 2> {log}
        ln -sf $(realpath {input.bai}) {output.bai} 2>> {log}
        """


rule link_provided_gvcf:
    """Symlink user-provided GVCF to results directory for GVCF input type samples."""
    wildcard_constraints:
        sample="|".join(GVCF_SAMPLES) if GVCF_SAMPLES else "NOMATCH",
    input:
        unpack(get_provided_gvcf),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    log:
        "logs/link_provided_gvcf/{sample}.txt",
    shell:
        """
        ln -sf $(realpath {input.gvcf}) {output.gvcf} 2> {log}
        ln -sf $(realpath {input.tbi}) {output.tbi} 2>> {log}
        """



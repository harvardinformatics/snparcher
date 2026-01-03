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
        prefetch --max-size 1T {wildcards.srr}
        prefetchExit=$?
        
        if [[ $prefetchExit -ne 0 ]]; then
            ffq --ftp {wildcards.srr} \
                | grep -Eo '"url": "[^"]*"' \
                | grep -o '"[^"]*"$' \
                | grep "fastq" \
                | xargs curl --remote-name-all --output-dir {params.outdir}
        else
            fasterq-dump {wildcards.srr} -O {params.outdir} -e {threads} -t {resources.tmpdir}
            pigz -p {threads} {params.outdir}{wildcards.srr}*.fastq
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
        json="results/summary_stats/{sample}/{unit}.fastp.out",
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

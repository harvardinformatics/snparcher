# The default mapping pipeline: fastp, bwa mem -M and per-library duplicate
# marking with sambamba, then the shared merges in rules/mapping/common.smk.

include: "fastp.smk"


rule bwa_mem:
    input:
        unpack(bwa_mem_input),
    output:
        bam=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam"),
    params:
        rg=get_read_group,
    threads: 8
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/bwa_mem/{sample}/{library}/{input_unit}.txt"
    log:
        "logs/bwa_mem/{sample}/{library}/{input_unit}.txt"
    shell:
        """
        bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} \
            | samtools sort -o {output.bam} - 2>> {log}
        """


rule markdup_library:
    input:
        unpack(dedup_library_input),
    output:
        bam=temp("results/bams/library_markdup/{sample}/{library}.bam"),
    threads: 4
    conda:
        "../../envs/sambamba.yaml"
    benchmark:
        "benchmarks/markdup/{sample}/{library}.txt"
    log:
        "logs/markdup/{sample}/{library}.txt"
    shell:
        """
        sambamba markdup -t {threads} {input.bam} {output.bam} 2> {log}
        rm -f {output.bam}.bai
        """

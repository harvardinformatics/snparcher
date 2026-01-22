rule picard_intervals:
    input:
        ref="results/{refGenome}/data/genome/{refGenome}.fna",
        fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
        dictf="results/{refGenome}/data/genome/{refGenome}.dict",
    output:
        intervals=temp("results/{refGenome}/intervals/picard_interval_list.list"),
    params:
        minNmer=int(config["minNmer"]),
    conda:
        "../envs/bam2vcf.yml"
    log:
        "logs/{refGenome}/picard_intervals/log.txt",
    benchmark:
        "benchmarks/{refGenome}/picard_intervals/benchmark.txt"
    shell:
        "picard ScatterIntervalsByNs -Xmx{resources.mem_mb_reduced}m REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} OUTPUT_TYPE=ACGT &> {log}\n"


rule format_interval_list:
    input:
        intervals="results/{refGenome}/intervals/picard_interval_list.list",
    output:
        intervals="results/{refGenome}/intervals/master_interval_list.list",
    run:
        with open(output.intervals, "w") as out:
            with open(input.intervals, "r") as inp:
                for line in inp:
                    if not line.startswith("@"):
                        line = line.strip().split("\t")
                        chrom, start, end = line[0], line[1], line[2]
                        print(f"{chrom}:{start}-{end}", file=out)


checkpoint create_db_intervals:
    input:
        ref="results/{refGenome}/data/genome/{refGenome}.fna",
        fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
        dictf="results/{refGenome}/data/genome/{refGenome}.dict",
        intervals="results/{refGenome}/intervals/master_interval_list.list",
    output:
        fof="results/{refGenome}/intervals/db_intervals/intervals.txt",
        out_dir=directory("results/{refGenome}/intervals/db_intervals"),
    params:
        max_intervals=get_db_interval_count,
    log:
        "logs/{refGenome}/db_intervals/log.txt",
    benchmark:
        "benchmarks/{refGenome}/db_intervals/benchmark.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        """
        gatk SplitIntervals --java-options '-Xmx{resources.mem_mb_reduced}m -Xms{resources.mem_mb_reduced}m' -L {input.intervals} \
        -O {output.out_dir} -R {input.ref} -scatter {params} \
        -mode INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof}
        """


checkpoint create_gvcf_intervals:
    input:
        ref="results/{refGenome}/data/genome/{refGenome}.fna",
        fai="results/{refGenome}/data/genome/{refGenome}.fna.fai",
        dictf="results/{refGenome}/data/genome/{refGenome}.dict",
        intervals="results/{refGenome}/intervals/master_interval_list.list",
    output:
        fof="results/{refGenome}/intervals/gvcf_intervals/intervals.txt",
        out_dir=directory("results/{refGenome}/intervals/gvcf_intervals"),
    params:
        max_intervals=config["num_gvcf_intervals"],
    log:
        "logs/{refGenome}/gvcf_intervals/log.txt",
    benchmark:
        "benchmarks/{refGenome}/gvcf_intervals/benchmark.txt"
    conda:
        "../envs/bam2vcf.yml"
    shell:
        """
        gatk SplitIntervals --java-options '-Xmx{resources.mem_mb_reduced}m -Xms{resources.mem_mb_reduced}m' -L {input.intervals} \
        -O {output.out_dir} -R {input.ref} -scatter {params} \
        -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY  &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof}
        """

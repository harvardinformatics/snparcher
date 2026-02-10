rule compute_d4:
    input:
        unpack(get_bams)
    output:
        dist = "results/{refGenome}/callable_sites/{sample}.mosdepth.global.dist.txt",
        d4="results/{refGenome}/callable_sites/{sample}.per-base.d4.gz",
        d4gzi ="results/{refGenome}/callable_sites/{sample}.per-base.d4.gz.gzi",
        summary="results/{refGenome}/callable_sites/{sample}.mosdepth.summary.txt"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{refGenome}/compute_d4/{sample}.txt"
    benchmark:
        "benchmarks/{refGenome}/compute_d4/{sample}.txt"
    params:
        prefix = subpath(output.summary, strip_suffix=".mosdepth.summary.txt"),
        d4 = subpath(output.d4, strip_suffix=".gz")
    shell:
        """
        mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}
        bgzip --index {params.d4}
        """



rule collect_covstats:
    input:
        unpack(get_input_covstats)
    output:
        "results/{refGenome}/summary_stats/all_cov_sumstats.txt"  
    run:
        covStats = collectCovStats(input.covStatFiles)
        with open(output[0], "w") as f:
            print("chrom\tmean_cov\tstdev_cov", file=f)
            for chrom in covStats:
                print(chrom, covStats[chrom]['mean'], covStats[chrom]['stdev'], sep="\t", file=f)

rule create_cov_thresholds:
    input:
        stats = "results/{refGenome}/summary_stats/all_cov_sumstats.txt",
    output:
        thresholds = "results/{refGenome}/callable_sites/{prefix}_callable_sites_thresholds.tsv"
    
    params:
        cov_threshold_stdev = config["cov_threshold_stdev"],
        cov_threshold_lower = config["cov_threshold_lower"],
        cov_threshold_upper = config["cov_threshold_upper"],
        cov_threshold_rel = config["cov_threshold_rel"]
    conda:
        "../envs/cov_filter.yml"
    script:
        "../scripts/create_coverage_thresholds.py"

rule clam_loci:
    input:
        unpack(get_input_for_coverage),
        thresholds = "results/{refGenome}/callable_sites/{prefix}_callable_sites_thresholds.tsv"
    output:
        cov = "results/{refGenome}/callable_sites/{prefix}/callable_sites.d4",
        bed = "results/{refGenome}/callable_sites/{prefix}/callable_sites.bed",
        tmp_bed = temp("results/{refGenome}/callable_sites/{prefix}/callable_sites_temp.bed") # temp fix until clam produces better bed files cm
    params:
        outdir = subpath(output.cov, parent=True)
    conda:
        "../envs/cov_filter.yml"
    log: 
        "logs/{refGenome}/covbed/{prefix}.txt"
    benchmark:
        "benchmarks/{refGenome}/covbed/{prefix}_benchmark.txt"
    shell:
        """
        clam loci -t {threads} --bed --thresholds-file {input.thresholds} -o {params.outdir} {input.d4} 2> {log}
        bedtk merge {output.bed} > {output.tmp_bed} 2>> {log}
        cp {output.tmp_bed} {output.bed} 2>> {log}
        """

rule callable_bed:
    input:
        cov = "results/{refGenome}/callable_sites/{prefix}/callable_sites.bed",
        map = "results/{refGenome}/callable_sites/{prefix}_callable_sites_map.bed"
    output:
        callable_sites = "results/{refGenome}/{prefix}_callable_sites.bed",
        tmp_cov = temp("results/{refGenome}/callable_sites/{prefix}_temp_cov.bed")
    conda:
        "../envs/cov_filter.yml"
    benchmark:
        "benchmarks/{refGenome}/callable_bed/{prefix}_benchmark.txt"
    params:
        merge = config['cov_merge']
    shell:
        """
        bedtools merge -d {params.merge} -i {input.cov} > {output.tmp_cov}
        bedtools intersect -a {output.tmp_cov} -b {input.map} | bedtools sort -i - | bedtools merge -i - > {output.callable_sites}
        """

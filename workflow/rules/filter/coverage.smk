"""
Coverage-based filtering rules for snpArcher v2.

Computes per-sample coverage using mosdepth and creates callable site BED files
based on coverage thresholds.

Rules:
- filter_compute_d4: Compute per-base coverage with mosdepth
- filter_collect_covstats: Aggregate coverage statistics across samples
- filter_create_thresholds: Calculate coverage thresholds
- filter_clam_loci: Identify callable loci based on coverage
- filter_callable_bed: Create final callable sites BED file
"""


rule filter_compute_d4:
    """Compute per-base coverage using mosdepth."""
    wildcard_constraints:
        sample="|".join(BAM_REQUIRED_SAMPLES) if BAM_REQUIRED_SAMPLES else "NOMATCH",
    input:
        unpack(get_final_bam),
    output:
        dist="results/callable_sites/{sample}.mosdepth.global.dist.txt",
        d4="results/callable_sites/{sample}.per-base.d4.gz",
        d4gzi="results/callable_sites/{sample}.per-base.d4.gz.gzi",
        summary="results/callable_sites/{sample}.mosdepth.summary.txt",
    conda:
        "../../envs/cov_filter.yml"
    log:
        "logs/filter_compute_d4/{sample}.txt"
    benchmark:
        "benchmarks/filter_compute_d4/{sample}.txt"
    threads: 4
    params:
        prefix="results/callable_sites/{sample}",
        d4="results/callable_sites/{sample}.per-base.d4",
    shell:
        """
        mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}
        bgzip --index {params.d4} &>> {log}
        """


rule filter_collect_covstats:
    """Aggregate coverage statistics across all samples."""
    input:
        cov_stat_files=get_mosdepth_summary_files(),
    output:
        "results/summary_stats/all_cov_sumstats.txt",
    run:
        covStats = collect_cov_stats(input.cov_stat_files)
        with open(output[0], "w") as f:
            print("chrom\tmean_cov\tstdev_cov", file=f)
            for chrom in covStats:
                print(chrom, covStats[chrom]["mean"], covStats[chrom]["stdev"], sep="\t", file=f)


rule filter_create_thresholds:
    """Calculate per-chromosome coverage thresholds for callable site detection."""
    input:
        stats="results/summary_stats/all_cov_sumstats.txt",
    output:
        thresholds=f"results/callable_sites/{REF_NAME}_callable_sites_thresholds.tsv",
    params:
        cov_threshold_stdev=config["callable_sites"]["coverage"]["stdev"],
    conda:
        "../../envs/cov_filter.yml"
    script:
        "../../scripts/create_coverage_thresholds.py"


rule filter_clam_loci:
    """Identify callable loci based on coverage thresholds using clam."""
    input:
        unpack(get_coverage_d4_files),
        thresholds=f"results/callable_sites/{REF_NAME}_callable_sites_thresholds.tsv",
    output:
        cov=f"results/callable_sites/{REF_NAME}/callable_sites.d4",
        bed=f"results/callable_sites/{REF_NAME}/callable_sites.bed",
        tmp_bed=temp(f"results/callable_sites/{REF_NAME}/callable_sites_temp.bed"),
    params:
        outdir=f"results/callable_sites/{REF_NAME}",
    conda:
        "../../envs/cov_filter.yml"
    log:
        "logs/filter_clam_loci.txt"
    benchmark:
        "benchmarks/filter_clam_loci.txt"
    threads: 4
    shell:
        """
        clam loci -t {threads} --bed --thresholds-file {input.thresholds} -o {params.outdir} {input.d4} 2> {log}
        bedtk merge {output.bed} > {output.tmp_bed} 2>> {log}
        cp {output.tmp_bed} {output.bed} 2>> {log}
        """


rule filter_callable_bed:
    """Create final callable sites BED by intersecting coverage and mappability."""
    input:
        cov=f"results/callable_sites/{REF_NAME}/callable_sites.bed",
        map=f"results/callable_sites/{REF_NAME}_callable_sites_map.bed",
    output:
        callable_sites=f"results/{REF_NAME}_callable_sites.bed",
        tmp_cov=temp(f"results/callable_sites/{REF_NAME}_temp_cov.bed"),
    conda:
        "../../envs/cov_filter.yml"
    log:
        "logs/filter_callable_bed.txt"
    benchmark:
        "benchmarks/filter_callable_bed.txt"
    params:
        merge=config["callable_sites"]["coverage"]["merge_distance"],
    shell:
        """
        bedtools merge -d {params.merge} -i {input.cov} > {output.tmp_cov} 2> {log}
        bedtools intersect -a {output.tmp_cov} -b {input.map} 2>> {log} \
            | bedtools sort -i - 2>> {log} \
            | bedtools merge -i - > {output.callable_sites} 2>> {log}
        """

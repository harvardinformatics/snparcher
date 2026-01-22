rule genmap:
    input:
        ref="results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        bg=temp("results/{refGenome}/genmap/{refGenome}.genmap.bedgraph"),
        sorted_bg="results/{refGenome}/genmap/sorted_mappability.bg",
    params:
        indir=os.path.join(DEFAULT_STORAGE_PREFIX, "results/{refGenome}/genmap_index"),
        outdir=os.path.join(DEFAULT_STORAGE_PREFIX, "results/{refGenome}/genmap"),
        kmer=config["mappability_k"],
    log:
        "logs/{refGenome}/genmap/log.txt",
    benchmark:
        "benchmarks/{refGenome}/genmap/benchmark.txt"
    conda:
        "../envs/mappability.yml"
    shell:
        # snakemake creates the output directory before the shell command, but genmap doesnt like this. so we remove the directory first.
        """
        rm -rf {params.indir} && genmap index -F {input.ref} -I {params.indir} &> {log}
        genmap map -K {params.kmer} -E 0 -I {params.indir} -O {params.outdir} -bg -T {threads} -v &> {log}
        sort -k1,1 -k2,2n {output.bg} > {output.sorted_bg} 2>> {log}
        """


rule mappability_bed:
    input:
        map="results/{refGenome}/genmap/sorted_mappability.bg",
    output:
        callable_sites=(
            "results/{refGenome}/callable_sites/{prefix}_callable_sites_map.bed"
            if config["cov_filter"]
            else "results/{refGenome}/{prefix}_callable_sites.bed"
        ),
        tmp_map=temp("results/{refGenome}/callable_sites/{prefix}_temp_map.bed"),
    conda:
        "../envs/mappability.yml"
    benchmark:
        "benchmarks/{refGenome}/mapbed/{prefix}_benchmark.txt"
    params:
        merge=config["mappability_merge"],
        mappability=config["mappability_min"],
    shell:
        """
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.mappability}) print $1,$2,$3 }}' {input.map} > {output.tmp_map}
        bedtools sort -i {output.tmp_map} | bedtools merge -d {params.merge} -i - > {output.callable_sites}
        """

"""
Mappability-based filtering rules for snpArcher v2.

Uses genmap to compute mappability and creates BED files of mappable regions.

Rules:
- filter_genmap: Compute mappability with genmap
- filter_mappability_bed: Create mappable regions BED file
"""


rule filter_genmap:
    """Compute genome mappability using genmap."""
    input:
        unpack(get_ref_bundle),
    output:
        bg=temp(f"results/genmap/{REF_NAME}.genmap.bedgraph"),
        sorted_bg=f"results/genmap/sorted_mappability.bg",
    params:
        indir=f"results/genmap_index",
        outdir=f"results/genmap",
        kmer=config["callable_sites"]["mappability"]["kmer"],
    conda:
        "../../envs/mappability.yml"
    log:
        "logs/filter_genmap.txt"
    benchmark:
        "benchmarks/filter_genmap.txt"
    threads: 8
    shell:
        """
        # genmap doesn't like pre-existing directories
        rm -rf {params.indir} && genmap index -F {input.ref} -I {params.indir} &> {log}
        genmap map -K {params.kmer} -E 0 -I {params.indir} -O {params.outdir} -bg -T {threads} -v &>> {log}
        sort -k1,1 -k2,2n {output.bg} > {output.sorted_bg} 2>> {log}
        """


rule filter_mappability_bed:
    """Create BED file of mappable regions based on mappability threshold."""
    input:
        map="results/genmap/sorted_mappability.bg",
    output:
        callable_sites=(
            f"results/callable_sites/{REF_NAME}_callable_sites_map.bed"
            if config["callable_sites"]["coverage"]["enabled"]
            else f"results/{REF_NAME}_callable_sites.bed"
        ),
        tmp_map=temp(f"results/callable_sites/{REF_NAME}_temp_map.bed"),
    conda:
        "../../envs/mappability.yml"
    log:
        "logs/filter_mappability_bed.txt"
    benchmark:
        "benchmarks/filter_mappability_bed.txt"
    params:
        merge=config["callable_sites"]["mappability"]["merge_distance"],
        mappability=config["callable_sites"]["mappability"]["min_score"],
    shell:
        """
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.mappability}) print $1,$2,$3 }}' {input.map} > {output.tmp_map} 2> {log}
        bedtools sort -i {output.tmp_map} 2>> {log} \
            | bedtools merge -d {params.merge} -i - > {output.callable_sites} 2>> {log}
        """

from pathlib import Path


GENMAP_SAMPLING_THRESHOLD_BYTES = 2 * 1024 * 1024 * 1024
GENMAP_SKEW_THRESHOLD_BYTES = 5 * 1024 * 1024 * 1024
GENMAP_SAMPLING_VALUE = 20


def get_callable_source_beds(_wildcards):
    beds = []
    if CALLABLE_COVERAGE_ENABLED:
        beds.append("results/callable_sites/coverage.bed")
    if CALLABLE_MAPPABILITY_ENABLED:
        beds.append("results/callable_sites/mappability.bed")
    return beds


rule mosdepth:
    """Compute per-base coverage with mosdepth."""
    input:
        bam=lambda wc: get_final_bam(wc.sample),
        bai=lambda wc: get_final_bam(wc.sample) + ".bai",
    output:
        d4=temp("results/callable_sites/depths/{sample}.per-base.d4"),
        summary="results/callable_sites/depths/{sample}.mosdepth.summary.txt",
    params:
        prefix="results/callable_sites/depths/{sample}",
    threads: 4
    conda:
        "../envs/mosdepth.yaml"
    benchmark:
        "benchmarks/mosdepth/{sample}.txt"
    log:
        "logs/mosdepth/{sample}.txt"
    shell:
        """
        mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}
        """


rule clam_collect:
    """Collect depth data from D4 files into Zarr store."""
    input:
        d4=expand(
            "results/callable_sites/depths/{sample}.per-base.d4",
            sample=SAMPLES_WITH_BAM,
        ),
    output:
        zarr=directory("results/callable_sites/depths.zarr"),
    threads: 8
    conda:
        "../envs/clam.yaml"
    benchmark:
        "benchmarks/clam_collect.txt"
    log:
        "logs/clam_collect.txt"
    shell:
        """
        clam collect \
            -o {output.zarr} \
            -t {threads} \
            {input.d4} \
            &> {log}
        """


rule callable_coverage_thresholds:
    """Calculate cohort-wide coverage thresholds for callable loci."""
    input:
        summaries=expand(
            "results/callable_sites/depths/{sample}.mosdepth.summary.txt",
            sample=SAMPLES_WITH_BAM,
        ),
    output:
        tsv="results/callable_sites/coverage_thresholds.tsv",
    params:
        min_coverage=config["callable_sites"]["coverage"]["min_coverage"],
        max_coverage=config["callable_sites"]["coverage"]["max_coverage"],
        script=Path(workflow.basedir) / "scripts" / "callable_coverage_thresholds.py",
    conda:
        "../envs/callable_sites.yaml"
    benchmark:
        "benchmarks/callable_coverage_thresholds.txt"
    log:
        "logs/callable_coverage_thresholds.txt"
    shell:
        """
        python {params.script} \
            --min-coverage {params.min_coverage} \
            --max-coverage {params.max_coverage} \
            --output {output.tsv} \
            {input.summaries} \
            &> {log}
        """


rule clam_loci:
    """Identify callable loci based on coverage."""
    input:
        zarr="results/callable_sites/depths.zarr",
        thresholds="results/callable_sites/coverage_thresholds.tsv",
    output:
        zarr=directory("results/callable_sites/callable_loci.zarr"),
    threads: 8
    conda:
        "../envs/clam.yaml"
    benchmark:
        "benchmarks/clam_loci.txt"
    log:
        "logs/clam_loci.txt"
    shell:
        """
        min_depth=$(awk -F '\t' 'NR==2 {{print $2}}' {input.thresholds})
        max_depth=$(awk -F '\t' 'NR==2 {{print $3}}' {input.thresholds})
        clam loci \
            -o {output.zarr} \
            -t {threads} \
            --per-sample \
            -m "${{min_depth}}" \
            -M "${{max_depth}}" \
            {input.zarr} \
            &> {log}
        """


rule coverage_bed:
    """Create BED of callable regions based on coverage."""
    input:
        zarr="results/callable_sites/callable_loci.zarr",
    output:
        bed="results/callable_sites/coverage.bed",
    params:
        fraction=config["callable_sites"]["coverage"]["fraction"],
        merge_distance=config["callable_sites"]["coverage"]["merge_distance"],
        script=Path(workflow.basedir) / "scripts" / "callable_zarr_to_bed.py",
    conda:
        "../envs/callable_sites.yaml"
    benchmark:
        "benchmarks/coverage_bed.txt"
    log:
        "logs/coverage_bed.txt"
    shell:
        """
        python {params.script} {input.zarr} /dev/stdout --fraction {params.fraction} 2> {log} \
            | bedtools sort -i - 2>> {log} \
            | bedtools merge -d {params.merge_distance} -i - > {output.bed} 2>> {log}
        """


rule genmap_index:
    """Create genmap index for mappability calculation."""
    input:
        ref=REF_FILES["ref"],
    output:
        idx=directory("results/callable_sites/genmap_index"),
    params:
        ref_decompressed=lambda wc, input: input.ref.replace(".gz", ""),
        sampling_threshold_bytes=GENMAP_SAMPLING_THRESHOLD_BYTES,
        skew_threshold_bytes=GENMAP_SKEW_THRESHOLD_BYTES,
        sampling_value=GENMAP_SAMPLING_VALUE,
    conda:
        "../envs/genmap.yaml"
    benchmark:
        "benchmarks/genmap_index.txt"
    log:
        "logs/genmap_index.txt"
    shell:
        """
        set -euo pipefail

        REF_FASTA="{params.ref_decompressed}"
        trap 'rm -f "$REF_FASTA"' EXIT

        gunzip -c "{input.ref}" > "$REF_FASTA"
        SIZE_BYTES=$(wc -c < "$REF_FASTA" | tr -d '[:space:]')

        : > "{log}"
        echo "Decompressed FASTA size (bytes): $SIZE_BYTES" >> "{log}"

        if (( SIZE_BYTES > {params.skew_threshold_bytes} )); then
            echo "GenMap index mode: skew (-A skew)" >> "{log}"
            genmap index -F "$REF_FASTA" -I "{output.idx}" -A skew >> "{log}" 2>&1
        elif (( SIZE_BYTES > {params.sampling_threshold_bytes} )); then
            echo "GenMap index mode: sampled (-S {params.sampling_value})" >> "{log}"
            genmap index -F "$REF_FASTA" -I "{output.idx}" -S {params.sampling_value} >> "{log}" 2>&1
        else
            echo "GenMap index mode: default (divsufsort)" >> "{log}"
            genmap index -F "$REF_FASTA" -I "{output.idx}" >> "{log}" 2>&1
        fi
        """


rule genmap_mappability:
    """Compute mappability scores with genmap."""
    input:
        idx="results/callable_sites/genmap_index",
    output:
        bg="results/callable_sites/mappability.bedgraph",
    params:
        kmer=config["callable_sites"]["mappability"]["kmer"],
        prefix="results/callable_sites/mappability",
    threads: 4
    conda:
        "../envs/genmap.yaml"
    benchmark:
        "benchmarks/genmap_mappability.txt"
    log:
        "logs/genmap_mappability.txt"
    shell:
        """
        genmap map \
            -K {params.kmer} \
            -E 2 \
            -I {input.idx} \
            -O {params.prefix} \
            -bg \
            -T {threads} \
            &> {log}
        """


rule mappability_bed:
    """Create BED of high-mappability regions."""
    input:
        bg="results/callable_sites/mappability.bedgraph",
    output:
        bed="results/callable_sites/mappability.bed",
    params:
        min_score=config["callable_sites"]["mappability"]["min_score"],
        merge_distance=config["callable_sites"]["mappability"]["merge_distance"],
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/mappability_bed.txt"
    log:
        "logs/mappability_bed.txt"
    shell:
        """
        awk -v min={params.min_score} '$4 >= min {{print $1"\\t"$2"\\t"$3}}' {input.bg} 2> {log} \
            | bedtools sort -i - 2>> {log} \
            | bedtools merge -d {params.merge_distance} -i - > {output.bed} 2>> {log}
        """


rule callable_sites_bed:
    """Combine enabled callable site sources into a final callable BED."""
    input:
        beds=get_callable_source_beds,
    output:
        bed="results/callable_sites/callable_sites.bed",
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/callable_sites_bed.txt"
    log:
        "logs/callable_sites_bed.txt"
    shell:
        """
        cat {input.beds} 2> {log} \
            | bedtools sort -i - 2>> {log} \
            | bedtools merge -i - > {output.bed} 2>> {log}
        """

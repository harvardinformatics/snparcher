from pathlib import Path


GENMAP_SAMPLING_THRESHOLD_BYTES = 2 * 1024 * 1024 * 1024
GENMAP_SKEW_THRESHOLD_BYTES = 5 * 1024 * 1024 * 1024
GENMAP_SAMPLING_VALUE = 20


def get_callable_source_beds(_wildcards):
    check_callable_sites_not_blocked()
    beds = []
    if CALLABLE_COVERAGE_BED_ENABLED:
        beds.append("results/callable_sites/coverage.bed")
    if CALLABLE_MAPPABILITY_ENABLED:
        beds.append("results/callable_sites/mappability.bed")
    return beds


def get_coverage_bed_zarr_input(_wildcards):
    check_callable_sites_not_blocked()
    if CALLABLE_COVERAGE_BED_ENABLED:
        return "results/callable_sites/callable_loci.zarr"

    if SAMPLES_WITH_GVCF:
        raise ValueError(
            "Cannot create results/callable_sites/coverage.bed when the sample sheet "
            "contains gVCF inputs, because coverage is unavailable for those samples. "
            "Coverage statistics can still be computed for BAM-backed samples."
        )

    raise ValueError(
        "Cannot create results/callable_sites/coverage.bed because "
        "callable_sites.coverage.enabled is false or no BAM-backed samples are available."
    )


def get_coverage_threshold_mappability(_wildcards):
    check_callable_sites_not_blocked()
    if CALLABLE_MAPPABILITY_ENABLED:
        return "results/callable_sites/mappability.bed"
    return []


# Samples excluded in sample_metadata are left out of the coverage mask, so
# it reflects the cohort that downstream filtering keeps.
CALLABLE_SAMPLES = [s for s in SAMPLES_WITH_BAM if s not in set(get_excluded_samples())]

if CALLABLE_COVERAGE_STATS_ENABLED and not CALLABLE_SAMPLES:
    raise WorkflowError(
        "callable_sites.coverage.enabled is true, but every BAM-backed sample is "
        "excluded in sample_metadata. Disable callable_sites.coverage or include "
        "at least one BAM-backed sample."
    )


rule mosdepth:
    """Compute per-base coverage with mosdepth."""
    input:
        bam=lambda wc: get_final_bam(wc.sample),
        bam_index=lambda wc: get_final_bam_index(wc.sample),
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
            sample=CALLABLE_SAMPLES,
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
    """Choose clam loci thresholds from the cohort's summed-depth distribution."""
    input:
        zarr="results/callable_sites/depths.zarr",
        mappability=get_coverage_threshold_mappability,
    output:
        tsv="results/callable_sites/coverage_thresholds.tsv",
    params:
        max_depth=config["callable_sites"]["coverage"]["max_depth"],
        max_depth_k=config["callable_sites"]["coverage"]["max_depth_k"],
        min_sample_depth=config["callable_sites"]["coverage"]["min_sample_depth"],
        fraction=config["callable_sites"]["coverage"]["fraction"],
        mappability_arg=lambda wc, input: (
            f"--mappability-bed {input.mappability}" if input.mappability else ""
        ),
        script=Path(workflow.basedir) / "scripts" / "callable_coverage_thresholds.py",
    conda:
        "../envs/callable_sites.yaml"
    benchmark:
        "benchmarks/callable_coverage_thresholds.txt"
    log:
        "logs/callable_coverage_thresholds.txt"
    shell:
        """
        python {params.script} {input.zarr} \
            {params.mappability_arg} \
            --max-depth {params.max_depth} \
            --max-depth-k {params.max_depth_k} \
            --min-sample-depth {params.min_sample_depth} \
            --fraction {params.fraction} \
            --output {output.tsv} \
            &> {log}
        """


rule clam_loci:
    """Count, per site, the samples that make it coverage-callable.

    Population-count mode: clam zeroes the count wherever the site fails a
    threshold, and the counts are what `clam stat` needs as denominators.
    """
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
        threshold() {{ awk -F '\t' -v key="$1" '$1 == key {{print $2}}' {input.thresholds}; }}
        clam loci \
            -o {output.zarr} \
            -t {threads} \
            -m "$(threshold min_sample_depth)" \
            -d "$(threshold clam_min_proportion)" \
            --max-mean-depth "$(threshold clam_max_mean_depth)" \
            {input.zarr} \
            &> {log}
        """


rule coverage_bed:
    """Create BED of callable regions based on coverage."""
    input:
        zarr=get_coverage_bed_zarr_input,
    output:
        bed="results/callable_sites/coverage.bed",
    params:
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
        python {params.script} {input.zarr} {output.bed} \
            --merge-distance {params.merge_distance} \
            &> {log}
        """


rule genmap_index:
    """Create genmap index for mappability calculation."""
    input:
        ref=REF_FILES["ref"],
    output:
        idx=directory("results/callable_sites/genmap_index"),
    params:
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

        # GenMap needs an uncompressed FASTA. Decompress it into a job-private
        # directory: results/reference/{{name}}.fa can be a FASTA that the user
        # or another rule owns, and this job must not overwrite or delete it.
        WORK=$(mktemp -d "{resources.tmpdir}/genmap_index_XXXXXX")
        trap 'rm -rf "$WORK"' EXIT
        REF_FASTA="$WORK/$(basename "{input.ref}" .gz)"

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
    """Intersect enabled callable site sources into a final callable BED."""
    input:
        beds=get_callable_source_beds,
    output:
        bed="results/callable_sites/callable_sites.bed",
    params:
        # A site is callable only if every enabled source passes it.
        intersect=lambda wc, input: " ".join(
            f"| bedtools intersect -a - -b {bed}" for bed in input.beds[1:]
        ),
    conda:
        "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/callable_sites_bed.txt"
    log:
        "logs/callable_sites_bed.txt"
    shell:
        """
        (cat {input.beds[0]} {params.intersect}) 2> {log} \
            | bedtools sort -i - 2>> {log} \
            | bedtools merge -i - > {output.bed} 2>> {log}
        """

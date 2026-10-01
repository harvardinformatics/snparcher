# Shared mapping code: external BAM staging, CSI indexing, the per-library
# merge rules used by the default and sentieon pipelines, and bam_stats.
# Pipeline-specific rules are in rules/mapping/<pipeline>.smk; see
# MAPPING_PIPELINES in rules/common.smk.

from pathlib import Path

def bwa_mem_input(wildcards):
    """Get input fastqs for alignment."""
    return {
        "r1": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}/{wildcards.input_unit}_1.fastq.gz",
        "r2": f"results/filtered_fastqs/{wildcards.sample}/{wildcards.library}/{wildcards.input_unit}_2.fastq.gz",
        **REF_FILES,
    }


def get_read_group(wildcards):
    """Generate read group string for BWA."""
    return (
        f"'@RG\\tID:{wildcards.library}.{wildcards.input_unit}"
        f"\\tSM:{wildcards.sample}\\tLB:{wildcards.library}\\tPL:ILLUMINA'"
    )


def get_library_rows(sample, library):
    """Get all rows for a sample/library pair."""
    sample_rows = get_sample_rows(sample)
    library_rows = sample_rows[sample_rows["library_id"] == library]
    if library_rows.empty:
        raise ValueError(f"No rows found for sample '{sample}', library '{library}'")
    return library_rows


def merge_library_bams_input(wildcards):
    """Get all per-input BAMs for a library."""
    library_rows = get_library_rows(wildcards.sample, wildcards.library)
    input_units = library_rows["input_unit"].tolist()
    return {
        "bams": [
            f"results/bams/raw/{wildcards.sample}/{wildcards.library}/{input_unit}.bam"
            for input_unit in input_units
        ],
    }


def dedup_library_input(wildcards):
    """Get input BAM for library-level duplicate marking."""
    bam = f"results/bams/library/{wildcards.sample}/{wildcards.library}.bam"
    return {
        "bam": bam,
    }


def merge_dedup_libraries_input(wildcards):
    """Get deduplicated library BAMs for sample-level merge."""
    if not get_sample_mark_duplicates(wildcards.sample):
        raise ValueError(
            f"Sample '{wildcards.sample}' has mark_duplicates=False, "
            "but merge_dedup_libraries was requested."
        )
    libraries = get_sample_libraries(wildcards.sample)
    return {
        "bams": [
            f"results/bams/library_markdup/{wildcards.sample}/{lib}.bam"
            for lib in libraries
        ],
    }


def merge_library_level_bams_input(wildcards):
    """Get library BAMs for sample-level merge when duplicate marking is disabled."""
    if get_sample_mark_duplicates(wildcards.sample):
        raise ValueError(
            f"Sample '{wildcards.sample}' has mark_duplicates=True, "
            "but merge_library_level_bams was requested."
        )
    libraries = get_sample_libraries(wildcards.sample)
    return {
        "bams": [f"results/bams/library/{wildcards.sample}/{lib}.bam" for lib in libraries],
    }


rule stage_external_bam:
    input:
        bam=lambda wc: get_external_bam(wc.sample),
    output:
        bam="results/bams/input/{sample}.bam",
    log:
        "logs/stage_external_bam/{sample}.txt"
    run:
        src = Path(input.bam).resolve()
        dest = Path(output.bam)
        dest.parent.mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)

        if not src.exists():
            raise FileNotFoundError(f"External BAM not found for sample {wildcards.sample}: {src}")

        if dest.exists() or dest.is_symlink():
            if dest.resolve() == src:
                with open(log[0], "w") as handle:
                    handle.write(f"External BAM already staged: {dest} -> {src}\n")
                return
            dest.unlink()

        dest.symlink_to(src)
        with open(log[0], "w") as handle:
            handle.write(f"Staged external BAM: {dest} -> {src}\n")


rule index_bam_csi:
    input:
        bam="{bam}.bam",
    output:
        csi="{bam}.bam.csi",
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/index_bam_csi/{bam}.txt"
    log:
        "logs/index_bam_csi/{bam}.txt"
    shell:
        """
        samtools index -c {input.bam} {output.csi} 2> {log}
        """


rule merge_library_bams:
    input:
        unpack(merge_library_bams_input),
    output:
        bam=temp("results/bams/library/{sample}/{library}.bam"),
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_library_bams/{sample}/{library}.txt"
    log:
        "logs/merge_library_bams/{sample}/{library}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        """


rule merge_dedup_libraries:
    input:
        unpack(merge_dedup_libraries_input),
    output:
        bam="results/bams/markdup/{sample}.bam",
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_dedup_libraries/{sample}.txt"
    log:
        "logs/merge_dedup_libraries/{sample}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        """


rule merge_library_level_bams:
    input:
        unpack(merge_library_level_bams_input),
    output:
        bam="results/bams/merged/{sample}.bam",
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/merge_bams/{sample}.txt"
    log:
        "logs/merge_bams/{sample}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        """


rule bam_stats:
    input:
        bam=lambda wc: get_final_bam(wc.sample),
    output:
        coverage=temp("results/qc_metrics/bam/{sample}_coverage.txt"),
        flagstat=temp("results/qc_metrics/bam/{sample}_flagstat.txt"),
    conda:
        "../../envs/samtools.yaml"
    benchmark:
        "benchmarks/bam_stats/{sample}.txt"
    log:
        "logs/bam_stats/{sample}.txt"
    shell:
        """
        samtools coverage {input.bam} -o {output.coverage} 2> {log}
        samtools flagstat -O tsv {input.bam} > {output.flagstat} 2>> {log}
        """

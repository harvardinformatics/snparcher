"""
snpArcher v2 - Common configuration, sample sheet parsing, and helper functions.
"""

import os
import statistics
from pathlib import Path

import pandas as pd
from snakemake.utils import min_version, validate

min_version("9.0")

validate(config, schema="../schemas/config.schema.json")

# =============================================================================
# Config Defaults
# =============================================================================
config.setdefault("variant_calling", {})
config["variant_calling"].setdefault("tool", "gatk")
config["variant_calling"].setdefault("gatk", {})
config["variant_calling"]["gatk"].setdefault("min_pruning", 1)
config["variant_calling"]["gatk"].setdefault("min_dangling", 1)
config["variant_calling"]["gatk"].setdefault("het_prior", 0.005)
config["variant_calling"]["gatk"].setdefault("ploidy", 2)
config["variant_calling"].setdefault("sentieon", {})
config["variant_calling"]["sentieon"].setdefault("license", "")

config.setdefault("intervals", {})
config["intervals"].setdefault("enabled", True)
config["intervals"].setdefault("min_nmer", 500)
config["intervals"].setdefault("num_gvcf_intervals", 50)
config["intervals"].setdefault("db_scatter_factor", 0.15)

config.setdefault("reads", {})
config["reads"].setdefault("mark_duplicates", True)
config["reads"].setdefault("sort", False)

config.setdefault("remote_reads", {})
config["remote_reads"].setdefault("enabled", False)
config["remote_reads"].setdefault("prefix", "")

config.setdefault("callable_sites", {})
config["callable_sites"].setdefault("coverage", {})
config["callable_sites"]["coverage"].setdefault("enabled", True)
config["callable_sites"]["coverage"].setdefault("stdev", 2)
config["callable_sites"]["coverage"].setdefault("merge_distance", 100)
config["callable_sites"].setdefault("mappability", {})
config["callable_sites"]["mappability"].setdefault("enabled", True)
config["callable_sites"]["mappability"].setdefault("min_score", 1)
config["callable_sites"]["mappability"].setdefault("kmer", 150)
config["callable_sites"]["mappability"].setdefault("merge_distance", 100)

config.setdefault("normalize_gvcfs", True)

REF_NAME = config["reference"]["name"]
REF_PATH = config["reference"].get("path")
REF_ACCESSION = config["reference"].get("accession")
REF_FILE = f"{REF_NAME}.fna.gz"  # Reference filename (bgzip-compressed)

samples = pd.read_csv(config["samples"])
validate(samples, schema="../schemas/samples.schema.json")


def _validate_library_ids_for_multi_row_samples(df: pd.DataFrame) -> None:
    """Error if a sample_id appears multiple times without explicit library_id."""
    sample_counts = df.groupby("sample_id").size()
    multi_row_samples = sample_counts[sample_counts > 1].index.tolist()

    if not multi_row_samples:
        return

    has_library_col = "library_id" in df.columns

    for sample_id in multi_row_samples:
        sample_rows = df[df["sample_id"] == sample_id]
        if has_library_col:
            missing_mask = sample_rows["library_id"].isna() | (sample_rows["library_id"] == "")
        else:
            missing_mask = pd.Series([True] * len(sample_rows), index=sample_rows.index)

        if missing_mask.any():
            missing_indices = missing_mask[missing_mask].index.tolist()
            raise ValueError(
                f"Sample '{sample_id}' has {len(sample_rows)} rows but is missing "
                f"library_id on rows {missing_indices}. "
                f"When a sample has multiple rows, library_id must be specified for all rows."
            )


def _validate_no_mixed_input_types(df: pd.DataFrame) -> None:
    """Error if a sample_id has rows with different input_types."""
    for sample_id, group in df.groupby("sample_id"):
        input_types = group["input_type"].unique()
        if len(input_types) > 1:
            raise ValueError(
                f"Sample '{sample_id}' has mixed input_types: {', '.join(input_types)}. "
                f"All rows for a sample must have the same input_type."
            )


def _validate_no_duplicate_rows(df: pd.DataFrame) -> None:
    """Error if exact duplicate rows exist (same sample_id + library_id + input)."""
    cols = ["sample_id", "input"]
    if "library_id" in df.columns:
        cols.append("library_id")

    duplicates = df[df.duplicated(subset=cols, keep=False)]
    if not duplicates.empty:
        first_dup = duplicates.iloc[0]
        lib_info = f", library_id='{first_dup['library_id']}'" if "library_id" in cols else ""
        raise ValueError(
            f"Duplicate row detected: sample_id='{first_dup['sample_id']}'{lib_info}, "
            f"input='{first_dup['input']}'"
        )


def _fill_missing_library_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing library_id values with sample_id (only for single-row samples)."""
    df = df.copy()
    if "library_id" not in df.columns:
        df["library_id"] = df["sample_id"]
    else:
        mask = df["library_id"].isna() | (df["library_id"] == "")
        df.loc[mask, "library_id"] = df.loc[mask, "sample_id"]
    return df


def _assign_unit_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Assign sequential unit IDs (0, 1, 2...) within each sample."""
    df = df.copy()
    df["_unit_id"] = df.groupby("sample_id").cumcount()
    return df


def _parse_fastq_paths(df: pd.DataFrame) -> pd.DataFrame:
    """Parse semicolon-separated FASTQ paths into _fq1 and _fq2 columns."""
    df = df.copy()
    fq1_list = []
    fq2_list = []

    for _, row in df.iterrows():
        if row["input_type"] == "fastq":
            parts = row["input"].split(";")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid fastq input for {row['sample_id']}: "
                    f"expected 'r1;r2' format, got '{row['input']}'"
                )
            fq1_list.append(parts[0].strip())
            fq2_list.append(parts[1].strip())
        else:
            fq1_list.append(None)
            fq2_list.append(None)

    df["_fq1"] = fq1_list
    df["_fq2"] = fq2_list
    return df


_validate_library_ids_for_multi_row_samples(samples)
_validate_no_mixed_input_types(samples)
_validate_no_duplicate_rows(samples)

samples = _fill_missing_library_ids(samples)
samples = _assign_unit_ids(samples)
samples = _parse_fastq_paths(samples)

SAMPLE_IDS = samples["sample_id"].unique().tolist()


# =============================================================================
# Sample Type Infrastructure
# =============================================================================
def get_sample_input_type(sample_id: str) -> str:
    """Get the input type for a sample (all rows must have same type)."""
    return samples.loc[samples["sample_id"] == sample_id, "input_type"].iloc[0]


def get_samples_by_type(input_type: str) -> list[str]:
    """Get list of sample IDs with a specific input type."""
    mask = samples["input_type"] == input_type
    return samples.loc[mask, "sample_id"].unique().tolist()


# Sample lists by input type for routing
FASTQ_SAMPLES = get_samples_by_type("fastq")
SRR_SAMPLES = get_samples_by_type("srr")
BAM_SAMPLES = get_samples_by_type("bam")
GVCF_SAMPLES = get_samples_by_type("gvcf")

# Samples requiring alignment (fastq + srr)
ALIGN_SAMPLES = FASTQ_SAMPLES + SRR_SAMPLES

# Samples with BAM files (generated or provided, not GVCF-only)
BAM_REQUIRED_SAMPLES = ALIGN_SAMPLES + BAM_SAMPLES


def get_provided_bam(wildcards) -> dict[str, str]:
    """Get user-provided BAM path for BAM input type samples."""
    row = samples.loc[samples["sample_id"] == wildcards.sample].iloc[0]
    if row["input_type"] != "bam":
        raise ValueError(f"Sample {wildcards.sample} is not BAM input type")
    bam_path = row["input"]
    return {"bam": bam_path, "bai": f"{bam_path}.bai"}


def get_provided_gvcf(wildcards) -> dict[str, str]:
    """Get user-provided GVCF path for GVCF input type samples."""
    row = samples.loc[samples["sample_id"] == wildcards.sample].iloc[0]
    if row["input_type"] != "gvcf":
        raise ValueError(f"Sample {wildcards.sample} is not GVCF input type")
    gvcf_path = row["input"]
    return {"gvcf": gvcf_path, "tbi": f"{gvcf_path}.tbi"}


wildcard_constraints:
    sample="[A-Za-z0-9_-]+",
    unit="[0-9]+",


def get_sample_units(sample_id: str) -> list[int]:
    """Get all unit IDs for a sample."""
    return samples.loc[samples["sample_id"] == sample_id, "_unit_id"].tolist()


def get_sample_libraries(sample_id: str) -> list[str]:
    """Get unique library IDs for a sample."""
    return samples.loc[samples["sample_id"] == sample_id, "library_id"].unique().tolist()


def get_units_for_library(sample_id: str, library_id: str) -> list[int]:
    """Get unit IDs for a specific sample+library combination."""
    mask = (samples["sample_id"] == sample_id) & (samples["library_id"] == library_id)
    return samples.loc[mask, "_unit_id"].tolist()


def get_unit_row(sample_id: str, unit_id: int) -> pd.Series:
    """Get the sample sheet row for a specific sample+unit."""
    mask = (samples["sample_id"] == sample_id) & (samples["_unit_id"] == unit_id)
    return samples.loc[mask].iloc[0]


def get_ref_path():
    """Get the reference genome path."""
    if REF_PATH:
        return REF_PATH
    elif REF_ACCESSION:
        return f"results/reference/{REF_FILE}"
    else:
        raise ValueError("No reference path or accession specified in config")


# =============================================================================
# Reference File Input Functions
# =============================================================================
def get_ref_bundle(wildcards=None) -> dict[str, str]:
    """Get all reference files as a bundle for rule inputs.

    Usage:
        input:
            unpack(get_ref_bundle),  # Unpacks to ref=, fai=, dictf=
    """
    return {
        "ref": f"results/reference/{REF_FILE}",
        "fai": f"results/reference/{REF_FILE}.fai",
        "dictf": f"results/reference/{REF_NAME}.dict",
    }


def get_ref_with_bwa(wildcards=None) -> dict[str, str | list[str]]:
    """Get reference files plus BWA indexes."""
    bundle = get_ref_bundle()
    bundle["bwa_idx"] = expand(
        f"results/reference/{REF_FILE}.{{ext}}",
        ext=["sa", "pac", "bwt", "ann", "amb"],
    )
    return bundle


def get_reads(wildcards):
    """Get read files for a sample unit."""
    row = get_unit_row(wildcards.sample, int(wildcards.unit))

    if row["input_type"] == "fastq":
        r1 = row["_fq1"]
        r2 = row["_fq2"]

        if config["remote_reads"]["enabled"]:
            return {"r1": storage(r1), "r2": storage(r2)}

        if not (os.path.exists(r1) and os.path.exists(r2)):
            raise WorkflowError(
                f"FASTQ files not found for {wildcards.sample} unit {wildcards.unit}: {r1}, {r2}"
            )
        return {"r1": r1, "r2": r2}

    elif row["input_type"] == "srr":
        srr = row["input"]
        return {
            "r1": f"results/fastq/{wildcards.sample}/{srr}_1.fastq.gz",
            "r2": f"results/fastq/{wildcards.sample}/{srr}_2.fastq.gz",
        }

    else:
        raise ValueError(f"Unsupported input_type: {row['input_type']}")


def get_reads_fastp(wildcards):
    """Get reads for fastp, optionally sorted."""
    if config["reads"]["sort"]:
        return {
            "r1": f"results/sorted_reads/{wildcards.sample}/{wildcards.unit}_1.fastq.gz",
            "r2": f"results/sorted_reads/{wildcards.sample}/{wildcards.unit}_2.fastq.gz",
        }
    else:
        return get_reads(wildcards)


def get_read_group(wildcards):
    """Generate read group string for BWA."""
    row = get_unit_row(wildcards.sample, int(wildcards.unit))
    lib = row["library_id"]
    return rf"'@RG\tID:{lib}_{wildcards.unit}\tSM:{wildcards.sample}\tLB:{lib}\tPL:ILLUMINA'"


def get_unit_bams(wildcards):
    """Get all pre-merge BAM files for a sample."""
    units = get_sample_units(wildcards.sample)
    return expand(
        "results/bams/preMerge/{sample}/{unit}.bam",
        sample=wildcards.sample,
        unit=units,
    )


def get_dedup_input(wildcards):
    """Get BAM input for deduplication."""
    units = get_sample_units(wildcards.sample)

    if len(units) == 1:
        return {
            "bam": f"results/bams/preMerge/{wildcards.sample}/{units[0]}.bam",
            "bai": f"results/bams/preMerge/{wildcards.sample}/{units[0]}.bam.bai",
        }
    else:
        return {
            "bam": f"results/bams/merged/{wildcards.sample}.bam",
            "bai": f"results/bams/merged/{wildcards.sample}.bam.bai",
        }


def get_final_bam(wildcards):
    """Get final BAM path for a sample.

    Routes based on input type:
    - BAM input: returns symlinked BAM path in results/bams/
    - fastq/srr: returns workflow-generated BAM

    Note: This function should not be called for GVCF samples.
    Use wildcard constraints to prevent GVCF samples from reaching rules that need BAMs.
    """
    input_type = get_sample_input_type(wildcards.sample)

    if input_type == "bam":
        # BAM samples use symlinked BAM in results/bams/
        return {
            "bam": f"results/bams/{wildcards.sample}.bam",
            "bai": f"results/bams/{wildcards.sample}.bam.bai",
        }

    # fastq/srr samples use workflow-generated BAM
    if config["reads"]["mark_duplicates"]:
        return {
            "bam": f"results/bams/{wildcards.sample}.bam",
            "bai": f"results/bams/{wildcards.sample}.bam.bai",
        }
    else:
        return get_dedup_input(wildcards)


def get_coverage_d4_files(wildcards=None):
    return {
        "d4": expand("results/callable_sites/{sample}.per-base.d4.gz", sample=SAMPLE_IDS),
        "d4gzi": expand("results/callable_sites/{sample}.per-base.d4.gz.gzi", sample=SAMPLE_IDS),
    }


def get_mosdepth_summary_files():
    return expand("results/callable_sites/{sample}.mosdepth.summary.txt", sample=SAMPLE_IDS)


def get_java_mem(wildcards, resources):
    """Calculate Java heap size from resources.mem_mb with 10% overhead for JVM.

    Use this as a params function to automatically calculate the Java heap size
    based on the mem_mb resource. This removes the need to manually specify
    mem_mb_reduced in rules.

    Args:
        wildcards: Snakemake wildcards object (unused, required for params signature)
        resources: Snakemake resources object

    Returns:
        Java heap size in MB as integer (90% of mem_mb)

    Example:
        rule example:
            resources:
                mem_mb=8000,
            params:
                java_mem=get_java_mem,
            shell:
                "gatk Tool --java-options '-Xmx{params.java_mem}m' ..."
    """
    import math

    mem_mb = getattr(resources, "mem_mb", 4096)
    return math.floor(mem_mb * 0.9)


def setup_curlrc():
    """Add -L flag to curlrc for pyd4 remote access."""
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"

    if curlrc_path.exists():
        with curlrc_path.open("r+") as f:
            if "-L" not in f.read():
                f.write(f"\n{entry}\n")
    else:
        with curlrc_path.open("a+") as f:
            f.write(f"{entry}\n")


def cleanup_curlrc():
    """Remove snpArcher-added -L flag from curlrc."""
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"

    if curlrc_path.exists():
        with curlrc_path.open("r") as f:
            lines = f.readlines()
        new_lines = [line for line in lines if line.strip() != entry.strip()]
        if len(new_lines) == 0:
            curlrc_path.unlink()
        else:
            with curlrc_path.open("w") as f:
                f.writelines(new_lines)


def collect_cov_stats(cov_sum_files):
    """Parse mosdepth coverage summary files."""
    sample_cov = {}

    for fn in cov_sum_files:
        with open(fn, "r") as f:
            for line in f:
                if "mean" in line:
                    continue
                fields = line.split()
                chrom = fields[0]
                cov = float(fields[3])

                if chrom in sample_cov:
                    sample_cov[chrom].append(cov)
                else:
                    sample_cov[chrom] = [cov]

    cov_stats = {}
    for chrom in sample_cov:
        mean_cov = statistics.mean(sample_cov[chrom])
        try:
            std_cov = statistics.stdev(sample_cov[chrom])
        except statistics.StatisticsError:
            std_cov = "NA"
        cov_stats[chrom] = {"mean": mean_cov, "stdev": std_cov}

    return cov_stats



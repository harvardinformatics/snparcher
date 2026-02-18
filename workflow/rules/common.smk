import os
from pathlib import Path
import pandas as pd
from snakemake.utils import validate


# --- Config defaults and validation ---

def set_defaults(cfg, defaults):
    """Recursively apply defaults to config dict."""
    for key, value in defaults.items():
        if key not in cfg:
            cfg[key] = value
        elif isinstance(value, dict) and isinstance(cfg[key], dict):
            set_defaults(cfg[key], value)


DEFAULTS = {
    "samples": "config/samples.csv",
    "variant_calling": {
        "expected_coverage": "low",
        "tool": "gatk",
        "ploidy": 2,
        "gatk": {
            "het_prior": 0.005,
        },
        "sentieon": {
            "license": "",
        },
    },
    "intervals": {
        "enabled": True,
        "min_nmer": 500,
        "num_gvcf_intervals": 50,
        "db_scatter_factor": 0.15,
    },
    "callable_sites": {
        "generate_bed_file": True,
        "coverage": {
            "enabled": True,
            "stdev": 2,
            "merge_distance": 100,
        },
        "mappability": {
            "enabled": True,
            "kmer": 150,
            "min_score": 1,
            "merge_distance": 100,
        },
    },
    "modules": {
        "qc": {
            "enabled": False,
            "clusters": 3,
            "min_depth": 2,
            "google_api_key": "",
        },
        "postprocess": {
            "enabled": False,
            "filtering": {
                "contig_size": 10000,
                "maf": 0.01,
                "missingness": 0.75,
                "exclude_scaffolds": "mtDNA,Y",
            },
        },
        "trackhub": {
            "enabled": False,
            "email": "example@email.com",
        },
    },
}

set_defaults(config, DEFAULTS)
validate(config, Path(workflow.basedir, "schemas/config.schema.yaml"))

USE_SENTIEON = config["variant_calling"]["tool"] == "sentieon"


# --- Sample sheet loading and validation ---

def _parse_mark_duplicates(values):
    truthy = {"true", "t", "yes", "y", "1"}
    falsy = {"false", "f", "no", "n", "0"}
    parsed = []
    invalid = []

    for idx, value in values.items():
        if pd.isna(value):
            parsed.append(True)
            continue

        if isinstance(value, bool):
            parsed.append(value)
            continue

        if isinstance(value, (int, float)):
            if value == 1:
                parsed.append(True)
                continue
            if value == 0:
                parsed.append(False)
                continue

        value_str = str(value).strip().lower()
        if value_str in truthy:
            parsed.append(True)
        elif value_str in falsy:
            parsed.append(False)
        else:
            invalid.append((idx, value))

    if invalid:
        idx, value = invalid[0]
        raise ValueError(
            f"Invalid mark_duplicates value at row {idx + 2}: {value!r}. "
            "Use true/false (or equivalent boolean values)."
        )

    return pd.Series(parsed, index=values.index, dtype=bool)


samples_df = pd.read_csv(config["samples"])

if "library_id" not in samples_df.columns:
    samples_df["library_id"] = samples_df["sample_id"]
else:
    samples_df["library_id"] = (
        samples_df["library_id"]
        .replace(r"^\s*$", pd.NA, regex=True)
        .fillna(samples_df["sample_id"])
    )

if "mark_duplicates" not in samples_df.columns:
    samples_df["mark_duplicates"] = True
else:
    samples_df["mark_duplicates"] = _parse_mark_duplicates(samples_df["mark_duplicates"])

validate(samples_df, Path(workflow.basedir, "schemas/samples.schema.yaml"))

for sample_id, group in samples_df.groupby("sample_id"):
    input_types = group["input_type"].dropna().unique().tolist()
    if len(input_types) > 1:
        raise ValueError(
            f"Sample '{sample_id}' has mixed input_type values: {input_types}"
        )

    mark_dups = group["mark_duplicates"].dropna().unique().tolist()
    if len(mark_dups) > 1:
        raise ValueError(
            f"Sample '{sample_id}' has mixed mark_duplicates values: {mark_dups}"
        )

    if len(group) > 1:
        library_ids = group["library_id"].tolist()
        if len(library_ids) != len(set(library_ids)):
            raise ValueError(
                f"Sample '{sample_id}' has duplicate library_id values: {library_ids}"
            )

samples_df = samples_df.set_index("sample_id", drop=False)


# --- Reference paths ---

REF_NAME = config["reference"]["name"]

REF_FILES = {
    "ref": f"results/reference/{REF_NAME}.fa.gz",
    "ref_fai": f"results/reference/{REF_NAME}.fa.gz.fai",
    "ref_dict": f"results/reference/{REF_NAME}.dict",
}

REF_BWA_IDX = multiext(
    f"results/reference/{REF_NAME}.fa.gz",
    ".sa", ".pac", ".bwt", ".ann", ".amb",
)


# --- Sample lists ---

SAMPLES_ALL = samples_df["sample_id"].unique().tolist()
SAMPLES_SRR = samples_df[samples_df["input_type"] == "srr"]["sample_id"].unique().tolist()
SAMPLES_NEED_ALIGNMENT = samples_df[
    samples_df["input_type"].isin(["srr", "fastq"])
]["sample_id"].unique().tolist()
SAMPLES_NEED_HAPLOTYPECALLER = samples_df[
    samples_df["input_type"] != "gvcf"
]["sample_id"].unique().tolist()
SAMPLES_WITH_BAM = samples_df[
    samples_df["input_type"].isin(["srr", "fastq", "bam"])
]["sample_id"].unique().tolist()
SAMPLES_WITH_FASTQ = samples_df[
    samples_df["input_type"].isin(["srr", "fastq"])
]["sample_id"].unique().tolist()


# --- Helper functions ---

def sample_has_multiple_libraries(sample):
    """Check if sample has multiple libraries."""
    return len(samples_df[samples_df["sample_id"] == sample]) > 1


def get_sample_libraries(sample):
    """Get all library IDs for a sample."""
    return samples_df[samples_df["sample_id"] == sample]["library_id"].tolist()


def get_final_bam(sample):
    """Get final BAM path for a sample."""
    sample_rows = samples_df[samples_df["sample_id"] == sample]
    input_type = sample_rows["input_type"].iloc[0]

    if input_type == "bam":
        return sample_rows["input"].iloc[0]

    mark_dups = sample_rows["mark_duplicates"].iloc[0]

    if mark_dups:
        return f"results/bams/markdup/{sample}.bam"
    elif sample_has_multiple_libraries(sample):
        return f"results/bams/merged/{sample}.bam"
    else:
        lib = get_sample_libraries(sample)[0]
        return f"results/bams/raw/{sample}/{lib}.bam"


def get_final_gvcf(sample):
    """Get final gVCF path for a sample."""
    sample_rows = samples_df[samples_df["sample_id"] == sample]
    input_type = sample_rows["input_type"].iloc[0]

    if input_type == "gvcf":
        return sample_rows["input"].iloc[0]

    result = f"results/gvcfs/{sample}.g.vcf.gz"
    return result

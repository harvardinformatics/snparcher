import os
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
    "sentieon": {
        "enabled": False,
        "license": "",
    },
    "variant_calling": {
        "expected_coverage": "low",
        "ploidy": 2,
        "gatk": {
            "het_prior": 0.005,
        },
    },
    "intervals": {
        "enabled": True,
        "min_nmer": 500,
        "num_gvcf_intervals": 50,
        "db_scatter_factor": 0.15,
    },
    "callable_sites": {
        "mappability": {
            "kmer": 150,
            "min_score": 1,
            "merge_distance": 100,
        },
    },
}

set_defaults(config, DEFAULTS)
validate(config, Path(workflow.basedir, "schemas/config.schema.yaml"))

USE_SENTIEON = config["sentieon"]["enabled"]


# --- Sample sheet loading and validation ---

samples_df = pd.read_csv(config["samples"])
validate(samples_df, Path(workflow.basedir, "schemas/samples.schema.yaml"))

print(f"After validate - samples_df:\n{samples_df}")
print(f"columns: {samples_df.columns.tolist()}")
print(f"input_type column: {samples_df['input_type'].tolist()}")

samples_df["library_id"] = samples_df["library_id"].fillna(samples_df["sample_id"])
if "mark_duplicates" not in samples_df.columns:
    samples_df["mark_duplicates"] = True
else:
    samples_df["mark_duplicates"] = samples_df["mark_duplicates"].fillna(True).astype(bool)

for sample_id, group in samples_df.groupby("sample_id"):
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

SAMPLES_ALL = samples_df.index.tolist()
SAMPLES_SRR = samples_df[samples_df["input_type"] == "srr"].index.tolist()
SAMPLES_NEED_ALIGNMENT = samples_df[samples_df["input_type"].isin(["srr", "fastq"])].index.tolist()
SAMPLES_NEED_HAPLOTYPECALLER = samples_df[samples_df["input_type"] != "gvcf"].index.tolist()
SAMPLES_WITH_BAM = samples_df[samples_df["input_type"].isin(["srr", "fastq", "bam"])].index.tolist()
SAMPLES_WITH_FASTQ = samples_df[samples_df["input_type"].isin(["srr", "fastq"])].index.tolist()


# --- Helper functions ---

def sample_has_multiple_libraries(sample):
    """Check if sample has multiple libraries."""
    return len(samples_df[samples_df["sample_id"] == sample]) > 1


def get_sample_libraries(sample):
    """Get all library IDs for a sample."""
    return samples_df[samples_df["sample_id"] == sample]["library_id"].tolist()


def get_final_bam(sample):
    """Get final BAM path for a sample."""
    input_type = samples_df.loc[sample, "input_type"]

    if input_type == "bam":
        return samples_df.loc[sample, "input"]

    mark_dups = samples_df.loc[sample, "mark_duplicates"]

    if mark_dups:
        return f"results/bams/markdup/{sample}.bam"
    elif sample_has_multiple_libraries(sample):
        return f"results/bams/merged/{sample}.bam"
    else:
        lib = get_sample_libraries(sample)[0]
        return f"results/bams/raw/{sample}/{lib}.bam"


def get_final_gvcf(sample):
    """Get final gVCF path for a sample."""
    input_type = samples_df.loc[sample, "input_type"]
    print(f"  inside get_final_gvcf: sample={sample}, input_type={input_type}, type={type(input_type)}")

    if input_type == "gvcf":
        return samples_df.loc[sample, "input"]

    result = f"results/gvcfs/{sample}.g.vcf.gz"
    print(f"  returning: {result}")
    return result
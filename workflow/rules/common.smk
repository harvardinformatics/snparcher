import os
from pathlib import Path
import pandas as pd
from snakemake.utils import validate
from snakemake.logging import logger


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
    "sample_metadata": "",
    "variant_calling": {
        "expected_coverage": "low",
        "tool": "gatk",
        "ploidy": 2,
        "gatk": {
            "het_prior": 0.005,
            "concat_batch_size": 250,
            "concat_max_rounds": 20,
        },
        "sentieon": {
            "license": "",
        },
        "bcftools": {
            "min_mapq": 20,
            "min_baseq": 20,
            "max_depth": 250,
        },
        "deepvariant": {
            "model_type": "WGS",
            "num_shards": 8,
        },
        "parabricks": {
            "container_image": "",
            "num_gpus": 1,
            "num_cpu_threads": 16,
            "extra_args": "",
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

VARIANT_TOOL = config["variant_calling"]["tool"]
USE_SENTIEON = VARIANT_TOOL == "sentieon"


# Canonical hard-filter definitions used by workflow/rules/variant_calling/hard_filters.smk.
GATK_HARD_FILTERS = [
    (
        "RPRS_filter",
        "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || "
        "((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || "
        "(vc.hasAttribute('QD') && QD < 2.0)",
    ),
    (
        "FS_SOR_filter",
        "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || "
        "((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))",
    ),
    (
        "MQ_filter",
        "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))",
    ),
    ("QUAL_filter", "QUAL < 30.0"),
]


def get_gatk_hard_filter_args():
    args = []
    for name, expr in GATK_HARD_FILTERS:
        args.extend(
            [
                f'--filter-name "{name}"',
                f'--filter-expression "{expr}"',
            ]
        )
    return " ".join(args)


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
    
    if input_types and input_types[0] in {"bam", "gvcf"} and len(group) > 1:
        raise ValueError(
            f"Sample '{sample_id}' has input_type '{input_types[0]}' but multiple rows. "
            "Only one row is supported for bam/gvcf inputs."
        )

if VARIANT_TOOL in {"bcftools", "deepvariant", "parabricks"}:
    gvcf_samples = (
        samples_df[samples_df["input_type"] == "gvcf"]["sample_id"]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    if gvcf_samples:
        raise ValueError(
            f"variant_calling.tool '{VARIANT_TOOL}' does not support samples with input_type='gvcf'. "
            f"Incompatible samples: {', '.join(gvcf_samples)}. "
            "Use 'gatk' or 'sentieon', or provide FASTQ/BAM inputs for all samples."
        )

if VARIANT_TOOL == "parabricks":
    image = config["variant_calling"]["parabricks"]["container_image"].strip()
    if not image:
        raise ValueError(
            "variant_calling.parabricks.container_image is required when variant_calling.tool='parabricks'."
        )

samples_df["input_unit"] = (
    samples_df.groupby(["sample_id", "library_id"]).cumcount().add(1).map(lambda x: f"u{x}")
)


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
    sample_rows = samples_df[samples_df["sample_id"] == sample]
    return sample_rows["library_id"].nunique() > 1


def get_sample_libraries(sample):
    """Get all library IDs for a sample."""
    sample_rows = samples_df[samples_df["sample_id"] == sample]
    return sample_rows["library_id"].drop_duplicates().tolist()


def get_sample_rows(sample):
    """Get all rows for a sample."""
    sample_rows = samples_df[samples_df["sample_id"] == sample]
    if sample_rows.empty:
        raise ValueError(f"Sample '{sample}' not found in sample sheet")
    return sample_rows


def get_sample_scalar(sample, column):
    """Get a scalar metadata value that must be constant within a sample."""
    sample_rows = get_sample_rows(sample)
    values = sample_rows[column].dropna().unique().tolist()
    if len(values) != 1:
        raise ValueError(
            f"Sample '{sample}' has non-unique values for '{column}': {values}"
        )
    return values[0]


def get_sample_input_type(sample):
    """Get sample input_type."""
    return get_sample_scalar(sample, "input_type")


def get_sample_mark_duplicates(sample):
    """Get sample mark_duplicates flag."""
    return bool(get_sample_scalar(sample, "mark_duplicates"))


def get_sample_inputs(sample):
    """Get sample rows as records."""
    return get_sample_rows(sample).to_dict("records")


def get_sample_srr_records(sample):
    """Get SRR input records for a sample."""
    rows = get_sample_rows(sample)
    rows = rows[rows["input_type"] == "srr"]
    return rows[["sample_id", "library_id", "input_unit", "input"]].to_dict("records")


def get_final_bam(sample):
    """Get final BAM path for a sample."""
    sample_rows = get_sample_rows(sample)
    input_type = get_sample_input_type(sample)

    if input_type == "bam":
        return sample_rows["input"].iloc[0]

    mark_dups = get_sample_mark_duplicates(sample)

    if mark_dups:
        return f"results/bams/markdup/{sample}.bam"
    else:
        return f"results/bams/merged/{sample}.bam"


def get_final_gvcf(sample):
    """Get final gVCF path for a sample."""
    sample_rows = get_sample_rows(sample)
    input_type = get_sample_input_type(sample)

    if input_type == "gvcf":
        return sample_rows["input"].iloc[0]

    result = f"results/gvcfs/{sample}.g.vcf.gz"
    return result


# --- Sample metadata loading (optional) ---

def _parse_bool_column(series, column_name):
    """Parse a boolean column from metadata CSV, tolerating string representations."""
    truthy = {"true", "t", "yes", "y", "1"}
    falsy = {"false", "f", "no", "n", "0", ""}
    parsed = []
    for idx, value in series.items():
        if pd.isna(value):
            parsed.append(False)
        elif isinstance(value, bool):
            parsed.append(value)
        elif isinstance(value, (int, float)):
            parsed.append(bool(value))
        else:
            v = str(value).strip().lower()
            if v in truthy:
                parsed.append(True)
            elif v in falsy:
                parsed.append(False)
            else:
                raise ValueError(
                    f"Invalid {column_name} value at row {idx + 2}: {value!r}. "
                    "Use true/false (or equivalent boolean values)."
                )
    return pd.Series(parsed, index=series.index, dtype=bool)


metadata_df = None

if config.get("sample_metadata"):
    metadata_df = pd.read_csv(config["sample_metadata"])
    validate(metadata_df, Path(workflow.basedir, "schemas/sample_metadata.schema.yaml"))

    # Validate sample_id values exist in the main sample sheet
    unknown = set(metadata_df["sample_id"]) - set(samples_df["sample_id"])
    if unknown:
        raise ValueError(
            f"sample_metadata contains sample_id(s) not in the samples CSV: {sorted(unknown)}"
        )

    # Parse boolean columns if present
    for col in ("exclude", "outgroup"):
        if col in metadata_df.columns:
            metadata_df[col] = _parse_bool_column(metadata_df[col], col)

    # Warn about samples without metadata
    missing = set(samples_df["sample_id"]) - set(metadata_df["sample_id"])
    if missing:
        logger.warning(
            f"Samples in sample sheet but not in sample_metadata: {sorted(missing)}. "
            "These samples will use default metadata values."
        )


def get_excluded_samples():
    """Return list of sample_ids where exclude=True. Empty list if no metadata."""
    if metadata_df is None or "exclude" not in metadata_df.columns:
        return []
    return metadata_df.loc[metadata_df["exclude"], "sample_id"].tolist()


def get_included_samples():
    """Return list of sample_ids NOT marked for exclusion."""
    excluded = set(get_excluded_samples())
    return [s for s in SAMPLES_ALL if s not in excluded]


def get_outgroup_samples():
    """Return list of sample_ids where outgroup=True. Empty list if no metadata."""
    if metadata_df is None or "outgroup" not in metadata_df.columns:
        return []
    return metadata_df.loc[metadata_df["outgroup"], "sample_id"].tolist()


def get_ingroup_samples():
    """Return list of sample_ids that are not outgroup and not excluded."""
    excluded = set(get_excluded_samples())
    outgroup = set(get_outgroup_samples())
    return [s for s in SAMPLES_ALL if s not in excluded and s not in outgroup]


def get_sample_coords():
    """Return DataFrame with sample_id, lat, long for samples that have coordinates.
    Returns empty DataFrame if no metadata or no lat/long columns."""
    if metadata_df is None:
        return pd.DataFrame(columns=["sample_id", "lat", "long"])
    if "lat" not in metadata_df.columns or "long" not in metadata_df.columns:
        return pd.DataFrame(columns=["sample_id", "lat", "long"])
    coords = metadata_df[["sample_id", "lat", "long"]].dropna(subset=["lat", "long"])
    return coords


def get_sample_metadata(sample_id, column, default=None):
    """Get a metadata value for a specific sample. Returns default if not available."""
    if metadata_df is None or column not in metadata_df.columns:
        return default
    rows = metadata_df.loc[metadata_df["sample_id"] == sample_id, column]
    if rows.empty:
        return default
    return rows.iloc[0]


def require_metadata(module_name, required_columns=None):
    """Validate that metadata is available. Called by modules that need it.
    Raises a clear error if metadata CSV is not configured."""
    if metadata_df is None:
        raise ValueError(
            f"Module '{module_name}' requires sample_metadata to be configured. "
            f"Set 'sample_metadata' in your config.yaml to a CSV file path."
        )
    if required_columns:
        missing = [c for c in required_columns if c not in metadata_df.columns]
        if missing:
            raise ValueError(
                f"Module '{module_name}' requires column(s) {missing} in sample_metadata CSV."
            )


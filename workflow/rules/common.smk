from importlib.metadata import PackageNotFoundError, version
import os
from pathlib import Path
import pandas as pd
import yaml
from jsonschema import Draft202012Validator
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate


# --- Config defaults and validation ---

def warn_if_old_snakemake_version():
    """Warn when running on Snakemake versions older than the supported v9 series."""
    try:
        raw_version = version("snakemake")
    except PackageNotFoundError:
        return

    parts = raw_version.split(".")
    try:
        major = int(parts[0])
    except (TypeError, ValueError, IndexError):
        return

    if major < 9:
        logger.warning(
            f"snpArcher v2 is tested with Snakemake >=9, but this run is using Snakemake {raw_version}. "
            "Checkpoint behavior and other workflow features may be unreliable on older versions."
        )

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
    "reads": {
        "mark_duplicates": True,
    },
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
        "mk": {
            "enabled": False,
            "gff": "",
        },
    },
}

V1_CONFIG_MARKERS = (
    "final_prefix",
    "resource_config",
    "intervals",
    "sentieon",
    "sentieon_lic",
    "remote_reads_prefix",
    "bigtmp",
    "cov_filter",
    "generate_trackhub",
    "trackhub_email",
    "mark_duplicates",
    "sort_reads",
    "refGenome",
    "refPath",
    "minNmer",
    "num_gvcf_intervals",
    "db_scatter_factor",
    "ploidy",
    "minP",
    "minD",
    "het_prior",
    "mappability_min",
    "mappability_k",
    "mappability_merge",
    "cov_merge",
    "cov_threshold",
    "cov_threshold_lower",
    "cov_threshold_upper",
    "cov_threshold_stdev",
    "cov_threshold_rel",
    "callable_merge",
    "nClusters",
    "GoogleAPIKey",
    "contig_size",
    "maf",
    "missingness",
    "scaffolds_to_exclude",
)


def _config_has_v2_blocks(cfg):
    return (
        isinstance(cfg.get("reference"), dict)
        and isinstance(cfg.get("variant_calling"), dict)
    )


def detect_v1_config_markers(cfg):
    """Return legacy v1-style keys when config does not have the v2 block structure."""
    found = [key for key in V1_CONFIG_MARKERS if key in cfg]

    if isinstance(cfg.get("remote_reads"), bool):
        found.append("remote_reads")

    if not found or _config_has_v2_blocks(cfg):
        return []

    return sorted(set(found))


def normalize_supported_config_aliases(cfg):
    """Normalize a few deprecated aliases for configs that are otherwise v2-shaped."""
    reference_cfg = cfg.get("reference")
    if isinstance(reference_cfg, dict) and "path" in reference_cfg:
        if "source" not in reference_cfg:
            logger.warning(
                "Config key 'reference.path' is deprecated; using it as 'reference.source'."
            )
            reference_cfg["source"] = reference_cfg["path"]
        else:
            logger.warning(
                "Ignoring deprecated config key 'reference.path' because 'reference.source' is set."
            )
        reference_cfg.pop("path", None)

    variant_cfg = cfg.get("variant_calling")
    if not isinstance(variant_cfg, dict):
        return

    gatk_cfg = variant_cfg.get("gatk")
    if isinstance(gatk_cfg, dict) and "ploidy" in gatk_cfg:
        if "ploidy" not in variant_cfg:
            logger.warning(
                "Config key 'variant_calling.gatk.ploidy' is deprecated; using it as 'variant_calling.ploidy'."
            )
            variant_cfg["ploidy"] = gatk_cfg["ploidy"]
        else:
            logger.warning(
                "Ignoring deprecated config key 'variant_calling.gatk.ploidy' because 'variant_calling.ploidy' is set."
            )
        gatk_cfg.pop("ploidy", None)


def _load_yaml(path):
    with open(path) as handle:
        return yaml.safe_load(handle)


def _schema_without_additional_properties(schema):
    if isinstance(schema, dict):
        relaxed = {}
        for key, value in schema.items():
            if key == "additionalProperties":
                continue
            relaxed[key] = _schema_without_additional_properties(value)
        return relaxed
    if isinstance(schema, list):
        return [_schema_without_additional_properties(item) for item in schema]
    return schema


def _collect_unknown_config_keys(cfg, schema, prefix=""):
    if not isinstance(cfg, dict) or not isinstance(schema, dict):
        return []

    properties = schema.get("properties", {})
    unknown = []
    for key, value in cfg.items():
        path = f"{prefix}.{key}" if prefix else key
        if key not in properties:
            unknown.append(path)
            continue
        unknown.extend(_collect_unknown_config_keys(value, properties[key], path))
    return unknown


def _format_validation_error(error):
    path = ".".join(str(part) for part in error.absolute_path) or "<root>"
    return f"- {path}: {error.message}"


def validate_config_with_warnings(cfg, schema_path):
    raw_schema = _load_yaml(schema_path)
    unknown_keys = sorted(set(_collect_unknown_config_keys(cfg, raw_schema)))
    if unknown_keys:
        logger.warning(
            f"Ignoring unsupported config key(s): {', '.join(unknown_keys)}"
        )

    relaxed_schema = _schema_without_additional_properties(raw_schema)
    validator = Draft202012Validator(relaxed_schema)
    errors = sorted(
        validator.iter_errors(cfg),
        key=lambda error: [str(part) for part in error.absolute_path],
    )
    if errors:
        formatted = [_format_validation_error(error) for error in errors[:10]]
        if len(errors) > 10:
            formatted.append(f"... and {len(errors) - 10} more validation error(s)")
        raise WorkflowError(
            "Error validating config file.\n"
            + "\n".join(formatted)
        )


warn_if_old_snakemake_version()

CONFIG_SCHEMA_PATH = Path(workflow.basedir, "schemas/config.schema.yaml")
v1_markers = detect_v1_config_markers(config)
if v1_markers:
    markers_text = ", ".join(v1_markers[:8])
    if len(v1_markers) > 8:
        markers_text += ", ..."
    raise WorkflowError(
        "Detected a v1-style snpArcher config. "
        "This workflow expects the v2 config structure and will not auto-migrate v1 configs.\n"
        f"Legacy keys detected: {markers_text}\n"
        "See docs/v2-migration.md for the v1 -> v2 mapping."
    )

normalize_supported_config_aliases(config)
set_defaults(config, DEFAULTS)
validate_config_with_warnings(config, CONFIG_SCHEMA_PATH)

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

def _parse_mark_duplicates(values, default):
    truthy = {"true", "t", "yes", "y", "1"}
    falsy = {"false", "f", "no", "n", "0"}
    parsed = []
    invalid = []

    for idx, value in values.items():
        if pd.isna(value):
            parsed.append(default)
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
global_mark_duplicates = bool(config["reads"]["mark_duplicates"])

if "library_id" not in samples_df.columns:
    samples_df["library_id"] = samples_df["sample_id"]
else:
    samples_df["library_id"] = (
        samples_df["library_id"]
        .replace(r"^\s*$", pd.NA, regex=True)
        .fillna(samples_df["sample_id"])
    )

if "mark_duplicates" not in samples_df.columns:
    samples_df["mark_duplicates"] = global_mark_duplicates
else:
    samples_df["mark_duplicates"] = _parse_mark_duplicates(
        samples_df["mark_duplicates"],
        default=global_mark_duplicates,
    )

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
SAMPLES_WITH_BAM = samples_df[
    samples_df["input_type"].isin(["srr", "fastq", "bam"])
]["sample_id"].unique().tolist()
SAMPLES_WITH_FASTQ = samples_df[
    samples_df["input_type"].isin(["srr", "fastq"])
]["sample_id"].unique().tolist()

if config["callable_sites"]["coverage"]["enabled"] and not SAMPLES_WITH_BAM:
    raise ValueError(
        "callable_sites.coverage.enabled requires at least one BAM-backed sample. "
        "Disable callable_sites.coverage.enabled for gVCF-only workflows."
    )


# --- Helper functions ---

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


def get_joint_gvcf_records():
    """Return sample IDs paired with their final gVCF paths."""
    return [(sample, get_final_gvcf(sample)) for sample in SAMPLES_ALL]


def get_joint_gvcf_paths():
    """Return final gVCF paths for joint genotyping."""
    return [gvcf for _, gvcf in get_joint_gvcf_records()]


def get_joint_gvcf_tbis():
    """Return final gVCF index paths for joint genotyping."""
    return [f"{gvcf}.tbi" for gvcf in get_joint_gvcf_paths()]


def write_joint_gvcf_mapfile(path):
    """Write a GenomicsDB sample-name map keyed by sample_id."""
    with open(path, "w") as handle:
        for sample, gvcf in get_joint_gvcf_records():
            print(sample, gvcf, sep="\t", file=handle)


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

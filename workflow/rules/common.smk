from importlib.metadata import PackageNotFoundError, version
import logging
import os
import platform
import pwd
import re
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import snakemake
import yaml
from jsonschema import Draft202012Validator
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate


# --- Config defaults and validation ---

def warn_if_old_snakemake_version():
    """Warn when running on Snakemake versions older than the supported v9 series."""
    raw_version = getattr(snakemake, "__version__", None)
    if raw_version is None:
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
            "fraction": 1.0,
            "min_coverage": "auto",
            "max_coverage": "auto",
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
            "max_sample_missingness": 0.49,
            "google_api_key": "",
            "exclude_scaffolds": "",
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
    },
}

REMOVED_MODULES = ("mk", "trackhub")

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


def remove_unsupported_module_blocks(cfg):
    """Warn about removed module blocks and drop them before validation."""
    modules_cfg = cfg.get("modules")
    if not isinstance(modules_cfg, dict):
        return

    for module_name in REMOVED_MODULES:
        if module_name in modules_cfg:
            logger.warning(
                f"Config block 'modules.{module_name}' was removed in snpArcher v2 and will be ignored."
            )
            modules_cfg.pop(module_name, None)


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
remove_unsupported_module_blocks(config)
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


SAMPLE_SHEET_DTYPES = {
    "sample_id": "string",
    "library_id": "string",
}

SAMPLE_METADATA_DTYPES = {
    "sample_id": "string",
}


samples_df = pd.read_csv(config["samples"], dtype=SAMPLE_SHEET_DTYPES)
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

MIXED_READ_INPUT_TYPES = {"fastq", "srr"}
SINGLE_ROW_INPUT_TYPES = {"bam", "gvcf"}

for sample_id, group in samples_df.groupby("sample_id"):
    input_types = sorted(group["input_type"].dropna().unique().tolist())
    input_type_set = set(input_types)

    if len(input_type_set) > 1 and input_type_set != MIXED_READ_INPUT_TYPES:
        raise ValueError(
            f"Sample '{sample_id}' has unsupported mixed input_type values: {input_types}. "
            "Only mixing 'srr' and 'fastq' is supported within one sample."
        )

    mark_dups = group["mark_duplicates"].dropna().unique().tolist()
    if len(mark_dups) > 1:
        raise ValueError(
            f"Sample '{sample_id}' has mixed mark_duplicates values: {mark_dups}"
        )

    for input_type in SINGLE_ROW_INPUT_TYPES:
        if input_type in input_type_set and len(group) > 1:
            raise ValueError(
                f"Sample '{sample_id}' has input_type '{input_type}' but multiple rows. "
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

CALLABLE_COVERAGE_ENABLED = bool(config["callable_sites"]["coverage"]["enabled"])
CALLABLE_MAPPABILITY_ENABLED = bool(config["callable_sites"]["mappability"]["enabled"])
CALLABLE_GENERATE_BED_FILE = bool(config["callable_sites"]["generate_bed_file"])

if CALLABLE_GENERATE_BED_FILE and not (
    CALLABLE_COVERAGE_ENABLED or CALLABLE_MAPPABILITY_ENABLED
):
    logger.warning(
        "callable_sites.generate_bed_file is true, but both "
        "callable_sites.coverage.enabled and callable_sites.mappability.enabled are false. "
        "Skipping results/callable_sites/callable_sites.bed generation."
    )
    CALLABLE_GENERATE_BED_FILE = False

CALLABLE_FINAL_BED_ENABLED = CALLABLE_GENERATE_BED_FILE and (
    CALLABLE_COVERAGE_ENABLED or CALLABLE_MAPPABILITY_ENABLED
)

POSTPROCESS_ENABLED = bool(config["modules"]["postprocess"]["enabled"])
if POSTPROCESS_ENABLED and not CALLABLE_FINAL_BED_ENABLED:
    logger.warning(
        "modules.postprocess.enabled is true, but results/callable_sites/callable_sites.bed "
        "will not be generated. Disabling postprocess."
    )
    POSTPROCESS_ENABLED = False

CALLABLE_TARGETS = (
    ["results/callable_sites/callable_sites.bed"]
    if CALLABLE_FINAL_BED_ENABLED
    else [
        *(
            ["results/callable_sites/callable_loci.zarr"]
            if CALLABLE_COVERAGE_ENABLED
            else []
        ),
        *(
            ["results/callable_sites/mappability.bed"]
            if CALLABLE_MAPPABILITY_ENABLED
            else []
        ),
    ]
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


def get_sample_input_types(sample):
    """Get unique input_type values for a sample."""
    sample_rows = get_sample_rows(sample)
    return sample_rows["input_type"].dropna().unique().tolist()


def sample_has_input_type(sample, input_type):
    """Return True when a sample contains at least one row of a given input_type."""
    return input_type in get_sample_input_types(sample)


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

    if sample_has_input_type(sample, "bam"):
        bam_rows = sample_rows[sample_rows["input_type"] == "bam"]
        return bam_rows["input"].iloc[0]

    mark_dups = get_sample_mark_duplicates(sample)

    if mark_dups:
        return f"results/bams/markdup/{sample}.bam"
    else:
        return f"results/bams/merged/{sample}.bam"


def get_final_gvcf(sample):
    """Get final gVCF path for a sample."""
    sample_rows = get_sample_rows(sample)

    if sample_has_input_type(sample, "gvcf"):
        gvcf_rows = sample_rows[sample_rows["input_type"] == "gvcf"]
        return gvcf_rows["input"].iloc[0]

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
    metadata_df = pd.read_csv(
        config["sample_metadata"],
        dtype=SAMPLE_METADATA_DTYPES,
    )
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


# --- Debug banner ---


def _get_snparcher_version():
    """Get snpArcher version from pyproject.toml, with optional git hash."""
    import tomllib

    ver = "unknown"
    toml_path = Path(workflow.basedir).parent / "pyproject.toml"
    try:
        with open(toml_path, "rb") as f:
            ver = tomllib.load(f)["project"]["version"]
    except Exception:
        pass
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True,
            cwd=Path(workflow.basedir).parent,
        )
        git_hash = result.stdout.strip()
        if git_hash:
            ver += "-" + git_hash
    except Exception:
        pass
    return ver


def print_debug_banner():
    """Print snpArcher startup banner with environment and config info."""
    indent = 24

    # Environment info
    username = pwd.getpwuid(os.getuid())[0]
    hostname = socket.gethostname()
    if platform.node() != hostname:
        hostname += "; " + platform.node()

    try:
        conda_proc = subprocess.run(
            ["conda", "--version"], capture_output=True, text=True
        )
        conda_ver = conda_proc.stdout.strip().removeprefix("conda ").strip()
        if not conda_ver:
            conda_ver = "n/a"
    except Exception:
        conda_ver = "n/a"

    conda_env_name = os.environ.get("CONDA_DEFAULT_ENV", "")
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_env_name or conda_prefix:
        conda_env = f"{conda_env_name} ({conda_prefix})"
    else:
        conda_env = "n/a"

    # Command line
    cmdline = sys.argv[0]
    for i in range(1, len(sys.argv)):
        if sys.argv[i].startswith("--"):
            cmdline += "\n" + (" " * indent) + sys.argv[i]
        else:
            cmdline += " " + sys.argv[i]

    # Profile (from --profile or --workflow-profile args)
    profile = "none"
    for i, arg in enumerate(sys.argv):
        if arg in ("--profile", "--workflow-profile") and i + 1 < len(sys.argv):
            profile = sys.argv[i + 1]
            break
        if arg.startswith("--profile="):
            profile = arg.split("=", 1)[1]
            break
        if arg.startswith("--workflow-profile="):
            profile = arg.split("=", 1)[1]
            break

    # Config files
    cfgfiles = [os.path.abspath(str(cfg)) for cfg in workflow.configfiles]
    cfgfiles_str = ("\n" + " " * indent).join(cfgfiles) if cfgfiles else "none"

    # Samples summary
    type_counts = samples_df.groupby("input_type")["sample_id"].nunique()
    samples_summary = f"{len(SAMPLES_ALL)} total"
    type_parts = [f"{count} {itype}" for itype, count in sorted(type_counts.items())]
    if type_parts:
        samples_summary += " (" + ", ".join(type_parts) + ")"

    # Thread overrides
    thread_overrides = getattr(
        workflow.resource_settings, "_parsed_overwrite_threads", {}
    )
    thread_str = (
        ", ".join(f"{k}={v}" for k, v in sorted(thread_overrides.items()))
        if thread_overrides
        else "none"
    )

    # Resource overrides
    resource_overrides = getattr(
        workflow.resource_settings, "_parsed_overwrite_resources", {}
    )
    resource_str = (
        ", ".join(
            f"{rule}: {', '.join(f'{k}={v}' for k, v in sorted(res.items()))}"
            for rule, res in sorted(resource_overrides.items())
        )
        if resource_overrides
        else "none"
    )

    logger.info("=" * 70)
    logger.info("    snpArcher " + _get_snparcher_version())
    logger.info("")
    logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("    Platform:           " + platform.platform())
    logger.info("    Host:               " + hostname)
    logger.info("    User:               " + username)
    logger.info("    Python:             " + sys.version.split(" ")[0])
    logger.info("    Snakemake:          " + str(snakemake.__version__))
    logger.info("    Conda:              " + conda_ver)
    logger.info("    Conda env:          " + conda_env)
    logger.info("")
    logger.info("    Base directory:     " + workflow.basedir)
    logger.info("    Working directory:  " + os.getcwd())
    logger.info("    Config file(s):     " + cfgfiles_str)
    logger.info("    Command:            " + cmdline)
    logger.info("    Profile:            " + profile)
    logger.info("")
    logger.info("    Reference:          " + REF_NAME + " (" + config["reference"]["source"] + ")")
    logger.info("    Variant tool:       " + VARIANT_TOOL)
    logger.info("    Ploidy:             " + str(config["variant_calling"]["ploidy"]))
    logger.info("    Expected coverage:  " + str(config["variant_calling"]["expected_coverage"]))
    logger.info("    Intervals:          " + str(config["intervals"]["enabled"]))
    logger.info("    Modules:            " + "qc=" + str(config["modules"]["qc"]["enabled"]) + ", postprocess=" + str(config["modules"]["postprocess"]["enabled"]))
    logger.info("    Samples:            " + samples_summary)
    logger.info("    Thread overrides:   " + thread_str)
    logger.info("    Resource overrides: " + resource_str)
    logger.info("")
    logger.info("=" * 70)
    logger.info("")


if workflow.is_main_process:
    print_debug_banner()

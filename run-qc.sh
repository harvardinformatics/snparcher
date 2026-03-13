#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run-qc.sh <snparcher_path> <samples_csv> <vcf_gz> <reference_fai> [--bam-stats <qc_report_tsv>] [--sample-metadata <metadata_csv>] [--cores <n>]

Description:
  Run the standalone snpArcher QC module against an existing VCF.

Required positional arguments:
  snparcher_path     Path to the snpArcher repository root
  samples_csv        Sample sheet CSV
  vcf_gz             Input bgzipped VCF
  reference_fai      Reference FASTA index (.fai) matching the VCF reference

Optional arguments:
  --bam-stats        QC metrics TSV for the mapping-rate panel
  --sample-metadata  Sample metadata CSV with optional exclude/outgroup/lat/long columns
  --cores            Number of cores to pass to Snakemake (default: 4)
  -h, --help         Show this help message
EOF
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
    usage
    exit 0
fi

if [[ $# -lt 4 ]]; then
    usage
    exit 1
fi

snparcher_path=$1
samples_csv=$2
vcf_gz=$3
reference_fai=$4
shift 4

bam_stats=""
sample_metadata=""
cores=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam-stats)
            [[ $# -ge 2 ]] || { echo "Missing value for --bam-stats" >&2; exit 1; }
            bam_stats=$2
            shift 2
            ;;
        --sample-metadata)
            [[ $# -ge 2 ]] || { echo "Missing value for --sample-metadata" >&2; exit 1; }
            sample_metadata=$2
            shift 2
            ;;
        --cores)
            [[ $# -ge 2 ]] || { echo "Missing value for --cores" >&2; exit 1; }
            cores=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

require_path() {
    local path=$1
    local label=$2
    if [[ ! -e "$path" ]]; then
        echo "$label not found: $path" >&2
        exit 1
    fi
}

require_path "$snparcher_path" "snpArcher path"
require_path "$samples_csv" "Sample sheet"
require_path "$vcf_gz" "VCF"
require_path "$reference_fai" "Reference .fai"

if [[ -n "$bam_stats" ]]; then
    require_path "$bam_stats" "BAM stats file"
fi

if [[ -n "$sample_metadata" ]]; then
    require_path "$sample_metadata" "Sample metadata file"
fi

snakefile="$snparcher_path/workflow/modules/qc/Snakefile"
configfile="$snparcher_path/workflow/modules/qc/config/config.yaml"
profile_dir="$snparcher_path/workflow-profiles/default"

require_path "$snakefile" "QC module Snakefile"
require_path "$configfile" "QC module config file"
require_path "$profile_dir" "Workflow profile"

cmd=(
    snakemake
    -s "$snakefile"
    --configfile "$configfile"
    --config
    "samples=$samples_csv"
    "vcf=$vcf_gz"
    "fai=$reference_fai"
    "sample_metadata=$sample_metadata"
    "qc_report=$bam_stats"
    --workflow-profile "$profile_dir"
    --cores "$cores"
    --use-conda
)

printf 'Running:'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

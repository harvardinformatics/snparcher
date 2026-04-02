#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run-postprocess.sh <snparcher_path> <samples_csv> <vcf> <reference_fai> <callable_sites_bed> [--sample-metadata <metadata_csv>] [--contig-size <bp>] [--maf <float>] [--missingness <float>] [--exclude-scaffolds <csv>] [--cores <n>]

Description:
  Run the standalone snpArcher postprocess module against an existing VCF.
  Output files are written to the current working directory.

Required positional arguments:
  snparcher_path      Path to the snpArcher repository root
  samples_csv         Sample sheet CSV
  vcf                 Input VCF (plain or bgzipped)
  reference_fai       Reference FASTA index (.fai) matching the VCF reference
  callable_sites_bed  Callable-sites BED used to restrict the filtered outputs

Optional arguments:
  --sample-metadata   Sample metadata CSV with optional exclude column
  --contig-size       Contig size threshold override for postprocess filtering
  --maf               MAF threshold override for postprocess filtering
  --missingness       Missingness threshold override for postprocess filtering
  --exclude-scaffolds Comma-separated scaffold list override (use "" to disable)
  --cores             Number of cores to pass to Snakemake (default: 4)
  -h, --help          Show this help message
EOF
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
    usage
    exit 0
fi

if [[ $# -lt 5 ]]; then
    usage
    exit 1
fi

snparcher_path=$1
samples_csv=$2
vcf=$3
reference_fai=$4
callable_sites_bed=$5
shift 5

sample_metadata=""
sample_metadata_set=false
contig_size=""
contig_size_set=false
maf=""
maf_set=false
missingness=""
missingness_set=false
exclude_scaffolds=""
exclude_scaffolds_set=false
cores=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-metadata)
            [[ $# -ge 2 ]] || { echo "Missing value for --sample-metadata" >&2; exit 1; }
            sample_metadata=$2
            sample_metadata_set=true
            shift 2
            ;;
        --contig-size)
            [[ $# -ge 2 ]] || { echo "Missing value for --contig-size" >&2; exit 1; }
            contig_size=$2
            contig_size_set=true
            shift 2
            ;;
        --maf)
            [[ $# -ge 2 ]] || { echo "Missing value for --maf" >&2; exit 1; }
            maf=$2
            maf_set=true
            shift 2
            ;;
        --missingness)
            [[ $# -ge 2 ]] || { echo "Missing value for --missingness" >&2; exit 1; }
            missingness=$2
            missingness_set=true
            shift 2
            ;;
        --exclude-scaffolds)
            [[ $# -ge 2 ]] || { echo "Missing value for --exclude-scaffolds" >&2; exit 1; }
            exclude_scaffolds=$2
            exclude_scaffolds_set=true
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
require_path "$vcf" "VCF"
require_path "$reference_fai" "Reference .fai"
require_path "$callable_sites_bed" "Callable-sites BED"

if [[ "$sample_metadata_set" == true ]]; then
    require_path "$sample_metadata" "Sample metadata file"
fi

snakefile="$snparcher_path/workflow/modules/postprocess/Snakefile"
configfile="$snparcher_path/workflow/modules/postprocess/config/config.yaml"
profile_dir="$snparcher_path/workflow-profiles/default"

require_path "$snakefile" "Postprocess module Snakefile"
require_path "$configfile" "Postprocess module config file"
require_path "$profile_dir" "Workflow profile"

config_args=(
    "samples=$samples_csv"
    "vcf=$vcf"
    "ref_fai=$reference_fai"
    "callable_sites_bed=$callable_sites_bed"
)

if [[ "$sample_metadata_set" == true ]]; then
    config_args+=("sample_metadata=$sample_metadata")
fi

if [[ "$contig_size_set" == true ]]; then
    config_args+=("contig_size=$contig_size")
fi

if [[ "$maf_set" == true ]]; then
    config_args+=("maf=$maf")
fi

if [[ "$missingness_set" == true ]]; then
    config_args+=("missingness=$missingness")
fi

if [[ "$exclude_scaffolds_set" == true ]]; then
    config_args+=("exclude_scaffolds=$exclude_scaffolds")
fi

cmd=(
    snakemake
    -s "$snakefile"
    --configfile "$configfile"
    --config
    "${config_args[@]}"
    --workflow-profile "$profile_dir"
    --cores "$cores"
    --use-conda
)

printf 'Running:'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

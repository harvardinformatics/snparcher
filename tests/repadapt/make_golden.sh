#!/usr/bin/env bash
# Produce the RepAdapt golden outputs used by snpArcher's RepAdapt tests.
#
# Runs RepAdapt's Nextflow pipeline, pinned at REPADAPT_COMMIT, twice on the
# simulated fixture in tests/data/repadapt/, checks that the two runs agree and
# that the fixture exercises both of RepAdapt's filters, then writes
# tests/data/repadapt/golden/.
#
# Usage:
#   NEXTFLOW=/path/to/nextflow tests/repadapt/make_golden.sh WORKDIR
#
# Requires Linux, apptainer or singularity, Nextflow <= 25.10.2 (the newest
# RepAdapt's README supports) running on Java 17+, curl, git and python3.
# Set JAVA_HOME if the default java is older than 17. WORKDIR holds the RepAdapt
# clone, the container images and both runs; the clone and images are reused
# when present.
set -euo pipefail

REPADAPT_URL=https://github.com/RepAdapt/nextflow_snp_calling_linux.git
REPADAPT_COMMIT=2077f6f311650c24ccc0df8c4088e616e624f5d8
DEPOT_URL=https://depot.galaxyproject.org/singularity
NEXTFLOW_MAX_VERSION=25.10.2
SAMPLES=(S1 S2 S3 S4)

# Nextflow process -> image tag, as in RepAdapt's nextflow.config. prepareDepth
# and joinDepth have no container and run on the host.
PROCESS_IMAGES=(
    "trimSequences fastp:0.20.1--h8b12597_0"
    "fastaIndex samtools:1.16.1--h6899075_0"
    "samtoolsSort samtools:1.16.1--h6899075_0"
    "samtoolsRealignedIndex samtools:1.16.1--h6899075_0"
    "samtoolsDedupIndex samtools:1.16.1--h6899075_0"
    "calculateDepth samtools:1.16.1--h6899075_0"
    "gatkIndex picard:2.26.3--hdfd78af_0"
    "addRG picard:2.26.3--hdfd78af_0"
    "dupRemoval picard:2.26.3--hdfd78af_0"
    "bwaIndex bwa:0.7.17--h5bf99c6_8"
    "bwaMap bwa:0.7.17--h5bf99c6_8"
    "realignIndel gatk:3.8--9"
    "calculateGenesDepth bedtools:2.27.1--0"
    "calculateWindowsDepth bedtools:2.27.1--0"
    "calculateWgDepth bedtools:2.27.1--0"
    "snpCalling bcftools:1.16--hfe4b78e_1"
    "concatVCFs bcftools:1.16--hfe4b78e_1"
)

die() {
    echo "make_golden.sh: $*" >&2
    exit 1
}

image_file() {
    echo "$WORK/images/${1/:/_}.sif"
}

[[ $# -eq 1 ]] || die "usage: $0 WORKDIR"
[[ $(uname -s) == Linux ]] || die "RepAdapt's images are linux-64 only"

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
FIX=$REPO/tests/data/repadapt
GOLDEN=$FIX/golden
HELPERS=$REPO/tests/repadapt
WORK=$(realpath -m "$1")
NEXTFLOW=${NEXTFLOW:-nextflow}
export NXF_HOME=${NXF_HOME:-$WORK/.nextflow}
export NXF_DISABLE_CHECK_LATEST=true

if command -v apptainer >/dev/null; then
    ENGINE=apptainer
elif command -v singularity >/dev/null; then
    ENGINE=singularity
else
    die "needs apptainer or singularity"
fi

NF_VERSION=$("$NEXTFLOW" -version 2>&1 | sed -n 's/^ *version \([0-9][0-9.]*\).*/\1/p' | head -n1)
[[ -n $NF_VERSION ]] || die "could not run '$NEXTFLOW -version'; it needs Java 17+ (set JAVA_HOME)"
[[ $(printf '%s\n' "$NF_VERSION" "$NEXTFLOW_MAX_VERSION" | sort -V | tail -n1) == "$NEXTFLOW_MAX_VERSION" ]] ||
    die "Nextflow $NF_VERSION is newer than $NEXTFLOW_MAX_VERSION, the newest RepAdapt supports"
JAVA_VERSION=$("${JAVA_HOME:+$JAVA_HOME/bin/}java" -version 2>&1 | head -n1)

for f in reference.fasta genes.gff; do
    [[ -s $FIX/$f ]] || die "missing $FIX/$f; run tests/data/repadapt/make_fixture.py"
done

mkdir -p "$WORK/images"

# --- RepAdapt, pinned -------------------------------------------------------
if [[ ! -d $WORK/repadapt/.git ]]; then
    git clone -q "$REPADAPT_URL" "$WORK/repadapt"
fi
git -C "$WORK/repadapt" checkout -q "$REPADAPT_COMMIT"
[[ $(git -C "$WORK/repadapt" rev-parse HEAD) == "$REPADAPT_COMMIT" ]] ||
    die "RepAdapt clone is not at $REPADAPT_COMMIT"

# --- Images -----------------------------------------------------------------
TAGS=$(printf '%s\n' "${PROCESS_IMAGES[@]}" | awk '{print $2}' | sort -u)
for tag in $TAGS; do
    img=$(image_file "$tag")
    if [[ ! -s $img ]]; then
        echo "Downloading $tag"
        curl -fsSL --retry 3 -o "$img.part" "$DEPOT_URL/$tag"
        mv "$img.part" "$img"
    fi
done

# --- Nextflow config --------------------------------------------------------
# Passed as -C, a top-level option that must come before 'run'. It replaces
# RepAdapt's nextflow.config, which only holds the author's local image paths
# and the container engine; CPUs are set in RepAdapt's process files.
{
    echo "process {"
    for entry in "${PROCESS_IMAGES[@]}"; do
        read -r process tag <<<"$entry"
        echo "    withName:$process { container = '$(image_file "$tag")' }"
    done
    echo "}"
    echo "$ENGINE {"
    echo "    enabled = true"
    echo "    autoMounts = true"
    echo "}"
} >"$WORK/golden.config"

NF_ARGS=(
    --ref_genome "$FIX/reference.fasta"
    --gff_file "$FIX/genes.gff"
    --reads "$FIX/fastq/*_{1,2}.fastq.gz"
)

run_repadapt() {
    local run=$WORK/run$1
    rm -rf "$run"
    mkdir -p "$run"
    echo "RepAdapt run $1 in $run"
    (
        cd "$run"
        LC_ALL=C "$NEXTFLOW" -C "$WORK/golden.config" run "$WORK/repadapt/main.nf" \
            -ansi-log false -work-dir "$run/work" "${NF_ARGS[@]}" --outdir "$run/out"
    ) >"$run/nextflow.out" 2>&1 || die "run $1 failed; see $run/nextflow.out"

    # RepAdapt uses errorStrategy 'ignore', so a failed sample is silently
    # dropped. Check every expected output exists.
    for s in "${SAMPLES[@]}"; do
        [[ -s $run/out/${s}_sorted_RG_dedup_realigned.bam ]] || die "run $1: no realigned BAM for $s"
    done
    [[ -s $run/out/final_variants.vcf.gz ]] || die "run $1: no final_variants.vcf.gz"

    # Pre-filter calls, one VCF per contig, from the snpCalling work dirs.
    local raw=()
    mapfile -t raw < <(find "$run/work" -name 'variants_chr_*.vcf.gz' | sort)
    [[ ${#raw[@]} -eq $(grep -c '^>' "$FIX/reference.fasta") ]] ||
        die "run $1: expected one variants_chr_*.vcf.gz per contig, found ${#raw[@]}"

    python3 "$HELPERS/normalize_vcf.py" --reference "$FIX/reference.fasta" "${raw[@]}" \
        -o "$run/raw_records.tsv"
    python3 "$HELPERS/normalize_vcf.py" --reference "$FIX/reference.fasta" \
        "$run/out/final_variants.vcf.gz" -o "$run/final_records.tsv"
    for s in "${SAMPLES[@]}"; do
        "$ENGINE" exec -B "$WORK" "$(image_file samtools:1.16.1--h6899075_0)" \
            samtools view "$run/out/${s}_sorted_RG_dedup_realigned.bam" >"$run/${s}.records.sam"
    done
}

run_repadapt 1
run_repadapt 2

# --- Determinism ------------------------------------------------------------
for f in raw_records.tsv final_records.tsv $(printf '%s.records.sam ' "${SAMPLES[@]}"); do
    cmp -s "$WORK/run1/$f" "$WORK/run2/$f" || die "RepAdapt is not deterministic on the fixture: $f differs"
done
echo "Runs 1 and 2 agree"

# --- Filter coverage --------------------------------------------------------
# The final records must be exactly the raw records that fail neither filter,
# with FILTER set to PASS (this also checks repadapt_failed_filters against
# bcftools 1.16), and the fixture must exercise each filter outcome.
FILTER_SUMMARY=$(
    python3 - "$HELPERS" "$WORK/run1" <<'EOF'
import sys
from collections import Counter

sys.path.insert(0, sys.argv[1])
from normalize_vcf import repadapt_failed_filters

run = sys.argv[2]
raw = [line.split("\t") for line in open(f"{run}/raw_records.tsv").read().splitlines()]
final = open(f"{run}/final_records.tsv").read().splitlines()

outcomes = Counter()
expected_final = []
for fields in raw:
    failed = repadapt_failed_filters(fields)
    outcomes[";".join(failed) or "PASS"] += 1
    if not failed:
        expected_final.append("\t".join(fields[:6] + ["PASS"] + fields[7:]))

if expected_final != final:
    sys.exit("final records are not the raw records that pass AC=AN || MQ < 30")
for outcome in ("PASS", "AllHomAlt", "LowMQ", "AllHomAlt;LowMQ"):
    if outcomes[outcome] == 0:
        sys.exit(f"fixture produced no raw records with outcome {outcome}; adjust make_fixture.py")
if outcomes["PASS"] < 20:
    sys.exit(f"only {outcomes['PASS']} PASS records; adjust make_fixture.py")
print(", ".join(f"{k}: {outcomes[k]}" for k in ("PASS", "AllHomAlt", "LowMQ", "AllHomAlt;LowMQ")))
EOF
) || die "filter coverage check failed"
echo "Raw record outcomes: $FILTER_SUMMARY"

# --- Golden outputs ---------------------------------------------------------
rm -rf "$GOLDEN"
mkdir -p "$GOLDEN"
cp "$WORK/run1/out/final_variants.vcf.gz" "$WORK/run1/raw_records.tsv" \
    "$WORK/run1/final_records.tsv" "$GOLDEN/"
for s in "${SAMPLES[@]}"; do
    cp "$WORK/run1/out/${s}_sorted_RG_dedup_realigned.bam" "$GOLDEN/"
done

{
    echo "# RepAdapt golden outputs"
    echo
    echo "Generated by \`tests/repadapt/make_golden.sh\`. Don't edit these files by hand; rerun the script."
    echo
    echo "- **Date:** $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "- **RepAdapt:** $REPADAPT_URL at \`$REPADAPT_COMMIT\`"
    echo "- **Nextflow:** $NF_VERSION ($JAVA_VERSION)"
    echo "- **Container engine:** $("$ENGINE" --version)"
    echo "- **Host:** $(uname -srm)"
    echo "- **Command** (run from an empty directory, twice; records identical across runs):"
    echo
    echo '  ```'
    echo "  LC_ALL=C nextflow -C golden.config run main.nf \\"
    echo "    --ref_genome tests/data/repadapt/reference.fasta --gff_file tests/data/repadapt/genes.gff \\"
    echo "    --reads 'tests/data/repadapt/fastq/*_{1,2}.fastq.gz' --outdir out"
    echo '  ```'
    echo
    echo "- **Raw record outcomes:** $FILTER_SUMMARY"
    echo
    echo "## Fixture inputs (sha256)"
    echo
    echo '```'
    (cd "$FIX" && sha256sum reference.fasta genes.gff fastq/*.fastq.gz)
    echo '```'
    echo
    echo "## Images"
    echo
    echo "| Image | sha256 |"
    echo "|---|---|"
    for tag in $TAGS; do
        echo "| $DEPOT_URL/$tag | \`$(sha256sum "$(image_file "$tag")" | cut -d' ' -f1)\` |"
    done
    echo
    echo "## Files"
    echo
    echo "- \`final_variants.vcf.gz\`: RepAdapt's published call set."
    echo "- \`final_records.tsv\`: its records, normalized by \`tests/repadapt/normalize_vcf.py\`."
    echo "- \`raw_records.tsv\`: RepAdapt's pre-filter calls (\`variants_chr_*.vcf.gz\` from the"
    echo "  snpCalling work directories), normalized the same way."
    echo "- \`{S}_sorted_RG_dedup_realigned.bam\`: RepAdapt's published realigned BAMs."
} >"$GOLDEN/PROVENANCE.md"

echo "Wrote $GOLDEN ($(du -sh --apparent-size "$GOLDEN" | cut -f1))"

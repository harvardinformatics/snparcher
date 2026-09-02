# --- SRA extraction tuning ---------------------------------------------------
#
# fasterq-dump is the right default: it is multi-threaded and beats fastq-dump
# by roughly 2x on ordinary runs. But for reference-aligned submissions (cSRA,
# loaded from BAM by bam-load) it re-pairs reads by building a spot-id lookup
# over the PRIMARY_ALIGNMENT table, sorted externally in its temp dir, and that
# phase does not scale. Measured against Pseudophryne corroboree data:
#
#   accession     aligned reads   reference    fasterq-dump    fastq-dump
#   SRR000001     0 (unaligned)   -                     3 s           6 s
#   SRR28349113   54.2M           0.69 Gbp          47 min        88 min
#   SRR28349114   1,108.7M        8.88 Gbp    ~577 h (est.)  ~3.3 h (est.)
#
# It is not simply that the third run has 20x more alignments: the cost *per
# alignment* is ~90x worse, so total cost is closer to alignments x f(reference
# size) than to alignments alone. Alignment count and reference size are
# confounded in the only dataset where this has been measured, so the shortcut
# below deliberately requires BOTH to be large -- it fires only inside the
# region actually observed to fail. Every other run tries fasterq-dump first and
# relies on the budget to catch a bad case, which costs one budget's wall time
# but cannot mis-classify.
#
# The two extractors' output is byte-identical given --split-3 --skip-technical:
# verified on SRR000001, and on all 307 GB of SRR28349113 (400,318,389 reads per
# mate, cmp-clean). fastq-dump defaults differ from fasterq-dump's on both flags,
# so neither may be dropped.
_SRA_CONFIG = config["reads"]["sra"]
SRA_FASTERQ_BUDGET_SECONDS = int(_SRA_CONFIG["fasterq_budget_minutes"]) * 60
SRA_DIRECT_MIN_ALIGNMENTS = int(_SRA_CONFIG["direct_fastq_dump_min_alignments"])
SRA_DIRECT_MIN_REFERENCE_BP = int(_SRA_CONFIG["direct_fastq_dump_min_reference_bp"])

# The SRA REFERENCE table stores the reference in fixed-size rows, so the row
# count reported by `vdb-dump --info` estimates the reference length.
SRA_REFERENCE_CHUNK_BP = 5000


def fastp_input(wildcards):
    """Get input fastqs for filtering."""
    sample_rows = get_sample_rows(wildcards.sample)
    row = sample_rows[
        (sample_rows["library_id"] == wildcards.library)
        & (sample_rows["input_unit"] == wildcards.input_unit)
    ]
    if len(row) != 1:
        raise ValueError(
            f"Expected one row for sample={wildcards.sample}, "
            f"library={wildcards.library}, input_unit={wildcards.input_unit}; "
            f"found {len(row)}."
        )
    row = row.iloc[0]
    
    input_type = row["input_type"]
    
    if input_type == "srr":
        accession = row["input"]
        return {
            "r1": (
                f"results/fastqs/{wildcards.sample}/{wildcards.library}/"
                f"{wildcards.input_unit}/{accession}_1.fastq.gz"
            ),
            "r2": (
                f"results/fastqs/{wildcards.sample}/{wildcards.library}/"
                f"{wildcards.input_unit}/{accession}_2.fastq.gz"
            ),
        }
    elif input_type == "fastq":
        r1, r2 = row["input"].split(";")
        return {"r1": r1, "r2": r2}
    else:
        raise ValueError(f"Cannot get fastqs for input_type '{input_type}'")


rule download_sra:
    output:
        r1=temp("results/fastqs/{sample}/{library}/{input_unit}/{accession}_1.fastq.gz"),
        r2=temp("results/fastqs/{sample}/{library}/{input_unit}/{accession}_2.fastq.gz"),
    params:
        outdir="results/fastqs/{sample}/{library}/{input_unit}",
        fasterq_budget=SRA_FASTERQ_BUDGET_SECONDS,
        direct_min_alignments=SRA_DIRECT_MIN_ALIGNMENTS,
        direct_min_reference_bp=SRA_DIRECT_MIN_REFERENCE_BP,
        reference_chunk_bp=SRA_REFERENCE_CHUNK_BP,
    threads: 4
    conda:
        "../envs/sra.yaml"
    benchmark:
        "benchmarks/download_sra/{sample}/{library}/{input_unit}/{accession}.txt"
    log:
        "logs/download_sra/{sample}/{library}/{input_unit}/{accession}.txt"
    shell:
        """
        rm -rf {wildcards.accession}
        mkdir -p {params.outdir}
        
        if prefetch --max-size 1T {wildcards.accession} 2>> {log}; then
            # Table row counts describing the run's shape. For aligned (cSRA)
            # submissions: SEQ = spots, PRIM = aligned reads, REF = reference
            # chunks. Logged on every download so the thresholds below can be
            # recalibrated from real data rather than re-derived by hand.
            SRA_INFO=$(vdb-dump --info {wildcards.accession} 2>/dev/null)
            echo "$SRA_INFO" | grep -E '^(SEQ|REF|PRIM|SEC) ' >> {log} || true

            PRIM=$(echo "$SRA_INFO" \
                | awk -F: '$1 ~ /^PRIM/ {{gsub(/[^0-9]/, "", $2); print $2; exit}}')
            REF_ROWS=$(echo "$SRA_INFO" \
                | awk -F: '$1 ~ /^REF/ {{gsub(/[^0-9]/, "", $2); print $2; exit}}')
            PRIM=${{PRIM:-0}}
            REF_BP=$(( ${{REF_ROWS:-0}} * {params.reference_chunk_bp} ))

            # Skip fasterq-dump only where its lookup is known to be hopeless:
            # many alignments AND a large reference. Either alone is not enough
            # evidence to route around the faster default.
            USE_FASTQ_DUMP=0
            if [ "{params.direct_min_alignments}" -gt 0 ] \
               && [ "$PRIM" -ge "{params.direct_min_alignments}" ] \
               && [ "$REF_BP" -ge "{params.direct_min_reference_bp}" ]; then
                USE_FASTQ_DUMP=1
                echo "$PRIM primary alignments against ~${{REF_BP}} bp of reference:" \
                     "skipping fasterq-dump, extracting with fastq-dump." >> {log}
            fi

            # Private temp dir so a timed-out attempt can be cleaned up whole.
            FASTERQ_TMP=$(mktemp -d "{resources.tmpdir}/fasterq_{wildcards.accession}_XXXXXX")

            if [ "$USE_FASTQ_DUMP" -eq 0 ]; then
                if timeout {params.fasterq_budget} fasterq-dump {wildcards.accession} \
                        -O {params.outdir} \
                        -e {threads} \
                        -t "$FASTERQ_TMP" \
                        >> {log} 2>&1; then
                    :
                else
                    # A killed fasterq-dump leaves a partial FASTQ behind. It
                    # must be deleted, not compressed: pigz would happily gzip a
                    # truncated file and hand silently corrupt reads to fastp.
                    echo "fasterq-dump did not finish within {params.fasterq_budget}s;" \
                         "discarding partial output and falling back to fastq-dump." >> {log}
                    rm -f {params.outdir}/{wildcards.accession}*.fastq
                    USE_FASTQ_DUMP=1
                fi
            fi

            if [ "$USE_FASTQ_DUMP" -eq 1 ]; then
                # --split-3 --skip-technical are required for parity with
                # fasterq-dump's defaults; see the note at the top of this file.
                fastq-dump {wildcards.accession} \
                    --split-3 \
                    --skip-technical \
                    -O {params.outdir} \
                    >> {log} 2>&1
            fi

            rm -rf "$FASTERQ_TMP"
            pigz -p {threads} {params.outdir}/{wildcards.accession}*.fastq
            rm -rf {wildcards.accession}
        else
            echo "Prefetch failed, trying ENA via ffq..." >> {log}
            ffq --ftp {wildcards.accession} 2>> {log} \
                | jq -r '.[].url' \
                | while IFS= read -r url; do
                    [[ -n "$url" ]] || continue
                    echo "Downloading: $url" >> {log}
                    curl -fSL "$url" \
                        -o {params.outdir}/$(basename "$url") \
                        >> {log} 2>&1
                done
        fi
        
        [[ -f {output.r1} ]] && [[ -f {output.r2} ]]
        """


rule fastp:
    input:
        unpack(fastp_input),
    output:
        r1="results/filtered_fastqs/{sample}/{library}/{input_unit}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{library}/{input_unit}_2.fastq.gz",
        json="results/fastp/{sample}/{library}/{input_unit}.json",
    threads: 4
    conda:
        "../envs/fastp.yaml"
    benchmark:
        "benchmarks/fastp/{sample}/{library}/{input_unit}.txt"
    log:
        "logs/fastp/{sample}/{library}/{input_unit}.txt"
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            -j {output.json} \
            -h /dev/null \
            &> {log}
        """

def collect_fastp_input(wildcards):
    """Get all fastp JSON files for a sample (one per row/input unit)."""
    records = get_sample_inputs(wildcards.sample)
    return {
        "jsons": [
            f"results/fastp/{wildcards.sample}/{record['library_id']}/{record['input_unit']}.json"
            for record in records
        ],
    }


rule collect_fastp_stats:
    input:
        unpack(collect_fastp_input),
    output:
        json="results/qc_metrics/fastp/{sample}.json",
    log:
        "logs/collect_fastp_stats/{sample}.txt"
    run:
        import json
        
        unfiltered = 0
        pass_filter = 0
        
        for fn in input.jsons:
            with open(fn) as f:
                data = json.load(f)
            unfiltered += data["summary"]["before_filtering"]["total_reads"]
            pass_filter += data["summary"]["after_filtering"]["total_reads"]
        
        out = {
            "summary": {
                "before_filtering": {"total_reads": unfiltered},
                "after_filtering": {"total_reads": pass_filter},
            }
        }
        
        with open(output.json, "w") as f:
            json.dump(out, f, indent=2)

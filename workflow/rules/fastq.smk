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
    threads: 4
    conda:
        "../envs/sra.yaml"
    benchmark:
        "benchmarks/download_sra/{sample}/{library}/{input_unit}/{accession}.txt"
    log:
        "logs/download_sra/{sample}/{library}/{input_unit}/{accession}.txt"
    shell:
        """
        # This rule pipes, and Snakemake does not enable pipefail by default.
        # Without it a failing sam-dump mid-stream still leaves samtools fastq
        # exiting 0, and the rule reports success on a truncated FASTQ.
        set -o pipefail

        mkdir -p {params.outdir}

        # Everything this rule downloads or spills lives in one job-private
        # directory, removed on exit however the job ends. prefetch previously
        # wrote into the working directory, so every concurrent download_sra
        # job created a sibling accession directory in the workflow root, and
        # the rule's own `rm -rf <accession>` could delete another job's
        # in-flight download whenever two sample rows shared an accession.
        WORK=$(mktemp -d "{resources.tmpdir}/sra_{wildcards.accession}_XXXXXX")
        trap 'rm -rf "$WORK"' EXIT

        SEQ=0
        ALIGNED=0

        if prefetch --max-size 1T -O "$WORK" {wildcards.accession} 2>> {log}; then
            # prefetch lays this out as <outdir>/<accession>/<accession>.sra,
            # with any reference object the run needs alongside it. fasterq-dump
            # accepts the file path, which keeps it independent of the working
            # directory, and names its output from the basename.
            SRA="$WORK/{wildcards.accession}/{wildcards.accession}.sra"
            [ -f "$SRA" ] || SRA="$WORK/{wildcards.accession}/{wildcards.accession}.sralite"

            # Belt and braces: every extractor below is given this path, so a
            # missing archive already fails rather than resolving the accession
            # over the network. Checking here names the cause instead of leaving
            # a tool-specific "item not found" in the log.
            if [ ! -f "$SRA" ]; then
                echo "ERROR: {wildcards.accession}: prefetch exited 0 but no .sra/.sralite" \
                     "is present under $WORK. Refusing to fall back to network streaming." >> {log}
                exit 1
            fi

            # Table row counts describing the run's shape: SEQ = spots,
            # PRIM = aligned reads, REF = reference in 5 kb chunks. Logged on
            # every download; PRIM selects the extractor below.
            SRA_INFO=$(vdb-dump --info "$SRA" 2>/dev/null)
            echo "$SRA_INFO" | grep -E '^(SEQ|REF|PRIM|SEC) ' >> {log} || true
            SEQ=$(echo "$SRA_INFO" \
                | awk -F: '$1 ~ /^SEQ/ {{gsub(/[^0-9]/, "", $2); print $2; exit}}')
            PRIM=$(echo "$SRA_INFO" \
                | awk -F: '$1 ~ /^PRIM/ {{gsub(/[^0-9]/, "", $2); print $2; exit}}')
            SEQ=${{SEQ:-0}}
            PRIM=${{PRIM:-0}}

            if [ "$PRIM" -gt 0 ]; then
                ALIGNED=1
                # Reference-aligned (cSRA) submission. Reads are stored in
                # reference order and reference-compressed, so a spot-ordered
                # extractor (fasterq-dump, fastq-dump) random-accesses the whole
                # reference once per read and degrades without bound as the
                # reference outgrows cache -- measured decaying from 5.4 GB/h to
                # 0.08 GiB/h on an 8.87 Gb reference, projecting ~52 days for a
                # single accession. sam-dump reads in reference order instead:
                # sequential, with no locality to lose. samtools collate then
                # restores name grouping and samtools fastq writes the pairs.
                #
                # One knob: sam-dump is single-threaded, and the remaining
                # threads are split between the two samtools stages.
                SAMTOOLS_THREADS=$(( ({threads} - 1) / 2 ))
                [ "$SAMTOOLS_THREADS" -lt 0 ] && SAMTOOLS_THREADS=0

                # collate gets -u. Without it, collate compresses its stdout and
                # samtools fastq immediately decompresses it again, burning CPU
                # on both sides of the pipe for nothing.
                sam-dump -u "$SRA" 2>> {log} \
                    | samtools collate -u -O -@ "$SAMTOOLS_THREADS" - "$WORK/collate" 2>> {log} \
                    | samtools fastq -F 0x900 -@ "$SAMTOOLS_THREADS" \
                        -1 {output.r1} \
                        -2 {output.r2} \
                        -0 "$WORK/unpaired.fastq.gz" \
                        -s "$WORK/singleton.fastq.gz" \
                        - 2>> {log}
            else
                # No alignment table: there is no reference to walk, so
                # fasterq-dump's multi-threaded read of the SEQUENCE table is
                # the right access pattern and sam-dump has no advantage.
                mkdir -p "$WORK/fasterq"
                fasterq-dump "$SRA" \
                    -O {params.outdir} \
                    -e {threads} \
                    -t "$WORK/fasterq" \
                    >> {log} 2>&1
                pigz -p {threads} {params.outdir}/{wildcards.accession}*.fastq
            fi
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
        
        # The previous test only asked whether the output files exist. That
        # passes for a pigz that compressed nothing and exited 0, for a
        # wall-killed extractor that left a partial FASTQ, and for silent read
        # loss from a mis-set flag -- all of which have happened on this rule.
        count_reads() {{
            [ -s "$1" ] || {{ echo 0; return; }}
            pigz -dc "$1" | wc -l | awk '{{print int($1 / 4)}}'
        }}

        N1=$(count_reads {output.r1})
        N2=$(count_reads {output.r2})
        echo "emitted: r1=$N1 r2=$N2" >> {log}

        if [ "$N1" -eq 0 ] || [ "$N2" -eq 0 ]; then
            echo "ERROR: {wildcards.accession}: empty FASTQ output (r1=$N1 r2=$N2)" >> {log}
            exit 1
        fi

        if [ "$N1" -ne "$N2" ]; then
            echo "ERROR: {wildcards.accession}: mate files disagree (r1=$N1 r2=$N2)" >> {log}
            exit 1
        fi

        if [ "$ALIGNED" -eq 1 ]; then
            # Every read in the archive must come out exactly once. The
            # unpaired and singleton files are counted here so that reads
            # landing there are not mistaken for reads lost; whether any read
            # landed there at all is asserted separately, because nothing
            # downstream consumes those files.
            NU=$(count_reads "$WORK/unpaired.fastq.gz")
            NS=$(count_reads "$WORK/singleton.fastq.gz")
            TOTAL=$(( N1 + N2 + NU + NS ))
            EXPECTED=$(( SEQ * 2 ))
            echo "completeness: emitted=$TOTAL expected=$EXPECTED unpaired=$NU singleton=$NS" >> {log}

            if [ "$TOTAL" -ne "$EXPECTED" ]; then
                echo "ERROR: {wildcards.accession}: emitted $TOTAL reads, archive holds" \
                     "$SEQ spots (expected $EXPECTED). Reads were lost or duplicated." >> {log}
                exit 1
            fi

            if [ "$NU" -ne 0 ] || [ "$NS" -ne 0 ]; then
                echo "ERROR: {wildcards.accession}: $NU unpaired and $NS singleton reads" \
                     "were extracted but nothing downstream consumes them." >> {log}
                exit 1
            fi
        fi
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

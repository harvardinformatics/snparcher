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
        EXPECTED=""
        EXPECTED_SOURCE=""
        ENA_PATH=0

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

                # sam-dump emits every stored read and -F 0x900 removes only
                # the secondary/supplementary copies, so the archive's own spot
                # count is the expectation. cSRA submissions are paired
                # alignments, two reads per spot.
                EXPECTED=$(( SEQ * 2 ))
                EXPECTED_SOURCE="2 x SEQ"
            else
                # No alignment table: there is no reference to walk, so
                # fasterq-dump's multi-threaded read of the SEQUENCE table is
                # the right access pattern and sam-dump has no advantage.
                mkdir -p "$WORK/fasterq"
                fasterq-dump "$SRA" \
                    -O {params.outdir} \
                    -e {threads} \
                    -t "$WORK/fasterq" \
                    > "$WORK/fasterq-dump.log" 2>&1
                cat "$WORK/fasterq-dump.log" >> {log}
                pigz -p {threads} {params.outdir}/{wildcards.accession}*.fastq

                # 2 x SEQ is not the right expectation here: fasterq-dump drops
                # technical and zero-length reads, and a spot need not hold two
                # biological reads. SRR000001 stores four reads per spot and
                # emits 707,026 of 1,883,940. Use the count the tool reports.
                EXPECTED=$(awk -F: '/reads written/ {{gsub(/[^0-9]/, "", $2); print $2; exit}}' \
                    "$WORK/fasterq-dump.log")
                EXPECTED_SOURCE="fasterq-dump reads-written"
            fi
        else
            # prefetch can fail for a single run while every other accession in
            # the same project succeeds -- NCBI returns "file unauthorized" for
            # some runs -- so this path is the only way to obtain otherwise
            # good data, and it must not be gated on something it cannot
            # produce. There is no extractor here to report a read count;
            # instead ENA publishes an md5 per file, which is a stronger
            # completeness check than a count because it catches truncation and
            # corruption exactly. That is what this path asserts.
            ENA_PATH=1
            echo "Prefetch failed, trying ENA via ffq..." >> {log}

            # Land the file list before looping. Feeding the loop by pipeline
            # runs it in a subshell, where a failure cannot stop the rule.
            ffq --ftp {wildcards.accession} 2>> {log} \
                | jq -r '.[] | .url + " " + (.md5 // "")' \
                > "$WORK/ena_files.txt"

            ENA_FILES=0
            ENA_NO_MD5=0
            while read -r url md5; do
                [[ -n "$url" ]] || continue
                dest={params.outdir}/$(basename "$url")
                echo "Downloading: $url" >> {log}
                if ! curl -fSL "$url" -o "$dest" >> {log} 2>&1; then
                    echo "ERROR: {wildcards.accession}: download failed for $url" >> {log}
                    exit 1
                fi
                ENA_FILES=$(( ENA_FILES + 1 ))

                if [ -n "$md5" ]; then
                    if echo "$md5  $dest" | md5sum -c - >> {log} 2>&1; then
                        echo "md5 verified: $dest" >> {log}
                    else
                        echo "ERROR: {wildcards.accession}: md5 mismatch for $dest;" \
                             "the download is corrupt or truncated." >> {log}
                        exit 1
                    fi
                else
                    ENA_NO_MD5=$(( ENA_NO_MD5 + 1 ))
                fi
            done < "$WORK/ena_files.txt"

            if [ "$ENA_FILES" -eq 0 ]; then
                echo "ERROR: {wildcards.accession}: prefetch failed and ffq returned" \
                     "no ENA files, so there is nothing to download." >> {log}
                exit 1
            fi

            if [ "$ENA_NO_MD5" -ne 0 ]; then
                echo "WARNING: {wildcards.accession}: $ENA_NO_MD5 of $ENA_FILES ENA files" \
                     "carried no md5, so their integrity is unverified." >> {log}
            fi
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

        # Reads that came out of the archive but are not part of the declared
        # pair. fasterq-dump's split-3 behaviour emits an orphan file for spots
        # carrying a single read, which pigz compresses alongside the pair;
        # samtools fastq writes the equivalent to -0/-s. They are counted so
        # that reads landing there are not mistaken for reads lost.
        NORPHAN=$(count_reads {params.outdir}/{wildcards.accession}.fastq.gz)
        NU=$(count_reads "$WORK/unpaired.fastq.gz")
        NS=$(count_reads "$WORK/singleton.fastq.gz")

        echo "emitted: r1=$N1 r2=$N2 orphan=$NORPHAN unpaired=$NU singleton=$NS" >> {log}

        if [ "$N1" -eq 0 ] || [ "$N2" -eq 0 ]; then
            echo "ERROR: {wildcards.accession}: empty FASTQ output (r1=$N1 r2=$N2)" >> {log}
            exit 1
        fi

        if [ "$N1" -ne "$N2" ]; then
            echo "ERROR: {wildcards.accession}: mate files disagree (r1=$N1 r2=$N2)" >> {log}
            exit 1
        fi

        # Every read the extractor produced must reach the output files. Both
        # paths check this; they differ only in what the expected count is,
        # because each extractor decides for itself which stored reads to emit.
        TOTAL=$(( N1 + N2 + NORPHAN + NU + NS ))
        echo "completeness: emitted=$TOTAL expected=${{EXPECTED:-unknown}} (${{EXPECTED_SOURCE:-none}})" >> {log}

        if [ -n "$EXPECTED" ]; then
            if [ "$TOTAL" -ne "$EXPECTED" ]; then
                echo "ERROR: {wildcards.accession}: emitted $TOTAL reads, expected" \
                     "$EXPECTED ($EXPECTED_SOURCE). Reads were lost or duplicated." >> {log}
                exit 1
            fi
        elif [ "$ENA_PATH" -eq 1 ]; then
            # Files came from ENA rather than an extractor, so there is no
            # read count to compare against. Completeness for this path was
            # established by the md5 check above, which already exited on
            # mismatch.
            echo "completeness: established by ENA md5, not by read count" >> {log}
        else
            echo "ERROR: {wildcards.accession}: could not determine the expected read" \
                 "count, so completeness cannot be checked. The extractor probably" \
                 "did not run to completion; see the log above." >> {log}
            exit 1
        fi

        # Reads the pipeline extracted but no downstream rule consumes. On the
        # aligned path these should never appear: a cSRA run's mates are both
        # present, so anything here means collate failed to pair them.
        if [ "$NU" -ne 0 ] || [ "$NS" -ne 0 ]; then
            echo "ERROR: {wildcards.accession}: $NU unpaired and $NS singleton reads" \
                 "were extracted but nothing downstream consumes them." >> {log}
            exit 1
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

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
        rm -rf {wildcards.accession}
        mkdir -p {params.outdir}
        
        if prefetch --max-size 1T {wildcards.accession} 2>> {log}; then
            fasterq-dump {wildcards.accession} \
                -O {params.outdir} \
                -e {threads} \
                -t {resources.tmpdir} \
                >> {log} 2>&1
            pigz -p {threads} {params.outdir}/{wildcards.accession}*.fastq
            rm -rf {wildcards.accession}
        else
            echo "Prefetch failed, trying ENA via ffq..." >> {log}
            ffq --ftp {wildcards.accession} 2>> {log} \
                | jq -r '.[].url' \
                | xargs -I {{}} curl -fSL -o {params.outdir}/$(basename {{}}) {{}} 2>> {log}
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

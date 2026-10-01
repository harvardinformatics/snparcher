# The repadapt mapping pipeline: RepAdapt's fastq -> BAM processing
# (https://github.com/RepAdapt/nextflow_snp_calling_linux, commit 2077f6f),
# with its pinned tools and result-affecting arguments, inside snpArcher's
# structure (rows mapped separately, read groups set at mapping, snpArcher's
# thread counts):
#
#   fastp 0.20.1, RepAdapt's defaults (no --detect_adapter_for_pe)
#   bwa 0.7.17 mem, no -M; -K 40000000 gives RepAdapt's -t 4 batching at any
#     thread count (bwa estimates insert sizes per batch)
#   samtools 1.16.1: view -q 10 | sort -n | fixmate -m | sort, RepAdapt's order
#   Picard 2.26.3 MarkDuplicates -REMOVE_DUPLICATES true, once per sample: one
#     library tag per sample, as RepAdapt's AddOrReplaceReadGroups sets
#   GATK 3.8-1 RealignerTargetCreator + IndelRealigner, unless
#     REPADAPT_INDEL_REALIGNMENT is off (long-contig mode or config)
#
# The rules named fastp and bwa_mem keep the default pipeline's names and
# output paths, so the profiles and collect_fastp_stats apply unchanged. The
# final BAM (see _repadapt_final_bam in common.smk) then goes through the
# shared CSI indexing, bam_stats, callable sites and calling.

REPADAPT_BAM_DIR = "results/bams/repadapt"
REPADAPT_QC_DIR = "results/qc_metrics/repadapt"
# GATK 3.8 can't read a bgzipped FASTA; it finds {REF_NAME}.dict (samtools
# dict, from index_reference) next to this uncompressed copy.
REPADAPT_REF_FASTA = f"results/reference/{REF_NAME}.fa"


def get_repadapt_read_group(wildcards):
    """Read group for one row: snpArcher's row-level ID, RepAdapt's single
    library tag per sample (Picard removes duplicates per library)."""
    return (
        f"'@RG\\tID:{wildcards.library}.{wildcards.input_unit}"
        f"\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}_LB\\tPL:ILLUMINA'"
    )


def _repadapt_rows(sample):
    return [(record["library_id"], record["input_unit"]) for record in get_sample_inputs(sample)]


def get_repadapt_row_bams(wildcards):
    """Every row's filtered, fixmate'd, coordinate-sorted BAM for a sample."""
    return [
        f"results/bams/raw/{wildcards.sample}/{library}/{unit}.bam"
        for library, unit in _repadapt_rows(wildcards.sample)
    ]


def get_repadapt_unrealigned_bam(sample):
    """The BAM that realignment reads: deduplicated, or merged when the sample
    has mark_duplicates: false."""
    if get_sample_mark_duplicates(sample):
        return f"{REPADAPT_BAM_DIR}/dedup/{sample}.bam"
    return f"{REPADAPT_BAM_DIR}/merged/{sample}.bam"


def _before_realignment(path):
    """Intermediate when realignment follows; otherwise the final BAM."""
    return temp(path) if REPADAPT_INDEL_REALIGNMENT else path


rule fastp:
    input:
        unpack(fastp_input),
    output:
        r1="results/filtered_fastqs/{sample}/{library}/{input_unit}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{library}/{input_unit}_2.fastq.gz",
        json="results/fastp/{sample}/{library}/{input_unit}.json",
    threads: 4
    conda:
        "../../envs/repadapt/fastp.yaml"
    benchmark:
        "benchmarks/fastp/{sample}/{library}/{input_unit}.txt"
    log:
        "logs/fastp/{sample}/{library}/{input_unit}.txt"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -w {threads} -j {output.json} -h /dev/null &> {log}
        """


rule bwa_mem:
    input:
        unpack(bwa_mem_input),
    output:
        bam=temp("results/bams/raw/{sample}/{library}/{input_unit}.bam"),
        # Flagstat of bwa's output, before the MAPQ filter, for mapping QC.
        flagstat=temp(REPADAPT_QC_DIR + "/flagstat/{sample}/{library}/{input_unit}.tsv"),
    params:
        rg=get_repadapt_read_group,
    threads: 8
    conda:
        "../../envs/repadapt/mapping.yaml"
    benchmark:
        "benchmarks/bwa_mem/{sample}/{library}/{input_unit}.txt"
    log:
        "logs/bwa_mem/{sample}/{library}/{input_unit}.txt"
    shell:
        # flagstat reads a copy of the stream through a FIFO and is waited on,
        # so its failure fails the job and its output is complete at the end.
        """
        tmp=$(mktemp -d {resources.tmpdir}/repadapt_bwa_mem.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        mkfifo "$tmp/sam"
        : > {log}
        samtools flagstat -O tsv "$tmp/sam" > {output.flagstat} 2>> {log} &
        flagstat_pid=$!
        bwa mem -K 40000000 -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2>> {log} \
            | tee "$tmp/sam" \
            | samtools view -u -q 10 - 2>> {log} \
            | samtools sort -n -u -T "$tmp/name" - 2>> {log} \
            | samtools fixmate -m - - 2>> {log} \
            | samtools sort -T "$tmp/coord" -o {output.bam} - 2>> {log}
        wait "$flagstat_pid"
        """


rule repadapt_remove_duplicates:
    """Picard reads all of a sample's rows (one -INPUT each) and removes
    duplicates per library, i.e. per sample."""
    input:
        bams=get_repadapt_row_bams,
    output:
        bam=_before_realignment(REPADAPT_BAM_DIR + "/dedup/{sample}.bam"),
        metrics=REPADAPT_QC_DIR + "/{sample}_duplicates.txt",
    params:
        inputs=lambda wc, input: " ".join(f"-INPUT {bam}" for bam in input.bams),
    threads: 1
    conda:
        "../../envs/repadapt/picard.yaml"
    benchmark:
        "benchmarks/repadapt_remove_duplicates/{sample}.txt"
    log:
        "logs/repadapt_remove_duplicates/{sample}.txt"
    shell:
        """
        picard -Xmx{resources.mem_mb_reduced}m MarkDuplicates {params.inputs} \
            -OUTPUT {output.bam} -METRICS_FILE {output.metrics} \
            -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT \
            --TMP_DIR {resources.tmpdir} &> {log}
        """


rule repadapt_merge_sample:
    """Rows of a sample with mark_duplicates: false, merged without duplicate
    removal."""
    input:
        bams=get_repadapt_row_bams,
    output:
        bam=_before_realignment(REPADAPT_BAM_DIR + "/merged/{sample}.bam"),
    conda:
        "../../envs/repadapt/mapping.yaml"
    benchmark:
        "benchmarks/repadapt_merge_sample/{sample}.txt"
    log:
        "logs/repadapt_merge_sample/{sample}.txt"
    shell:
        """
        samtools merge {output.bam} {input.bams} 2> {log}
        """


if REPADAPT_INDEL_REALIGNMENT:

    def repadapt_realignment_input(wildcards):
        bam = get_repadapt_unrealigned_bam(wildcards.sample)
        return {
            "bam": bam,
            "bai": f"{bam}.bai",
            "ref": REPADAPT_REF_FASTA,
            "ref_fai": f"{REPADAPT_REF_FASTA}.fai",
            "ref_dict": REF_FILES["ref_dict"],
        }

    rule repadapt_reference_fasta:
        input:
            ref=REF_FILES["ref"],
        output:
            fasta=temp(REPADAPT_REF_FASTA),
            fai=temp(f"{REPADAPT_REF_FASTA}.fai"),
        conda:
            "../../envs/repadapt/mapping.yaml"
        log:
            "logs/repadapt_reference_fasta.txt"
        shell:
            """
            bgzip -dc {input.ref} > {output.fasta} 2> {log}
            samtools faidx {output.fasta} 2>> {log}
            """

    rule repadapt_index_bai:
        """GATK3 needs a BAI; snpArcher's indexes are otherwise CSI."""
        input:
            bam=REPADAPT_BAM_DIR + "/{stage}/{sample}.bam",
        output:
            bai=temp(REPADAPT_BAM_DIR + "/{stage}/{sample}.bam.bai"),
        wildcard_constraints:
            stage="dedup|merged",
        conda:
            "../../envs/repadapt/mapping.yaml"
        log:
            "logs/repadapt_index_bai/{stage}/{sample}.txt"
        shell:
            """
            samtools index -b {input.bam} {output.bai} 2> {log}
            """

    rule repadapt_realigner_targets:
        input:
            unpack(repadapt_realignment_input),
        output:
            intervals=temp(REPADAPT_BAM_DIR + "/realign/{sample}.intervals"),
        threads: 1
        conda:
            "../../envs/repadapt/gatk3.yaml"
        benchmark:
            "benchmarks/repadapt_realigner_targets/{sample}.txt"
        log:
            "logs/repadapt_realigner_targets/{sample}.txt"
        shell:
            """
            gatk3 -Xmx{resources.mem_mb_reduced}m -Djava.io.tmpdir={resources.tmpdir} \
                -T RealignerTargetCreator -R {input.ref} -I {input.bam} \
                -o {output.intervals} &> {log}
            """

    rule repadapt_indel_realigner:
        """Single-threaded, as RepAdapt runs it; with its fixed random seed the
        output is deterministic."""
        input:
            unpack(repadapt_realignment_input),
            intervals=REPADAPT_BAM_DIR + "/realign/{sample}.intervals",
        output:
            bam=REPADAPT_BAM_DIR + "/realigned/{sample}.bam",
            # GATK writes this index itself; the final index is CSI.
            bai=temp(REPADAPT_BAM_DIR + "/realigned/{sample}.bai"),
        threads: 1
        conda:
            "../../envs/repadapt/gatk3.yaml"
        benchmark:
            "benchmarks/repadapt_indel_realigner/{sample}.txt"
        log:
            "logs/repadapt_indel_realigner/{sample}.txt"
        shell:
            """
            gatk3 -Xmx{resources.mem_mb_reduced}m -Djava.io.tmpdir={resources.tmpdir} \
                -T IndelRealigner -R {input.ref} -I {input.bam} \
                -targetIntervals {input.intervals} --consensusDeterminationModel USE_READS \
                -o {output.bam} &> {log}
            """


def repadapt_mapping_qc_input(wildcards):
    sample = wildcards.sample
    inputs = {
        "flagstats": [
            f"{REPADAPT_QC_DIR}/flagstat/{sample}/{library}/{unit}.tsv"
            for library, unit in _repadapt_rows(sample)
        ],
        "coverage": f"results/qc_metrics/bam/{sample}_coverage.txt",
    }
    if get_sample_mark_duplicates(sample):
        inputs["duplicates"] = f"{REPADAPT_QC_DIR}/{sample}_duplicates.txt"
    return inputs


rule repadapt_mapping_qc:
    """Per-sample mapping QC with parse_bam_stats' keys. Mapping counts come
    from bwa's output before the MAPQ filter (summed over rows), duplicates
    from Picard's metrics (Picard removes them, so the final BAM has none), and
    depth from the final BAM."""
    input:
        unpack(repadapt_mapping_qc_input),
    output:
        json=REPADAPT_QC_DIR + "/{sample}.json",
    log:
        "logs/repadapt_mapping_qc/{sample}.txt"
    run:
        import json

        # samtools flagstat -O tsv lines: 0 total, 6 mapped,
        # 10 paired in sequencing, 13 properly paired.
        counts = {"total": 0, "mapped": 0, "paired": 0, "proper": 0}
        for path in input.flagstats:
            with open(path) as handle:
                lines = handle.read().splitlines()
            for key, index in (("total", 0), ("mapped", 6), ("paired", 10), ("proper", 13)):
                counts[key] += int(lines[index].split("\t")[0])

        num_dups = 0
        pct_dups = 0.0
        if "duplicates" in input.keys():
            rows = read_picard_metrics(input.duplicates)
            if len(rows) != 1:
                raise ValueError(f"Expected one library in {input.duplicates}, found {len(rows)}")
            row = rows[0]
            num_dups = int(row["UNPAIRED_READ_DUPLICATES"]) + 2 * int(row["READ_PAIR_DUPLICATES"])
            pct_dups = float(row["PERCENT_DUPLICATION"]) * 100 if row["PERCENT_DUPLICATION"] != "?" else 0.0

        mean_depth, covered_bases = parse_samtools_coverage(input.coverage)
        total = counts["total"]
        out = {
            "sample": wildcards.sample,
            "total_reads": total,
            "num_mapped": counts["mapped"],
            "percent_mapped": counts["mapped"] / total * 100 if total else 0,
            "num_duplicates": num_dups,
            "percent_duplicates": pct_dups,
            "percent_properly_paired": counts["proper"] / counts["paired"] * 100 if counts["paired"] else 0,
            "mean_depth": mean_depth,
            "covered_bases": covered_bases,
        }
        with open(output.json, "w") as handle:
            json.dump(out, handle, indent=2)

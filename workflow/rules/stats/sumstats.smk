"""
Summary statistics rules for snpArcher v2.

Collects BAM statistics, fastp output, and generates summary reports.

Rules:
- stats_sentieon_bam: Compute additional Sentieon-specific BAM stats
- stats_parse_fastp: Parse fastp JSON for a single sample
- stats_parse_flagstat: Parse samtools flagstat TSV for a single sample
- stats_parse_coverage: Parse samtools coverage for a single sample
- stats_parse_inserts: Parse Sentieon insert size metrics for a single sample
- stats_collect_sumstats: Aggregate all parsed stats into final summary file
"""


rule stats_sentieon_bam:
    """Compute Sentieon-specific BAM statistics."""
    input:
        unpack(get_final_bam),
        unpack(get_ref_bundle),
    output:
        insert_file="results/summary_stats/{sample}_insert_metrics.txt",
        qd="results/summary_stats/{sample}_qd_metrics.txt",
        gc="results/summary_stats/{sample}_gc_metrics.txt",
        gc_summary="results/summary_stats/{sample}_gc_summary.txt",
        mq="results/summary_stats/{sample}_mq_metrics.txt",
    params:
        lic=config["variant_calling"]["sentieon"]["license"]
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/stats_sentieon_bam/{sample}.txt"
    benchmark:
        "benchmarks/stats_sentieon_bam/{sample}.txt"
    threads: 4
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} \
            -t {threads} -i {input.bam} \
            --algo MeanQualityByCycle {output.mq} \
            --algo QualDistribution {output.qd} \
            --algo GCBias --summary {output.gc_summary} {output.gc} \
            --algo InsertSizeMetricAlgo {output.insert_file} \
            2> {log}
        """


rule stats_parse_fastp:
    """Parse fastp JSON for a single sample."""
    input:
        fastp="results/fastp/{sample}.fastp.json",
    output:
        json="results/summary_stats/{sample}.fastp_parsed.json",
    run:
        import json

        with open(input.fastp) as f:
            data = json.load(f)

        unfiltered = data["summary"]["before_filtering"]["total_reads"]
        pass_filter = data["summary"]["after_filtering"]["total_reads"]
        result = {
            "fraction_pass": float(pass_filter / unfiltered) if unfiltered else 0,
            "num_pass": pass_filter,
        }

        with open(output.json, "w") as f:
            json.dump(result, f)


rule stats_parse_flagstat:
    """Parse samtools flagstat TSV for a single sample."""
    input:
        flagstat="results/bams/stats/{sample}.flagstat.tsv",
    output:
        json="results/summary_stats/{sample}.flagstat_parsed.json",
    run:
        import json

        result = {}
        with open(input.flagstat) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    count = int(fields[0])
                    metric = fields[2]
                    if metric.startswith("total"):
                        result["total"] = count
                    elif metric == "mapped":
                        result["mapped"] = count
                    elif metric == "duplicates":
                        result["duplicates"] = count
                    elif metric == "properly paired":
                        result["properly_paired"] = count

        with open(output.json, "w") as f:
            json.dump(result, f)


rule stats_parse_coverage:
    """Parse samtools coverage for a single sample."""
    input:
        coverage="results/bams/stats/{sample}.coverage.txt",
    output:
        json="results/summary_stats/{sample}.coverage_parsed.json",
    run:
        import json

        num_sites = []
        cov_bases = 0
        depths = []

        with open(input.coverage) as f:
            for line in f:
                if not line.startswith("#rname"):
                    fields = line.split()
                    num_sites.append(int(fields[2]) - int(fields[1]) + 1)
                    depths.append(float(fields[6]))
                    cov_bases += float(fields[4])

        total = sum(num_sites)
        depths_mean = sum(d * n / total for d, n in zip(depths, num_sites)) if total else 0
        result = {"depth_mean": depths_mean, "covered_bases": cov_bases}

        with open(output.json, "w") as f:
            json.dump(result, f)


rule stats_parse_inserts:
    """Parse Sentieon insert size metrics for a single sample."""
    input:
        inserts="results/summary_stats/{sample}_insert_metrics.txt",
    output:
        json="results/summary_stats/{sample}.inserts_parsed.json",
    run:
        import json

        result = {"median": None, "std": None}
        with open(input.inserts) as f:
            for i, line in enumerate(f):
                if i == 2:
                    fields = line.strip().split()
                    if len(fields) >= 2:
                        result["median"] = fields[0]
                        result["std"] = fields[1]

        with open(output.json, "w") as f:
            json.dump(result, f)


def _get_sumstats_jsons(wildcards):
    """Get all parsed JSON files for aggregation."""
    out = {}

    if BAM_REQUIRED_SAMPLES:
        out["flagstat"] = expand(
            "results/summary_stats/{sample}.flagstat_parsed.json",
            sample=BAM_REQUIRED_SAMPLES,
        )
        out["coverage"] = expand(
            "results/summary_stats/{sample}.coverage_parsed.json",
            sample=BAM_REQUIRED_SAMPLES,
        )

    if ALIGN_SAMPLES:
        out["fastp"] = expand(
            "results/summary_stats/{sample}.fastp_parsed.json",
            sample=ALIGN_SAMPLES,
        )

    if config["sentieon"] and BAM_REQUIRED_SAMPLES:
        out["inserts"] = expand(
            "results/summary_stats/{sample}.inserts_parsed.json",
            sample=BAM_REQUIRED_SAMPLES,
        )

    return out


rule stats_collect_sumstats:
    """Aggregate all parsed stats into final summary file."""
    input:
        unpack(_get_sumstats_jsons),
    output:
        f"results/summary_stats/{REF_NAME}_bam_sumstats.txt",
    run:
        import json
        import os

        # Load flagstat data
        aln_metrics = {}
        if hasattr(input, "flagstat"):
            for fn in input.flagstat:
                sample = os.path.basename(fn).replace(".flagstat_parsed.json", "")
                with open(fn) as f:
                    data = json.load(f)
                total = data.get("total", 1)
                mapped = data.get("mapped", 0)
                duplicates = data.get("duplicates", 0)
                properly_paired = data.get("properly_paired", 0)
                aln_metrics[sample] = {
                    "total": total,
                    "percent_mapped": f"{100 * mapped / total:.2f}" if total else "0.00",
                    "duplicates": duplicates,
                    "percent_properly_paired": f"{100 * properly_paired / total:.2f}" if total else "0.00",
                }

        # Load coverage data
        seq_depths = {}
        covered_bases = {}
        if hasattr(input, "coverage"):
            for fn in input.coverage:
                sample = os.path.basename(fn).replace(".coverage_parsed.json", "")
                with open(fn) as f:
                    data = json.load(f)
                seq_depths[sample] = data.get("depth_mean", 0)
                covered_bases[sample] = data.get("covered_bases", 0)

        # Load fastp data
        fraction_pass = {}
        num_pass = {}
        if hasattr(input, "fastp"):
            for fn in input.fastp:
                sample = os.path.basename(fn).replace(".fastp_parsed.json", "")
                with open(fn) as f:
                    data = json.load(f)
                fraction_pass[sample] = data.get("fraction_pass", 0)
                num_pass[sample] = data.get("num_pass", 0)

        # Load insert size data (Sentieon only)
        med_inserts = {}
        med_insert_std = {}
        if hasattr(input, "inserts"):
            for fn in input.inserts:
                sample = os.path.basename(fn).replace(".inserts_parsed.json", "")
                with open(fn) as f:
                    data = json.load(f)
                med_inserts[sample] = data.get("median")
                med_insert_std[sample] = data.get("std")

        # Determine sample list from available data
        sample_list = list(aln_metrics.keys()) if aln_metrics else []

        # Build header
        header = [
            "Sample",
            "Total_Reads",
            "Percent_mapped",
            "Num_duplicates",
            "Percent_properly_paired",
            "Fraction_reads_pass_filter",
            "NumReadsPassingFilters",
        ]

        if med_inserts:
            header.extend(["MedianInsertSize", "MedianAbsDev_InsertSize"])

        # Write output
        with open(output[0], "w") as f:
            print("\t".join(header), file=f)

            for samp in sample_list:
                row = [
                    samp,
                    str(aln_metrics.get(samp, {}).get("total", "NA")),
                    str(aln_metrics.get(samp, {}).get("percent_mapped", "NA")),
                    str(aln_metrics.get(samp, {}).get("duplicates", "NA")),
                    str(aln_metrics.get(samp, {}).get("percent_properly_paired", "NA")),
                    str(fraction_pass.get(samp, "NA")),
                    str(num_pass.get(samp, "NA")),
                ]

                if med_inserts:
                    row.extend([
                        str(med_inserts.get(samp, "NA")),
                        str(med_insert_std.get(samp, "NA")),
                    ])

                print("\t".join(row), file=f)

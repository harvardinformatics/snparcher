"""
snpArcher v2 - Common configuration, sample sheet parsing, and helper functions.
"""

import json
import os
import random
import statistics
import string
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from snakemake.utils import min_version, validate

min_version("9.0")

validate(config, schema="../schemas/config.schema.json")

REF_NAME = config["reference"]["name"]
REF_PATH = config["reference"].get("path")
REF_ACCESSION = config["reference"].get("accession")
REF_FILE = f"{REF_NAME}.fna.gz"  # Reference filename (bgzip-compressed)

samples = pd.read_csv(config["samples"])
validate(samples, schema="../schemas/samples.schema.json")


def _validate_library_ids_for_multi_row_samples(df: pd.DataFrame) -> None:
    """Error if a sample_id appears multiple times without explicit library_id."""
    sample_counts = df.groupby("sample_id").size()
    multi_row_samples = sample_counts[sample_counts > 1].index.tolist()

    if not multi_row_samples:
        return

    has_library_col = "library_id" in df.columns

    for sample_id in multi_row_samples:
        sample_rows = df[df["sample_id"] == sample_id]
        if has_library_col:
            missing_mask = sample_rows["library_id"].isna() | (sample_rows["library_id"] == "")
        else:
            missing_mask = pd.Series([True] * len(sample_rows), index=sample_rows.index)

        if missing_mask.any():
            missing_indices = missing_mask[missing_mask].index.tolist()
            raise ValueError(
                f"Sample '{sample_id}' has {len(sample_rows)} rows but is missing "
                f"library_id on rows {missing_indices}. "
                f"When a sample has multiple rows, library_id must be specified for all rows."
            )


def _validate_no_mixed_input_types(df: pd.DataFrame) -> None:
    """Error if a sample_id has rows with different input_types."""
    for sample_id, group in df.groupby("sample_id"):
        input_types = group["input_type"].unique()
        if len(input_types) > 1:
            raise ValueError(
                f"Sample '{sample_id}' has mixed input_types: {', '.join(input_types)}. "
                f"All rows for a sample must have the same input_type."
            )


def _validate_no_duplicate_rows(df: pd.DataFrame) -> None:
    """Error if exact duplicate rows exist (same sample_id + library_id + input)."""
    cols = ["sample_id", "input"]
    if "library_id" in df.columns:
        cols.append("library_id")

    duplicates = df[df.duplicated(subset=cols, keep=False)]
    if not duplicates.empty:
        first_dup = duplicates.iloc[0]
        lib_info = f", library_id='{first_dup['library_id']}'" if "library_id" in cols else ""
        raise ValueError(
            f"Duplicate row detected: sample_id='{first_dup['sample_id']}'{lib_info}, "
            f"input='{first_dup['input']}'"
        )


def _fill_missing_library_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing library_id values with sample_id (only for single-row samples)."""
    df = df.copy()
    if "library_id" not in df.columns:
        df["library_id"] = df["sample_id"]
    else:
        mask = df["library_id"].isna() | (df["library_id"] == "")
        df.loc[mask, "library_id"] = df.loc[mask, "sample_id"]
    return df


def _assign_unit_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Assign sequential unit IDs (0, 1, 2...) within each sample."""
    df = df.copy()
    df["_unit_id"] = df.groupby("sample_id").cumcount()
    return df


def _parse_fastq_paths(df: pd.DataFrame) -> pd.DataFrame:
    """Parse semicolon-separated FASTQ paths into _fq1 and _fq2 columns."""
    df = df.copy()
    fq1_list = []
    fq2_list = []

    for _, row in df.iterrows():
        if row["input_type"] == "fastq":
            parts = row["input"].split(";")
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid fastq input for {row['sample_id']}: "
                    f"expected 'r1;r2' format, got '{row['input']}'"
                )
            fq1_list.append(parts[0].strip())
            fq2_list.append(parts[1].strip())
        else:
            fq1_list.append(None)
            fq2_list.append(None)

    df["_fq1"] = fq1_list
    df["_fq2"] = fq2_list
    return df


_validate_library_ids_for_multi_row_samples(samples)
_validate_no_mixed_input_types(samples)
_validate_no_duplicate_rows(samples)

samples = _fill_missing_library_ids(samples)
samples = _assign_unit_ids(samples)
samples = _parse_fastq_paths(samples)

SAMPLE_IDS = samples["sample_id"].unique().tolist()

wildcard_constraints:
    sample="[A-Za-z0-9_-]+",
    unit="[0-9]+",


def get_sample_units(sample_id: str) -> list[int]:
    """Get all unit IDs for a sample."""
    return samples.loc[samples["sample_id"] == sample_id, "_unit_id"].tolist()


def get_sample_libraries(sample_id: str) -> list[str]:
    """Get unique library IDs for a sample."""
    return samples.loc[samples["sample_id"] == sample_id, "library_id"].unique().tolist()


def get_units_for_library(sample_id: str, library_id: str) -> list[int]:
    """Get unit IDs for a specific sample+library combination."""
    mask = (samples["sample_id"] == sample_id) & (samples["library_id"] == library_id)
    return samples.loc[mask, "_unit_id"].tolist()


def get_unit_row(sample_id: str, unit_id: int) -> pd.Series:
    """Get the sample sheet row for a specific sample+unit."""
    mask = (samples["sample_id"] == sample_id) & (samples["_unit_id"] == unit_id)
    return samples.loc[mask].iloc[0]


def get_ref_path():
    """Get the reference genome path."""
    if REF_PATH:
        return REF_PATH
    elif REF_ACCESSION:
        return f"results/reference/{REF_FILE}"
    else:
        raise ValueError("No reference path or accession specified in config")


def get_reads(wildcards):
    """Get read files for a sample unit."""
    row = get_unit_row(wildcards.sample, int(wildcards.unit))

    if row["input_type"] == "fastq":
        r1 = row["_fq1"]
        r2 = row["_fq2"]

        if config.get("remote_reads", False):
            return {"r1": storage(r1), "r2": storage(r2)}

        if not (os.path.exists(r1) and os.path.exists(r2)):
            raise WorkflowError(
                f"FASTQ files not found for {wildcards.sample} unit {wildcards.unit}: {r1}, {r2}"
            )
        return {"r1": r1, "r2": r2}

    elif row["input_type"] == "srr":
        srr = row["input"]
        return {
            "r1": f"results/fastq/{wildcards.sample}/{srr}_1.fastq.gz",
            "r2": f"results/fastq/{wildcards.sample}/{srr}_2.fastq.gz",
        }

    else:
        raise ValueError(f"Unsupported input_type: {row['input_type']}")


def get_reads_fastp(wildcards):
    """Get reads for fastp, optionally sorted."""
    if config.get("sort_reads", False):
        return {
            "r1": f"results/sorted_reads/{wildcards.sample}/{wildcards.unit}_1.fastq.gz",
            "r2": f"results/sorted_reads/{wildcards.sample}/{wildcards.unit}_2.fastq.gz",
        }
    else:
        return get_reads(wildcards)


def get_read_group(wildcards):
    """Generate read group string for BWA."""
    row = get_unit_row(wildcards.sample, int(wildcards.unit))
    lib = row["library_id"]
    return rf"'@RG\tID:{lib}_{wildcards.unit}\tSM:{wildcards.sample}\tLB:{lib}\tPL:ILLUMINA'"


def get_unit_bams(wildcards):
    """Get all pre-merge BAM files for a sample."""
    units = get_sample_units(wildcards.sample)
    return expand(
        "results/bams/preMerge/{sample}/{unit}.bam",
        sample=wildcards.sample,
        unit=units,
    )


def get_dedup_input(wildcards):
    """Get BAM input for deduplication."""
    units = get_sample_units(wildcards.sample)

    if len(units) == 1:
        return {
            "bam": f"results/bams/preMerge/{wildcards.sample}/{units[0]}.bam",
            "bai": f"results/bams/preMerge/{wildcards.sample}/{units[0]}.bam.bai",
        }
    else:
        return {
            "bam": f"results/bams/merged/{wildcards.sample}.bam",
            "bai": f"results/bams/merged/{wildcards.sample}.bam.bai",
        }


def get_final_bam(wildcards):
    """Get final BAM path for a sample."""
    if config.get("mark_duplicates", True):
        return {
            "bam": f"results/bams/{wildcards.sample}.bam",
            "bai": f"results/bams/{wildcards.sample}.bam.bai",
        }
    else:
        return get_dedup_input(wildcards)


def get_coverage_d4_files(wildcards=None):
    return {
        "d4": expand("results/callable_sites/{sample}.per-base.d4.gz", sample=SAMPLE_IDS),
        "d4gzi": expand("results/callable_sites/{sample}.per-base.d4.gz.gzi", sample=SAMPLE_IDS),
    }


def get_mosdepth_summary_files():
    return expand("results/callable_sites/{sample}.mosdepth.summary.txt", sample=SAMPLE_IDS)


def get_sumstats_input(wildcards):
    aln = expand("results/summary_stats/{sample}_AlnSumMets.txt", sample=SAMPLE_IDS)
    cov = expand("results/summary_stats/{sample}_coverage.txt", sample=SAMPLE_IDS)
    fastp = expand("results/summary_stats/{sample}_fastp.out", sample=SAMPLE_IDS)

    out = {
        "alnSumMetsFiles": aln,
        "fastpFiles": fastp,
        "coverageFiles": cov,
    }

    if config.get("sentieon", False):
        out["insert_files"] = expand(
            "results/summary_stats/{sample}_insert_metrics.txt", sample=SAMPLE_IDS
        )
        out["qc_files"] = expand(
            "results/summary_stats/{sample}_qd_metrics.txt", sample=SAMPLE_IDS
        )
        out["mq_files"] = expand(
            "results/summary_stats/{sample}_mq_metrics.txt", sample=SAMPLE_IDS
        )
        out["gc_files"] = expand(
            "results/summary_stats/{sample}_gc_metrics.txt", sample=SAMPLE_IDS
        )
        out["gc_summary"] = expand(
            "results/summary_stats/{sample}_gc_summary.txt", sample=SAMPLE_IDS
        )

    return out


def collect_fastp_stats_input(wildcards):
    units = get_sample_units(wildcards.sample)
    return expand(
        "results/summary_stats/{sample}/{unit}.fastp.out",
        sample=wildcards.sample,
        unit=units,
    )


def get_big_temp(wildcards):
    """Get temp directory for rules needing large temp space."""
    if config.get("bigtmp"):
        base = config["bigtmp"].rstrip("/")
        suffix = "".join(random.choices(string.ascii_uppercase, k=12))
        return f"{base}/{suffix}/"
    else:
        return tempfile.gettempdir()


def setup_curlrc():
    """Add -L flag to curlrc for pyd4 remote access."""
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"

    if curlrc_path.exists():
        with curlrc_path.open("r+") as f:
            if "-L" not in f.read():
                f.write(f"\n{entry}\n")
    else:
        with curlrc_path.open("a+") as f:
            f.write(f"{entry}\n")


def cleanup_curlrc():
    """Remove snpArcher-added -L flag from curlrc."""
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"

    if curlrc_path.exists():
        with curlrc_path.open("r") as f:
            lines = f.readlines()
        new_lines = [line for line in lines if line.strip() != entry.strip()]
        if len(new_lines) == 0:
            curlrc_path.unlink()
        else:
            with curlrc_path.open("w") as f:
                f.writelines(new_lines)


def collectCovStats(covSumFiles):
    """Parse mosdepth coverage summary files."""
    sampleCov = {}

    for fn in covSumFiles:
        with open(fn, "r") as f:
            for line in f:
                if "mean" in line:
                    continue
                fields = line.split()
                chrom = fields[0]
                cov = float(fields[3])

                if chrom in sampleCov:
                    sampleCov[chrom].append(cov)
                else:
                    sampleCov[chrom] = [cov]

    covStats = {}
    for chrom in sampleCov:
        mean_cov = statistics.mean(sampleCov[chrom])
        try:
            std_cov = statistics.stdev(sampleCov[chrom])
        except statistics.StatisticsError:
            std_cov = "NA"
        covStats[chrom] = {"mean": mean_cov, "stdev": std_cov}

    return covStats


def collectFastpOutput(fastpFiles):
    """Parse fastp JSON output files."""
    FractionReadsPassFilter = defaultdict(float)
    NumReadsPassFilter = defaultdict(int)

    for fn in fastpFiles:
        sample = os.path.basename(fn).replace("_fastp.out", "")

        with open(fn, "r") as f:
            data = json.load(f)

        unfiltered = data["summary"]["before_filtering"]["total_reads"]
        pass_filter = data["summary"]["after_filtering"]["total_reads"]
        FractionReadsPassFilter[sample] = float(pass_filter / unfiltered) if unfiltered else 0
        NumReadsPassFilter[sample] = pass_filter

    return (FractionReadsPassFilter, NumReadsPassFilter)


def combine_fastp_files(fastpFiles, outputfile):
    """Merge multiple fastp outputs into single summary."""
    unfiltered = 0
    pass_filter = 0

    for fn in fastpFiles:
        with open(fn, "r") as f:
            data = json.load(f)
        unfiltered += data["summary"]["before_filtering"]["total_reads"]
        pass_filter += data["summary"]["after_filtering"]["total_reads"]

    out = {
        "summary": {
            "before_filtering": {"total_reads": unfiltered},
            "after_filtering": {"total_reads": pass_filter},
        }
    }

    with open(outputfile[0], "w") as f:
        json.dump(out, f)


def collectAlnSumMets(alnSumMetsFiles):
    """Parse alignment summary metrics files."""
    aln_metrics = defaultdict(dict)

    for fn in alnSumMetsFiles:
        sample = os.path.basename(fn).replace("_AlnSumMets.txt", "")

        with open(fn, "r") as f:
            lines = f.readlines()

        total_aligns = int(lines[0].split()[0])
        num_dups = int(lines[4].split()[0])
        percent_mapped = lines[7].split()[0].strip("%")
        percent_proper_paired = lines[14].split()[0].strip("%")

        aln_metrics[sample]["Total alignments"] = total_aligns
        aln_metrics[sample]["Percent Mapped"] = percent_mapped
        aln_metrics[sample]["Num Duplicates"] = num_dups
        aln_metrics[sample]["Percent Properly Paired"] = percent_proper_paired

    return aln_metrics


def collectCoverageMetrics(coverageFiles):
    """Parse samtools coverage output files."""
    SeqDepths = defaultdict(float)
    CoveredBases = defaultdict(float)

    for fn in coverageFiles:
        sample = os.path.basename(fn).replace("_coverage.txt", "")
        numSites = []
        covbases = 0
        depths = []

        with open(fn, "r") as f:
            for line in f:
                if not line.startswith("#rname"):
                    fields = line.split()
                    numSites.append(int(fields[2]) - int(fields[1]) + 1)
                    depths.append(float(fields[6]))
                    covbases += float(fields[4])

        total = sum(numSites)
        depthsMean = sum(d * n / total for d, n in zip(depths, numSites)) if total else 0
        SeqDepths[sample] = depthsMean
        CoveredBases[sample] = covbases

    return (SeqDepths, CoveredBases)


def collect_inserts(files):
    """Parse insert size metrics files."""
    med_inserts = defaultdict(float)
    med_insert_std = defaultdict(float)

    for file in files:
        sample = os.path.basename(file).replace("_insert_metrics.txt", "")

        with open(file, "r") as f:
            for i, line in enumerate(f):
                if i == 2:
                    fields = line.strip().split()
                    try:
                        med_inserts[sample] = fields[0]
                        med_insert_std[sample] = fields[1]
                    except IndexError:
                        continue

    return med_inserts, med_insert_std


def printBamSumStats(
    depths,
    covered_bases,
    aln_metrics,
    FractionReadsPassFilter,
    NumReadsPassFilter,
    out_file,
    med_insert_sizes=None,
    med_abs_insert_std=None,
):
    """Write final summary statistics TSV file."""
    sample_list = list(depths.keys())

    header = [
        "Sample",
        "Total_Reads",
        "Percent_mapped",
        "Num_duplicates",
        "Percent_properly_paired",
        "Fraction_reads_pass_filter",
        "NumReadsPassingFilters",
    ]

    if med_insert_sizes is not None:
        header.extend(["MedianInsertSize", "MedianAbsDev_InsertSize"])

    with open(out_file, "w") as f:
        print("\t".join(header), file=f)

        for samp in sample_list:
            row = [
                samp,
                str(aln_metrics[samp]["Total alignments"]),
                str(aln_metrics[samp]["Percent Mapped"]),
                str(aln_metrics[samp]["Num Duplicates"]),
                str(aln_metrics[samp]["Percent Properly Paired"]),
                str(FractionReadsPassFilter[samp]),
                str(NumReadsPassFilter[samp]),
            ]

            if med_insert_sizes is not None:
                row.extend([str(med_insert_sizes[samp]), str(med_abs_insert_std[samp])])

            print("\t".join(row), file=f)

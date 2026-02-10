import sys
import os
import tempfile
import random
import string
import statistics
import json
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
import snparcher_utils
from yaml import safe_load

# Can't be less than 7 cuz of min version in snakefile.
# Avoid pkg_resources/setuptools dependency by extracting the major version.
_snakemake_version = str(snakemake.__version__)
_snakemake_major_match = re.search(r"\d+", _snakemake_version)
if _snakemake_major_match is None:
    raise RuntimeError(f"Unable to parse snakemake version: {_snakemake_version}")
SNAKEMAKE_MAJOR = int(_snakemake_major_match.group(0))
SNAKEMAKE_VERSION = 8 if SNAKEMAKE_MAJOR >= 8 else 7
logger.warning(f"snpArcher: Using Snakemake {snakemake.__version__}")
if SNAKEMAKE_VERSION >= 8:
    DEFAULT_STORAGE_PREFIX = StorageSettings.default_storage_prefix if StorageSettings.default_storage_prefix is not None else ""
else:
    # backwards compatability w/ snakemake <= 7
    DEFAULT_STORAGE_PREFIX = workflow.default_remote_prefix
    if config["remote_reads"]:
        from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
        GS = GSRemoteProvider()
        GS_READS_PREFIX = config['remote_reads_prefix']

samples = snparcher_utils.parse_sample_sheet(config)


def get_output():
    
    if config["final_prefix"] == "":
        raise (WorkflowError("'final_prefix' is not set in config."))
    out = []
    genomes = samples["refGenome"].unique().tolist()
    sample_counts = samples.drop_duplicates(
        subset=["BioSample", "refGenome"]
    ).value_counts(
        subset=["refGenome"]
    )  # get BioSample for each refGenome
    out.extend
    if config["final_prefix"] == "":
        raise (WorkflowError("'final_prefix' is not set in config."))
    out = []
    
    sample_counts = samples.drop_duplicates(
        subset=["BioSample", "refGenome"]
    ).value_counts(
        subset=["refGenome"]
    )  # get BioSample for each refGenome
    
    for ref in genomes:
        # Workaround for Snakemake issue 2762. There is problem with running nested checkpoints in snakemake8. Adding mapfile in rule all forces gvcf interval checkpoint to run.
        # This is actually kind of a good thing to do since it makes dryrun more clear (shows all bam>gvcf jobs now).
        if not config["sentieon"]:
            out.extend(expand("results/{refGenome}/genomics_db_import/DB_mapfile.txt", refGenome=ref))
        out.extend(
            expand( "results/{refGenome}/{prefix}_raw.vcf.gz",refGenome=ref, prefix=config["final_prefix"]))
        out.extend(
            expand( "results/{refGenome}/summary_stats/{prefix}_bam_sumstats.txt", refGenome=ref, prefix=config["final_prefix"]))
        out.extend(
            expand("results/{refGenome}/{prefix}_callable_sites.bed", refGenome=ref, prefix=config["final_prefix"]))
        if sample_counts[ref] > 2:
            out.append(rules.qc_all.input)
        if "SampleType" in samples.columns:
            out.append(rules.postprocess_all.input)
            if all(
                i in samples["SampleType"].tolist() for i in ["ingroup", "outgroup"]
            ):
                out.append(rules.mk_all.input)
            if config["generate_trackhub"]:
                if not config["trackhub_email"]:
                    raise (
                        WorkflowError(
                            "If generating trackhub, you must provide an email in the config file."
                        )
                    )
                out.append(rules.trackhub_all.input)
    return out


def merge_bams_input(wc):
    return expand(
        "results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam",
        run=samples.loc[samples["BioSample"] == wc.sample]["Run"].tolist(),
    )

def setup_curlrc():
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"
    if curlrc_path.exists():
        with curlrc_path.open("r+") as f:
            if "-L" not in f.read():
                f.write(f"\n{entry}\n")
                logger.info(f"Added -L to {curlrc_path} for pyd4")
            
    else:
        with curlrc_path.open("a+") as f:
            f.write(f"{entry}\n")

def cleanup_curlrc():
    curlrc_path = Path("~/.curlrc").expanduser()
    marker = "# Added by snpArcher"
    entry = f"-L {marker}\n"
    logger.info(f"Removing -L we added from {curlrc_path}...")
    if curlrc_path.exists():
        with curlrc_path.open("r") as f:
            lines = f.readlines()
            # remove entry if its there
            new_lines = [line for line in lines if line.strip() != entry.strip()]
            if len(new_lines) == 0:
                # our entry was only thing there, we can delete .curlrc
                curlrc_path.unlink()
            else:
            # write back any options that were there.
                with curlrc_path.open("w") as f:
                    f.writelines(new_lines)

def get_ref(wildcards):
    
    if "refPath" in samples.columns:
        _refs = (
            samples.loc[(samples["refGenome"] == wildcards.refGenome)]["refPath"]
            .dropna()
            .unique()
            .tolist()
        )             
        if _refs:
            return _refs
    # if not user-specified refpath, force MissingInputError in copy_ref with dummyfile, which allows download_ref to run b/c of ruleorder.
    logger.warning(f"snpArcher: refPath specified in sample sheet header, but no path provided for refGenome '{wildcards.refGenome}'\n" + 
                    f"Will try to download '{wildcards.refGenome}' from NCBI. If this is a genome accession, you can ignore this warning.")
    return []

def get_bams(wc):
    out = {"bam": None, "bai": None}
    if config["mark_duplicates"]:
        out["bam"] = "results/{refGenome}/bams/{sample}_final.bam"
        out["bai"] = "results/{refGenome}/bams/{sample}_final.bam.bai"
        return out
    else:
        return dedup_input(wc)

def sentieon_combine_gvcf_cmd_line(wc):
    gvcfs = sentieon_combine_gvcf_input(wc)["gvcfs"]
    return " ".join(["-v " + gvcf for gvcf in gvcfs])


def get_interval_gvcfs(wc):
    checkpoint_output = checkpoints.create_gvcf_intervals.get(**wc).output[0]
    with checkpoint_output.open() as f:
        lines = [l.strip() for l in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]
    list_numbers = [f.replace("-scattered.interval_list", "") for f in list_files]
    gvcfs = expand(
        "results/{{refGenome}}/interval_gvcfs/{{sample}}/{l}.raw.g.vcf.gz",
        l=list_numbers,
    )

    return gvcfs


def get_interval_gvcfs_idx(wc):
    tbis = [f + ".tbi" for f in get_interval_gvcfs(wc)]
    return tbis


def get_db_interval_count(wc):
    _samples = (
        samples.loc[(samples["refGenome"] == wc.refGenome)]["BioSample"]
        .unique()
        .tolist()
    )
    out = max(
        int((config["db_scatter_factor"]) * len(_samples) * config["num_gvcf_intervals"]),1)
    return out


def get_interval_vcfs(wc):
    checkpoint_output = checkpoints.create_db_intervals.get(**wc).output[0]
    with checkpoint_output.open() as f:
        lines = [l.strip() for l in f.readlines()]
    list_files = [os.path.basename(x) for x in lines]

    list_numbers = [f.replace("-scattered.interval_list", "") for f in list_files]
    vcfs = expand("results/{{refGenome}}/vcfs/intervals/filtered_L{l}.vcf.gz", l=list_numbers)

    return vcfs


def get_interval_vcfs_idx(wc):
    tbis = [f + ".tbi" for f in get_interval_vcfs(wc)]
    return tbis


def get_gvcfs_db(wc):
    _samples = samples.loc[(samples["refGenome"] == wc.refGenome)]["BioSample"].unique().tolist()
    gvcfs = expand("results/{{refGenome}}/gvcfs/{sample}.g.vcf.gz", sample=_samples)
    tbis = expand("results/{{refGenome}}/gvcfs/{sample}.g.vcf.gz.tbi", sample=_samples)
    return {"gvcfs": gvcfs, "tbis": tbis}


def dedup_input(wc):
    runs = samples.loc[samples["BioSample"] == wc.sample]["Run"].tolist()

    if len(runs) == 1:
        bam = expand("results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam", run=runs)
        bai = expand("results/{{refGenome}}/bams/preMerge/{{sample}}/{run}.bam.bai", run=runs)
    else:
        bam = "results/{refGenome}/bams/postMerge/{sample}.bam"
        bai = "results/{refGenome}/bams/postMerge/{sample}.bam.bai"
    return {"bam": bam, "bai": bai}


def sentieon_combine_gvcf_input(wc):
    _samples = samples["BioSample"].unique().tolist()
    gvcfs = expand("results/{{refGenome}}/gvcfs/{sample}.g.vcf.gz", sample=_samples)
    tbis = expand("results/{{refGenome}}/gvcfs/{sample}.g.vcf.gz.tbi", sample=_samples)
    return {"gvcfs": gvcfs, "tbis": tbis}


def get_reads(wc):
    """Returns local read files if present. Defaults to SRR if no local reads in sample sheet."""
    row = samples.loc[samples["Run"] == wc.run]
    r1 = f"results/data/fastq/{wc.refGenome}/{wc.sample}/{wc.run}_1.fastq.gz"
    r2 = f"results/data/fastq/{wc.refGenome}/{wc.sample}/{wc.run}_2.fastq.gz"
    if "fq1" in samples.columns and "fq2" in samples.columns:
        if row["fq1"].notnull().any() and row["fq2"].notnull().any():
            r1 = row.fq1.item()
            r2 = row.fq2.item()
            if config["remote_reads"]:
                if SNAKEMAKE_VERSION>=8:
                    # remote read path must have full remote prefix, eg: gs://reads_bucket/sample1/...
                    # depends on snakemake>8 to figure out proper remote provider from prefix using storage()
                    return {"r1": storage(r1), "r2": storage(r2)}
                else:
                    return get_remote_reads(wc)
            if os.path.exists(row.fq1.item()) and os.path.exists(row.fq2.item()):
                return {"r1": r1, "r2": r2}
            else:
                raise WorkflowError(
                    f"fq1 and fq2 specified for {wc.sample}, but files were not found."
                )
        else:
            # this allows mixed srr and user-specified paths for reads
            return {"r1": r1, "r2": r2}
    else:
        return {"r1": r1, "r2": r2}

def get_reads_fastp(wc):
    if config.get("sort_reads", False):
        return {"r1":"results/{refGenome}/sorted_reads/{sample}/{run}_1.fastq.gz", "r2":"results/{refGenome}/sorted_reads/{sample}/{run}_2.fastq.gz"}
    else:
        return get_reads(wc)

def get_remote_reads(wildcards):
    """Use this for reads on a different remote bucket than the default. For backwards compatibility."""
    # print(wildcards)
    row = samples.loc[samples["Run"] == wildcards.run]
    r1 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq1.item()))
    r2 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq2.item()))
    return {"r1": r1, "r2": r2}

def collect_fastp_stats_input(wc):
    return expand(
        "results/{{refGenome}}/summary_stats/{{sample}}/{run}.fastp.out",
        run=samples.loc[samples["BioSample"] == wc.sample]["Run"].tolist(),
    )


def get_read_group(wc):
    """Denote sample name and library_id in read group."""
    libname = samples.loc[samples["Run"] == wc.run]["LibraryName"].tolist()[0]
    return r"'@RG\tID:{lib}\tSM:{sample}\tLB:{lib}\tPL:ILLUMINA'".format(
        sample=wc.sample, lib=libname
    )


def get_input_sumstats(wildcards):
    _samples = samples.loc[(samples["refGenome"] == wildcards.refGenome)]["BioSample"].unique().tolist()
    aln = expand("results/{{refGenome}}/summary_stats/{sample}_AlnSumMets.txt", sample=_samples)
    cov = expand("results/{{refGenome}}/summary_stats/{sample}_coverage.txt", sample=_samples)
    fastp = expand("results/{{refGenome}}/summary_stats/{sample}_fastp.out", sample=_samples)
    insert = expand("results/{{refGenome}}/summary_stats/{sample}_insert_metrics.txt",sample=_samples)
    qd = expand("results/{{refGenome}}/summary_stats/{sample}_qd_metrics.txt", sample=_samples)
    mq = expand("results/{{refGenome}}/summary_stats/{sample}_mq_metrics.txt", sample=_samples)
    gc = expand("results/{{refGenome}}/summary_stats/{sample}_gc_metrics.txt", sample=_samples)
    gc_summary = expand("results/{{refGenome}}/summary_stats/{sample}_gc_summary.txt", sample=_samples)
    if config["sentieon"]:
        out = {
            "alnSumMetsFiles": aln,
            "fastpFiles": fastp,
            "coverageFiles": cov,
            "insert_files": insert,
            "qc_files": qd,
            "mq_files": mq,
            "gc_files": gc,
            "gc_summary": gc_summary,
        }
        return out
    else:
        out = {
            "alnSumMetsFiles": aln,
            "fastpFiles": fastp,
            "coverageFiles": cov,
        }
        return out


def get_input_for_mapfile(wildcards):
    sample_names = samples.loc[(samples["refGenome"] == wildcards.refGenome)]["BioSample"].unique().tolist()
    if config["intervals"]:
        return expand("results/{{refGenome}}/gvcfs_norm/{sample}.g.vcf.gz", sample=sample_names)
    else:
        return expand("results/{{refGenome}}/gvcfs/{sample}.g.vcf.gz", sample=sample_names)


def get_input_for_coverage(wildcards):
    # Gets the correct sample given the organism and reference genome for the bedgraph merge step
    _samples = samples.loc[(samples["refGenome"] == wildcards.refGenome)]["BioSample"].unique().tolist()
    
    d4files = expand("results/{{refGenome}}/callable_sites/{sample}.per-base.d4.gz", sample=_samples)
    d4gzifiles = expand("results/{{refGenome}}/callable_sites/{sample}.per-base.d4.gz.gzi", sample=_samples)
    return {"d4":d4files, "d4gzi":d4gzifiles}


def get_input_covstats(wildcards):
    # Gets the correct sample given the organism and reference genome for the bedgraph merge step
    _samples = samples.loc[(samples["refGenome"] == wildcards.refGenome)]["BioSample"].unique().tolist()
    
    covstats = expand("results/{{refGenome}}/callable_sites/{sample}.mosdepth.summary.txt",sample=_samples)
    return {"covStatFiles": covstats}


def get_bedgraphs(wildcards):
    """Snakemake seems to struggle with unpack() and default_remote_prefix. So we have to do this one by one."""
    _samples = (
        samples.loc[
            (samples["Organism"] == wildcards.Organism)
            & (samples["refGenome"] == wildcards.refGenome)
        ]["BioSample"]
        .unique()
        .tolist()
    )
    bedgraphFiles = expand(
        config["output"]
        + "{{Organism}}/{{refGenome}}/"
        + config["bamDir"]
        + "preMerge/{sample}"
        + ".sorted.bg",
        sample=_samples,
    )
    return bedgraphFiles


def get_big_temp(wildcards):
    """Sets a temp dir for rules that need more temp space that is typical on some cluster environments. Defaults to system temp dir."""
    if config["bigtmp"]:
        if config["bigtmp"].endswith("/"):
            return (
                config["bigtmp"]
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
        else:
            return (
                config["bigtmp"]
                + "/"
                + "".join(random.choices(string.ascii_uppercase, k=12))
                + "/"
            )
    else:
        return tempfile.gettempdir()


def collectCovStats(covSumFiles):
    covStats = {}
    sampleCov = {}

    for fn in covSumFiles:
        f = open(fn, "r")
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

    for chr in sampleCov:
        mean_cov = statistics.mean(sampleCov[chr])
        try:
            std_cov = statistics.stdev(sampleCov[chr])
        except:
            std_cov = "NA"
        covStats[chr] = {"mean": mean_cov, "stdev": std_cov}

    return covStats


def collectFastpOutput(fastpFiles):
    FractionReadsPassFilter = defaultdict(float)
    NumReadsPassFilter = defaultdict(int)

    for fn in fastpFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_fastp.out", "")
        
        
        with open(fn, "r") as f:
            data = json.load(f)
        unfiltered = data["summary"]["before_filtering"]["total_reads"]
        pass_filter = data["summary"]["after_filtering"]["total_reads"]
        FractionReadsPassFilter[sample] = float(pass_filter / unfiltered)
        NumReadsPassFilter[sample] = pass_filter

    return (FractionReadsPassFilter, NumReadsPassFilter)

def combine_fastp_files(fastpFiles, outputfile):
    unfiltered = 0
    pass_filter = 0
    for fn in fastpFiles:
        with open(fn, "r") as f:
            data = json.load(f)
        unfiltered += data["summary"]["before_filtering"]["total_reads"]
        pass_filter += data["summary"]["after_filtering"]["total_reads"]
    out = {"summary": {"before_filtering":{"total_reads": unfiltered}, "after_filtering":{"total_reads": pass_filter}}}
    with open(outputfile[0], "w") as f:
        json.dump(out, f)

def collectAlnSumMets(alnSumMetsFiles):
    aln_metrics = defaultdict(dict)
    for fn in alnSumMetsFiles:
        sample = os.path.basename(fn)
        sample = sample.replace("_AlnSumMets.txt", "")
        with open(fn, "r") as f:
            lines = f.readlines()
            total_aligns = int(lines[0].split()[0])
            num_dups = int(lines[4].split()[0])
            num_mapped = int(lines[6].split()[0])
            percent_mapped = lines[7].split()[0].strip("%")
            percent_proper_paired = lines[14].split()[0].strip("%")

            percent_dups = num_dups / total_aligns if total_aligns != 0 else 0
            aln_metrics[sample]["Total alignments"] = total_aligns
            aln_metrics[sample]["Percent Mapped"] = percent_mapped
            aln_metrics[sample]["Num Duplicates"] = num_dups
            aln_metrics[sample]["Percent Properly Paired"] = percent_proper_paired

    return aln_metrics


def collectCoverageMetrics(coverageFiles):
    SeqDepths = defaultdict(float)
    CoveredBases = defaultdict(float)

    for fn in coverageFiles:
        # these files contain coverage data by scaffold; take weighted average
        sample = os.path.basename(fn)
        sample = sample.replace("_coverage.txt", "")
        numSites = []
        covbases = 0
        depths = []  # samtools coverage function prints depth per scaffold
        f = open(fn, "r")
        for line in f:
            if not line.startswith("#rname"):
                line = line.split()
                numSites.append(int(line[2]) - int(line[1]) + 1)
                depths.append(float(line[6]))
                covbases += float(line[4])
        f.close()
        total = sum(numSites)
        depthsMean = 0
        for i in range(len(depths)):
            depthsMean += depths[i] * numSites[i] / (total)
        SeqDepths[sample] = depthsMean
        CoveredBases[sample] = covbases
    return (SeqDepths, CoveredBases)


def collect_inserts(files):
    med_inserts = defaultdict(float)
    med_insert_std = defaultdict(float)

    for file in files:
        sample = os.path.basename(file)
        sample = sample.replace("_insert_metrics.txt", "")
        with open(file, "r") as f:
            for i, line in enumerate(f):
                if i == 2:
                    line = line.strip().split()
                    try:
                        med_inserts[sample] = line[0]
                        med_insert_std[sample] = line[1]
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
    samples = depths.keys()
    if med_insert_sizes is None:
        with open(out_file, "w") as f:
            print(
                "Sample\tTotal_Reads\tPercent_mapped\tNum_duplicates\tPercent_properly_paired\tFraction_reads_pass_filter\tNumReadsPassingFilters",
                file=f,
            )
            for samp in samples:
                print(
                    samp,
                    "\t",
                    aln_metrics[samp]["Total alignments"],
                    "\t",
                    aln_metrics[samp]["Percent Mapped"],
                    "\t",
                    aln_metrics[samp]["Num Duplicates"],
                    "\t",
                    aln_metrics[samp]["Percent Properly Paired"],
                    "\t",
                    FractionReadsPassFilter[samp],
                    "\t",
                    NumReadsPassFilter[samp],
                    file=f,
                )
    else:
        with open(out_file, "w") as f:
            print(
                "Sample\tTotal_Reads\tPercent_mapped\tNum_duplicates\tPercent_properly_paired\tFraction_reads_pass_filter\tNumReadsPassingFilters\tMedianInsertSize\tMedianAbsDev_InsertSize",
                file=f,
            )
            for samp in samples:
                print(
                    samp,
                    "\t",
                    aln_metrics[samp]["Total alignments"],
                    "\t",
                    aln_metrics[samp]["Percent Mapped"],
                    "\t",
                    aln_metrics[samp]["Num Duplicates"],
                    "\t",
                    aln_metrics[samp]["Percent Properly Paired"],
                    "\t",
                    FractionReadsPassFilter[samp],
                    "\t",
                    NumReadsPassFilter[samp],
                    "\t",
                    med_insert_sizes[samp],
                    "\t",
                    med_abs_insert_std[samp],
                    file=f,
                )

import os
import pandas as pd

USE_SENTIEON = config["sentieon"]["enabled"]

# Load and prepare samples dataframe
samples_df = pd.read_csv(config["samples"])

# Set defaults
samples_df["library_id"] = samples_df["library_id"].fillna(samples_df["sample_id"])
if "mark_duplicates" not in samples_df.columns:
    samples_df["mark_duplicates"] = True
else:
    samples_df["mark_duplicates"] = samples_df["mark_duplicates"].fillna(True).astype(bool)

# Validate: if sample_id appears multiple times, library_id must be unique
for sample_id, group in samples_df.groupby("sample_id"):
    if len(group) > 1:
        library_ids = group["library_id"].tolist()
        if len(library_ids) != len(set(library_ids)):
            raise ValueError(
                f"Sample '{sample_id}' has multiple rows with duplicate library_id values: {library_ids}. "
                f"Each row for a sample must have a unique library_id."
            )

samples_df = samples_df.set_index("sample_id", drop=False)

REF_NAME = config["reference"]["name"]

REF_FILES = {
    "ref": f"results/reference/{REF_NAME}.fa.gz",
    "ref_fai": f"results/reference/{REF_NAME}.fa.gz.fai",
    "ref_dict": f"results/reference/{REF_NAME}.dict",
}

REF_BWA_IDX = multiext(
    f"results/reference/{REF_NAME}.fa.gz",
    ".sa", ".pac", ".bwt", ".ann", ".amb",
)

# Sample lists by input type
# Sample lists by input type
SAMPLES_ALL = samples_df.index.tolist()
SAMPLES_NEED_ALIGNMENT = samples_df[samples_df["input_type"].isin(["srr", "fastq"])].index.tolist()
SAMPLES_NEED_HAPLOTYPECALLER = samples_df[samples_df["input_type"] != "gvcf"].index.tolist()
SAMPLES_WITH_BAM = samples_df[samples_df["input_type"].isin(["srr", "fastq", "bam"])].index.tolist()
SAMPLES_WITH_FASTQ = samples_df[samples_df["input_type"].isin(["srr", "fastq"])].index.tolist()
SAMPLES_SRR = samples_df[samples_df["input_type"] == "srr"].index.tolist()



def sample_has_multiple_libraries(sample):
    """Check if sample has multiple libraries (multiple rows in sample sheet)."""
    return len(samples_df[samples_df["sample_id"] == sample]) > 1


def get_sample_libraries(sample):
    """Get all library IDs for a sample."""
    return samples_df[samples_df["sample_id"] == sample]["library_id"].tolist()


def get_final_bam(sample):
    """Get final BAM path for a sample."""
    input_type = samples_df.loc[sample, "input_type"]
    
    if input_type == "bam":
        return samples_df.loc[sample, "input"]
    
    mark_dups = samples_df.loc[sample, "mark_duplicates"]
    
    if mark_dups:
        return f"results/bams/markdup/{sample}.bam"
    elif sample_has_multiple_libraries(sample):
        return f"results/bams/merged/{sample}.bam"
    else:
        lib = get_sample_libraries(sample)[0]
        return f"results/bams/raw/{sample}/{lib}.bam"


def get_final_gvcf(sample):
    """Get final gVCF path for a sample based on input type."""
    input_type = samples_df.loc[sample, "input_type"]
    
    # User provided gVCF directly
    if input_type == "gvcf":
        return samples_df.loc[sample, "input"]
    
    return f"results/gvcfs/{sample}.g.vcf.gz"


def java_mem_overhead(wildcards, resources):
    mem_mb = resources.mem_mb
    heap_mb = int(mem_mb * 0.9)
    return f"-Xmx{heap_mb}m"
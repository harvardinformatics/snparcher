# Sample Sheet Specification

The sample sheet is a CSV file that defines the input samples for snpArcher.

## Format

CSV with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Sample identifier. Used in final VCF and output files. |
| `input_type` | Yes | Type of input data: `fastq`, `srr`, `bam`, `gvcf` |
| `input` | Yes | Input path(s) or SRR accession |
| `library_id` | No | Library identifier for duplicate marking |

## Column Details

### sample_id

The sample identifier that will appear in the final VCF and all output files.

**Validation**: Must match pattern `^[A-Za-z0-9_-]+$` (alphanumeric, underscores, hyphens only).

```csv
# Valid
sample_id
sample1
my-sample
Sample_001

# Invalid
sample id      # spaces not allowed
sample.1       # periods not allowed
sample@home    # special characters not allowed
```

### input_type

Specifies the type of input data.

| Value | Description | Status |
|-------|-------------|--------|
| `fastq` | Paired-end FASTQ files | Supported |
| `srr` | NCBI SRA accession | Supported |
| `bam` | Aligned BAM file | Planned |
| `gvcf` | GATK GVCF file | Planned |

### input

The input path(s) or accession, depending on `input_type`.

**For `fastq`**: Semicolon-separated paths to R1 and R2 files:

```csv
input_type,input
fastq,/path/to/sample_R1.fq.gz;/path/to/sample_R2.fq.gz
```

**For `srr`**: SRA run accession:

```csv
input_type,input
srr,SRR12345678
```

### library_id

Library identifier used for read group assignment and duplicate marking.

**Purpose**: The `library_id` is written into the BAM read group (`LB` tag) and used by duplicate marking tools to identify which reads came from the same library preparation.

**Behavior**:

- If not specified, defaults to `sample_id`
- Reads with the same `library_id` are considered the same library for deduplication
- Different `library_id` values = different library preparations = deduplicated separately

**Validation**: Must match pattern `^[A-Za-z0-9_-]*$` (alphanumeric, underscores, hyphens only).

## Handling Multiple Sequencing Runs

When a sample has been sequenced multiple times, each sequencing run should be a separate row in the sample sheet. snpArcher will:

1. **Map each row independently** to produce per-run BAM files
2. **Merge all BAMs** for the same `sample_id` into a single BAM
3. **Mark duplicates** using the `library_id` to identify which reads share a library

!!! note "Internal unit tracking"
    snpArcher automatically assigns an internal unit identifier to each row for intermediate file paths. You do not need to specify run IDs manually.

### Same library, multiple sequencing runs

When the same library was sequenced multiple times (e.g., to increase depth):

```csv
sample_id,input_type,input,library_id
sample_A,fastq,/data/sample_A_lane1_R1.fq.gz;/data/sample_A_lane1_R2.fq.gz,lib_A
sample_A,fastq,/data/sample_A_lane2_R1.fq.gz;/data/sample_A_lane2_R2.fq.gz,lib_A
```

**Result**: Both runs are merged, and duplicates are marked together (same library = duplicates can span runs).

### Different libraries, same sample

When multiple library preparations were made from the same sample:

```csv
sample_id,input_type,input,library_id
sample_A,fastq,/data/sample_A_libX_R1.fq.gz;/data/sample_A_libX_R2.fq.gz,lib_X
sample_A,fastq,/data/sample_A_libY_R1.fq.gz;/data/sample_A_libY_R2.fq.gz,lib_Y
```

**Result**: Both runs are merged into one BAM, but duplicates are marked separately per library (different libraries = independent duplicate marking).

### Complex example

Consider 2 samples with different sequencing scenarios:

| sample_id | library_id | Description |
|-----------|------------|-------------|
| sample_A | lib_A_1 | First library, lane 1 |
| sample_A | lib_A_1 | First library, lane 2 (same lib, re-sequenced) |
| sample_A | lib_A_2 | Second library |
| sample_B | (default) | Single run, single library |

```csv
sample_id,input_type,input,library_id
sample_A,fastq,/data/A_lib1_lane1_R1.fq.gz;/data/A_lib1_lane1_R2.fq.gz,lib_A_1
sample_A,fastq,/data/A_lib1_lane2_R1.fq.gz;/data/A_lib1_lane2_R2.fq.gz,lib_A_1
sample_A,fastq,/data/A_lib2_R1.fq.gz;/data/A_lib2_R2.fq.gz,lib_A_2
sample_B,fastq,/data/B_R1.fq.gz;/data/B_R2.fq.gz,
```

**Processing**:

- `sample_A`: 3 BAMs merged → duplicates marked (lib_A_1 together, lib_A_2 separate) → final BAM
- `sample_B`: 1 BAM → duplicates marked → final BAM

## Metadata Columns

Additional columns are allowed and will be passed through to the QC module if applicable:

```csv
sample_id,input_type,input,lat,long,population
sample1,fastq,/data/s1_R1.fq.gz;/data/s1_R2.fq.gz,42.3601,-71.0589,Boston
sample2,fastq,/data/s2_R1.fq.gz;/data/s2_R2.fq.gz,40.7128,-74.0060,NYC
```

## Validation

The sample sheet is validated at workflow start. Validation includes schema checks and semantic rules.

### Schema Validation

Basic structure validation using JSON Schema:

- Required columns are present (`sample_id`, `input_type`, `input`)
- `sample_id` matches allowed pattern (alphanumeric, underscores, hyphens)
- `input_type` is a valid value (`fastq`, `srr`, `bam`, `gvcf`)
- `input` is not empty

### Semantic Validation Rules

**1. Multi-row samples require explicit library_id**

If a `sample_id` appears on multiple rows, ALL rows for that sample must have an explicit `library_id`. This prevents incorrect duplicate marking.

```csv
# ERROR - sample_A has 2 rows but no library_id
sample_id,input_type,input
sample_A,fastq,/data/run1_R1.fq.gz;/data/run1_R2.fq.gz
sample_A,fastq,/data/run2_R1.fq.gz;/data/run2_R2.fq.gz

# OK - single-row sample, library_id defaults to sample_id
sample_id,input_type,input
sample_B,fastq,/data/B_R1.fq.gz;/data/B_R2.fq.gz

# OK - multi-row sample with explicit library_id on all rows
sample_id,input_type,input,library_id
sample_A,fastq,/data/run1_R1.fq.gz;/data/run1_R2.fq.gz,lib1
sample_A,fastq,/data/run2_R1.fq.gz;/data/run2_R2.fq.gz,lib1
```

**2. No mixed input types per sample**

All rows for a `sample_id` must have the same `input_type`. You cannot mix FASTQ and BAM inputs for the same sample.

```csv
# ERROR - sample_A has mixed input_types (fastq and bam)
sample_id,input_type,input,library_id
sample_A,fastq,/data/run1_R1.fq.gz;/data/run1_R2.fq.gz,lib1
sample_A,bam,/data/existing.bam,lib2
```

**3. No duplicate rows**

Exact duplicate rows (same `sample_id`, `library_id`, and `input`) are not allowed.

### Viewing Validation Errors

To see validation errors, run a dry-run:

```bash
snakemake -n --configfile config/config.yaml
```

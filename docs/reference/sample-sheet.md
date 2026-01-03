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

Optional library identifier used for duplicate marking. 

- If a `sample_id` appears on multiple rows (multiple sequencing runs), `library_id` distinguishes them
- If not specified, defaults to `sample_id`
- Samples with the same `library_id` are considered the same library for deduplication

## Examples

### Basic usage

```csv
sample_id,input_type,input
sample1,fastq,/data/sample1_R1.fq.gz;/data/sample1_R2.fq.gz
sample2,fastq,/data/sample2_R1.fq.gz;/data/sample2_R2.fq.gz
sample3,srr,SRR12345678
```

### Multiple sequencing runs per sample

When a sample was sequenced multiple times:

```csv
sample_id,input_type,input,library_id
sampleA,fastq,/data/sampleA_run1_R1.fq.gz;/data/sampleA_run1_R2.fq.gz,libA
sampleA,fastq,/data/sampleA_run2_R1.fq.gz;/data/sampleA_run2_R2.fq.gz,libA
sampleB,fastq,/data/sampleB_R1.fq.gz;/data/sampleB_R2.fq.gz,
```

In this example:
- `sampleA` has two sequencing runs from the same library (`libA`)
- Both runs will be merged and deduplicated together
- `sampleB` has one run; `library_id` defaults to `sampleB`

### Multiple libraries per sample

When a sample has multiple library preparations:

```csv
sample_id,input_type,input,library_id
sampleA,fastq,/data/sampleA_lib1_R1.fq.gz;/data/sampleA_lib1_R2.fq.gz,lib1
sampleA,fastq,/data/sampleA_lib2_R1.fq.gz;/data/sampleA_lib2_R2.fq.gz,lib2
```

In this example:
- `sampleA` has two different library preparations
- Each library is deduplicated separately, then merged

## Metadata Columns

Additional columns are allowed and will be passed through to the QC module if applicable:

```csv
sample_id,input_type,input,lat,long,population
sample1,fastq,/data/s1_R1.fq.gz;/data/s1_R2.fq.gz,42.3601,-71.0589,Boston
sample2,fastq,/data/s2_R1.fq.gz;/data/s2_R2.fq.gz,40.7128,-74.0060,NYC
```

## Validation

The sample sheet is validated at workflow start using JSON Schema. Validation checks:

- Required columns are present
- `sample_id` matches allowed pattern
- `input_type` is a valid value
- `input` is not empty

To see validation errors, run:

```bash
snakemake -n --configfile config/config.yaml
```

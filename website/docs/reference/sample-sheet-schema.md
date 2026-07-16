# Sample sheet schema

CSV file specifying samples and their input data.
Validated at startup against `workflow/schemas/samples.schema.yaml`.

## Required columns

| Column | Type | Constraints | Description |
|--------|------|-------------|-------------|
| `sample_id` | string | 1-80 characters, pattern `^[A-Za-z0-9._-]+$` | Unique sample identifier. Used as the key across the entire workflow. |
| `input_type` | string | `srr`, `fastq`, `bam`, `gvcf` | Type of input data for this row. |
| `input` | string | Non-empty. For `fastq`, must contain a semicolon. | Input path, SRA accession, or semicolon-separated FASTQ pair. Interpretation depends on `input_type` (see below). |

## Optional columns

| Column | Type | Default | Constraints | Description |
|--------|------|---------|-------------|-------------|
| `library_id` | string | Same as `sample_id` | 1-80 characters, pattern `^[A-Za-z0-9._-]+$` | Library identifier. Used to group rows within a sample for duplicate marking and read group assignment. |
| `mark_duplicates` | boolean | Value of `reads.mark_duplicates` in config | See [boolean parsing](#boolean-parsing) | Whether to mark duplicates for this sample. Must be consistent across all rows with the same `sample_id`. |

## Input type semantics

The meaning of the `input` column depends on `input_type`:

| `input_type` | `input` interpretation | Multiple rows per sample | Notes |
|--------------|----------------------|:------------------------:|-------|
| `srr` | SRA run accession (e.g., `SRR12345678`) | yes | Each row represents one SRA run. FASTQ files are downloaded automatically. |
| `fastq` | Semicolon-separated pair of FASTQ paths: `R1;R2` | yes | Paths can be absolute or relative to the working directory. Files must be gzipped (`.fastq.gz`). |
| `bam` | Path to an aligned BAM file | no | Exactly one row per sample. The BAM is used directly (no alignment step). |
| `gvcf` | Path to a gVCF file | no | Exactly one row per sample. Skips alignment and per-sample variant calling. |

!!! warning "gVCF compatibility"
    `gvcf` input is only supported with the `gatk` and `sentieon` variant calling tools.
    It is rejected at startup when `variant_calling.tool` is `bcftools`, `deepvariant`, or `parabricks`.

## Multiple rows per sample

A single `sample_id` can appear on multiple rows to represent multiple sequencing runs or libraries.
The rules are:

- **`srr` and `fastq`:** Multiple rows are allowed. Each row represents one run or lane.
- **`bam` and `gvcf`:** Exactly one row per `sample_id`.
- **`input_type` consistency:** All rows for a given `sample_id` must have the same `input_type`.
- **`mark_duplicates` consistency:** All rows for a given `sample_id` must have the same `mark_duplicates` value (or leave it blank on all rows).
- **`library_id`:** Rows with the same `library_id` are treated as the same library (e.g., multiple lanes). Rows with different `library_id` values are treated as different libraries. Within each `(sample_id, library_id)` group, rows are assigned ordinal input units (`u1`, `u2`, ...) internally.

## Boolean parsing

The `mark_duplicates` column accepts the following values (case-insensitive):

| True values | False values |
|-------------|--------------|
| `true`, `t`, `yes`, `y`, `1` | `false`, `f`, `no`, `n`, `0` |

Missing or `NA` values fall back to the global `reads.mark_duplicates` config setting.
Invalid values cause an error with the row number reported.

## Optional metadata sheet

A separate CSV file can be provided via the `sample_metadata` config key.
This sheet is validated against `workflow/schemas/sample_metadata.schema.yaml`.
Additional columns beyond those listed below are permitted (the schema allows `additionalProperties`).

| Column | Type | Required | Default | Constraints | Used by |
|--------|------|:--------:|---------|-------------|---------|
| `sample_id` | string | yes | | Must match a `sample_id` in the main samples CSV. Pattern `^[A-Za-z0-9._-]+$`. | All modules |
| `exclude` | boolean | no | `false` | Boolean-like (same parsing as `mark_duplicates`). | Postprocess: excluded samples are removed from filtered VCFs. |
| `lat` | number | no | | -90 to 90 | QC: latitude for geographic map panels. |
| `long` | number | no | | -180 to 180 | QC: longitude for geographic map panels. |

## Examples

### Minimal sample sheet (SRA accessions)

```csv
sample_id,input_type,input
bird_A,srr,SRR12345678
bird_B,srr,SRR12345679
bird_C,srr,SRR12345680
```

### Local FASTQ files with library IDs

```csv
sample_id,input_type,input,library_id,mark_duplicates
bird_A,fastq,/data/bird_A_lane1_R1.fq.gz;/data/bird_A_lane1_R2.fq.gz,lib1,true
bird_A,fastq,/data/bird_A_lane2_R1.fq.gz;/data/bird_A_lane2_R2.fq.gz,lib1,true
bird_B,fastq,/data/bird_B_R1.fq.gz;/data/bird_B_R2.fq.gz,lib1,true
```

### Mixed input types

```csv
sample_id,input_type,input
bird_A,srr,SRR12345678
bird_A,srr,SRR12345679
bird_B,fastq,/data/bird_B_R1.fq.gz;/data/bird_B_R2.fq.gz
bird_C,bam,/data/bird_C.bam
bird_D,gvcf,/data/bird_D.g.vcf.gz
```

### Metadata sheet

```csv
sample_id,exclude,lat,long
bird_A,false,42.38,-71.12
bird_B,true,42.36,-71.06
bird_C,false,43.07,-70.76
bird_D,false,,
```

## v1 to v2 migration

!!! warning "v1 sample sheets are not compatible with v2"
    snpArcher v2 uses a different sample sheet schema.
    If you are migrating from v1, you must convert your sample sheet.

| v1 column | v2 column | Notes |
|-----------|-----------|-------|
| `BioSample` | `sample_id` | Must match pattern `^[A-Za-z0-9._-]+$`. |
| `Run` (SRR/ERR/DRR) | `input` | Set `input_type` to `srr`. |
| `fq1`, `fq2` | `input` | Join as `fq1;fq2` and set `input_type` to `fastq`. |
| `LibraryName` | `library_id` | Optional; defaults to `sample_id`. |
| `refGenome`, `refPath` | _(removed)_ | Reference genome is now specified in `config.yaml` under `reference.name` and `reference.source`. |
| `SampleType` | `exclude` (in metadata sheet) | `exclude: true` in the metadata sheet replaces `SampleType: exclude`. |
| `lat`, `long` | `lat`, `long` (in metadata sheet) | Moved from the sample sheet to the optional metadata sheet. |

See the [Changelog](changelog.md) for the full list of v1-to-v2 breaking changes.

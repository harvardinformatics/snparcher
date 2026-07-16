# How to create a sample sheet

The sample sheet is a CSV file describing your samples and their data locations.
Each row represents one sequencing run for one sample.

For the full column specification, see the [sample sheet schema reference](../reference/sample-sheet-schema.md).

## Minimal example

Every sample sheet needs three columns: `sample_id`, `input_type`, and `input`.
Here is a minimal sheet for two samples from the NCBI Sequence Read Archive:

```csv
sample_id,input_type,input
sample_A,srr,SRR12345678
sample_B,srr,SRR12345679
```

snpArcher will download the reads automatically.

## From SRA accessions

To create a sample sheet from an existing NCBI BioProject:

1. Go to the BioProject page on the SRA.
2. In the **Project Data** table, click the link in the **Number of Links** column on the **SRA Experiments** row.
3. Click **Send results to Run Selector** near the top of the search results.
4. Select or deselect samples using the checkboxes.
5. Click the **Metadata** button in the **Download** column to download `SraRunTable.txt`.
6. Open the file and create a new CSV with three columns:

```csv
sample_id,input_type,input
SAMN12345678,srr,SRR12345678
SAMN12345679,srr,SRR12345679
SAMN12345680,srr,SRR12345680
```

Use the `BioSample` column from the SRA table as `sample_id` and the `Run` column as `input`.
Set `input_type` to `srr` for every row.

!!! tip
    The `sample_id` value will appear as the sample name in your final VCF.
    Choose something meaningful and consistent with your metadata.

## From local FASTQ files

For local paired-end FASTQ files, set `input_type` to `fastq` and join the forward and reverse read paths with a semicolon:

```csv
sample_id,input_type,input
bird_1,fastq,/data/reads/bird_1_R1.fq.gz;/data/reads/bird_1_R2.fq.gz
bird_2,fastq,/data/reads/bird_2_R1.fq.gz;/data/reads/bird_2_R2.fq.gz
bird_3,fastq,/data/reads/bird_3_R1.fq.gz;/data/reads/bird_3_R2.fq.gz
```

!!! warning
    Use absolute paths for FASTQ files.
    Relative paths are resolved from the project directory (the `--directory` argument to Snakemake), which can lead to confusing errors.

## From BAM files

To start from aligned BAM files (skipping trimming and alignment):

```csv
sample_id,input_type,input
sample_A,bam,/data/bams/sample_A.bam
sample_B,bam,/data/bams/sample_B.bam
```

Each sample must have exactly one BAM row.

## From gVCF files

To start from per-sample gVCF files (skipping everything up to joint genotyping):

```csv
sample_id,input_type,input
sample_A,gvcf,/data/gvcfs/sample_A.g.vcf.gz
sample_B,gvcf,/data/gvcfs/sample_B.g.vcf.gz
```

Each sample must have exactly one gVCF row.

!!! note
    gVCF input is only supported with the `gatk` and `sentieon` variant callers.
    If you set `variant_calling.tool` to `bcftools`, `deepvariant`, or `parabricks`, gVCF rows will cause a validation error.

## Mixing input types

You can mix input types in a single sample sheet.
snpArcher will start each sample at the appropriate pipeline stage:

```csv
sample_id,input_type,input,library_id
sample_A,fastq,/data/sample_A_R1.fq.gz;/data/sample_A_R2.fq.gz,lib_A
sample_B,srr,SRR12345678,
sample_C,bam,/data/sample_C.bam,
```

## Multiple runs per sample

When a sample has been sequenced across multiple lanes or runs, create one row per run with the same `sample_id`.
snpArcher will align each run independently, then merge the BAMs before variant calling.

The `library_id` column controls duplicate-marking behavior during the merge:

- Rows with the **same** `library_id` are treated as technical replicates of one library (e.g., the same library sequenced on multiple lanes). Duplicates are marked across the merged set.
- Rows with **different** `library_id` values are treated as independent libraries. Duplicates are marked within each library separately.

```csv
sample_id,input_type,input,library_id
sample_A,fastq,/data/sample_A_lane1_R1.fq.gz;/data/sample_A_lane1_R2.fq.gz,lib1
sample_A,fastq,/data/sample_A_lane2_R1.fq.gz;/data/sample_A_lane2_R2.fq.gz,lib1
sample_A,fastq,/data/sample_A_newlib_R1.fq.gz;/data/sample_A_newlib_R2.fq.gz,lib2
sample_B,fastq,/data/sample_B_R1.fq.gz;/data/sample_B_R2.fq.gz,lib3
```

In this example, `sample_A` has three sequencing runs: two from the same library (`lib1`) and one from a different library (`lib2`).

!!! note
    If you omit `library_id`, it defaults to the `sample_id`.
    This means all runs for a sample are treated as coming from one library.

## Disabling duplicate marking

If duplicate marking is inappropriate for your protocol (e.g., amplicon sequencing), set `mark_duplicates` to `false`:

```csv
sample_id,input_type,input,mark_duplicates
sample_A,fastq,/data/sample_A_R1.fq.gz;/data/sample_A_R2.fq.gz,false
sample_B,fastq,/data/sample_B_R1.fq.gz;/data/sample_B_R2.fq.gz,false
```

You can also set the global default in `config.yaml` under `reads.mark_duplicates` and omit the column from the sample sheet.

## Column reference

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Unique sample identifier. Alphanumeric, periods, underscores, and hyphens only. |
| `input_type` | Yes | One of `srr`, `fastq`, `bam`, `gvcf`. |
| `input` | Yes | SRA accession, semicolon-separated FASTQ paths, BAM path, or gVCF path. |
| `library_id` | No | Library identifier. Defaults to `sample_id`. |
| `mark_duplicates` | No | Boolean (`true`/`false`). Defaults to `reads.mark_duplicates` in config. |

For the full specification and validation rules, see the [sample sheet schema reference](../reference/sample-sheet-schema.md).

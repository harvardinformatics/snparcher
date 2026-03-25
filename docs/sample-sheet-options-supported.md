# snpArcher Sample Sheet Options Parsed in Code

This document summarizes sample sheet behavior currently implemented in:

- `workflow/schemas/samples.schema.yaml`
- `workflow/rules/common.smk`
- `workflow/rules/fastq.smk`

## 1) Primary sample sheet (`samples` config key)

Expected format: CSV with these columns.

| Column | Required | Type | Allowed values / constraints | Default behavior |
|---|---:|---|---|---|
| `sample_id` | yes | string | `^[A-Za-z0-9._-]+$` | Used as sample key across workflow |
| `input_type` | yes | string | `srr`, `fastq`, `bam`, `gvcf` | Must be consistent within each `sample_id` |
| `input` | yes | string | Non-empty. For `fastq`, must match `^.+;.+$` | Semantics depend on `input_type` |
| `library_id` | no | string | `^[A-Za-z0-9._-]+$` | If missing/blank, set to `sample_id` |
| `mark_duplicates` | no | boolean-like | Parsed as bool per rules below | Defaults to `reads.mark_duplicates` |

## 2) `input_type` semantics

| `input_type` | How `input` is interpreted | Row count support per sample |
|---|---|---|
| `srr` | SRA accession ID | Multiple rows allowed |
| `fastq` | `R1;R2` paths (semicolon-separated pair) | Multiple rows allowed |
| `bam` | BAM path | Exactly one row per sample |
| `gvcf` | gVCF path | Exactly one row per sample |

Additional runtime compatibility:

- `gvcf` samples are rejected when `variant_calling.tool` is `bcftools`, `deepvariant`, or `parabricks`.
- `gvcf` samples are supported with `gatk` and `sentieon`.

## 3) `mark_duplicates` parsing behavior

Accepted true values: `true`, `t`, `yes`, `y`, `1` (case-insensitive), numeric `1`, boolean `True`.

Accepted false values: `false`, `f`, `no`, `n`, `0` (case-insensitive), numeric `0`, boolean `False`.

Missing (`NA`) uses the global default from `reads.mark_duplicates`.

Invalid values raise an error with row number context.

Within each `sample_id`, `mark_duplicates` must be consistent across rows.

## 4) Derived/internal fields

The workflow derives these internally (not user-provided columns):

- `input_unit`: per `(sample_id, library_id)` ordinal `u1`, `u2`, ...

## 5) Optional secondary metadata sheet (`sample_metadata` config key)

If `sample_metadata` is provided, it is validated against `workflow/schemas/sample_metadata.schema.yaml` and parsed for module behavior.

Parsed columns used by code include:

- `sample_id` (must exist in primary sample sheet)
- `exclude` (boolean-like parser)
- `outgroup` (boolean-like parser)
- `lat`, `long` (used by QC map logic if present and non-null)


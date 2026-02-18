# v1 to v2 Migration Guide

This guide maps legacy snpArcher v1-style inputs to the v2 contracts used by the current workflow.

## Sample sheet migration

v2 uses a minimal manifest:
- required: `sample_id`, `input_type`, `input`
- optional: `library_id`, `mark_duplicates`

### Column mapping

| v1 column | v2 column | Notes |
| --- | --- | --- |
| `BioSample` | `sample_id` | Must be alphanumeric plus `._-` only. |
| `Run` (SRR/ERR/DRR) | `input` | Use with `input_type: srr`. |
| `fq1`,`fq2` | `input` | Join as `fq1;fq2` and set `input_type: fastq`. |
| `LibraryName` | `library_id` | Optional; defaults to `sample_id`. |
| `mark_duplicates` (global or per-row) | `mark_duplicates` | Optional; defaults to `true`. |
| `refGenome`,`refPath` | config `reference` | v2 keeps reference in config, not sample sheet. |

### Example: v1 fastq row -> v2 row

v1:

```csv
BioSample,LibraryName,Run,fq1,fq2
bird_1,bird_1_lib,1,/data/bird_1_R1.fastq.gz,/data/bird_1_R2.fastq.gz
```

v2:

```csv
sample_id,input_type,input,library_id,mark_duplicates
bird_1,fastq,/data/bird_1_R1.fastq.gz;/data/bird_1_R2.fastq.gz,bird_1_lib,true
```

### Example: v1 SRA row -> v2 row

v1:

```csv
BioSample,Run
bird_2,SRR12345678
```

v2:

```csv
sample_id,input_type,input
bird_2,srr,SRR12345678
```

## Config migration

v2 uses nested config keys. Core keys are:
- `samples`
- `reference.name`
- `reference.source`
- `variant_calling.*`
- `intervals.*`
- `callable_sites.*`
- `modules.*`

### Common legacy -> v2 config mapping

| v1 key | v2 key | Notes |
| --- | --- | --- |
| `refGenome` | `reference.name` | Required in v2. |
| `refPath` | `reference.source` | Can be local path, URL, or accession. |
| `sentieon` | `variant_calling.tool` | Set to `sentieon` (else `gatk`). |
| `sentieon_lic` | `variant_calling.sentieon.license` | String path or license value. |
| `minNmer` | `intervals.min_nmer` | Integer. |
| `cov_filter` | `callable_sites.coverage.enabled` | Boolean. |
| `mappability_k` | `callable_sites.mappability.kmer` | Integer. |
| `mappability_min` | `callable_sites.mappability.min_score` | Number in [0,1]. |
| `mappability_merge` | `callable_sites.mappability.merge_distance` | Integer. |

## Notes

- v2 no longer expects per-sample reference columns in the manifest.
- Optional metadata columns from v1 (for example `lat`, `long`, `SampleType`) are not part of the v2 core manifest contract.
- If a sample appears in multiple rows, `library_id` values must be unique within that sample.

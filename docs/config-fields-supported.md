# snpArcher Config Fields Supported in Code

This document is based on the current `main` branch code and schema definitions in:

- `workflow/schemas/config.schema.yaml`
- `workflow/rules/common.smk`
- `workflow/Snakefile`
- `workflow/modules/*/Snakefile`

## 1) Canonical v2 config fields (main workflow)

These are the supported v2 config keys used by the main workflow. Unknown keys are warned and ignored.

| Field | Type | Required | Default | Allowed values / constraints |
|---|---|---:|---|---|
| `samples` | string | yes | `config/samples.csv` | Path to samples CSV |
| `sample_metadata` | string | no | `""` | Optional metadata CSV path |
| `reads.mark_duplicates` | boolean | no | `true` | Global default for sample-level `mark_duplicates` |
| `reference.name` | string | yes | none | Reference name label for outputs |
| `reference.source` | string | yes | none | Path, URL, or NCBI accession |
| `variant_calling.expected_coverage` | string | no | `low` | `low`, `high` |
| `variant_calling.tool` | string | no | `gatk` | `gatk`, `sentieon`, `bcftools`, `deepvariant`, `parabricks` |
| `variant_calling.ploidy` | integer | no | `2` | `>= 1` |
| `variant_calling.gatk.het_prior` | number | no | `0.005` | `0..1` |
| `variant_calling.gatk.concat_batch_size` | integer | no | `250` | `>= 2` |
| `variant_calling.gatk.concat_max_rounds` | integer | no | `20` | `>= 1` |
| `variant_calling.sentieon.license` | string | no | `""` | Sentieon license value/path |
| `variant_calling.bcftools.min_mapq` | integer | no | `20` | `>= 0` |
| `variant_calling.bcftools.min_baseq` | integer | no | `20` | `>= 0` |
| `variant_calling.bcftools.max_depth` | integer | no | `250` | `>= 1` |
| `variant_calling.deepvariant.model_type` | string | no | `WGS` | `WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA` |
| `variant_calling.deepvariant.num_shards` | integer | no | `8` | `>= 1` |
| `variant_calling.parabricks.container_image` | string | no | `""` | Required at runtime if tool is `parabricks` |
| `variant_calling.parabricks.num_gpus` | integer | no | `1` | `>= 1` |
| `variant_calling.parabricks.num_cpu_threads` | integer | no | `16` | `>= 1` |
| `variant_calling.parabricks.extra_args` | string | no | `""` | Extra CLI args |
| `intervals.enabled` | boolean | no | `true` | Enables interval generation for GATK paths |
| `intervals.min_nmer` | integer | no | `500` | `>= 1` |
| `intervals.num_gvcf_intervals` | integer | no | `50` | `>= 1` |
| `intervals.db_scatter_factor` | number | no | `0.15` | `>= 0` |
| `callable_sites.generate_bed_file` | boolean | no | `true` | Controls final callable BED target |
| `callable_sites.coverage.enabled` | boolean | no | `true` | Requires BAM-backed samples if true |
| `callable_sites.coverage.fraction` | number | no | `1.0` | `0..1` |
| `callable_sites.coverage.min_coverage` | number or string | no | `auto` | `>= 0` or `auto` |
| `callable_sites.coverage.max_coverage` | number or string | no | `auto` | `>= 0` or `auto` |
| `callable_sites.coverage.merge_distance` | integer | no | `100` | `>= 0` |
| `callable_sites.mappability.enabled` | boolean | no | `true` | Enables mappability branch |
| `callable_sites.mappability.kmer` | integer | no | `150` | `>= 1` |
| `callable_sites.mappability.min_score` | number | no | `1` | `0..1` |
| `callable_sites.mappability.merge_distance` | integer | no | `100` | `>= 0` |
| `modules.qc.enabled` | boolean | no | `false` | Enable QC module |
| `modules.qc.clusters` | integer | no | `3` | `>= 1` |
| `modules.qc.min_depth` | number | no | `2` | `>= 0` |
| `modules.qc.google_api_key` | string | no | `""` | Optional for map panel |
| `modules.qc.exclude_scaffolds` | string | no | `""` | Comma-separated scaffold list |
| `modules.postprocess.enabled` | boolean | no | `false` | Enable postprocess module |
| `modules.postprocess.filtering.contig_size` | integer | no | `10000` | `>= 0` |
| `modules.postprocess.filtering.maf` | number | no | `0.01` | `0..1` |
| `modules.postprocess.filtering.missingness` | number | no | `0.75` | `0..1` |
| `modules.postprocess.filtering.exclude_scaffolds` | string | no | `mtDNA,Y` | Comma-separated scaffold list |

## 2) Aliases and compatibility parsing

These keys are explicitly parsed by code (not canonical schema fields):

| Key | Status in code | Behavior |
|---|---|---|
| `reference.path` | deprecated alias | If `reference.source` missing, copied to `reference.source`; then removed |
| `variant_calling.gatk.ploidy` | deprecated alias | If `variant_calling.ploidy` missing, copied there; then removed |

Legacy v1 markers (for example `final_prefix`, `trackhub_email`, `refGenome`, etc.) are detected and cause a hard error if config is not v2-shaped.

Removed module blocks `modules.mk` and `modules.trackhub` are warned about and ignored when present in an otherwise valid v2 config.

## 3) Standalone module flat config keys

When running module Snakefiles directly (outside the main workflow module-import mapping), these flat keys are parsed:

### QC module (`workflow/modules/qc/Snakefile`)

`samples`, `sample_metadata`, `vcf`, `fai`, `qc_report`, `clusters`, `min_depth`, `google_api_key`, `exclude_scaffolds`

### Postprocess module (`workflow/modules/postprocess/Snakefile`)

`samples`, `sample_metadata`, `vcf`, `ref_fai`, `callable_sites_bed`, `contig_size`, `maf`, `missingness`, `exclude_scaffolds`

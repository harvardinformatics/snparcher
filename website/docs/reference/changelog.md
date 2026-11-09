# Changelog

## v2.3 (unreleased)

### Breaking changes

- **Callable sites coverage mask redesigned.**
  `callable_sites.coverage.min_coverage` and `max_coverage` are replaced by `max_depth`, `max_depth_k` and `min_sample_depth`, and the default `fraction` changes from 1.0 to 0.25.
  A site is now coverage-callable when its depth summed across samples is not an outlier above the genome-wide median (or, with a numeric `max_depth`, not above that multiple of the median), and at least `fraction` of samples have `min_sample_depth` reads there.
  The old per-sample thresholds kept only a small fraction of the genome in realistic cohorts.
  Configs that still set the old keys skip callable sites, and postprocess with it, with a warning; the rest of the workflow runs.
- **`callable_loci.zarr` now holds clam population counts** (samples with data per site) instead of per-sample masks, with all samples as one population.
  Use it with `clam stat` for a single population; for statistics by population, build a store with your populations from `depths.zarr` (see [Filtering](../explanation/filtering.md#diversity-and-divergence-with-clam)).

### Fixes

- `coverage.bed` no longer misplaces callable sites. clam's per-sample masks were decoded in the wrong bit order whenever the number of samples was not a multiple of 8; the masks are no longer used.
- `fraction` now converts to an exact sample count: 0.55 of 200 samples requires 110, not 111.
- GenMap no longer overwrites or deletes `results/reference/{name}.fa`.
- Samples marked `exclude` in the sample metadata are left out of the coverage mask.
- Sample IDs containing a dot keep their full name in `depths.zarr`. clam used to cut each name at the first dot, so `ind.1` was stored as `ind`, and `ind.1` with `ind.2` failed `clam_collect`.

## v2.0

v2.0 is a major rewrite of snpArcher with breaking changes to the sample sheet schema, config structure, and module interfaces.

### Breaking changes

- **Sample sheet schema redesigned.**
  The v1 columns (`BioSample`, `Run`, `LibraryName`, `fq1`, `fq2`, `refGenome`, `refPath`, `SampleType`, `lat`, `long`) are removed.
  v2 uses `sample_id`, `input_type`, `input`, with optional `library_id` and `mark_duplicates`.
- **Reference genome moved to config.**
  `refGenome` and `refPath` are no longer sample sheet columns.
  Use `reference.name` and `reference.source` in `config.yaml`.
- **Sample metadata split into a separate file.**
  Columns like `exclude`, `lat`, and `long` are now in an optional metadata CSV specified by the `sample_metadata` config key.
- **Config structure is now hierarchical.**
  Flat v1 keys (e.g., `minNmer`, `cov_filter`, `mappability_k`) are replaced by nested keys (e.g., `intervals.min_nmer`, `callable_sites.coverage.enabled`, `callable_sites.mappability.kmer`).
- **Module config is nested under `modules`.**
  QC and postprocess options are now under `modules.qc.*`, `modules.postprocess.*`, etc.
- **v1 config keys trigger a hard error.**
  If legacy v1 keys (e.g., `final_prefix`, `refGenome`) are detected, the workflow exits with an error message.

### Sample sheet column mapping (v1 to v2)

| v1 column | v2 column | Notes |
|-----------|-----------|-------|
| `BioSample` | `sample_id` | Must match `^[A-Za-z0-9._-]+$`. |
| `Run` | `input` | Set `input_type: srr`. |
| `fq1`, `fq2` | `input` | Join as `fq1;fq2`, set `input_type: fastq`. |
| `LibraryName` | `library_id` | Optional; defaults to `sample_id`. |
| `refGenome` | `reference.name` (config) | Moved to config. |
| `refPath` | `reference.source` (config) | Moved to config. |
| `SampleType` | `exclude` (metadata sheet) | `exclude: true` in metadata replaces `SampleType: exclude`. |
| `lat` | `lat` (metadata sheet) | Moved to optional metadata CSV. |
| `long` | `long` (metadata sheet) | Moved to optional metadata CSV. |
| `mark_duplicates` | `mark_duplicates` | Now a per-row sample sheet column or global config default. |

### Config key mapping (v1 to v2)

| v1 key | v2 key |
|--------|--------|
| `refGenome` | `reference.name` |
| `refPath` | `reference.source` |
| `sentieon` | `variant_calling.tool: sentieon` |
| `sentieon_lic` | `variant_calling.sentieon.license` |
| `minNmer` | `intervals.min_nmer` |
| `cov_filter` | `callable_sites.coverage.enabled` |
| `mappability_k` | `callable_sites.mappability.kmer` |
| `mappability_min` | `callable_sites.mappability.min_score` |
| `mappability_merge` | `callable_sites.mappability.merge_distance` |

### New features

- **New input types:** `bam` and `gvcf` input types allow entering the pipeline at later stages.
- **Callable sites overhaul:** Coverage analysis uses mosdepth D4 and Zarr-based processing. Auto-thresholding computes min/max coverage from the depth distribution.
- **Interval-based parallelization for GATK:** Hierarchical gVCF/VCF concatenation with configurable batch sizes.
- **Parabricks support:** GPU-accelerated variant calling via NVIDIA Parabricks.
- **DeepVariant support:** Deep learning variant caller with GLnexus joint genotyping.
- **bcftools pathway:** Region-parallelized bcftools mpileup/call as an alternative to GATK.

### Deprecated aliases

These keys are accepted with a warning but should be updated:

| Deprecated key | Replacement |
|----------------|-------------|
| `reference.path` | `reference.source` |
| `variant_calling.gatk.ploidy` | `variant_calling.ploidy` |

## v1.x

The original snpArcher release described in:

> Mirchandani et al. (2024). "A fast, reproducible, high-throughput variant calling workflow for population genomics." *Molecular Biology and Evolution*, 41(1), msad270. [doi:10.1093/molbev/msad270](https://doi.org/10.1093/molbev/msad270)

v1 supported GATK and Sentieon variant calling with a flat config structure and a single sample sheet combining sample data, reference paths, and metadata.
See the [MBE paper](https://doi.org/10.1093/molbev/msad270) for benchmarking results and design rationale.

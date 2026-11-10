# Changelog

## v2.3.0

### Upgrading from 2.2

- **Callable sites coverage mask redesigned.**
  `callable_sites.coverage.min_coverage` and `max_coverage` are replaced by `max_depth`, `max_depth_k` and `min_sample_depth`, and the default `fraction` changes from 1.0 to 0.25.
  A site is now coverage-callable when its depth summed across samples is not an outlier above the genome-wide median (or, with a numeric `max_depth`, not above that multiple of the median), and at least `fraction` of samples have `min_sample_depth` reads there.
  The old per-sample thresholds kept only a small fraction of the genome in realistic cohorts.
  With coverage enabled, configs that still set the old keys skip callable sites, and postprocess with it, with a warning; the rest of the workflow runs.
- **DeepVariant needs Apptainer.** Add `--use-apptainer`, or `use-apptainer: True` in your workflow profile. Without it, the run now stops while the job graph is built instead of failing after mapping.
- **Stricter checks at startup.** These used to fail later, or silently do the wrong thing:
  - `ploidy` above 2 with bcftools or Parabricks, or other than 2 with DeepVariant;
  - an empty Sentieon licence;
  - a FASTQ row that isn't exactly one `R1;R2` pair, or a FASTQ file or SRA accession listed twice for one sample;
  - a reference source that isn't an existing file, a URL, or an NCBI accession with its version (`GCF_000002315.6`);
  - a metadata `exclude` value other than true/false, yes/no or 1/0;
  - a postprocess `maf` above 0.5.

### Using v2.3 in an existing project

v2.3 changes many rules, and Snakemake's default rerun triggers rerun any job whose rule changed, including jobs whose temporary outputs were deleted long ago, and everything downstream of them.
In a project started with 2.2, this can mean redoing mapping and variant calling, and re-downloading SRA samples.
Dry-run first (`-n`) to see what would run, then choose:

- **Recompute everything the new code touches:** run with the default triggers.
- **Keep the existing results:** run with `--rerun-triggers mtime`.
- **Recompute only some outputs:** delete them, then run with `--rerun-triggers mtime`. The outputs with v2.3 fixes worth rebuilding are `results/callable_sites/` (new coverage mask), `results/postprocess/`, `results/qc_metrics/` and, if you use the QC module, `results/qc/`.

Long-contig GATK projects started with 2.2 should delete `results/intervals/db/` before their first v2.3 run.
Those lists can genotype part of a chromosome twice, and long-contig mode kept the duplicate records, so also check VCFs from those runs for duplicate records at the same position.

### Changed outputs

- **GATK runs with hard filtering keep `filtered.vcf.gz` as the call set**, and no longer keep `raw.vcf.gz`. `filtered.vcf.gz` has every raw record, with FILTER set.
- **`callable_loci.zarr` holds clam population counts** (samples with data per site) instead of per-sample masks, with all samples as one population.
  Use it with `clam stat` for a single population; for statistics by population, build a store with your populations from `depths.zarr` (see [Filtering](../explanation/filtering.md#diversity-and-divergence-with-clam)).
- **bcftools and DeepVariant name samples by `sample_id`**, not the read-group `SM` of the sample's BAM. With bcftools, BAMs that shared an `SM` were merged into one sample. Sentieon still uses `SM`, so with Sentieon an external BAM's `SM` must equal its `sample_id`.
- **QC:** site filters now see only the samples kept by `min_depth`, and the MAF filter is two-sided, so QC's SNP set, PCA and ADMIXTURE change. Reruns now give the same PCA and ADMIXTURE. Intermediate files use `qcctgN` contig names.
- **Postprocess:** indels whose first base isn't callable are no longer kept, and `callable_sites_filtered.bed` no longer contains excluded scaffolds.
- **`qc_report.tsv`:** `total_reads` and the mapping and duplicate rates count primary records only.
- **Sample metadata:** the `outgroup` column is no longer read. It was never used; files that have it still work.

### Fixes

**Variant calling**
- GATK with `intervals.enabled: false` runs; GenomicsDBImport used to reject its interval argument after every HaplotypeCaller job had finished.
- bcftools, DeepVariant, Sentieon and GATK runs without intervals deliver an index with the call set; long-contig GATK runs no longer re-genotype every interval on a later `all`.
- bcftools handles contig names with shell characters (such as NCBI `gi|…|` names) and assemblies with more contigs than the open-file limit.
- `snakemake gvcfs` with bcftools explains that bcftools writes no per-sample gVCFs.
- Long-contig GATK gathers (`picard SortVcf`) no longer run out of memory above about 30 samples.
- DeepVariant deletes its intermediate files when each sample finishes.
- `deepvariant.num_shards` and `parabricks.num_cpu_threads` take effect; the default profile overrode them.
- `variant_calling.tool: sentieon` works again: it failed to parse under Snakemake 9 and could not find its conda environment.

**Callable sites**
- `coverage.bed` no longer misplaces callable sites. clam's per-sample masks were decoded in the wrong bit order whenever the number of samples was not a multiple of 8; the masks are no longer used.
- `fraction` and a numeric `max_depth` convert to exact thresholds: 0.55 of 200 samples requires 110, not 111.
- Samples marked `exclude` in the sample metadata are left out of the coverage mask.
- Sample IDs containing a dot keep their full name in `depths.zarr`. clam used to cut each name at the first dot, so `ind.1` was stored as `ind`, and `ind.1` with `ind.2` failed `clam_collect`.
- GenMap no longer overwrites or deletes `results/reference/{name}.fa`.

**Postprocess**
- Editing `exclude` in the sample metadata reruns postprocess with the new samples.
- `exclude_scaffolds` values with spaces work, and the excluded scaffolds are removed from the filtered BED as well as the VCF.
- Records left with only the `*` allele are dropped; records with a real allele beside `*` are kept.
- A run with no SNPs left, or no callable region left, finishes or stops with an explanation instead of a bcftools error.
- The standalone module runs outside a project directory, and `run-postprocess.sh --long-contig-mode` writes CSI indexes.
- A stray `modules.postprocess.filtering.vcf` key no longer replaces the input VCF.

**QC**
- Contigs named `chrX`, `X`, `Y`, `MT` or `chrM` no longer crash PLINK, and cohorts of 3–10 samples no longer crash the PCA or the dashboard.
- The dashboard no longer counts `PASS` records as filtered out, handles more than 10 clusters, and renders without writing into the install directory.

**Samples, config and downloads**
- Building the job graph is no longer quadratic in the number of samples: a 200-sample dry run went from 272 s to 44 s.
- A `library_id` or `sample_id` of `NA`, `null` or `None` is kept as written, instead of being read as missing.
- Metadata `exclude` accepts yes/no and 1/0, as documented.
- A v2 config with an `intervals:` block is no longer mistaken for a v1 config.
- SRA re-downloads no longer fail on the previous attempt's files.
- An unpaired external BAM no longer crashes QC metrics.

**Documentation**
- `--workflow-profile` is required; several pages said it defaulted to `workflow-profiles/default`, and their example commands failed without it.
- The postprocess missingness default is described correctly everywhere: 0.75 keeps sites where at least 25% of samples are called.

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

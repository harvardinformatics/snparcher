# Config schema

Values are validated at startup against `workflow/schemas/config.schema.yaml`.
Unknown keys produce a warning and are ignored.

## Required settings

| Key | Type | Description |
|-----|------|-------------|
| `samples` | string | Path to the samples CSV file. |
| `reference.name` | string | Reference genome name, used in output filenames. |
| `reference.source` | string | Path, URL, or NCBI accession for the reference genome FASTA. |

An optional key is also accepted at the top level:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `sample_metadata` | string | _(none)_ | Path to an optional sample metadata CSV. Must contain a `sample_id` column matching entries in the samples CSV. Used by modules for `exclude`, `lat`, `long`, etc. See [Sample sheet schema](sample-sheet-schema.md#optional-metadata-sheet). |

!!! note "Deprecated alias: `reference.path`"
    If `reference.source` is not set but `reference.path` is present, the value is copied to `reference.source` automatically.
    Use `reference.source` in new configurations.

## `reference`

| Key | Type | Required | Default | Description |
|-----|------|:--------:|---------|-------------|
| `reference.name` | string | yes | _(none)_ | Reference genome name label for output files. |
| `reference.source` | string | yes | _(none)_ | Local path, URL, or NCBI RefSeq/GenBank accession for the reference FASTA. Compressed (`.gz`) files are accepted. |

## `reads`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `reads.mark_duplicates` | boolean | `true` | Global default for duplicate marking. Can be overridden per row in the sample sheet. |

## `variant_calling`

| Key | Type | Default | Allowed values | Description |
|-----|------|---------|----------------|-------------|
| `variant_calling.expected_coverage` | string | `"low"` | `low`, `high` | Expected sequencing coverage regime. |
| `variant_calling.tool` | string | `"gatk"` | `gatk`, `sentieon`, `bcftools`, `deepvariant`, `parabricks` | Variant caller to use. |
| `variant_calling.ploidy` | integer | `2` | >= 1 | Ploidy of the organism. |

!!! note "Deprecated alias: `variant_calling.gatk.ploidy`"
    If `variant_calling.ploidy` is not set but `variant_calling.gatk.ploidy` is present, the value is copied to `variant_calling.ploidy`.
    Use `variant_calling.ploidy` in new configurations.

### `variant_calling.gatk`

Applies when `variant_calling.tool` is `gatk`.

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `variant_calling.gatk.het_prior` | number | `0.005` | 0 to 1 | Prior probability of heterozygosity, passed to `HaplotypeCaller --heterozygosity` and `GenotypeGVCFs --heterozygosity`. |
| `variant_calling.gatk.concat_batch_size` | integer | `250` | >= 2 | Number of interval gVCFs/VCFs concatenated per batch during hierarchical merging. |
| `variant_calling.gatk.concat_max_rounds` | integer | `20` | >= 1 | Maximum number of hierarchical concatenation rounds. |

See [Variant calling](../explanation/variant-calling.md) for guidance on setting this value.

### `variant_calling.sentieon`

Applies when `variant_calling.tool` is `sentieon`.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `variant_calling.sentieon.license` | string | `""` | Sentieon license server address or license file path. Required at runtime when `tool` is `sentieon`. |

### `variant_calling.bcftools`

Applies when `variant_calling.tool` is `bcftools`.

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `variant_calling.bcftools.min_mapq` | integer | `20` | >= 0 | Minimum mapping quality for reads included in pileup. |
| `variant_calling.bcftools.min_baseq` | integer | `20` | >= 0 | Minimum base quality for bases included in pileup. |
| `variant_calling.bcftools.max_depth` | integer | `250` | >= 1 | Maximum per-sample read depth in pileup. |

### `variant_calling.deepvariant`

Applies when `variant_calling.tool` is `deepvariant`.

| Key | Type | Default | Allowed values | Description |
|-----|------|---------|----------------|-------------|
| `variant_calling.deepvariant.model_type` | string | `"WGS"` | `WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA` | DeepVariant model type matching the sequencing technology. |
| `variant_calling.deepvariant.num_shards` | integer | `8` | >= 1 | Number of shards for parallelizing DeepVariant inference. |

### `variant_calling.parabricks`

Applies when `variant_calling.tool` is `parabricks`.

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `variant_calling.parabricks.container_image` | string | `""` | | Docker/Singularity image for Parabricks. Required at runtime when `tool` is `parabricks`. |
| `variant_calling.parabricks.num_gpus` | integer | `1` | >= 1 | Number of GPUs to allocate per Parabricks job. |
| `variant_calling.parabricks.num_cpu_threads` | integer | `16` | >= 1 | Number of CPU threads per Parabricks job. |
| `variant_calling.parabricks.extra_args` | string | `""` | | Additional CLI arguments passed to the Parabricks command. |

## `intervals`

Interval-based parallelization settings.
See [Parallelization](../explanation/parallelization.md) for background.

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `intervals.enabled` | boolean | `true` | | Enable interval-based parallelization for GATK. When `false`, GATK runs whole-genome HaplotypeCaller and single-pass GenomicsDBImport. |
| `intervals.min_nmer` | integer | `500` | >= 1 | Minimum length of N-mer run (assembly gap) used to split the genome into intervals. |
| `intervals.num_gvcf_intervals` | integer | `50` | >= 1 | Target number of intervals for per-sample gVCF calling. |
| `intervals.db_scatter_factor` | number | `0.15` | >= 0 | Scaling factor for GenomicsDB import parallelism. The number of database intervals is `db_scatter_factor * num_samples * num_gvcf_intervals` (rounded, minimum 1). Scales with both cohort size and interval count. |

## `callable_sites`

Callable-sites BED generation settings.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `callable_sites.generate_bed_file` | boolean | `true` | When `true`, produce the final intersected callable sites BED file. When `false`, individual coverage and mappability outputs are still produced if their sub-sections are enabled, but the final intersection is skipped. |

### `callable_sites.coverage`

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `callable_sites.coverage.enabled` | boolean | `true` | | Enable coverage-based callable site detection. |
| `callable_sites.coverage.fraction` | number | `1.0` | 0 to 1 | Fraction of samples that must meet coverage thresholds at a site for it to be called "callable". |
| `callable_sites.coverage.min_coverage` | number or `"auto"` | `"auto"` | >= 0 or `"auto"` | Minimum per-sample depth threshold. `"auto"` computes the threshold from the depth distribution. |
| `callable_sites.coverage.max_coverage` | number or `"auto"` | `"auto"` | >= 0 or `"auto"` | Maximum per-sample depth threshold. `"auto"` computes the threshold from the depth distribution. |
| `callable_sites.coverage.merge_distance` | integer | `100` | >= 0 | Merge callable regions within this many bases of each other. |

### `callable_sites.mappability`

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `callable_sites.mappability.enabled` | boolean | `true` | | Enable mappability-based callable site detection using GenMap. |
| `callable_sites.mappability.kmer` | integer | `150` | >= 1 | K-mer length for GenMap mappability computation. Should match read length. |
| `callable_sites.mappability.min_score` | number | `1` | 0 to 1 | Minimum mappability score. `1` means only uniquely mappable regions are retained. |
| `callable_sites.mappability.merge_distance` | integer | `100` | >= 0 | Merge mappable regions within this many bases of each other. |

## `modules`

See the [Modules](modules.md) reference for full details.

### `modules.qc`

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `modules.qc.enabled` | boolean | `false` | | Enable the QC dashboard module. |
| `modules.qc.clusters` | integer | `3` | >= 1 | Number of clusters for PCA-based clustering. |
| `modules.qc.min_depth` | number | `2` | >= 0 | Samples with mean depth below this value are excluded from QC analyses. |
| `modules.qc.google_api_key` | string | `""` | | Google Maps API key for geographic map panels. Optional. |
| `modules.qc.exclude_scaffolds` | string | `""` | | Comma-separated list of scaffold/contig names to exclude from QC analyses. |

### `modules.postprocess`

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `modules.postprocess.enabled` | boolean | `false` | | Enable the postprocessing/filtering module. |
| `modules.postprocess.filtering.contig_size` | integer | `10000` | >= 0 | Exclude SNPs on contigs with length equal to or smaller than this value (bp). |
| `modules.postprocess.filtering.maf` | number | `0.01` | 0 to 1 | Minimum minor allele frequency. Sites below this threshold are excluded. |
| `modules.postprocess.filtering.missingness` | number | `0.75` | 0 to 1 | Minimum genotyping rate. Sites with a called-genotype fraction below this threshold are excluded. |
| `modules.postprocess.filtering.exclude_scaffolds` | string | `"mtDNA,Y"` | | Comma-separated list of scaffold/contig names to exclude from filtered output. |

## Minimal example

```yaml
samples: "config/samples.csv"

reference:
  name: "my_organism"
  source: "/data/reference/genome.fna.gz"
```

All other keys use their documented defaults.

## Full example

```yaml
samples: "config/samples.csv"
sample_metadata: "config/sample_metadata.csv"

reference:
  name: "my_organism"
  source: "/data/reference/genome.fna.gz"

reads:
  mark_duplicates: true

variant_calling:
  expected_coverage: "low"
  tool: "gatk"
  ploidy: 2
  gatk:
    het_prior: 0.005

intervals:
  enabled: true
  min_nmer: 500
  num_gvcf_intervals: 50
  db_scatter_factor: 0.15

callable_sites:
  generate_bed_file: true
  coverage:
    enabled: true
    fraction: 1.0
    min_coverage: auto
    max_coverage: auto
    merge_distance: 100
  mappability:
    enabled: true
    kmer: 150
    min_score: 1
    merge_distance: 100

modules:
  qc:
    enabled: true
    clusters: 3
    min_depth: 2
    google_api_key: ""
    exclude_scaffolds: ""
  postprocess:
    enabled: true
    filtering:
      contig_size: 10000
      maf: 0.01
      missingness: 0.75
      exclude_scaffolds: "mtDNA,Y"
```

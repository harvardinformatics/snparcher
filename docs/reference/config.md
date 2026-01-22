# Configuration Reference

snpArcher is configured via `config/config.yaml`. This page documents all available options.

## Required Options

### samples

Path to sample sheet CSV file.

```yaml
samples: "config/samples.csv"
```

### reference

Reference genome configuration. Must specify `name` and either `path` or `accession`.

```yaml
reference:
  name: "my_genome"               # Identifier used in output naming
  path: "/path/to/ref.fna.gz"     # Local FASTA file (bgzip-compressed)
```

Or using an NCBI accession:

```yaml
reference:
  name: "ecoli_k12"
  accession: "GCF_000005845.2"
```

#### Reference format requirements

When providing a local reference via `path`:

- **Bgzip-compressed FASTA** (`.fna.gz`, `.fa.gz`) is **required**
- The file is **symlinked** to results directory (no copying)
- Regular gzip is **not supported** and will cause indexing to fail

If your reference is uncompressed or gzip-compressed, convert it to bgzip format:

```bash
# From uncompressed FASTA
bgzip genome.fna               # Creates genome.fna.gz

# From regular gzip
gunzip genome.fna.gz           # Decompress
bgzip genome.fna               # Recompress with bgzip
```

Bgzip (block gzip) allows random access, which is required by `samtools faidx` and 
other indexing tools. The `bgzip` command is part of [htslib](http://www.htslib.org/).

When using an NCBI `accession`, snpArcher downloads and bgzip-compresses 
reference automatically.

## Variant Calling

### variant_calling.tool

Variant calling tool to use: `"gatk"` or `"sentieon"`.

```yaml
variant_calling:
  tool: "gatk"      # or "sentieon"
```

### variant_calling.gatk

GATK variant calling parameters.

```yaml
variant_calling:
  gatk:
    min_pruning: 1              # HaplotypeCaller --min-pruning (default: 1)
    min_dangling: 1             # HaplotypeCaller --min-dangling-branch-length (default: 1)
    het_prior: 0.005            # Prior probability of heterozygous site (default: 0.005)
    ploidy: 2                    # Ploidy for variant calling (default: 2)
```

**Coverage tuning:** Use low-coverage settings for <10x, high-coverage for >10x.

- **Low coverage (<10x):** `min_pruning: 1`, `min_dangling: 1`
- **High coverage (>10x):** `min_pruning: 2`, `min_dangling: 4`

### variant_calling.sentieon

Sentieon variant calling parameters.

```yaml
variant_calling:
  sentieon:
    license: "/path/to/sentieon/license"    # Path to Sentieon license file
```

## Intervals

### intervals.enabled

Use interval-based approach for variant calling. Recommended for most use cases.

```yaml
intervals:
  enabled: true  # default: true
```

### Interval Options

Only applicable when `intervals.enabled: true`.

```yaml
intervals:
  min_nmer: 500              # Minimum N-mer for splitting (default: 500, min: 50)
  num_gvcf_intervals: 50    # Maximum GVCF intervals to create (default: 50)
  db_scatter_factor: 0.15   # GenomicsDB scatter factor (default: 0.15)
```

## Read Processing

### reads.mark_duplicates

Mark optical duplicates before variant calling.

```yaml
reads:
  mark_duplicates: true  # default: true
```

### reads.sort

Sort reads by name before adapter trimming.

```yaml
reads:
  sort: false  # default: false
```

## Remote Reads

### remote_reads.enabled

Enable cloud execution with reads in remote storage.

```yaml
remote_reads:
  enabled: false  # default: false
```

### remote_reads.prefix

Remote storage prefix for reads.

```yaml
remote_reads:
  prefix: "gs://my-bucket/"  # Remote storage prefix
```

## Callable Sites

### callable_sites.coverage.enabled

Enable coverage-based filtering of callable sites.

```yaml
callable_sites:
  coverage:
    enabled: true  # default: true
```

### callable_sites.coverage.stdev

Filter regions outside mean ± N standard deviations (assumes Poisson distribution).

```yaml
callable_sites:
  coverage:
    stdev: 2  # Filter mean ± 2 standard deviations (default: 2)
```

### callable_sites.coverage.merge_distance

Merge passing coverage regions separated by this many bp into a single region.

```yaml
callable_sites:
  coverage:
    merge_distance: 100  # default: 100
```

### callable_sites.mappability.enabled

Enable mappability-based filtering.

```yaml
callable_sites:
  mappability:
    enabled: true  # default: true
```

### callable_sites.mappability.min_score

Minimum mappability score for callable sites.

```yaml
callable_sites:
  mappability:
    min_score: 1  # default: 1
```

### callable_sites.mappability.kmer

K-mer size for mappability calculation.

```yaml
callable_sites:
  mappability:
    kmer: 150  # default: 150
```

### callable_sites.mappability.merge_distance

Merge passing mappability regions separated by this many bp into a single region.

```yaml
callable_sites:
  mappability:
    merge_distance: 100  # default: 100
```

## Modules

### modules.qc.enabled

Enable quality control module.

```yaml
modules:
  qc:
    enabled: false  # default: false
```

### modules.qc.clusters

Number of clusters for PCA analysis.

```yaml
modules:
  qc:
    clusters: 3  # default: 3
```

### modules.qc.min_depth

Minimum average depth for QC inclusion.

```yaml
modules:
  qc:
    min_depth: 2  # default: 2
```

### modules.qc.google_api_key

Google Maps API key for geographic plots (optional).

```yaml
modules:
  qc:
    google_api_key: ""  # Optional
```

### modules.postprocess.enabled

Enable postprocessing module for clean VCF generation.

```yaml
modules:
  postprocess:
    enabled: false  # default: false
```

### modules.postprocess.filtering

VCF filtering parameters for clean VCF generation.

```yaml
modules:
  postprocess:
    filtering:
      contig_size: 10000           # Exclude SNPs on contigs smaller than this (default: 10000)
      maf: 0.01                    # Minimum minor allele frequency (default: 0.01)
      missingness: 0.75            # Maximum missingness threshold (default: 0.75)
      exclude_scaffolds: "mtDNA,Y"  # Comma-separated scaffolds to exclude (default: "mtDNA,Y")
```

### modules.trackhub.enabled

Generate UCSC Genome Browser trackhub with population genomic statistics.

**Important:** trackhub module requires `modules.postprocess.enabled: true`.

```yaml
modules:
  trackhub:
    enabled: false  # default: false
```

### modules.trackhub.email

Email address for trackhub (required by UCSC if enabled).

```yaml
modules:
  trackhub:
    email: "email@example.com"  # Required if trackhub enabled
```

## Migration Guide

If you have a v1 config file, here's how to migrate to v2:

| Old Config | New Config |
|-------------|-------------|
| `sentieon: True` | `variant_calling: {tool: "sentieon"}` |
| `sentieon_lic: "/path/to/license"` | `variant_calling: {sentieon: {license: "/path/to/license"}}` |
| `minP: 1` | `variant_calling: {gatk: {min_pruning: 1}}` |
| `minD: 1` | `variant_calling: {gatk: {min_dangling: 1}}` |
| `het_prior: 0.005` | `variant_calling: {gatk: {het_prior: 0.005}}` |
| `ploidy: 2` | `variant_calling: {gatk: {ploidy: 2}}` |
| `intervals: True` | `intervals: {enabled: true}` |
| `minNmer: 500` | `intervals: {min_nmer: 500}` |
| `num_gvcf_intervals: 50` | `intervals: {num_gvcf_intervals: 50}` |
| `db_scatter_factor: 0.15` | `intervals: {db_scatter_factor: 0.15}` |
| `cov_filter: True` | `callable_sites: {coverage: {enabled: true}}` |
| `cov_threshold_stdev: 2` | `callable_sites: {coverage: {stdev: 2}}` |
| `cov_merge: 100` | `callable_sites: {coverage: {merge_distance: 100}}` |
| `mappability_min: 1` | `callable_sites: {mappability: {min_score: 1}}` |
| `mappability_k: 150` | `callable_sites: {mappability: {kmer: 150}}` |
| `mappability_merge: 100` | `callable_sites: {mappability: {merge_distance: 100}}` |
| `sort_reads: False` | `reads: {sort: false}` |
| `mark_duplicates: True` | `reads: {mark_duplicates: true}` |
| `nClusters: 3` | `modules: {qc: {clusters: 3}}` |
| `min_depth: 2` | `modules: {qc: {min_depth: 2}}` |
| `GoogleAPIKey: ""` | `modules: {qc: {google_api_key: ""}}` |

**Removed options (simplified coverage filtering):**
- `cov_threshold_lower` - Use `callable_sites.coverage.stdev` instead
- `cov_threshold_upper` - Use `callable_sites.coverage.stdev` instead
- `cov_threshold_rel` - Use `callable_sites.coverage.stdev` instead

## Example Configuration

```yaml
# Required
samples: "config/samples.csv"
reference:
  name: "my_organism"
  path: "/data/reference/genome.fna.gz"  # Must be bgzip-compressed

# Variant calling
variant_calling:
  tool: "gatk"
  gatk:
    min_pruning: 1
    min_dangling: 1
    het_prior: 0.005
    ploidy: 2

# Intervals
intervals:
  enabled: true
  min_nmer: 500
  num_gvcf_intervals: 50
  db_scatter_factor: 0.15

# Read processing
reads:
  mark_duplicates: true
  sort: false

# Callable sites
callable_sites:
  coverage:
    enabled: true
    stdev: 2
    merge_distance: 100
  mappability:
    enabled: true
    min_score: 1
    kmer: 150
    merge_distance: 100

# Modules
modules:
  qc:
    enabled: false
  postprocess:
    enabled: false
    trackhub:
    enabled: false
```

# Configuration Reference

snpArcher is configured via `config/config.yaml`. This page documents all available options.

## Required Options

### samples

Path to the sample sheet CSV file.

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
- The file is **symlinked** to the results directory (no copying)
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

When using an NCBI `accession`, snpArcher downloads and bgzip-compresses the 
reference automatically.

## Variant Calling Options

### intervals

Use interval-based approach for variant calling. Recommended for most use cases.

```yaml
intervals: true  # default: true
```

### sentieon

Use Sentieon tools instead of GATK. Requires valid license.

```yaml
sentieon: false      # default: false
sentieon_lic: ""     # Path to Sentieon license file
```

### ploidy

Ploidy for variant calling.

```yaml
ploidy: 2  # default: 2
```

### Interval Options

Only applicable when `intervals: true`.

```yaml
minNmer: 500              # Minimum N-mer for splitting (default: 500, min: 50)
num_gvcf_intervals: 50    # Maximum GVCF intervals (default: 50)
db_scatter_factor: 0.15   # GenomicsDB scatter factor (default: 0.15)
```

### HaplotypeCaller Options

Tuning for coverage depth. Use low-coverage settings for <10x, high-coverage for >10x.

```yaml
# Low coverage (<10x) - default
minP: 1  # --min-pruning
minD: 1  # --min-dangling-branch-length

# High coverage (>10x)
# minP: 2
# minD: 4

het_prior: 0.005  # Prior probability of heterozygous site
```

## Read Processing Options

### mark_duplicates

Mark optical duplicates before variant calling.

```yaml
mark_duplicates: true  # default: true
```

### sort_reads

Sort reads by name before adapter trimming.

```yaml
sort_reads: false  # default: false
```

## Callable Sites Options

### Mappability

```yaml
mappability_min: 1      # Minimum mappability score (default: 1)
mappability_k: 150      # K-mer size for mappability (default: 150)
mappability_merge: 100  # Merge regions within this distance (default: 100)
```

### Coverage Filtering

Enable coverage-based filtering of callable sites.

```yaml
cov_filter: true  # default: true
cov_merge: 100    # Merge regions within this distance (default: 100)
```

**Coverage threshold options** (use only one approach):

Option 1 - Absolute thresholds:
```yaml
cov_threshold_lower: 1
cov_threshold_upper: 10000
```

Option 2 - Standard deviations (assumes Poisson distribution):
```yaml
cov_threshold_stdev: 2  # Remove regions outside X standard deviations
```

Option 3 - Relative scaling:
```yaml
cov_threshold_rel: 2  # Remove regions outside (mean/X) to (mean*X)
```

## QC Options

```yaml
nClusters: 3       # Number of clusters for PCA (default: 3)
min_depth: 2       # Minimum average depth for QC inclusion (default: 2)
GoogleAPIKey: ""   # Google Maps API key for geographic plots (optional)
```

## Filtering Options

Applied during postprocessing to create "clean" VCF.

```yaml
contig_size: 10000           # Exclude SNPs on contigs smaller than this (default: 10000)
maf: 0.01                    # Minimum minor allele frequency (default: 0.01)
missingness: 0.75            # Maximum missingness (default: 0.75)
scaffolds_to_exclude: "mtDNA,Y"  # Comma-separated scaffolds to exclude
```

## Output Options

### generate_trackhub

Generate UCSC Genome Browser trackhub with population genomic statistics.

```yaml
generate_trackhub: true   # default: true
trackhub_email: ""        # Required if trackhub enabled
```

## Advanced Options

### Remote Reads

For cloud execution with reads in remote storage.

```yaml
remote_reads: false       # default: false
remote_reads_prefix: ""   # Remote storage prefix
```

## Example Configuration

```yaml
# Required
samples: "config/samples.csv"
reference:
  name: "my_organism"
  path: "/data/reference/genome.fna.gz"  # Must be bgzip-compressed

# Variant calling
intervals: true
sentieon: false
ploidy: 2

# Read processing
mark_duplicates: true
sort_reads: false

# Callable sites
cov_filter: true
mappability_min: 1
cov_threshold_lower: 1
cov_threshold_upper: 10000

# QC
nClusters: 3
min_depth: 2

# Filtering
contig_size: 10000
maf: 0.01
missingness: 0.75

# Output
generate_trackhub: false
```

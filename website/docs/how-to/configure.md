# How to configure your run

snpArcher is configured through a YAML file at `config/config.yaml` in your project directory.
The settings below are the ones you are most likely to change.
For the full specification of every option, see the [config schema reference](../reference/config-schema.md).

## Minimal configuration

You need to set three things: the path to your sample sheet, a reference genome name, and a reference genome source.

```yaml
samples: "config/samples.csv"  # <-- change this

reference:
  name: "my_organism"  # <-- change this (used in output filenames)
  source: "/data/reference/genome.fna.gz"  # <-- change this
```

The `reference.source` field accepts:

- A **local file path** to a FASTA (optionally gzipped)
- A **URL** to a remote FASTA
- A **GenBank or RefSeq accession** (e.g., `GCF_000002315.6`). snpArcher will download and index it automatically.

## Choose a variant caller

Set `variant_calling.tool` to select the caller.
The default is `gatk`.

```yaml
variant_calling:
  tool: "gatk"  # <-- change this
```

| Tool | Description | Notes |
|---|---|---|
| `gatk` | GATK HaplotypeCaller + GenomicsDB + GenotypeGVCFs | Default. Well-tested for population genomics. |
| `sentieon` | Sentieon DNAseq (GATK-compatible, faster) | Requires a valid license (`variant_calling.sentieon.license`). |
| `bcftools` | bcftools mpileup + call | No per-sample gVCF step. Cannot accept `gvcf` input. |
| `deepvariant` | Google DeepVariant + GLnexus | ML-based caller. Cannot accept `gvcf` input. |
| `parabricks` | NVIDIA Parabricks | GPU-accelerated. Requires Apptainer image and NVIDIA GPUs. Cannot accept `gvcf` input. |

!!! note
    You only need to configure the sub-block for the tool you selected.
    Settings for other callers are ignored.

## Set the heterozygosity prior (GATK)

If you are using `gatk` and working with a non-model organism, you should adjust the heterozygosity prior.
The GATK default (0.001) is calibrated for humans.
Many non-model organisms have higher diversity, and raising this value improves variant sensitivity.

```yaml
variant_calling:
  tool: "gatk"
  gatk:
    het_prior: 0.005  # <-- adjust for your species
```

0.005 is a reasonable starting point for species with moderate to high nucleotide diversity.
For high-diversity species (e.g., some invertebrates), try 0.01 or higher.

For more background, see [Variant calling approaches](../explanation/variant-calling.md).

## Set ploidy and coverage

```yaml
variant_calling:
  ploidy: 2  # <-- change if not diploid
  expected_coverage: "low"  # "low" for <30x, "high" for >30x
```

The `expected_coverage` setting controls internal caller tuning.
Use `low` for most resequencing projects (5-15x per sample) and `high` for datasets with consistently deep coverage (>30x).

## Enable modules

Optional modules are off by default.
Enable them under the `modules` block:

```yaml
modules:
  qc:
    enabled: true  # <-- interactive QC dashboard
    clusters: 3  # number of k-means clusters for PCA
    min_depth: 2  # exclude samples below this depth from QC
  postprocess:
    enabled: true  # <-- sample exclusion and site filtering
```

!!! tip
    We recommend enabling the QC module for every run.
    It costs little extra computation and helps catch problems early.

## Configure callable sites

Callable sites produce a BED file of regions that pass coverage and mappability filters.
This is important for analyses that need to distinguish "no variant" from "no data" (e.g., estimating nucleotide diversity).

Callable sites are enabled by default.
To customize the coverage thresholds:

```yaml
callable_sites:
  generate_bed_file: true
  coverage:
    enabled: true
    fraction: 1.0  # fraction of samples that must be callable
    min_coverage: auto  # "auto" = floor(mean_depth / 2)
    max_coverage: auto  # "auto" = ceil(mean_depth * 2)
    merge_distance: 100
  mappability:
    enabled: true
    kmer: 150
    min_score: 1  # 1 = uniquely mappable only
    merge_distance: 100
```

To disable callable sites entirely:

```yaml
callable_sites:
  generate_bed_file: false
  coverage:
    enabled: false
  mappability:
    enabled: false
```

## Configure postprocessing filters

When the postprocess module is enabled, you can control site-level filtering:

```yaml
modules:
  postprocess:
    enabled: true
    filtering:
      contig_size: 10000  # exclude contigs smaller than this (bp)
      maf: 0.01  # minimum minor allele frequency
      missingness: 0.75  # maximum fraction of missing genotypes per site
      exclude_scaffolds: "mtDNA,Y"  # <-- comma-separated, no spaces
```

Sample exclusions are controlled through a separate metadata file.
See [How to filter and postprocess](postprocess.md) for details.

## Full realistic example

Here is a complete `config.yaml` for a bird resequencing project using GATK with the QC and postprocess modules enabled:

```yaml
samples: "config/samples.csv"
sample_metadata: "config/sample_metadata.csv"

reference:
  name: "bTaeGut2"
  source: "GCF_003957565.2"

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

callable_sites:
  generate_bed_file: true
  coverage:
    enabled: true
    fraction: 1.0
    min_coverage: auto
    max_coverage: auto
  mappability:
    enabled: true
    kmer: 150
    min_score: 1

modules:
  qc:
    enabled: true
    clusters: 3
    min_depth: 2
  postprocess:
    enabled: true
    filtering:
      contig_size: 10000
      maf: 0.01
      missingness: 0.75
      exclude_scaffolds: "Z,W,mtDNA"
```

## Next steps

- [Run locally](run-local.md) on a workstation
- [Run on HPC](run-hpc.md) on a SLURM cluster
- See the [config schema reference](../reference/config-schema.md) for every available option

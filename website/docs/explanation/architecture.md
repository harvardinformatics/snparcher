# Pipeline architecture

snpArcher is a Snakemake workflow that takes raw sequencing reads (or pre-processed BAMs or gVCFs) and produces a joint-called, hard-filtered VCF.
This page explains how the pipeline is structured, why key design decisions were made, and how the pieces fit together.

## The core pipeline

The core pipeline follows the widely adopted GATK Best Practices workflow for germline short-variant discovery, adapted for non-model organisms.
The main stages are:

1. **Data acquisition**: If input data are SRA accessions, snpArcher downloads reads automatically.
   If the reference genome is specified as a URL or NCBI accession, it is fetched and indexed.
2. **Read preprocessing**: Adapter trimming and quality filtering with fastp.
3. **Alignment**: Reads are mapped to the reference genome with BWA-MEM.
4. **Duplicate marking**: PCR and optical duplicates are flagged (can be disabled per sample for amplicon protocols).
5. **Per-sample variant calling**: Each sample is genotyped independently, producing a gVCF (genomic VCF) that records genotype likelihoods at every site, not just variant positions.
6. **Database import**: Per-sample gVCFs are combined into a GenomicsDB datastore, which is an efficient on-disk format for storing multi-sample genotype data.
7. **Joint genotyping**: All samples are genotyped simultaneously from the GenomicsDB, producing a raw multi-sample VCF.
8. **Hard filtering**: GATK VariantFiltration applies quality filters (strand bias, mapping quality, read position bias, quality by depth) and tags sites as passing or failing.

The result is a filtered VCF with all samples joint-called together.
Every site in the VCF retains its filter tag, so you can always recover the unfiltered data.

!!! note
    snpArcher does not remove filtered sites from the VCF. It tags them.
    Downstream tools can then include or exclude filtered sites as appropriate.

## Why Snakemake?

snpArcher is built on Snakemake because:

- **Reproducibility.** Snakemake tracks inputs, outputs, and parameters for every step.
  If you change a configuration value, only the affected downstream steps are re-run.
  If a run is interrupted, it resumes from where it left off.
- **Conda integration.** Each rule specifies its own conda environment.
  Users never need to install BWA, GATK, samtools, or any other bioinformatics tool manually. Snakemake creates isolated environments automatically on first execution and caches them for reuse.
- **Cluster and cloud support.** Snakemake's executor plugin system handles job submission to SLURM, PBS, LSF, Google Cloud, and other backends.
  The same workflow runs on a laptop and on a 10,000-core cluster without modification to the pipeline logic.
- **Transparency.** The entire workflow is defined in human-readable Snakemake rules.
  There is no compiled binary or opaque orchestration layer; you can read exactly what commands are being run and why.

## Module system

snpArcher separates the core variant-calling pipeline from optional downstream analyses using Snakemake's module system.
The core pipeline always runs.
Optional modules are toggled on or off in the configuration file and extend the pipeline with additional rules.

The current modules are:

| Module | Purpose | Dependencies |
|--------|---------|--------------|
| **QC** | Interactive HTML dashboard with PCA, relatedness, depth, admixture, and geographic visualizations | Core pipeline |
| **Postprocess** | Sample exclusion and site-level filtering (MAF, missingness, callable sites, scaffold exclusion) | Core pipeline + callable sites |

Users who only need a VCF do not pay the computational cost of QC or postprocessing, and each module can be developed, tested, and versioned independently.
The [module contribution guidelines](https://github.com/harvardinformatics/snparcher) define criteria for adding new modules.

!!! tip "Enabling modules"
    Modules are controlled under the `modules` block in `config.yaml`.
    For example, to enable QC and postprocessing:

    ```yaml
    modules:
      qc:
        enabled: true
      postprocess:
        enabled: true
    ```

    See the [configuration how-to](../how-to/configure.md) for full details.

## Input flexibility

Not every project starts from raw FASTQ files, and not every sample comes from the same source.
snpArcher supports four input types that enter the pipeline at different points:

| Input type | Format | Enters pipeline at |
|------------|--------|--------------------|
| `fastq` | Paired-end FASTQ files | Read preprocessing (beginning) |
| `srr` | SRA run accessions | Data acquisition (snpArcher downloads reads) |
| `bam` | Aligned BAM files | Per-sample variant calling (skips alignment) |
| `gvcf` | Per-sample gVCF files | Database import (skips individual calling) |

These input types can be mixed freely within a single sample sheet.
This is particularly useful when combining newly sequenced samples with published data from the SRA, or when adding new samples to an existing dataset without re-processing the old ones.

For example, suppose you sequenced 20 new samples locally and want to combine them with 15 published samples from the SRA and 5 samples that were already processed through HaplotypeCaller in a previous study.
Your sample sheet might look like:

```text
sample_id,input_type,input
new_sample_01,fastq,/data/new_01_R1.fq.gz;/data/new_01_R2.fq.gz
...
sra_sample_01,srr,SRR12345678
...
old_sample_01,gvcf,/data/old_01.g.vcf.gz
...
```

snpArcher routes each sample through the appropriate pipeline stages and merges them all at the joint genotyping step.

## One reference genome per run

In version 2, snpArcher requires that all samples in a single run share one reference genome.
This is a deliberate simplification from version 1, which allowed multiple reference genomes per sample sheet (grouping samples by their `refGenome` column).

Joint genotyping requires all samples to be called against the same coordinate system.
While v1 technically supported multiple references by processing each group independently, this created complexity in the sample sheet, confusion about which samples were being jointly genotyped, and edge cases in module behavior.
In v2, the reference is specified once in the configuration file, and every sample in the sample sheet is processed against it.

If you need to call variants against two different reference assemblies, run snpArcher twice with separate configuration files and sample sheets.

## Workflow profiles: decoupling resources from logic

A recurring pain point in bioinformatics pipelines is the mixing of computational resource settings (how many cores, how much memory, which partition) with analytical logic (what tools to run, what parameters to use).
snpArcher separates these concerns using Snakemake's workflow profile system.

The **configuration file** (`config/config.yaml`) controls analytical decisions:

- Which variant caller to use
- What heterozygosity prior to set
- Whether to enable the QC module
- Interval generation parameters

The **workflow profile** (`workflow-profiles/default/config.yaml`) controls resource allocation:

- Threads per rule
- Memory per rule
- SLURM partition, account, and runtime limits
- Temporary file locations

When you share a configuration with a collaborator, they can adjust resources for their own cluster without changing the scientific parameters of the analysis.

snpArcher ships with a default workflow profile that provides reasonable resource settings for most datasets.
Many rules use Snakemake's `attempt` multiplier for memory, so if a job fails with an out-of-memory error, the retry will automatically request more memory:

```yaml
set-resources:
  gatk_haplotypecaller:
    mem_mb: attempt * 4000
```

The first attempt gets 4 GB, the second attempt gets 8 GB, and so on.
Combined with the `--retries` flag, this provides automatic recovery from transient memory issues without manual intervention.

## How the DAG unfolds

When you launch snpArcher, Snakemake resolves the full directed acyclic graph (DAG) of jobs before executing anything.
For a typical 20-sample dataset with 50 calling intervals, the DAG might contain several thousand jobs.
The rough shape is:

1. **Fan-out by sample**: Each sample gets its own trimming, alignment, and duplicate-marking jobs.
2. **Fan-out by interval**: Per-sample variant calling is parallelized across genomic intervals (see [parallelization](parallelization.md)).
   With 20 samples and 50 intervals, that is 1,000 HaplotypeCaller jobs.
3. **Fan-out by database interval**: GenomicsDB import is parallelized across a potentially different number of intervals controlled by `db_scatter_factor`.
4. **Fan-in**: Joint genotyping, hard filtering, and VCF concatenation merge results back into a single output.
5. **Module branches**: QC and postprocess modules branch off from the core outputs.

The parallelism in steps 2 and 3 is what makes snpArcher fast on clusters: hundreds of independent HaplotypeCaller jobs can run simultaneously, each processing a small genomic region.
See [parallelization](parallelization.md) for a detailed discussion of how intervals are generated and why the scatter-by-Ns strategy matters.

## Callable sites: a parallel track

The callable sites subsystem runs alongside the core pipeline and produces a BED file of genomic regions that pass both coverage and mappability filters.
This BED file is not used by the core pipeline itself, but it is critical for downstream analyses that need to distinguish "no variant called" from "no data at this site."

The callable sites track has two components:

- **Coverage filtering** uses aligned BAMs to identify regions where depth falls within acceptable bounds (by default, between half and twice the global mean coverage across all BAM-backed samples).
  Regions outside these bounds are masked.
- **Mappability filtering** uses GenMap to identify regions where short reads cannot map uniquely.
  Regions with ambiguous mappability (score below the threshold, default 1.0 = perfectly unique) are masked.

The final callable sites BED is the intersection of the coverage-passing and mappability-passing regions.
This BED file is then used by the postprocess module to restrict the filtered VCF to reliably genotyped sites.

Understanding what the callable sites BED represents is essential for analyses that depend on the denominator: the total number of sites examined, not just the number of variant sites.
Estimates of nucleotide diversity, Tajima's D, and demographic parameters from the SFS all require knowing how many sites were genotyped, not just how many were polymorphic.
See [filtering philosophy](filtering.md) for a deeper discussion.

## Design philosophy

snpArcher ships with defaults that work for most non-model population genomics projects (GATK HaplotypeCaller, ploidy 2, heterozygosity prior of 0.005, scatter-by-Ns parallelization), but every default can be overridden.
A first-time user gets a reasonable result with minimal configuration; an experienced user can tune every parameter.

Hard filtering rather than VQSR is a deliberate choice.
VQSR requires a truth set of known variants, which does not exist for non-model organisms.
Hard filtering works for every species, does not require training data, and produces reproducible results.
See [variant calling](variant-calling.md) for more on this decision.

Every step of the pipeline is a Snakemake rule with explicit inputs, outputs, and shell commands.
If something goes wrong, you can read the rule, check the log, and understand exactly what happened.

Finally, standardization enables apples-to-apples comparisons across studies.
When two labs use the same pipeline with the same configuration, differences in the results reflect biology rather than methodology.
The California Conservation Genomics Project, for example, used snpArcher to process datasets from approximately 150 species through one analytical framework.

## Further reading

- [Parallelization](parallelization.md): How scatter-by-Ns works and why it matters for non-model genomes.
- [Variant calling](variant-calling.md): The five supported callers and when to use each.
- [Non-model organisms](non-model.md): How snpArcher's design addresses the specific challenges of non-model species.
- [Configuration how-to](../how-to/configure.md): Practical guide to setting up `config.yaml`.

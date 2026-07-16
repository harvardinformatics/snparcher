# Parallelization: the scatter-by-Ns strategy

Variant calling is computationally expensive.
A single sample at 10x coverage against a 1 Gb genome can take many hours of CPU time for HaplotypeCaller alone, and a 50-sample dataset multiplies that proportionally.
snpArcher splits the genome into intervals and processes each interval independently.

This page explains how interval generation works, why snpArcher uses runs of N bases as split points, and how to tune the parallelization parameters for your dataset.

## The basic idea

Variant calling operates on reads aligned to a reference genome.
Because reads at position 1,000,000 are independent of reads at position 50,000,000, there is no reason to process them sequentially.
By dividing the genome into non-overlapping intervals, snpArcher can run HaplotypeCaller (or another caller) on each interval as a separate job.
On a cluster, these jobs run simultaneously, turning a 20-hour single-job run into a 20-minute wall-clock run with enough nodes.

This divide-and-conquer approach applies to three stages of the pipeline:

1. **Per-sample variant calling** (HaplotypeCaller): Each sample is called independently on each interval.
2. **GenomicsDB import**: Per-sample gVCFs are combined into the database, interval by interval.
3. **Joint genotyping** (GenotypeGVCFs): Joint calls are made on each interval, then concatenated.

## Why split at N-mers?

The traditional approach (used by GATK's own tools and many other pipelines) is to split the genome into fixed-size intervals (for example, every 1 Mb) or to split by chromosome.
Both approaches have problems for non-model genomes.

**Fixed-size splitting** can break the genome in the middle of a gene, a repeat, or a region of active variation.
While variant callers handle interval boundaries correctly in most cases, the approach ignores any structural information about the assembly.
More importantly, fixed-size intervals distribute work unevenly: a 1 Mb interval in a gene desert requires much less computation than a 1 Mb interval in a highly polymorphic region.

**Splitting by chromosome** works well for model organisms with chromosome-level assemblies (human has 24 chromosomes, providing 24 intervals).
But non-model genomes are often fragmented assemblies with hundreds or thousands of scaffolds.
A scaffold-level avian genome might have 1,000 scaffolds ranging from 100 kb to 100 Mb.
Using each scaffold as an interval produces wildly uneven job sizes: the largest scaffold takes hours while the smallest finishes in seconds, wasting cluster resources on scheduling overhead.

**Scatter-by-Ns** is snpArcher's alternative.
It splits the genome at runs of N characters in the reference assembly.
Runs of Ns represent gaps in the assembly: places where the assembler could not resolve the sequence, typically at the boundaries between contigs within a scaffold.
These gaps are natural split points because:

- No real sequence data exists in the gap, so nothing is lost by splitting there.
- The gaps often correspond to biological boundaries (centromeres, heterochromatic regions, repetitive elements) where variant calling would be unreliable anyway.
- For fragmented assemblies, there are many gaps, providing ample split points.

The result is intervals that follow the gap structure of the assembly and distribute work more evenly than chromosome-based splitting.

## How interval generation works

snpArcher's interval generation is controlled by three parameters in the `intervals` section of the configuration file:

### `min_nmer` (default: 500)

This parameter sets the minimum length of an N-run that qualifies as a split point.
The genome is scanned for consecutive runs of N characters, and only runs of at least `min_nmer` bases are used to define interval boundaries.

A lower value (e.g., 50-100) splits at smaller gaps, producing more intervals and finer-grained parallelism; useful for highly fragmented assemblies where large gaps are rare.
A higher value (e.g., 500-1000) splits only at major gaps, producing fewer, larger intervals appropriate for chromosome-level assemblies.

The default of 500 works well for most vertebrate genomes, which typically have gaps of several hundred to several thousand Ns between contigs within a scaffold.

### `num_gvcf_intervals` (default: 50)

After the genome is split at N-runs, the resulting segments are grouped into a target number of intervals for gVCF calling.
This parameter controls that target.
Segments are combined (in order along the genome) until the desired number of intervals is reached, with the goal of producing intervals of roughly equal genomic length.

Increasing this (e.g., 100-200) gives more parallelism at the cost of more scheduling overhead and more intermediate files (each interval produces a separate gVCF per sample).
Decreasing it (e.g., 10-20) reduces overhead but increases per-job runtime.

For a typical vertebrate genome (1-3 Gb), 50 intervals provide a good balance: each interval is roughly 20-60 Mb, taking 10-30 minutes per sample on HaplotypeCaller with moderate resources.

!!! tip "Tuning for your cluster"
    If your cluster has many nodes but jobs take a long time to start (long queue wait times), fewer larger intervals may be more efficient.
    If your cluster starts jobs quickly and you have many available slots, more intervals will reduce wall-clock time.

### `db_scatter_factor` (default: 0.15)

GenomicsDB import and joint genotyping use a different number of intervals than gVCF calling.
The number of database intervals is calculated as:

```text
num_db_intervals = db_scatter_factor * num_samples * num_gvcf_intervals
```

This formula scales the parallelism of the joint genotyping step with the number of samples, because GenomicsDB import processes all samples at each interval and memory requirements grow with sample count.

For example, with 20 samples and 50 gVCF intervals:

```text
num_db_intervals = 0.15 * 20 * 50 = 150
```

This produces 150 intervals for GenomicsDB import and joint genotyping, each covering a smaller region of the genome and requiring less memory per job.

For large cohorts (100+ samples) where GenomicsDB import is memory-constrained, raise this value to produce more, smaller database intervals.
For small cohorts where memory is not a concern, the default or a lower value keeps overhead minimal.

## Comparison with fixed-size interval splitting

The benchmarks in Mirchandani et al. (2024) directly compared scatter-by-Ns against GATK's standard approach of splitting by chromosome or fixed-size intervals.
For non-model organisms with fragmented assemblies, scatter-by-Ns consistently produced better load balancing and shorter wall-clock times.

The advantage is most pronounced for assemblies with many scaffolds.
Consider a genome with 500 scaffolds:

- **By chromosome/scaffold**: You get 500 intervals, but most are tiny and finish instantly while a few large scaffolds dominate runtime.
- **By fixed size (1 Mb)**: For a 1 Gb genome, you get ~1,000 intervals.
  Some split in the middle of problematic regions; load balancing is better but not optimal.
- **Scatter-by-Ns (50 intervals)**: The genome is split at natural gap boundaries and then grouped into 50 roughly equal intervals.
  Each interval takes approximately the same amount of time.

For chromosome-level assemblies with few gaps, scatter-by-Ns produces fewer split points, and the intervals may be less balanced.
In that case, increasing `num_gvcf_intervals` or decreasing `min_nmer` can help.

## When scatter-by-Ns does not work well

There are a few scenarios where the default scatter-by-Ns strategy may need adjustment:

### Highly contiguous genomes

A telomere-to-telomere (T2T) assembly may have very few runs of Ns.
If the genome has only 5 gaps, scatter-by-Ns can only produce at most 6 intervals, regardless of the `num_gvcf_intervals` setting.
In this case, you may want to disable interval generation entirely and let GATK process each chromosome as a single unit, or consider supplementing with artificial split points.

Most current non-model genome assemblies still have enough gaps for scatter-by-Ns to work well; even high-quality assemblies from the Vertebrate Genomes Project typically retain gaps at centromeres and in heterochromatic regions.

### Very small genomes

For compact genomes (e.g., bacterial genomes at 5 Mb, or some insect genomes under 200 Mb), the overhead of creating many intervals may outweigh the parallelism benefit.
Reducing `num_gvcf_intervals` to 5-10 or disabling intervals entirely (`intervals.enabled: false`) can simplify the workflow without meaningful performance loss.

### Extremely fragmented assemblies

Assemblies with tens of thousands of very short contigs (e.g., a draft assembly with N50 of 10 kb) present the opposite problem: there are too many potential split points and many contigs are shorter than a single interval.
In this case, the grouping step (combining segments to reach `num_gvcf_intervals`) handles the situation naturally, but the assembly quality itself may be a bigger concern for variant calling accuracy than the parallelization strategy.

## Performance trade-offs

More intervals generally means faster wall-clock time but greater total overhead:

| Factor | More intervals | Fewer intervals |
|--------|---------------|-----------------|
| Wall-clock time | Lower (more parallelism) | Higher (less parallelism) |
| Total CPU time | Slightly higher (overhead per job) | Slightly lower |
| Intermediate files | More (one gVCF per sample per interval) | Fewer |
| Disk space | Higher | Lower |
| Scheduler pressure | Higher (more jobs submitted) | Lower |
| Concatenation time | Longer (more files to merge) | Shorter |

For most users, the default of 50 gVCF intervals is a reasonable default.
If you find that your cluster's job scheduler is a bottleneck (e.g., jobs spend more time in the queue than running), reduce the number of intervals.
If individual jobs are running for many hours and your cluster has idle capacity, increase the number of intervals.

## Disabling interval parallelization

If you prefer not to use scatter-by-Ns (for example, because you are running bcftools, which has its own parallelization), you can disable it:

```yaml
intervals:
  enabled: false
```

When intervals are disabled, the entire genome is processed as a single unit for variant calling.
This is simpler but much slower for large genomes with GATK.

!!! note
    Parabricks always uses interval-based parallelization for the joint genotyping step (GenomicsDB import and GenotypeGVCFs), regardless of the `intervals.enabled` setting.
    The `intervals.enabled` setting controls only the per-sample HaplotypeCaller step for the GATK backend.

## gVCF concatenation

After per-interval variant calling and joint genotyping, the interval-level VCFs and gVCFs must be concatenated back into whole-genome files.
snpArcher uses a staged concatenation strategy controlled by two parameters:

- **`concat_batch_size`** (default: 250): The maximum number of interval files merged per concatenation job.
  Lower values reduce per-job file handle pressure; higher values reduce the number of merge rounds.
- **`concat_max_rounds`** (default: 20): A safety limit on the number of staged merging rounds.
  If the number of intervals requires more rounds than this, snpArcher raises a configuration error rather than silently producing an extremely deep merge tree.

For most datasets with the default 50 gVCF intervals, concatenation completes in a single round.
Large datasets with many database intervals may require two or three rounds.
This staged approach avoids the operating system limit on open file descriptors that can cause problems when trying to merge thousands of files simultaneously.

## A worked example

Suppose you have a scaffold-level bird genome (1.2 Gb, ~800 scaffolds, N50 = 20 Mb) and 30 samples at 10x coverage.

With the default settings (`min_nmer: 500`, `num_gvcf_intervals: 50`, `db_scatter_factor: 0.15`):

1. The genome is scanned for N-runs of at least 500 bases.
   A typical avian genome might have ~200 such gaps, splitting the genome into ~200 segments.
2. These segments are grouped into 50 intervals of roughly 24 Mb each.
3. HaplotypeCaller runs on 30 samples x 50 intervals = 1,500 jobs.
   Each job processes ~24 Mb of genome for one sample, taking roughly 15-30 minutes.
4. GenomicsDB import uses `0.15 * 30 * 50 = 225` intervals.
5. Joint genotyping runs on 225 intervals, each covering ~5 Mb of genome across all 30 samples.

On a cluster with 100 available cores, the HaplotypeCaller step completes in roughly 15 batches of ~100 jobs each, for a wall-clock time of about 4-8 hours.
Without interval parallelization, the same step would take 30 x 4-8 hours = 120-240 hours of serial processing.

## Further reading

- [Pipeline architecture](architecture.md): How interval parallelization fits into the overall DAG.
- [Variant calling](variant-calling.md): How the callers use intervals.
- [Configuration reference](../reference/config-schema.md): Full specification of interval parameters.
- [Running on HPC](../how-to/run-hpc.md): Practical guidance for cluster execution.

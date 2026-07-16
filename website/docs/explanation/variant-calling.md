# Variant calling

snpArcher supports five variant calling backends.
This page explains what each one does, when to choose it, and how snpArcher's approach to joint genotyping and hard filtering works.

## The five supported callers

| Caller | Type | License | Hardware | Joint genotyping path |
|--------|------|---------|----------|----------------------|
| **GATK HaplotypeCaller** | Local reassembly + PairHMM | Open source | CPU | gVCF -> GenomicsDB -> GenotypeGVCFs |
| **bcftools** | Pileup-based (mpileup + call) | Open source | CPU | Direct multi-sample calling |
| **DeepVariant** | Deep learning (CNN) | Open source | CPU (GPU optional) | gVCF -> GLnexus |
| **Sentieon** | GATK-compatible, optimized | Commercial | CPU | gVCF -> GenomicsDB -> GenotypeGVCFs |
| **Parabricks** | GPU-accelerated GATK | Commercial (NVIDIA) | GPU | gVCF -> GenomicsDB -> GenotypeGVCFs |

You choose a caller by setting `variant_calling.tool` in your configuration file.
Only one caller is active per run.

## GATK HaplotypeCaller (default)

GATK HaplotypeCaller is the default and most thoroughly tested caller in snpArcher.
It works by performing local de novo assembly of haplotypes in regions where there is evidence of variation, then evaluating each read against each candidate haplotype using a pair hidden Markov model (PairHMM).
This approach is more sensitive than simple pileup-based methods, particularly for indels and in repetitive regions.

HaplotypeCaller produces gVCFs (genomic VCFs), which record genotype likelihoods at every position in the genome, not just variant sites.
This is essential for joint genotyping: when samples are combined later, the genotype likelihoods at non-variant sites in one sample inform the joint call across all samples.

**When to use GATK:**
GATK is the best choice for most projects.
It is well-validated across a wide range of organisms, has extensive documentation, and is the basis for GATK Best Practices.
The benchmarks in Mirchandani et al. (2024) showed that GATK produces high-quality calls across diverse vertebrate taxa when paired with appropriate hard filters.

**Key configuration:**

```yaml
variant_calling:
  tool: "gatk"
  ploidy: 2
  expected_coverage: "low"
  gatk:
    het_prior: 0.005
```

### The heterozygosity prior

The `het_prior` parameter deserves special attention because it is one of the few settings where the default may need adjustment for your organism.

GATK's genotype caller uses a Bayesian framework.
The heterozygosity prior is the expected probability that any given site in the genome is heterozygous.
GATK's own default (0.001) was calibrated for humans, who have relatively low nucleotide diversity (~0.001 per bp).
snpArcher sets a default of 0.005, which is more appropriate for many non-model organisms, but you may need to adjust it further.

**Why does it matter?**
A prior that is too low makes the caller conservative about calling heterozygous genotypes.
In a species with high diversity (say, a marine invertebrate with per-site heterozygosity of 0.02), a prior of 0.001 will miss real heterozygous sites because the model does not expect them to be common.
Conversely, a prior that is too high can inflate false-positive heterozygous calls in low-diversity species.

**How to choose a value:**
If you have a published estimate of nucleotide diversity (pi or theta) for your species or a close relative, use that.
If not, the snpArcher default of 0.005 is a reasonable middle ground for many vertebrates.
For highly diverse invertebrates or marine organisms, values of 0.01-0.02 may be appropriate.
For species with very low diversity (island populations, recent bottlenecks), the GATK default of 0.001 may actually be correct.

!!! tip "Checking your prior"
    After a first run, examine the ratio of heterozygous to homozygous variant calls in your VCF.
    If the het/hom ratio seems unusually low for your species, the prior may be too conservative.
    You can re-run joint genotyping (from the existing gVCFs) with an adjusted prior without re-running the entire pipeline.

See [non-model organisms](non-model.md) for a broader discussion of why priors calibrated for humans often fail for other species.

## bcftools

bcftools uses a pileup-based approach: it counts alleles at each position across all reads and applies a statistical model to call variants.
This is conceptually simpler than GATK's local reassembly and substantially faster.

In snpArcher, bcftools calling uses `bcftools mpileup` followed by `bcftools call`.
The key parameters are minimum mapping quality (`min_mapq`, default 20), minimum base quality (`min_baseq`, default 20), and maximum per-file depth (`max_depth`, default 250).

**When to use bcftools:**
bcftools is a good choice for preliminary analyses where speed matters more than maximum sensitivity, or for organisms where GATK's local reassembly adds little value (e.g., when you are only interested in SNPs, not indels).
It is also useful for very large datasets where the computational cost of GATK is prohibitive.

**Limitations:**
bcftools does not produce gVCFs, so it cannot take advantage of the gVCF -> GenomicsDB -> joint genotyping workflow that GATK uses.
Instead, it performs multi-sample calling directly.
You cannot incrementally add new samples to an existing callset without re-running the entire calling step.
Additionally, bcftools is generally less sensitive for indels and in low-complexity regions.

!!! warning
    When `variant_calling.tool` is set to `bcftools`, sample rows with `input_type: gvcf` are not supported.
    bcftools works from BAM files, not gVCFs.

## DeepVariant

DeepVariant uses a convolutional neural network (CNN) to call variants.
It converts pileup data into image-like tensors and classifies each candidate site as homozygous reference, heterozygous, or homozygous variant.
In PrecisionFDA Truth Challenge benchmarks, DeepVariant has consistently ranked among the highest-accuracy open-source callers for both SNPs and indels.

In snpArcher, DeepVariant produces gVCFs that are then merged with GLnexus for joint genotyping.
The key configuration parameter is `model_type` (default `WGS` for whole-genome sequencing data; other options include `WES`, `PACBIO`, and `ONT_R104` for different data types).

**When to use DeepVariant:**
DeepVariant is the best choice when genotyping accuracy is the top priority and computational resources are available.
It is particularly strong for indel calling and performs well even with relatively low coverage.
However, it is slower than bcftools and can be slower than GATK depending on the parallelization setup.

**Considerations for non-model organisms:**
DeepVariant's CNN was trained primarily on human data.
It generalizes well to many vertebrates, but its performance on genomes with high repeat content, extreme GC bias, or polyploidy lacks validation.
For highly divergent taxa, validate a subset of calls against an independent method.

## Sentieon

Sentieon provides a commercially licensed, performance-optimized implementation of the GATK algorithms.
It produces near-identical results to GATK but runs significantly faster through better CPU utilization.

Sentieon is a drop-in replacement if your institution already has a license.
The workflow is identical to the GATK path (gVCF -> GenomicsDB -> GenotypeGVCFs), and the outputs are compatible.

A license must be specified in the configuration:

```yaml
variant_calling:
  tool: "sentieon"
  sentieon:
    license: "/path/to/license/or/server:port"
```

## Parabricks

NVIDIA Parabricks provides GPU-accelerated implementations of GATK algorithms.
HaplotypeCaller on GPU can be 10-30x faster than on CPU, but Parabricks requires NVIDIA GPU hardware and a container image.

Parabricks makes sense when you have access to GPU nodes (increasingly common on modern HPC systems) and need to process large datasets quickly.

**Configuration:**

```yaml
variant_calling:
  tool: "parabricks"
  parabricks:
    container_image: "/path/to/parabricks.sif"
    num_gpus: 1
    num_cpu_threads: 16
```

Parabricks uses the same joint genotyping path as GATK (GenomicsDB -> GenotypeGVCFs) and uses interval-based parallelization for the joint genotyping step regardless of the `intervals.enabled` setting.

## Joint genotyping: why it matters

All five callers feed into joint genotyping, which is fundamental to how snpArcher works.

**The problem with single-sample calling:**
If you call variants in each sample independently and then merge the VCFs, you face two issues.
First, a site that is variant in sample A but reference in sample B will simply be absent from sample B's VCF, so you cannot tell whether sample B is truly homozygous reference or simply was not called at that site.
Second, you lose statistical power: the evidence for a variant at a given site accumulates across samples, and joint genotyping uses that shared evidence.

**The gVCF + GenomicsDB approach (GATK, Sentieon, Parabricks):**
These callers produce gVCFs that record genotype likelihoods at every site, including non-variant positions.
The gVCFs are then loaded into a GenomicsDB datastore, which provides efficient columnar access to the multi-sample data.
GenotypeGVCFs then jointly genotypes all samples simultaneously, using the evidence from all samples at each site to make the final call.

GenomicsDB can handle thousands of samples, and the gVCF format means that adding new samples to an existing dataset requires only running HaplotypeCaller on the new samples and re-importing them into the database. Existing gVCFs do not need to be recomputed.

**Scalability considerations:**
GenomicsDB import is parallelized across genomic intervals (controlled by `db_scatter_factor` in the snpArcher configuration).
Memory requirements scale with the number of samples and the number of variants per interval.
For large cohorts (hundreds of samples), GenomicsDB import can be the most memory-intensive step in the pipeline.
See [parallelization](parallelization.md) for details on how `db_scatter_factor` controls this.

## Hard filtering vs. VQSR

After joint genotyping, variant sites need to be filtered to remove artifacts.
There are two main approaches: hard filtering and Variant Quality Score Recalibration (VQSR).

### VQSR

VQSR is a machine-learning approach that trains a Gaussian mixture model on a set of known-true variants (a "truth set") to learn which combinations of quality annotations distinguish real variants from artifacts.
It then applies this learned model to the full callset, assigning each variant a quality score.

VQSR produces excellent results when a high-quality truth set is available.
For humans, this means resources like dbSNP, HapMap, and the 1000 Genomes Project.
For a handful of well-studied model organisms (mouse, *Drosophila*, *Arabidopsis*), comparable truth sets exist.

### Hard filtering

Hard filtering applies fixed thresholds to variant quality annotations.
snpArcher uses the GATK-recommended hard filters:

- **QD** (Quality by Depth): Variant quality normalized by depth.
  Low values suggest the variant call is driven by few reads.
- **FS** (Fisher Strand): Strand bias estimated by Fisher's exact test.
  High values indicate reads supporting the variant come disproportionately from one strand.
- **SOR** (Strand Odds Ratio): Another strand bias metric, more robust to high-depth sites.
- **MQ** (Mapping Quality): Average mapping quality across all reads at the site.
  Low values indicate reads that do not map uniquely.
- **MQRankSum**: Comparison of mapping quality between reads supporting reference and variant alleles.
- **ReadPosRankSum**: Whether variant-supporting reads are concentrated at the ends of reads (suggesting alignment artifacts).

### Why snpArcher uses hard filters

For non-model organisms, there is no truth set.
You cannot train VQSR without one, and this constraint applies to any pipeline, not just snpArcher.

Hard filtering works for every species: no training data, deterministic output, and thresholds well-characterized from the GATK literature.
Some real variants will be filtered, and some artifacts will pass, but the defaults are a reasonable baseline that you can refine with additional filters in the postprocessing step.

!!! note "Can I use VQSR instead?"
    If you have a truth set for your organism, you can run VQSR on the snpArcher output VCF outside the pipeline.
    snpArcher does not include VQSR as a built-in option because it would only work for a small number of species, and misapplying VQSR (e.g., with an inadequate truth set) can produce worse results than hard filtering.

See [filtering philosophy](filtering.md) for a detailed discussion of how to evaluate and refine filtering using the site frequency spectrum, and [non-model organisms](non-model.md) for more on why the absence of truth sets shapes so many design decisions in non-model genomics.

## Choosing a caller: practical guidance

| Situation | Caller | Why |
|-----------|--------|-----|
| Most projects | GATK | Well-validated, widely used, supports gVCF-based incremental analysis |
| Maximum per-site accuracy | DeepVariant | Strongest in current benchmarks with GLnexus joint genotyping, but higher compute cost |
| Preliminary exploration or very large datasets | bcftools | Fast; good for a first look before committing to a full GATK run |
| Sentieon license available | Sentieon | Equivalent results to GATK, faster. Purely a performance decision. |
| GPU nodes available | Parabricks | 10-30x speedup on GPU; useful when CPU queue times are long |

In practice, most snpArcher users run GATK.
The heterozygosity prior is the single most impactful configuration decision for call quality in non-model organisms. Spend more time thinking about `het_prior` than about which caller to use.

## Further reading

- [Pipeline architecture](architecture.md): How the callers fit into the overall pipeline.
- [Parallelization](parallelization.md): How variant calling is parallelized across genomic intervals.
- [Non-model organisms](non-model.md): Why het_prior, hard filtering, and caller choice matter more for non-model species.
- [Filtering philosophy](filtering.md): What happens after variant calling, evaluating and refining filters.
- [Configuration reference](../reference/config-schema.md): Full specification of all variant calling parameters.

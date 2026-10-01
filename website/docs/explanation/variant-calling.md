# Variant calling

snpArcher supports six variant calling backends.
This page explains what each one does, when to choose it, and how snpArcher's approach to joint genotyping and hard filtering works.

## The supported callers

| Caller | Type | License | Hardware | Joint genotyping path |
|--------|------|---------|----------|----------------------|
| **GATK HaplotypeCaller** | Local reassembly + PairHMM | Open source | CPU | gVCF -> GenomicsDB -> GenotypeGVCFs |
| **bcftools** | Pileup-based (mpileup + call) | Open source | CPU | Direct multi-sample calling |
| **RepAdapt model** | Pileup-based (bcftools 1.16, per-sample `call -G -`) | Open source | CPU | Direct multi-sample calling |
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

## RepAdapt calling model

`tool: repadapt` runs the calling model of [RepAdapt](https://github.com/RepAdapt/nextflow_snp_calling_linux)'s SNP-calling pipeline on snpArcher's BAMs.
Use it to produce call sets that can be combined or compared with RepAdapt datasets.

Per contig, snpArcher runs RepAdapt's command:

```
bcftools mpileup -Ou -f REF -r CONTIG BAM... -q 10 -I -a FMT/AD,FMT/DP \
  | bcftools call -G - -f GQ -mv --ploidy PLOIDY
```

and then applies RepAdapt's filter, `AC=AN || MQ < 30`, as two soft filters in `results/vcfs/filtered.vcf.gz` (see [filtering](filtering.md)).

**How it differs from bcftools calling.**
`call -G -` calls each sample independently: a site is kept if any sample supports it, and each sample's genotype prior comes from its own reads rather than from a cohort-wide allele frequency under Hardy-Weinberg equilibrium.
This avoids assuming one randomly mating population, which suits structured sample sets.
The costs are more false-positive sites at low coverage, and genotypes that lean homozygous at low depth, because a heterozygous call needs reads from both alleles.
`-I` skips indels, so the call set is SNPs only.
Base quality and depth use bcftools' defaults (minimum base quality 1, maximum depth 250 per file), and the `variant_calling.bcftools` settings have no effect.

**How it differs from RepAdapt's literal command.**
RepAdapt filters its BAMs to MAPQ >= 10 before calling and runs mpileup with `-q 5`; snpArcher uses `-q 10` in mpileup instead, which selects the same reads.
snpArcher also passes `--ploidy`, which gives identical calls for diploids.

**Sample order at extreme depth.**
At sites with extreme depth, such as collapsed repeats far above mpileup's 250-reads-per-file cap, bcftools' genotype likelihoods depend on the order in which BAMs are given.
snpArcher passes BAMs in sample-sheet order, so its output is reproducible.
RepAdapt passes them in the order its Nextflow tasks finish, which can change between runs.
On 10 *Candida albicans* isolates called from RepAdapt's own BAMs, snpArcher matched RepAdapt's output exactly except for 4 of 296,311 records, all in the collapsed rDNA, where one sample's PL and GQ differed.
With RepAdapt's sample order, all records were identical.

**Pinned bcftools.**
The caller uses bcftools 1.16, RepAdapt's version, with the exact package builds from RepAdapt's container image.
bcftools 1.16 stores INFO/MQ as an integer and newer versions as a float, which can move sites across the `MQ < 30` filter.
Those builds exist only for Linux, so full runs don't work on macOS (dry runs do).
Newer bcftools versions used elsewhere in the pipeline warn that `MQ should be declared as Type=Float` when reading these VCFs; the warning is harmless.

**Upstream differences.**
With snpArcher's default mapping, the BAMs differ from RepAdapt's: snpArcher uses `bwa mem -M`, does not filter BAMs by MAPQ, marks duplicates with sambamba, and does not realign indels.
The calling model is RepAdapt's, but the calls are "RepAdapt-adjacent" rather than identical.

!!! warning
    `tool: repadapt` supports `ploidy` 1 or 2 only (bcftools 1.16 accepts only its predefined ploidy aliases), and does not support samples with `input_type: gvcf`.
    Because it calls SNPs only, `modules.postprocess.filtering.split_by_type` produces an empty `clean_indels.vcf.gz`; snpArcher warns about this.

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

All the callers feed into joint genotyping, which is fundamental to how snpArcher works.

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
| Comparing or combining with RepAdapt datasets | RepAdapt model | RepAdapt's bcftools settings and filter, on snpArcher's BAMs |
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

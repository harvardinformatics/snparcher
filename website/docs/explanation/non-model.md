# Non-model organisms

snpArcher was designed specifically for population genomics in non-model organisms.
Many of its design decisions (hard filtering instead of VQSR, configurable heterozygosity priors, scatter-by-Ns parallelization, the callable sites module) exist because the standard approaches developed for human genomics do not transfer cleanly to other species.

This page explains what makes non-model organisms different and how snpArcher addresses each challenge.

## What makes an organism "non-model"?

In genomics, "model organism" is less about the biology of the species and more about the infrastructure available for studying it.
A model organism typically has:

- A high-quality, chromosome-level reference genome (often telomere-to-telomere)
- Curated gene annotations
- Databases of known variants (e.g., dbSNP for humans)
- Well-characterized population structure and demographic history
- Published best-practice pipelines tuned for that species
- Large community of researchers cross-validating results

A non-model organism lacks some or all of these.
The species you are studying might have a scaffold-level reference assembly published last year, no variant database, and a handful of published population genomic datasets.
You might be the first person to call variants in this species at this scale.

Many standard genomics tools were designed with model-organism assumptions baked in.
GATK's default heterozygosity prior assumes human-level diversity.
VQSR assumes you have a truth set of known variants.
Interval-splitting assumes a chromosome-level assembly.
Quality metrics like mapping rate have been calibrated against human reference genomes with decades of curation.

## Reference genome quality

The quality of your reference genome is the single biggest factor determining the quality of your variant calls, because every downstream step (read mapping, variant calling, genotyping, filtering) depends on the assembly being a faithful representation of the true genome.

### Scaffold-level vs. chromosome-level assemblies

Most non-model genomes available today are scaffold-level assemblies, with gaps (runs of Ns) where the assembler could not resolve the sequence.

A typical scaffold-level avian genome might have:

- Total assembly size: 1.0-1.3 Gb
- Number of scaffolds: 500-5,000
- N50: 5-50 Mb
- Total gap length: 50-200 Mb

Chromosome-level assemblies (from the Vertebrate Genomes Project, Darwin Tree of Life, or similar efforts) are increasingly available and substantially reduce the gaps and misassemblies that cause variant-calling artifacts.
When a chromosome-level assembly is available for your species or a close relative, use it.

### How reference quality affects variant calling

**Collapsed duplications.**
When the assembler cannot resolve two paralogous copies of a sequence, it may collapse them into one.
Reads from both copies then map to the collapsed region, creating:

- Artificially high depth (reads from two copies piling up)
- False heterozygous calls (fixed differences between copies appear as heterozygosity)
- A characteristic peak in the SFS at frequency 1/2

Collapsed duplications are more common in scaffold-level assemblies, but they can occur even in chromosome-level assemblies, particularly in recently duplicated regions.

**Missing sequence.**
Gaps in the assembly (runs of Ns) represent sequence that is absent from the reference.
Reads originating from these missing regions either fail to map (increasing the apparent unmapped fraction) or mismap to similar sequences elsewhere in the genome (creating false variants).
The scatter-by-Ns parallelization strategy (see [parallelization](parallelization.md)) uses these gaps as natural breakpoints, which also means that the regions immediately flanking gaps are not split in the middle of a calling interval.

**Assembly errors.**
Misassemblies (inversions, translocations, or chimeric scaffolds) create regions where the reference sequence does not match the true genome.
Reads spanning a misassembly breakpoint will have split or discordant alignments, leading to poor mapping quality and unreliable variant calls.
The hard filters for mapping quality (MQ) and mapping quality rank sum (MQRankSum) help catch some of these artifacts.

## Heterozygosity priors

This is one of the most consequential settings in snpArcher for non-model organisms.

### The problem with GATK's default

GATK HaplotypeCaller uses a Bayesian framework for genotyping.
The heterozygosity prior (`--heterozygosity` in GATK, `variant_calling.gatk.het_prior` in snpArcher) is the prior probability that any given site is heterozygous.
This prior influences the genotype likelihoods: a higher prior makes the model more willing to call heterozygous genotypes.

GATK's built-in default is 0.001, calibrated for humans.
Human nucleotide diversity is approximately 0.001 per base pair, so this prior reflects the expected rate of polymorphism in a typical human genome.

The problem is that many non-model organisms have much higher diversity:

| Taxon | Typical nucleotide diversity (pi) |
|-------|----------------------------------|
| Humans | ~0.001 |
| Great apes | ~0.001-0.002 |
| Most mammals | ~0.001-0.005 |
| Most birds | ~0.002-0.008 |
| Many fish | ~0.003-0.015 |
| Many insects | ~0.005-0.030 |
| Marine invertebrates | ~0.010-0.050+ |

For a marine invertebrate with pi = 0.03, using a prior of 0.001 means the model expects heterozygous sites to be 30 times rarer than they actually are, causing it to miss real heterozygous calls and underestimate diversity.

snpArcher sets a default of 0.005, which is more appropriate than 0.001 for many non-model vertebrates but still a compromise.
Use a published diversity estimate for your species if one exists; for highly diverse invertebrates, values of 0.01-0.02 are appropriate; for species with very low diversity, the GATK default of 0.001 may be correct.
See [variant calling](variant-calling.md) for details on choosing and validating this parameter.

## Hard filtering: the only universal option

VQSR requires a truth set of known variants (dbSNP, HapMap, etc.), and for non-model organisms these resources do not exist.
Applying VQSR with an inadequate truth set can produce worse results than not filtering at all.

Hard filtering works for every species because it uses statistical properties of the calls themselves rather than external training data.
The thresholds are not as precise as a well-trained VQSR model, but they are universal and reproducible.
If you do have truth data for your organism, you can always apply VQSR to the snpArcher output VCF as a post-hoc step.
See [variant calling](variant-calling.md) for the full rationale and filter specifications.

## Paralogy and collapsed duplications

Paralogy is one of the most persistent sources of error in non-model genomics.
Duplicated regions that are not properly resolved in the reference assembly create systematic false variants that can dominate the dataset if not filtered.

When two paralogous regions are collapsed into one in the reference assembly, reads from both copies map to the same location, inflating depth and creating false heterozygous calls.
The hallmark is a spike at frequency 0.5 in the SFS.

snpArcher addresses this through coverage filtering (the callable sites module flags regions with depth > 2x the genome-wide mean), the QD hard filter, and SFS diagnostics.
See [filtering philosophy](filtering.md) for the full explanation of how these filters interact and how to assess whether they are sufficient.

For organisms with recent whole-genome duplications (e.g., salmonids, many plant species), paralogy may require specialized approaches beyond what snpArcher provides.

## Sample size challenges

Non-model population genomics often works with small sample sizes.
A typical project might have 5-30 individuals, far fewer than the thousands of samples common in human genomics.

Small sample sizes affect several aspects of the analysis:

### Allele frequency estimation

With small samples, allele frequencies are estimated imprecisely.
A variant present in 1 of 10 diploid individuals could have a true population frequency anywhere from ~1% to ~20%.
This uncertainty propagates to all frequency-dependent analyses: SFS-based demographic inference, selection scans, and diversity estimates.

### Rare variant detection

Rare variants (those present in only one or a few individuals) are informative for demographic inference but are also the variants most likely to be genotyping errors.
In small samples, you have less power to distinguish real rare variants from errors, and MAF filtering removes a larger proportion of real variants.

### Statistical power for QC

The QC dashboard's ability to detect batch effects, population structure, and relatedness depends on sample size.
With 5 samples, PCA has only 4 possible dimensions and k-means clustering is nearly meaningless.
ADMIXTURE is unreliable with very small samples.
Relatedness estimation is less precise.

snpArcher's QC module requires a minimum of 3 samples, but realistically, the visualizations become useful at around 10+ samples and informative at 20+.

## Mixed data sources

Non-model genomics projects frequently combine data from multiple sources:

- Newly sequenced samples from your own collection
- Published raw reads from the SRA
- Pre-existing BAMs from a collaborator
- gVCFs from a previous snpArcher run

snpArcher's input flexibility (four input types that can be mixed in a single sample sheet) is designed for exactly this situation.
But combining data from different sources introduces its own challenges:

### Batch effects

Samples sequenced on different platforms, at different facilities, or with different library preparation protocols may have systematic quality differences.
These batch effects can mimic population structure in PCA and bias allele frequency estimates.

The QC dashboard's depth-PC correlation plot is your primary tool for detecting batch effects.
If PC1 separates samples by sequencing facility rather than geography, you have a batch effect problem.

### Variable coverage

Combining deeply sequenced new samples with shallowly sequenced archival data creates heterogeneous coverage.
This is the setting where the depth-PC correlation is most problematic: depth variation across samples can drive the top principal components.

The protocol paper recommends considering whether depth filtering (removing the lowest-coverage samples) resolves the batch effect.
If it does, the cost of losing a few samples may be worth the gain in data quality.

### Reference genome mismatch

If some samples were originally aligned to a different reference genome, they need to be re-processed from FASTQ (not from the existing BAMs) unless the same reference is used.
snpArcher v2 requires all samples in a run to share one reference genome, which prevents the subtle errors that arise from mixing coordinate systems.

## Why standardized pipelines matter

Using a standardized pipeline rather than a custom script enables comparability across studies.

When two labs independently develop variant-calling pipelines, they inevitably make different choices: different read trimmers, different alignment parameters, different hard-filter thresholds, different approaches to interval splitting.
Each choice is defensible in isolation, but the cumulative effect is that the two labs produce callsets that differ not only because of biological differences between their study organisms but also because of methodological differences.

This makes meta-analysis difficult.
If study A reports theta = 0.005 for species X and study B reports theta = 0.003 for species Y, how much of the difference is biological and how much is methodological?

When both studies use snpArcher with documented configurations, the methodological differences are minimized and clearly recorded.
This is especially valuable for large collaborative projects.
The California Conservation Genomics Project, for example, used snpArcher to process datasets from approximately 150 species through the same analytical framework, enabling meaningful cross-species comparisons of diversity, population structure, and demographic history.

## snpArcher design decisions through the non-model lens

To summarize, here is how each major snpArcher design decision relates to the non-model context:

| Design decision | Non-model rationale |
|----------------|-------------------|
| Hard filtering (not VQSR) | No truth sets exist for most species |
| Configurable het_prior (default 0.005) | GATK default (0.001) is wrong for most non-humans |
| Scatter-by-Ns parallelization | Fragmented assemblies have many gaps; chromosome-level splitting does not work |
| Callable sites module | Collapsed duplications and unmappable regions are common in draft assemblies |
| Multiple input types | Non-model projects routinely combine SRA, local FASTQ, and pre-existing BAMs |
| Single reference per run | Avoids subtle errors from mixing coordinate systems |
| QC dashboard | No established QC baselines exist; researchers need to assess each dataset on its own terms |
| Postprocess module | Filtering choices depend on species biology and data quality, not a one-size-fits-all recipe |
| Workflow profiles | The same analysis may run on a laptop, a university cluster, or a cloud instance |

Each of these decisions represents a trade-off.
Hard filtering is less precise than VQSR when a truth set is available.
A default het_prior of 0.005 is wrong for some species.
Scatter-by-Ns does not help with telomere-to-telomere assemblies.
But the combination of these defaults covers the vast majority of non-model population genomics projects out of the box.

## The big picture

Working with non-model organisms means working with less information.
You do not know the heterozygosity rate before you measure it.
You do not have a truth set to validate against.
Your reference genome has gaps and errors you have not characterized yet.
Your sample size is limited by specimen availability, not sequencing budget.

snpArcher is built for this situation.
Rather than assuming you know the answer before you start (as VQSR does), it provides:

- Defaults that err on the side of caution
- The QC dashboard and SFS analysis for assessing your specific dataset
- Every parameter exposed for adjustment as you learn more about your organism
- Full transparency: every step is a readable Snakemake rule with explicit shell commands

The goal is to handle the engineering (dependency management, parallelization, file routing, resource allocation) so that you can focus on the biology.

## Further reading

- [Pipeline architecture](architecture.md): How the overall design accommodates diverse organisms and assemblies.
- [Variant calling](variant-calling.md): Why het_prior matters and how to choose a caller.
- [Parallelization](parallelization.md): The scatter-by-Ns strategy for fragmented assemblies.
- [Filtering philosophy](filtering.md): SFS-guided filtering and callable sites for non-model genomes.
- [QC metrics explained](qc-metrics.md): Interpreting QC visualizations without species-specific baselines.

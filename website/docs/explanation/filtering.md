# Filtering philosophy

Variant calling does not end with a VCF.
The raw output of any caller (GATK, bcftools, DeepVariant) contains a mixture of real variants and artifacts.
Mapping errors, paralogous sequences collapsed in the assembly, low-complexity regions, and stochastic genotyping failures all contribute false or misleading calls.

This page explains snpArcher's approach to filtering, the reasoning behind each filter type, and how to use the site frequency spectrum (SFS) as a diagnostic tool.

## Why filtering matters

A raw VCF from a 30-sample dataset against a scaffold-level reference contains a mixture of:

- True SNPs from population-level variation
- False heterozygous calls from duplicated regions collapsed in the reference assembly
- Low-quality calls in regions of extreme depth (too high or too low)
- Genotyping errors at sites where only one or two reads support the variant allele
- Calls in unmappable regions where reads from multiple genomic locations pile up
- Artifacts from systematic biases in library preparation or sequencing

These artifacts do not affect all analyses equally.
A simple population assignment based on thousands of SNPs may be robust to moderate error rates.
But analyses that depend on allele frequencies (nucleotide diversity estimates, demographic inference from the SFS, selection scans) can be severely distorted by even a small proportion of artifactual calls.

## The two levels of filtering

snpArcher applies filtering at two levels, corresponding to two stages of the pipeline:

### 1. Hard filtering (core pipeline)

The core pipeline applies GATK hard filters immediately after joint genotyping.
These filters tag individual variant sites based on quality annotations:

| Filter | Metric | What it detects |
|--------|--------|-----------------|
| QD < 2.0 | Quality by Depth | Low-confidence calls driven by few reads |
| FS > 60.0 | Fisher Strand bias | Variant supported by reads from one strand only |
| SOR > 3.0 | Strand Odds Ratio | Strand bias (more robust at high depth) |
| MQ < 40.0 | Mapping Quality | Reads that do not map uniquely |
| MQRankSum < -12.5 | Mapping Quality Rank Sum | Reference allele reads map better than variant allele reads |
| ReadPosRankSum < -8.0 | Read Position Rank Sum | Variant-supporting reads concentrated at read ends |

These are the GATK-recommended thresholds for hard filtering when VQSR is not available.
They remove the most egregious artifacts while retaining most real variants.

!!! note
    Hard filtering tags sites in the FILTER column but does not remove them from the VCF.
    Downstream tools can choose whether to include or exclude filtered sites.

### 2. Postprocessing (postprocess module)

The postprocess module applies additional filters after you have reviewed the QC dashboard and made decisions about sample inclusion.
These filters operate on both samples and sites:

**Sample-level filtering:**

- Exclude specific samples based on QC review (contamination, mislabeling, low coverage).
  Samples are marked in the sample metadata with an `exclude` designation.

**Site-level filtering:**

- **Minor allele frequency (MAF)**: Remove sites below a frequency threshold (default 0.01).
- **Missingness**: Remove sites where more than a fraction of samples have missing genotypes (default 0.75 = at most 75% of genotypes called).
- **Scaffold/contig exclusion**: Remove sites on specific scaffolds (e.g., sex chromosomes, mitochondria, very small contigs).
- **Minimum contig size**: Remove sites on scaffolds below a length threshold (default 10,000 bp).
- **Callable sites intersection**: Restrict to regions that pass coverage and mappability filters.

## The site frequency spectrum as a diagnostic

The site frequency spectrum (SFS) is the distribution of allele frequencies across variant sites.
For a sample of diploid individuals, the folded SFS counts how many SNPs have a minor allele observed in 1 copy, 2 copies, 3 copies, and so on, up to half the total number of haploid genomes.

The SFS is useful for diagnosing filtering problems because different types of artifacts leave characteristic signatures:

### Expected SFS shape

Under standard neutral models the folded SFS declines roughly as 1/frequency.
Real populations deviate from this due to demography and selection, but those deviations are typically smooth.

### Artifact signatures

Artifacts, by contrast, create sharp peaks or unusual features in the SFS:

- **Peak at frequency 1/2 of haploid sample size**: A hallmark of collapsed duplications.
  When two paralogous regions are collapsed into one in the reference assembly, reads from both copies map to the same location.
  Fixed differences between the paralogs appear as heterozygous sites in every sample, producing a spike at allele frequency 0.5 (or N/2 in terms of allele counts).

- **Excess of singletons**: While rare variants are expected, an extreme excess of singletons (variants observed in only one individual) can indicate genotyping errors, especially in low-coverage samples.

- **Peak at frequency 1/4 of haploid sample size**: In species with heterogametic sex determination (ZW birds, XY mammals), hemizygous sex chromosome variants appear at characteristic frequencies.
  If sex chromosomes are not excluded, they create a peak in the SFS.

### Using the SFS to evaluate filters

Compute the SFS after each filtering step and observe how the distribution changes (the protocol paper demonstrates this with a burrowing owl dataset):

1. **Raw SFS**: Often shows peaks at 1/2 and possibly 1/4 frequency, plus an excess of singletons.
2. **After removing low-coverage samples**: The singleton excess typically decreases because the worst-quality genotype calls are removed.
   The total number of segregating sites may paradoxically increase as the removal of high-missingness samples reveals variants that were previously uncalled.
3. **After excluding sex chromosomes**: The 1/4-frequency peak should diminish or disappear.
4. **After callable sites filtering** (removing high-heterozygosity and low-coverage regions): The 1/2-frequency peak should diminish or disappear, indicating that collapsed duplications have been masked.

If the SFS still shows artifactual peaks after all filters, the filters may need to be more aggressive.
If the SFS looks smooth and follows the expected declining shape (possibly with demographic deviations), the filtering is likely adequate.

## Coverage filtering: the callable sites approach

Coverage-based filtering addresses two classes of problems:

### Low coverage

Regions with very low read depth produce unreliable genotype calls.
A site covered by only 1-2 reads has a high probability of missing the heterozygous allele entirely, leading to a false homozygous reference call.
Systematically, low-coverage regions create missing data and underestimate heterozygosity.

### Excessive coverage

Regions with much higher depth than expected are suspicious.
Common causes include:

- **Collapsed duplications**: Two or more paralogous regions in the true genome are represented by a single sequence in the reference assembly, so reads from all copies map to one location.
- **Repetitive elements**: Highly similar repeats can attract reads from across the genome.
- **Copy number variants**: Regions that are duplicated in some individuals.

These high-depth regions produce an excess of heterozygous calls (from fixed differences between paralogs) and distort allele frequency estimates.

### How snpArcher handles coverage

The callable sites module computes per-site depth across all BAM-backed samples and identifies regions that fall within acceptable bounds.
By default, the bounds are:

- **Minimum coverage**: `auto` = `max(1, floor(global_mean_coverage / 2))`
- **Maximum coverage**: `auto` = `ceil(global_mean_coverage * 2)`

The `fraction` parameter (default 1.0) controls what proportion of BAM-backed samples must be callable at a site for it to be included.
A value of 1.0 means every sample must be within the depth bounds; lower values allow some samples to be outside bounds.

Regions passing the coverage filter are merged (with a configurable merge distance, default 100 bp) into contiguous intervals, producing a BED file of coverage-callable regions.

## Mappability filtering

Regions with low mappability produce false variants because reads from multiple genomic locations pile up at the same reference position.

snpArcher uses GenMap to compute mappability scores.
For each position in the genome, GenMap reports a score between 0 and 1, where 1.0 means a read of the specified length (default `kmer: 150`, matching typical Illumina read length) can only map to that single location.

The default `min_score: 1` means only perfectly uniquely mappable regions are retained.
This is conservative but appropriate for most analyses.
Regions with ambiguous mappability are masked.

## The callable sites BED: putting it together

The final callable sites BED file is the intersection of the coverage-passing and mappability-passing regions.
It represents the portion of the genome where:

1. Read depth is within acceptable bounds across samples (not too low, not too high)
2. Reads can be uniquely mapped

This BED file serves two purposes:

- **Site filtering**: The postprocess module intersects the VCF with the callable sites BED, retaining only variants in callable regions.
- **Denominator for rate calculations**: Any analysis that computes a rate (variants per site, nucleotide diversity, etc.) needs to know the denominator: how many sites were examined.
  The callable sites BED defines that denominator.

!!! warning "Effective genome length"
    When you filter your VCF to callable sites, the effective length of the genome changes.
    If you started with a 1 Gb genome and 800 Mb passes callable sites filters, your denominator for diversity calculations is 800 Mb, not 1 Gb.
    Forgetting to adjust the denominator leads to underestimates of nucleotide diversity and biases in demographic inference.

## Sample filtering: QC-guided exclusion

The QC dashboard (see [QC metrics explained](qc-metrics.md)) helps identify problematic samples.
Common reasons to exclude a sample:

- **Very low coverage**: High missingness, unreliable genotype calls.
  A good threshold depends on your dataset, but samples below ~5x typically contribute more noise than signal.
- **Contamination**: Strongly negative inbreeding coefficient (F << 0), excess heterozygosity genome-wide.
- **Mislabeling**: Sample clusters with the wrong population in PCA, or has unexpected kinship.
- **Close relatives**: Kinship coefficient > 0.354 with another sample.
  Including close relatives biases estimates of effective population size and nucleotide diversity.
- **Extreme depth outliers**: Samples with much higher or lower depth than the cohort may drive batch-effect-like artifacts.

The protocol paper demonstrates this with the burrowing owl dataset: removing 18 samples with average depth below 8x (out of 137 total) resulted in a tenfold increase in the number of segregating sites, because the high-missingness samples had been suppressing variant detection.

## MAF and missingness: trade-offs

### Minor allele frequency

Filtering by MAF removes rare variants.
The default threshold (0.01) removes singletons and very rare variants, which are enriched for genotyping errors.

**Trade-offs:**

- **Stringent MAF filter (e.g., 0.05)**: Removes more noise, but also removes real rare variants.
  Appropriate for analyses that focus on common variation (PCA, GWAS).
- **Lenient MAF filter (e.g., 0.01 or none)**: Retains more of the frequency spectrum.
  Appropriate for demographic inference from the SFS, where rare variants carry the most information.

### Missingness

The missingness filter removes sites where too many samples lack a genotype call.
The default threshold (0.75) means at least 75% of samples must have a called genotype.

**Trade-offs:**

- **Stringent missingness filter (e.g., 0.95)**: Retains only well-genotyped sites, but may drastically reduce the number of variants, especially if depth varies across samples.
- **Lenient missingness filter (e.g., 0.50)**: Retains more sites, but some will have genotypes from only half the samples.
  This can bias allele frequency estimates.

The right balance depends on your analysis.
For the SFS, you generally want low missingness (high completeness).
For PCA, moderate missingness is tolerable because the dimensionality reduction is robust to missing data.

## Scaffold and contig exclusion

Several types of genomic regions are routinely excluded from population genomic analyses:

### Sex chromosomes

In species with heterogametic sex determination, the sex chromosomes (Z/W in birds, X/Y in mammals) have different ploidy in males and females.
This creates systematic genotyping differences: hemizygous regions in the heterogametic sex appear as homozygous, with half the expected depth.
Mixing autosomal and sex-linked variants in a single analysis confounds allele frequency estimation and can create artificial peaks in the SFS.

The default `exclude_scaffolds` setting in the postprocess module is `"mtDNA,Y"`.
For birds, you should add the Z chromosome scaffolds.
For other organisms, adjust accordingly.

### Mitochondrial genome

The mitochondrial genome is haploid, maternally inherited, and present in hundreds of copies per cell (leading to extremely high depth).
It violates most of the assumptions of diploid population genomic methods.
Exclude it.

### Small scaffolds

Very short scaffolds (below the `contig_size` threshold, default 10,000 bp) are often unplaced fragments of uncertain genomic origin.
They may represent poorly assembled regions, contamination, or organellar sequences.
They contribute few variants and disproportionate noise.

## Summary principles

Start with stringent filters and relax if needed.
You can always re-run the postprocess module with different thresholds; you cannot un-propagate artifacts.

Use the SFS as your guide: if the filtered SFS is smooth with no sharp peaks at characteristic frequencies, your filtering is probably adequate.
If peaks persist, investigate.

Every filter that removes genomic regions changes the effective genome length.
Any rate statistic (diversity, divergence, recombination rate) must be calculated relative to the filtered genome length, not the full assembly size.
Report site and sample counts before and after each filtering step so that readers can evaluate the impact of your filtering choices.

## A worked example: burrowing owl

The protocol paper illustrates filtering with a dataset of 137 burrowing owl (*Athene cunicularia*) genomes mapped to a scaffold-level reference.
The key observations:

1. **Raw VCF**: SFS shows clear peaks at 1/2 and 1/4 of the folded frequency spectrum.
2. **After removing 18 low-coverage samples** (depth < 8x): Total segregating sites increase tenfold (high-missingness samples were suppressing variant counts).
   The SFS peaks persist.
3. **After excluding Z chromosome scaffolds**: The 1/4-frequency peak diminishes.
   This peak was caused by hemizygous Z-linked variants in females.
4. **After filtering high-heterozygosity regions**: Both 1/2 and 1/4 peaks are removed.
   The filtered SFS is smooth and consistent with neutral demographic expectations.

Demographic inference (using blockbuster) on the raw SFS produced spurious bottleneck signals that disappeared after proper filtering.
The final filtered SFS supported a scenario of population size increase, consistent with the species' conservation status ("Least Concern").

Inadequate filtering can lead to qualitatively wrong biological conclusions.

## Further reading

- [Postprocessing how-to](../how-to/postprocess.md): Step-by-step guide to running the postprocess module.
- [QC metrics explained](qc-metrics.md): How to interpret the QC dashboard to inform filtering decisions.
- [Non-model organisms](non-model.md): Why filtering is especially challenging (and important) for non-model species.
- [Configuration reference](../reference/config-schema.md): Full specification of callable sites and postprocess parameters.

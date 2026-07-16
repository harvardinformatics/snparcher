# QC metrics explained

snpArcher's QC module produces an interactive HTML dashboard with ten visualizations designed to surface common problems in population genomic datasets.
This page explains what each metric captures, how to interpret it, and what red flags to watch for.

For step-by-step instructions on generating and navigating the dashboard, see the [QC report how-to](../how-to/qc-report.md).

## What the QC report is (and is not)

The QC report is a diagnostic tool for identifying problems in your dataset before committing to downstream analyses.
It is based on a random, LD-pruned subset of approximately 100,000 biallelic SNPs, after removing indels, sites with MAF < 0.01, sites with >75% missing data, and samples below the minimum depth threshold.

This subset is sufficient to capture the major patterns in your data (population structure, batch effects, relatedness, contamination) without the computational cost of analyzing millions of variants.
However, it is not a final analysis.
The exact values of PCA loadings, admixture proportions, or nucleotide diversity reported in the QC dashboard should not be cited in publications.
They are approximations meant to guide decisions about sample inclusion, filtering thresholds, and whether the data are suitable for your intended analyses.

!!! note
    The QC module requires at least 3 samples.
    With fewer samples, most of the metrics (PCA, relatedness, admixture) are not meaningful.

## The header: summary statistics

The top of the report displays four summary values:

- **Number of individuals**: How many samples passed the minimum depth filter.
- **Number of SNPs**: Total biallelic SNPs in the filtered VCF (before the 100K subsample).
- **Nucleotide diversity** (Watterson's theta): An estimate of population-level diversity, calculated from the number of segregating sites scaled by sample size and genome length.
  This is a rough estimate; different filtering choices or different estimators will give different values.
- **Mean depth**: Average genome-wide sequencing depth across all samples, based on read depth from mosdepth.
  This differs from "variant depth" (depth at polymorphic sites only), which does not account for monomorphic positions.

## Figure 1: Principal Component Analysis (PCA)

**What it shows:**
A scatterplot of the first two principal components from a PCA of the genotype matrix, with samples colored by k-means cluster assignments.
PCA is calculated using PLINK2 on the 100K SNP subset.

**Limitations:**
PCA is a visualization of genetic distance, not a model of population history.
Two populations can overlap in PCA space despite having distinct evolutionary histories, and a single continuous population can appear as multiple clusters if sampling is uneven.

**Interpreting clusters:**
The k-means coloring (default k=3) is a visual convenience, not an assertion about population number.
Set `clusters` in the configuration to a value that matches your expected population structure, or just use it to track samples across panels.

**What to look for:**

- **Clear clusters separated along PC1 or PC2**: Likely real population structure.
  Check whether the separation corresponds to geography, collection site, or known biological groupings.
- **A single sample far from all others**: Possible contamination, mislabeling, or a sample from a different species.
- **Separation that corresponds to sequencing batch or library prep**: Batch effects.
  See the depth-PC correlation figure below.

## Figure 2: Depth vs. principal components

**What it shows:**
Sequencing depth plotted against each of the first 10 principal components.

**Why it matters:**
A common artifact in moderate-coverage datasets is that sequencing depth variation among samples creates apparent population structure.
Low-coverage samples have more missing genotypes and more genotyping errors (particularly, fewer heterozygous calls) than high-coverage samples.
These systematic differences can be captured by PCA, creating clusters that reflect sequencing effort rather than biology.

**What to look for:**

- **Strong correlation between depth and a top PC**: This indicates that depth variation is a major source of apparent genetic variation in your dataset.
  If PC1 correlates with depth, the dominant pattern in your PCA is likely technical, not biological.
- **No correlation**: Depth variation is not driving the PCA. No action needed.
- **Correlation with a lower PC (e.g., PC5) but not PC1**: The depth effect exists but is subordinate to real biological structure.
  This is common and usually acceptable for most analyses.

!!! tip
    If depth correlates strongly with PC1, consider whether additional depth filtering (removing the lowest-coverage samples) might resolve the issue.
    The [filtering page](filtering.md) discusses how to make these decisions.

## Figure 3: Depth vs. missingness

**What it shows:**
For each sample, mean SNP depth is plotted against SNP missingness (the fraction of genotype calls that are missing).
Samples are colored by PCA cluster.

**What to look for:**

- **Expected pattern**: Higher depth = lower missingness.
  Samples with very low depth will have high missingness.
- **Cluster-specific patterns**: If one PCA cluster has systematically higher missingness at the same depth as another cluster, this may indicate reference bias (samples more divergent from the reference genome have more unmapped reads and more missing genotypes).
- **Outliers**: A sample with unusually high missingness for its depth level may have data quality issues (degraded DNA, contamination, wrong species).
- **Threshold decisions**: This plot helps you decide on a minimum depth cutoff for sample inclusion.
  You can visually identify the depth below which missingness increases sharply.

## Figure 4: Mapping rate vs. depth

**What it shows:**
The proportion of reads that successfully align to the reference genome, plotted against sequencing depth.

**What to look for:**

- **Mapping rate below ~80%**: Suggests that a substantial fraction of reads come from a source other than the target genome.
  Common causes: bacterial contamination, sample from a different (or very divergent) species, poor reference genome quality.
- **Cluster of samples with low mapping rate**: If several samples share low mapping rates, they may be from a divergent population or subspecies.
  Cross-reference with PCA clusters.
- **Single sample with very low mapping rate**: Likely a mislabeled sample or a sample from the wrong species.
- **Mapping rate variation with no depth dependence**: Expected, since mapping rate depends on genetic distance from the reference and DNA quality, not sequencing depth.

## Figure 5: Inbreeding coefficient (F)

**What it shows:**
Per-sample inbreeding coefficient, estimated by PLINK.
F measures the departure from Hardy-Weinberg expected heterozygosity.

**How to interpret F values:**

| F value | Interpretation |
|---------|---------------|
| Near 0 | Consistent with Hardy-Weinberg expectations (random mating) |
| Positive (toward +1) | Excess homozygosity: inbreeding, population subdivision, or low coverage underestimating heterozygosity |
| Negative (toward -1) | Excess heterozygosity: contamination (two individuals mixed into one library), outbreeding, or paralogous mapping |

**What to look for:**

- **Strongly negative F (e.g., below -0.1)**: A classic signature of contamination.
  When DNA from two individuals is mixed in one library, the resulting reads appear to come from a single very heterozygous individual.
  Check these samples carefully.
- **Outliers with strongly positive F**: Could be truly inbred (island populations, captive individuals) or could indicate very low coverage causing heterozygous sites to be missed.
- **F varies systematically by PCA cluster**: Expected if there is real population structure.
  Subdivision inflates F within populations and creates differences between them.

!!! warning
    F is sensitive to population structure.
    If your dataset contains multiple distinct populations, the overall F values will be shifted upward for all samples because the calculation assumes a single panmictic population.
    Interpret F in the context of the PCA clusters.

## Figure 6: Neighbor-joining tree

**What it shows:**
An unrooted neighbor-joining (NJ) tree constructed from a pairwise genetic distance matrix across all samples.
Terminal branches are colored by PCA cluster assignment.

**What it captures:**
The NJ tree preserves all pairwise distances rather than projecting onto two axes, so it can reveal relationships that PCA misses.

**What to look for:**

- **Concordance with PCA clusters**: Samples from the same PCA cluster should generally group together on the tree.
  Discordance may indicate that PCA is collapsing structure that the tree preserves, or vice versa.
- **Long terminal branches**: A sample on a very long branch is genetically distinct from all others.
  This could be real (a divergent individual) or artifactual (contamination, mislabeling).
- **Samples that cluster on the tree but not in PCA**: The tree uses all pairwise distances, not just the top two PCs, so it may reveal structure that PCA misses.

## Figure 7: Relatedness heatmap (KING kinship)

**What it shows:**
A heatmap of the KING-robust kinship coefficient matrix, calculated by PLINK.
Each cell shows the kinship between a pair of samples.

**How to interpret KING kinship values:**

| Kinship | Relationship |
|---------|-------------|
| ~0.5 | Identical samples or monozygotic twins |
| ~0.25 | First-degree relatives (parent-offspring, full siblings) |
| ~0.125 | Second-degree relatives (half-siblings, grandparent-grandchild, avuncular) |
| ~0.0625 | First cousins |
| ~0 | Unrelated |

**What to look for:**

- **Kinship near 0.5**: Two samples that are effectively identical.
  This usually means a duplicate sample (same individual sequenced twice, or a labeling error).
  One should be removed.
- **Kinship above 0.354**: The standard threshold for "close relative."
  Including close relatives can bias population genomic estimates (effective population size, nucleotide diversity, selection scans).
  Consider removing one individual from each close-relative pair.
- **Block structure**: Groups of samples with uniformly high kinship may be from a family group or a highly inbred population.
- **Unexpected patterns**: If two samples from distant collection sites have high kinship, suspect mislabeling.

## Figures 8 and 9: Geographic maps

**What they show:**
If latitude and longitude coordinates are provided in the sample metadata, Figure 8 plots sample locations colored by PCA cluster.
If a Google Maps API key is also provided, Figure 9 renders the same points on a satellite basemap.

**What to look for:**

- **Geographic concordance with PCA clusters**: Populations often correspond to geographic regions.
  If PCA clusters align with geography, that supports a biological interpretation.
- **Isolation by distance**: A gradual shift in cluster composition across space suggests continuous population structure rather than discrete populations.
- **Samples in unexpected locations**: A collection coordinate error, or an indication that the sample was mislabeled.

## Figure 10: ADMIXTURE

**What it shows:**
ADMIXTURE barplots for K=2 and K=3 ancestral populations.
Each bar represents a sample, and the colors indicate the estimated proportion of ancestry from each hypothetical ancestral population.

**What ADMIXTURE captures:**
ADMIXTURE explicitly estimates ancestry proportions under a model of K discrete ancestral populations, which makes it useful for identifying admixed individuals.

**Important caveats:**

- **K=2 and K=3 are arbitrary.**
  The QC report runs ADMIXTURE at two fixed K values for quick visualization, with no cross-validation for optimal K.
  Do not cite these results as evidence for a specific number of ancestral populations.
- **ADMIXTURE coloring is independent of PCA coloring.**
  The colors in the ADMIXTURE plot represent ancestry components, not PCA clusters.
  Do not assume that "blue" in ADMIXTURE corresponds to "blue" in PCA.
- **ADMIXTURE can be misleading for continuous populations.**
  When population structure is clinal (gradual change across space), ADMIXTURE will still assign individuals to discrete ancestral populations, producing a barplot that suggests abrupt boundaries where none exist.

**What to look for:**

- **Clear separation at K=2**: A major division in your dataset, which should correspond to a pattern visible in PCA.
- **Admixed individuals (mixed colors)**: Samples with intermediate ancestry proportions, suggesting gene flow between populations.
- **Discordance with PCA**: If ADMIXTURE groups samples differently than PCA clusters, this may reveal aspects of population structure that PCA misses (or vice versa).

## Why 100,000 SNPs?

The QC report uses a random, LD-pruned subset of 100,000 SNPs rather than the full variant set.
This choice balances three considerations:

1. **Speed**: PCA, ADMIXTURE, and distance matrix calculations on 100K SNPs take seconds to minutes.
   On millions of SNPs, they could take hours.
2. **LD-pruning**: Population structure analyses assume approximately independent markers.
   Using all SNPs without LD-pruning would give disproportionate weight to regions of the genome in strong linkage disequilibrium, distorting PCA and admixture estimates.
3. **Signal retention**: 100,000 LD-pruned SNPs are more than sufficient to resolve population structure, detect relatedness, and identify outlier samples in datasets typical of non-model population genomics (10-200 samples).
   Benchmark comparisons show that PCA results stabilize well below this number.

The LD-pruning is done by randomly retaining one SNP per genomic window, where the window size is calculated to yield approximately 100,000 SNPs for the given genome size.

## Common QC red flags

Here is a summary of the most common problems the QC report can detect, and how they manifest:

### Batch effects

PCA clusters correspond to sequencing lane, library prep date, or extraction batch rather than biological populations, and the depth-PC correlation is strong.
The underlying cause is systematic differences in sequencing quality or depth across batches.

If batches correspond perfectly to populations, the effect cannot be separated from biology, and there is no clean fix.
If batches are distributed across populations, the effect is usually minor.
In either case, document the batch structure and consider batch-aware methods (e.g., including batch as a covariate in downstream models).

### Contamination

Samples with strongly negative inbreeding coefficients (F well below 0) and possibly lower mapping rates are likely contaminated: DNA from two or more individuals mixed into a single library, creating apparent heterozygosity genome-wide.
These samples typically appear as outliers in PCA.
Remove them.
If contamination is widespread (more than one or two samples), the problem is likely in the library preparation protocol.

### Sample mislabeling

A sample that clusters with an unexpected population in PCA, has high kinship with samples it should be unrelated to, or whose geographic coordinates do not match its genetic placement is probably mislabeled.
Cross-reference with collection metadata; if the true identity can be determined, relabel the sample, otherwise remove it.

### Coverage outliers

One or more samples have much higher or lower depth than the rest, driving a correlation between depth and PC1.
This is common when samples are pooled with unequal DNA concentrations.
Remove extreme low-depth outliers (or set a minimum depth threshold in the postprocess module), and consider downsampling extreme high-depth samples.
The [filtering page](filtering.md) discusses depth-based sample exclusion in detail.

### Duplicate samples

Two samples with kinship coefficient near 0.5 that overlap perfectly in PCA are the same individual sequenced twice (or a labeling error).
Remove one.

### Reference bias

One PCA cluster has systematically higher missingness than another at comparable depth levels, because reads from samples more divergent from the reference fail to map.
This is inherent to reference-based mapping.
If severe, consider using a reference genome more closely related to the study populations.
See [non-model organisms](non-model.md) for more discussion.

## Further reading

- [QC report how-to](../how-to/qc-report.md): Step-by-step instructions for generating and navigating the dashboard.
- [Filtering philosophy](filtering.md): How QC findings inform site-level and sample-level filtering decisions.
- [Non-model organisms](non-model.md): Why QC is especially important for non-model datasets.

# How to interpret QC reports

The snpArcher QC module produces an interactive HTML dashboard that summarizes individual-level and population-level quality metrics.

For deeper discussion of what each metric means biologically, see [QC metrics explained](../explanation/qc-metrics.md).

## Open the report

After a run with `modules.qc.enabled: true`, find the dashboard at:

```text
results/qc/qc_dashboard.html
```

Open it in a web browser.
All plots are interactive (built with Plotly): zoom, pan, and hover over points to see sample IDs.

## Report header

The top of the report shows a summary of:

- **Number of individuals** processed
- **Number of SNPs** detected
- **Nucleotide diversity** (Watterson's theta estimate)
- **Average genome-wide depth** (based on read depth)

!!! note
    The QC report is based on a random, LD-pruned subset of approximately 100,000 biallelic SNPs, after removing indels, sites with MAF < 0.01, sites with > 75% missing data, and samples below the `min_depth` threshold.
    This subset provides a good approximation of dataset-wide patterns but is not a substitute for final analysis.

## Panels and what to look for

### PCA (Figure 1)

A principal component analysis colored by k-means cluster assignments (default k = 3).

**What to look for:**

- Distinct clusters may reflect real population structure or batch effects.
- Outlier samples far from any cluster warrant investigation: they could indicate contamination, sample swaps, or a different species.
- The number of clusters (`modules.qc.clusters` in config) can be adjusted if the default does not match your expected population structure.

### Depth vs. principal components (Figure 2)

Sequencing depth plotted against PCs 1 through 10.

**What to look for:**

- A strong correlation between depth and a top PC suggests that depth variation is driving apparent structure.
- Low-to-moderate coverage datasets are especially susceptible, since low-depth samples have reduced heterozygote detection.

**Red flag:** if depth explains a large fraction of variance in PC1, population structure results should be interpreted cautiously.

### Depth vs. missingness (Figure 3)

Mean SNP depth plotted against SNP missingness for each individual.

**What to look for:**

- The overall distribution of depth and missingness across samples.
- Whether one PCA cluster shows systematically higher missingness at a given depth, which could indicate reference bias.
- How many samples would be excluded under different depth thresholds (useful for deciding on a depth cutoff).

### Mapping rate vs. depth (Figure 4)

Proportion of reads that successfully aligned to the reference, plotted against depth.

**What to look for:**

- Most samples should cluster at high mapping rates (> 90%).
- Samples with mapping rates below 80% may indicate contamination, sample misidentification, or a divergent species.

**Red flag:** a sample with very low mapping rate but typical depth could contain substantial non-target DNA (e.g., bacterial contamination).

### Inbreeding coefficient (Figure 5)

Per-individual estimates of the inbreeding coefficient (F).

**What to look for:**

- Values near 0 indicate Hardy-Weinberg proportions.
- Positive values indicate excess homozygosity (inbreeding or population structure).
- Strongly negative values (excess heterozygosity) can indicate cross-contamination between samples.

**Red flag:** samples that are PCA outliers with strongly negative F values are likely contaminated.

### Neighbor-joining tree (Figure 6)

An unrooted tree from pairwise genetic distances, with branches colored by PCA cluster.

**What to look for:**

- Concordance between tree structure and PCA clusters.
- Very long terminal branches, which may indicate highly divergent or problematic samples.

### Relatedness heatmap (Figure 7)

A heatmap of KING kinship coefficients between all pairs of samples.

**What to look for:**

| Kinship coefficient | Relationship |
|---|---|
| ~0.5 | Identical / duplicate |
| ~0.25 | First-degree relatives (parent-offspring, full siblings) |
| ~0.125 | Second-degree relatives (half-siblings, grandparent-grandchild) |
| < 0.0625 | Unrelated |

**Red flag:** kinship coefficients near 0.5 suggest duplicate samples or monozygotic twins that should be investigated.
Consider removing samples with kinship > 0.354 to avoid biasing population genetic estimates.

### Geographic map (Figures 8 and 9)

If you provided latitude and longitude coordinates in the sample metadata file, the report includes:

- **Figure 8:** sample locations colored by PCA cluster on a simple map.
- **Figure 9:** the same locations on a satellite basemap (requires a Google Maps API key in config).

**What to look for:**

- Whether PCA clusters correspond to geographic regions (expected for real population structure).
- Geographically isolated samples that cluster unexpectedly, which may indicate sample labeling errors.

### ADMIXTURE (Figure 10)

Model-based ancestry estimates from ADMIXTURE for k = 2 and k = 3.

**What to look for:**

- General concordance with PCA clusters.
- Samples with mixed ancestry that might represent admixed populations or zones of secondary contact.

!!! note
    The k values (2 and 3) are chosen for exploratory purposes.
    No cross-validation is performed.
    Use formal ADMIXTURE analysis with cross-validation for any publication.

## Deciding which samples to exclude

After reviewing the QC report, exclude samples that show:

- Very low sequencing depth or high missingness
- Evidence of contamination (low mapping rate, negative inbreeding coefficient)
- Duplicate samples (high kinship)
- Sample misidentification (wrong species or unexpected placement in PCA)

To exclude samples from downstream analysis, see [How to filter and postprocess](postprocess.md).

## Next steps

- [Filter and postprocess](postprocess.md) to remove problematic samples and apply site-level filters
- [QC metrics explained](../explanation/qc-metrics.md) for deeper background on each metric

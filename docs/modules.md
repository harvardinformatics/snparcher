# Modules
A key goal in the design of the snpArcher pipeline is to allow seamless extensibility with downstream processing. We implement this using Snakemake modules, which allow additional rules to easily extend the main pipeline. We present several modular extensions of snpArcher here, but we hope also that user-developed modules will grow the set of tools linked to snpArcher in order to facilitate diverse analysis.
## Module Contribution Guidelines
We developed a set of criteria for including additional user-developed modules into snpArcher. This project is designed to be modular and easily extensible as we and workflow users develop additional features and downstream analysis tools. To ensure that contributed modules are reproducible and easily implemented, we propose the following evaluation criteria:

1. Each module must include Snakemake workflow that defines necessary environments using Conda. 
2. The module must be freely distributed via Github with sufficient documentation that users can adapt it to their needs
3. The module must provide a unit test based on either existing test datasets available from the main snpArcher workflow or via a module-specific minimal test dataset
4. Each module should be registered within the main project page to enhance discoverability and ensure the above criteria are met.

If you are interested in developing a module please reach out via email or Github, we'd love to know and chat about it. 
## Quality Control
The quality control module aggregates various statistics from the workflow and produces preliminary analyses and plots in an interactive HTML file, offering visualizations of summary statistics related to population structure, batch effects, sequencing depth, genetic relatedness, geography, and admixture. Most summaries are based on a random sample of 100,000 SNPs, while others provide high-level summaries of the full variant dataset. These visualizations help identify outliers, potential biases, and sequencing artifacts, such as cryptic genetic variation, batch effects, and reference bias. Additionally, an interactive heatmap aids in quickly identifying close relatives within the dataset, and spatial maps provide a visualization of PCA clusters in space.
### Config Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`modules.qc.clusters`| Number of clusters for PCA visualization.| `int`|
|`modules.qc.google_api_key`| Google Maps API key for the terrain panel (optional).| `str`|
|`modules.qc.min_depth`| Samples with average depth below this will be excluded for QC analysis.| `int`|
|`modules.qc.exclude_scaffolds`| Comma-separated scaffolds to exclude from QC SNP sampling.| `str`|

```{note}
To generate the QC dashboard, you must have at least 3 samples specified in your sample sheet.
```
```{note}
The output of the QC module should not be considered a final analysis and is solely intended to direct quality control of the dataset.
```
## Postprocessing
The postprocessing module is designed to be run after the main workflow once you have decided whether any samples should be excluded from downstream analyses. To exclude samples, provide a `sample_metadata` file with an `exclude` column; samples with `exclude=true` are removed from the postprocessed outputs.

This module produces a filtered VCF by removing excluded samples, restricting to callable regions, excluding small contigs, and applying user-defined SNP/indel filters.
### Config Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`modules.postprocess.filtering.contig_size`| Variants on contigs this size or smaller are excluded from clean outputs.| `int`|
|`modules.postprocess.filtering.maf`| Variants with MAF below this value are excluded.| `float`|
|`modules.postprocess.filtering.missingness`| Variants with missingness above this value are excluded.| `float`|
|`modules.postprocess.filtering.exclude_scaffolds` | Comma-separated scaffolds/contigs to exclude from clean outputs.| `str`|

```{hint}
If you want to keep all samples in postprocess, omit the `exclude` column or set every row to `false`.
```

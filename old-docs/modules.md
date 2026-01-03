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
|`nClusters`| Number of clusters for PCA| `int`|
|`GoogleAPIKey`| Google Maps API key (optional).| `str`|
|`min_depth`| Samples with average depth below this will be excluded for QC analysis| `int`|

```{note}
To generate the QC dashboard, you must have at least 3 samples specified in your sample sheet.
```
```{note}
The output of the QC module should not be considered a final analysis and is solely intended to direct quality control of the dataset.
```
## Postprocessing
The postprocessing module is designed to be run after snpArcher has intially been run and you have determined if there are samples that you would like to exclude from downstream analyses. In order to trigger this module, you must add the `SampleType` column to your sample sheet, and mark samples for inclusion with the value `include` and exclusion with the value `exclude`. 

This module produces a filtered VCF by filtering excluded samples as well as sites not passing default and user defined filters.
### Config Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`contig_size`| SNPs on contigs this size or smaller will be excluded from 'clean' VCF | `int`|
|`maf`| SNPs with MAF below this will be excluded from clean VCF| `float`|
|`missingness`| SNPs with missingness below this will be excluded from clean VCF| `float`|
|`scaffolds_to_exclude` | Comma separated, no spaces list of scaffolds/contigs to exclude from clean VCF|

```{hint}
If you'd like to run the postprocessing module by default, you can add the `SampleType` column in your sample sheet, and mark all samples as `include`.
```
## Trackhubs
The trackhub module generates UCSC Genome Browser track files to explore population variation data from the VCF produced by snpArcher. This module computes and generates genome browser tracks for traditional population genomic summary statistics such as windowed estimates of Tajima’s D, SNP density, Pi, Minor Allele Frequency, SNP depth. To trigger this module, you must set the [config](./setup.md#core-configuration) option to `True` and supply a email (a requirement for tracks displayed on the UCSC Genome Browser).

```{warning}
The Trackhubs module is dependent on the postprocessing module.
```

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
|`modules.qc.max_sample_missingness`| Samples with >49% missing genotypes in the pruned QC SNP set are excluded before PLINK PCA/GRM.| `float`|
|`modules.qc.exclude_scaffolds`| Comma-separated scaffolds to exclude from QC SNP sampling.| `str`|

```{note}
To generate the QC dashboard, you must have at least 3 samples specified in your sample sheet.
```
```{note}
The output of the QC module should not be considered a final analysis and is solely intended to direct quality control of the dataset.
```
## Hard filtering (GATK-lineage callers)

For GATK-lineage callers (`gatk`, `sentieon`, `parabricks`), the main workflow annotates
the raw calls with the GATK hard-filter FILTER column (`gatk VariantFiltration`), writing
`results/vcfs/filtered.vcf.gz`. This file has the same records as `results/vcfs/raw.vcf.gz`
and differs only in the FILTER column — sites that fail a hard filter are *tagged*, not
removed, so downstream tools can choose to include or exclude them.

`bcftools` and `deepvariant` do not emit the annotations these filters use
(`QD`/`FS`/`SOR`/`MQ`/`MQRankSum`/`ReadPosRankSum`), so hard filtering is skipped for them
and their raw VCF is the final call set.

`variant_calling.generate_filtered_vcf` (default `true`) controls whether
`results/vcfs/filtered.vcf.gz` is produced as a default output alongside `raw.vcf.gz`
(reproducing the v1 raw + filtered output set). It only applies to GATK-lineage callers; for
`bcftools`/`deepvariant` it is ignored with a warning (they have no GATK hard-filter step).
Set it to `false` for a raw-only default. When `false`, the filtered VCF is still built on
demand when postprocess/QC run or you request the `call_variants` target.

## Postprocessing
The postprocessing module is designed to be run after the main workflow once you have decided whether any samples should be excluded from downstream analyses. To exclude samples, provide a `sample_metadata` file with an `exclude` column; samples with `exclude=true` are removed from the postprocessed outputs.

The module consumes the final call set (the hard-filtered `results/vcfs/filtered.vcf.gz` for
GATK-lineage callers, otherwise `results/vcfs/raw.vcf.gz`) and produces the strictly-filtered
`results/postprocess/filtered.vcf.gz` by removing excluded samples and hard-filter-failing
sites, restricting to callable regions, excluding small contigs, and applying MAF/missingness
filters. Note this is a different file from `results/vcfs/filtered.vcf.gz` (the main-workflow
hard-filtered VCF, all sites retained). By default it also emits SNP-only and indel-only
subsets (`clean_snps.vcf.gz`, `clean_indels.vcf.gz`).

For standalone runs against an existing VCF, you can use `run-postprocess.sh` as a thin wrapper around `workflow/modules/postprocess/Snakefile`.
### Config Options
| Option | Description | Type |
| ---- | -------------| ------ |
|`modules.postprocess.filtering.contig_size`| Variants on contigs this size or smaller are excluded from clean outputs.| `int`|
|`modules.postprocess.filtering.maf`| Variants with MAF below this value are excluded.| `float`|
|`modules.postprocess.filtering.missingness`| Variants with missingness above this value are excluded.| `float`|
|`modules.postprocess.filtering.exclude_scaffolds` | Comma-separated scaffolds/contigs to exclude from clean outputs.| `str`|
|`modules.postprocess.filtering.split_by_type`| Also emit `clean_snps.vcf.gz` and `clean_indels.vcf.gz` (the strict-filtered VCF split by variant type). Default `true`.| `bool`|
|`modules.postprocess.filtering.keep_basic_filter`| Retain the intermediate basic-filter VCF (`results/postprocess/basic.vcf.gz`) instead of discarding it. Default `false`.| `bool`|

```{hint}
If you want to keep all samples in postprocess, omit the `exclude` column or set every row to `false`.
```

# 🚀 snpArcher

<img src="./img/logo.png" alt="snpArcher logo" height="300"/>

Have resequencing data as fastq files and a reference genome? Want a VCF file of genotypes? Use snpRarcher as your one-stop shop to quickly and efficiently produce an analysis-ready dataset. No need to create a hand-tailored workflow cobbled together by tape and error-ridden chatGPT code — use snpArcher for all your variant calling needs on your laptop, on your server, or up in the clouds. 

snpArcher is a reproducible workflow optimized for nonmodel organisms and comparisons across datasets, built on the [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html#) workflow management system. It provides a streamlined approach to dataset acquisition, variant calling, quality control, and downstream analysis.

Snakemake makes it easy to bundle together the many steps involved in running a bioinformatics pipeline. The workflow involves mapping reads to a reference genome, calling SNPs using GATK's haplotypecaller, and calling variants at the population level using GATK's combineGVCFs. Each of these steps can be slow and tiresome to run on their own and the workflow has been carefully designed and tested to maximize efficiency. We use intervals to break up jobs into smaller chunks so that time and memory-hungry steps like haplotypecaller run quickly. If you have access to a Sentieon license for accelerated variant calling, we include options for this. 

Finally, the pipeline makes it easy to evaluate how the data looks. Review the HTML file in the QC folder at the end of a run to see how your samples relate to each other and also a number of metrics for evaluating variant-calling quality. 

Remember to examine the config.yaml file to edit options for each step. We have carefully chosen default options that should work for most users, but these can be tweaked here. 

## Requirements
- Illumina paired-end fastq files for one or more individuals
- A reference genome
- A sample sheet with sample names matched to the read names
- Snakemake and Mamba installed on your system
- If using Google Cloud, you will need to have set up an account on the GCP console

## Using snpArcher
- To get started quickly, check out the quick start tutorial!
- Otherwise start [here](./setup.md).

## Citing
- Please cite our preprint [here](https://www.biorxiv.org/content/10.1101/2023.06.22.546168v1)
- Also, make sure to cite the tools you used within snpArcher.

## Contributing to snpArcher
- If you encounter a bug or want to request a feature, please open a issue on our [github page](https://github.com/harvardinformatics/snpArcher).
- If you'd like to contribute a module, check out our [module contribution guidelines](./modules.md#module-contribution-guidelines).

```{toctree}
:hidden: True
./setup.md
./executing.md
./examples.md
./modules.md
./datasets.md
```

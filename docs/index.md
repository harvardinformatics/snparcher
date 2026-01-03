# snpArcher

A reproducible Snakemake workflow for population genomic variant calling, optimized for non-model organisms.

## What is snpArcher?

snpArcher streamlines the process of going from raw sequencing data (FASTQ files) to analysis-ready variant calls (VCF files). Built on Snakemake, it provides:

- **Reproducible variant calling** using GATK or Sentieon
- **Automated quality control** with interactive reports
- **Flexible input handling** for local files or SRA accessions
- **Scalable execution** from laptops to cloud clusters

## Quick Start

```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup
git clone https://github.com/harvardinformatics/snpArcher.git
cd snpArcher
pixi install

# Run test dataset
pixi run test-dry-run
```

## Documentation

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } **Tutorials**

    ---

    Step-by-step guides to get you started with snpArcher

    [:octicons-arrow-right-24: Get started](tutorials/index.md)

-   :material-tools:{ .lg .middle } **How-to Guides**

    ---

    Practical recipes for common tasks

    [:octicons-arrow-right-24: How-to guides](how-to/index.md)

-   :material-book-open-variant:{ .lg .middle } **Reference**

    ---

    Technical specifications and configuration options

    [:octicons-arrow-right-24: Reference](reference/index.md)

-   :material-lightbulb-on:{ .lg .middle } **Explanation**

    ---

    Background concepts and design decisions

    [:octicons-arrow-right-24: Explanation](explanation/index.md)

</div>

## Requirements

- Illumina paired-end FASTQ files
- A reference genome (local FASTA or NCBI accession)
- Sample sheet mapping samples to reads

## Citing

Please cite our paper when using snpArcher:

> [snpArcher preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.06.22.546168v1)

Also cite the underlying tools used in your analysis.

## Contributing

- Report bugs or request features on [GitHub](https://github.com/harvardinformatics/snpArcher/issues)
- See [AGENTS.md](https://github.com/harvardinformatics/snpArcher/blob/main/AGENTS.md) for development guidelines

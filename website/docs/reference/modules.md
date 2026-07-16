# Modules

Two optional modules extend the core pipeline, each enabled via `modules.<name>.enabled: true` in `config.yaml`.

## QC

The QC module produces an interactive HTML dashboard with visualizations for quality control of population genomic datasets.

### Requirements

- At least 3 samples in the sample sheet.
- The core pipeline must produce `results/vcfs/raw.vcf.gz` and `results/qc_metrics/qc_report.tsv`.

### Config options

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `modules.qc.enabled` | boolean | `false` | | Enable the module. |
| `modules.qc.clusters` | integer | `3` | >= 1 | Number of clusters for PCA-based clustering in the dashboard. |
| `modules.qc.min_depth` | number | `2` | >= 0 | Samples with mean depth below this value are excluded from QC analyses. |
| `modules.qc.google_api_key` | string | `""` | | Google Maps JavaScript API key for geographic map panels. If empty, map panels are omitted. |
| `modules.qc.exclude_scaffolds` | string | `""` | | Comma-separated list of scaffold/contig names to exclude from QC analyses (e.g., `"scaffold_mt,scaffold_Z"`). |

### Dashboard panels

The HTML dashboard (`results/qc/qc_dashboard.html`) contains:

- PCA plots (PC1 vs PC2, PC1 vs PC3) with cluster assignments
- Per-sample depth vs. PC correlation plots
- Relatedness heatmap (kinship coefficients)
- Inbreeding coefficient (F) distribution
- Neighbor-joining tree
- ADMIXTURE bar plots (K=2 and K=3)
- Per-sample missingness and heterozygosity
- Geographic map of PCA clusters (requires `google_api_key` and `lat`/`long` in sample metadata)
- Filter summary counts

### Outputs

| Path | Format | Description |
|------|--------|-------------|
| `results/qc/qc_dashboard.html` | HTML | Interactive QC dashboard. |
| `results/qc/pruned.vcf.gz` | VCF | LD-pruned, subsampled VCF used for analyses. |
| `results/qc/plink.eigenvec` | Text | PCA eigenvectors. |
| `results/qc/plink.eigenval` | Text | PCA eigenvalues. |
| `results/qc/plink.2.Q` | Text | ADMIXTURE Q matrix (K=2). |
| `results/qc/plink.3.Q` | Text | ADMIXTURE Q matrix (K=3). |
| `results/qc/individuals.idepth` | Text | Per-individual depth. |
| `results/qc/individuals.imiss` | Text | Per-individual missingness. |
| `results/qc/individuals.het` | Text | Per-individual heterozygosity. |

## Postprocess

The postprocess module filters the raw VCF to produce clean SNP and indel call sets by removing excluded samples, applying site-level filters, and intersecting with callable sites.

Sample exclusion is controlled via the `exclude` column in the [sample metadata sheet](sample-sheet-schema.md#optional-metadata-sheet).

### Requirements

- The core pipeline must produce `results/vcfs/raw.vcf.gz`.
- Callable sites BED (`results/callable_sites/callable_sites.bed`) is used if available.

### Config options

| Key | Type | Default | Constraints | Description |
|-----|------|---------|-------------|-------------|
| `modules.postprocess.enabled` | boolean | `false` | | Enable the module. |
| `modules.postprocess.filtering.contig_size` | integer | `10000` | >= 0 | Exclude SNPs on contigs with length equal to or smaller than this value (bp). Set to `0` to disable. |
| `modules.postprocess.filtering.maf` | number | `0.01` | 0 to 1 | Minimum minor allele frequency threshold. |
| `modules.postprocess.filtering.missingness` | number | `0.75` | 0 to 1 | Minimum genotyping rate (fraction of samples with a called genotype). |
| `modules.postprocess.filtering.exclude_scaffolds` | string | `"mtDNA,Y"` | | Comma-separated scaffold/contig names to exclude. |

### Outputs

| Path | Format | Description |
|------|--------|-------------|
| `results/postprocess/filtered.vcf.gz` | VCF | VCF with excluded samples removed and default GATK hard filters applied. |
| `results/postprocess/clean_snps.vcf.gz` | VCF | Final clean biallelic SNP VCF after all site-level filters. |
| `results/postprocess/clean_indels.vcf.gz` | VCF | Final clean indel VCF after all site-level filters. |
| `results/postprocess/callable_sites_filtered.bed` | BED | Callable sites BED after small contig and scaffold exclusion. |
| `results/postprocess/include_samples.txt` | Text | List of retained sample IDs. |


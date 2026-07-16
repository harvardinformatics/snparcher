# Outputs

All output files are written under the `results/` directory relative to the Snakemake working directory.
Paths use the placeholders `{ref}` for `reference.name` and `{sample}` for `sample_id` unless otherwise noted.

## Core pipeline outputs

These files are always produced by the default `all` target.

| Path | Format | Description |
|------|--------|-------------|
| `results/vcfs/raw.vcf.gz` | VCF (bgzipped) | Joint-genotyped, unfiltered VCF containing all samples and all sites called by the chosen variant caller. |
| `results/vcfs/raw.vcf.gz.tbi` | Tabix index | Index for the raw VCF. |
| `results/vcfs/filtered.vcf.gz` | VCF (bgzipped) | Hard-filtered VCF. GATK `VariantFiltration` tags are applied; sites are annotated but not removed. |
| `results/vcfs/filtered.vcf.gz.tbi` | Tabix index | Index for the filtered VCF. |
| `results/qc_metrics/qc_report.tsv` | TSV | Per-sample QC metrics table aggregated from fastp, BAM stats, and coverage summaries. |

## Reference files

| Path | Format | Description |
|------|--------|-------------|
| `results/reference/{ref}.fa.gz` | FASTA (bgzipped) | Reference genome, downloaded or copied from `reference.source`. |
| `results/reference/{ref}.fa.gz.fai` | FAI index | Samtools FASTA index. |
| `results/reference/{ref}.fa.gz.gzi` | GZI index | Bgzip index. |
| `results/reference/{ref}.dict` | Sequence dictionary | Picard/GATK sequence dictionary. |
| `results/reference/{ref}.fa.gz.sa` | BWA index | BWA SA index (plus `.pac`, `.bwt`, `.ann`, `.amb`). |

## Per-sample intermediate files

Many intermediate files are marked as temporary and removed after downstream rules complete.
The tables below list retained per-sample outputs alongside key temporary files.

### FASTQ processing

| Path | Format | Temporary | Description |
|------|--------|:---------:|-------------|
| `results/fastqs/{sample}/{library}/{unit}/{accession}_1.fastq.gz` | FASTQ | yes | Downloaded R1 reads (SRA samples only). |
| `results/fastqs/{sample}/{library}/{unit}/{accession}_2.fastq.gz` | FASTQ | yes | Downloaded R2 reads (SRA samples only). |
| `results/filtered_fastqs/{sample}/{library}/{unit}_1.fastq.gz` | FASTQ | | Adapter-trimmed R1 reads (fastp output). |
| `results/filtered_fastqs/{sample}/{library}/{unit}_2.fastq.gz` | FASTQ | | Adapter-trimmed R2 reads (fastp output). |
| `results/fastp/{sample}/{library}/{unit}.json` | JSON | | Per-unit fastp QC report. |
| `results/qc_metrics/fastp/{sample}.json` | JSON | | Aggregated fastp metrics for the sample. |

### Alignment

| Path | Format | Temporary | Description |
|------|--------|:---------:|-------------|
| `results/bams/raw/{sample}/{library}/{unit}.bam` | BAM | yes | Per-unit aligned BAM (BWA-MEM or Sentieon). |
| `results/bams/library_markdup/{sample}/{library}.bam` | BAM | yes | Per-library BAM after duplicate marking (when `mark_duplicates` is true). |
| `results/bams/library/{sample}/{library}.bam` | BAM | yes | Per-library merged BAM (units merged, duplicates optionally marked). |
| `results/bams/markdup/{sample}.bam` | BAM | | Final merged, deduplicated BAM for the sample (multi-library merge). |
| `results/bams/markdup/{sample}.bam.bai` | BAI index | | Index for the final BAM. |
| `results/bams/merged/{sample}.bam` | BAM | | Final merged BAM when duplicate marking is disabled. |
| `results/bams/merged/{sample}.bam.bai` | BAI index | | Index for the merged BAM. |

### BAM QC metrics

| Path | Format | Temporary | Description |
|------|--------|:---------:|-------------|
| `results/qc_metrics/bam/{sample}_coverage.txt` | Text | yes | Per-sample genome-wide coverage (samtools coverage). |
| `results/qc_metrics/bam/{sample}_flagstat.txt` | Text | yes | Per-sample flagstat output. |
| `results/qc_metrics/bam/{sample}.json` | JSON | | Parsed BAM QC metrics. |

### Sentieon QC metrics

Produced when `variant_calling.tool` is `sentieon`.

| Path | Format | Description |
|------|--------|-------------|
| `results/qc_metrics/sentieon/{sample}_insert_metrics.txt` | Text | Insert size metrics. |
| `results/qc_metrics/sentieon/{sample}_qd_metrics.txt` | Text | Quality distribution metrics. |
| `results/qc_metrics/sentieon/{sample}_gc_metrics.txt` | Text | GC bias metrics. |
| `results/qc_metrics/sentieon/{sample}_gc_summary.txt` | Text | GC bias summary. |
| `results/qc_metrics/sentieon/{sample}_mq_metrics.txt` | Text | Mapping quality metrics. |
| `results/qc_metrics/sentieon/{sample}.json` | JSON | Parsed Sentieon QC metrics. |

### Variant calling (per-sample)

| Path | Format | Temporary | Description |
|------|--------|:---------:|-------------|
| `results/gvcfs/{sample}.g.vcf.gz` | gVCF (bgzipped) | | Per-sample gVCF (GATK, Sentieon, DeepVariant, Parabricks). |
| `results/gvcfs/{sample}.g.vcf.gz.tbi` | Tabix index | | Index for the per-sample gVCF. |

## Interval files

Produced when `intervals.enabled` is `true` (GATK pathway) or when using the Parabricks pathway.

| Path | Format | Description |
|------|--------|-------------|
| `results/intervals/gvcf/intervals.txt` | Text | File-of-filenames listing gVCF interval BED files. |
| `results/intervals/gvcf/` | Directory | Directory of interval list files for per-sample calling. |
| `results/intervals/db/intervals.txt` | Text | File-of-filenames listing GenomicsDB interval BED files. |
| `results/intervals/db/` | Directory | Directory of interval list files for GenomicsDB import. |

## Callable sites outputs

Produced when `callable_sites.coverage.enabled` and/or `callable_sites.mappability.enabled` is `true`.

| Path | Format | Condition | Description |
|------|--------|-----------|-------------|
| `results/callable_sites/depths/{sample}.per-base.d4` | D4 | coverage enabled | Per-sample per-base depth (mosdepth D4). Temporary. |
| `results/callable_sites/depths/{sample}.mosdepth.summary.txt` | Text | coverage enabled | Per-sample mosdepth summary. |
| `results/callable_sites/depths.zarr` | Zarr | coverage enabled | Combined depth data in Zarr format. |
| `results/callable_sites/coverage_thresholds.tsv` | TSV | coverage enabled | Per-sample min/max coverage thresholds (computed or user-specified). |
| `results/callable_sites/callable_loci.zarr` | Zarr | coverage enabled | Callable loci mask in Zarr format. |
| `results/callable_sites/coverage.bed` | BED | coverage enabled | Coverage-based callable sites BED. |
| `results/callable_sites/genmap_index/` | Directory | mappability enabled | GenMap index of the reference genome. |
| `results/callable_sites/mappability.bedgraph` | BedGraph | mappability enabled | Raw mappability scores from GenMap. |
| `results/callable_sites/mappability.bed` | BED | mappability enabled | Mappability-filtered callable sites BED. |
| `results/callable_sites/callable_sites.bed` | BED | `generate_bed_file: true` | Final callable sites BED: intersection of coverage and mappability beds. |

## QC module outputs

Produced when `modules.qc.enabled` is `true`.
Requires at least 3 samples.

| Path | Format | Description |
|------|--------|-------------|
| `results/qc/qc_dashboard.html` | HTML | Interactive QC dashboard with PCA, relatedness, depth, inbreeding, NJ tree, ADMIXTURE, and optional geographic map panels. |
| `results/qc/qc_report.tsv` | TSV | Copy of the pipeline QC metrics report used as dashboard input. |
| `results/qc/pruned.vcf.gz` | VCF | LD-pruned, subsampled VCF used for QC analyses. |
| `results/qc/snpqc.txt` | Text | SNP QC statistics. |
| `results/qc/individuals.idepth` | Text | Per-individual depth statistics (vcftools). |
| `results/qc/individuals.imiss` | Text | Per-individual missingness (vcftools). |
| `results/qc/individuals.het` | Text | Per-individual heterozygosity (vcftools). |
| `results/qc/individuals.FILTER.summary` | Text | Filter summary counts. |
| `results/qc/plink.eigenvec` | Text | PCA eigenvectors (PLINK). |
| `results/qc/plink.eigenval` | Text | PCA eigenvalues (PLINK). |
| `results/qc/plink.bed` | Binary | PLINK binary genotype file. |
| `results/qc/plink.bim` | Text | PLINK variant file. |
| `results/qc/plink.fam` | Text | PLINK sample file. |
| `results/qc/plink.2.Q` | Text | ADMIXTURE Q matrix (K=2). |
| `results/qc/plink.3.Q` | Text | ADMIXTURE Q matrix (K=3). |
| `results/qc/coords.txt` | Text | Sample coordinates (lat/long) for geographic plots. |

## Postprocess module outputs

Produced when `modules.postprocess.enabled` is `true`.

| Path | Format | Description |
|------|--------|-------------|
| `results/postprocess/filtered.vcf.gz` | VCF (bgzipped) | VCF with excluded samples removed and default filters applied. |
| `results/postprocess/filtered.vcf.gz.csi` | CSI index | Index for the filtered VCF. |
| `results/postprocess/callable_sites_filtered.bed` | BED | Callable sites BED after removing small contigs and excluded scaffolds. |
| `results/postprocess/clean_snps.vcf.gz` | VCF (bgzipped) | Final clean SNP-only VCF after MAF, missingness, callable-sites, and scaffold filters. |
| `results/postprocess/clean_snps.vcf.gz.tbi` | Tabix index | Index for the clean SNPs VCF. |
| `results/postprocess/clean_indels.vcf.gz` | VCF (bgzipped) | Final clean indel-only VCF after the same filtering criteria. |
| `results/postprocess/clean_indels.vcf.gz.tbi` | Tabix index | Index for the clean indels VCF. |
| `results/postprocess/include_samples.txt` | Text | List of sample IDs retained after exclusion filtering. |


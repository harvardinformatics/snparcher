# CLI reference

snpArcher is invoked through Snakemake.

## Basic invocation

```bash
snakemake --cores 8 --use-conda --workflow-profile workflow-profiles/default
```

## Named targets

In addition to the default `all` target, snpArcher defines several named targets that run subsets of the pipeline:

| Target | Description |
|--------|-------------|
| `all` | Full pipeline: alignment, variant calling, hard filtering, QC metrics, callable sites, and all enabled modules. |
| `setup` | Prepare reference genome, indices, and interval files. Useful for resolving checkpoints before a full run. |
| `download_reads` | Download FASTQs for all SRA samples. |
| `map_samples` | Align all samples that need alignment. |
| `call_variants` | Run variant calling through to the hard-filtered VCF. |
| `qc_report` | Generate the QC metrics TSV report. |
| `callable_sites` | Generate callable sites BED from coverage and mappability. |
| `gvcfs` | Generate per-sample gVCFs for all samples. |

## Commonly used Snakemake flags

| Flag | Short | Description |
|------|-------|-------------|
| `--cores N` | `-c N` | Maximum number of local CPU cores to use. Required for local execution. |
| `--dry-run` | `-n` | Show what would be executed without running anything. |
| `--configfile PATH` | | Override the default `config/config.yaml` with a custom config file. |
| `--directory PATH` | `-d PATH` | Set the working directory. All relative paths in config and sample sheets resolve from here. |
| `--workflow-profile PATH` | | Path to a directory containing a Snakemake profile YAML. Controls resource defaults and thread counts. |
| `--executor NAME` | `-e NAME` | Snakemake executor plugin (e.g., `slurm`, `googlebatch`). Default is local execution. |
| `--jobs N` | `-j N` | Maximum number of jobs submitted simultaneously to a cluster scheduler. |
| `--latency-wait N` | `-w N` | Seconds to wait for output files to appear (compensates for filesystem latency on network storage). |
| `--retries N` | | Number of times to retry a failed job. |
| `--use-conda` | | Enable conda environment creation for rules that specify a `conda:` directive. Required for snpArcher. |
| `--rerun-triggers {mtime,params,input,software-env,code}` | | Conditions under which to rerun a rule. Snakemake default is `mtime`. |
| `--printshellcmds` | `-p` | Print the shell command of each job before execution. |
| `--snakefile PATH` | `-s PATH` | Explicit path to the Snakefile. Needed when running from outside the snpArcher directory. |

## Profiles

snpArcher uses two types of Snakemake profiles:

### Workflow profile (`--workflow-profile`)

Controls thread counts and resource defaults for snpArcher rules.
The bundled default profile is at `workflow-profiles/default/config.yaml`.

### Execution profile (`--profile`)

Controls how Snakemake submits jobs to a scheduler.
A SLURM template is provided at `workflow-profiles/slurm/config.yaml`.
For cluster execution profiles, see [Run on HPC](../how-to/run-hpc.md).

## Default resource allocation

The bundled workflow profile (`workflow-profiles/default/config.yaml`) sets these defaults:

### Default resources (all rules)

| Resource | Value | Description |
|----------|-------|-------------|
| `mem_mb` | `attempt * 16000` | Memory in MB, scaling with retry attempt. |
| `mem_mb_reduced` | `(attempt * 16000) * 0.9` | Memory allocated to Java processes (GATK), 90% of total to prevent OOM. |
| `tmpdir` | `system_tmpdir` | Temporary directory. Replace with an absolute path to force a custom location. |

### Thread allocations per rule

| Rule | Threads | Stage |
|------|--------:|-------|
| `genmap` | 1 | Mappability |
| `get_fastq_pe` | 6 | FASTQ download |
| `fastp` | 6 | Read trimming |
| `bwa_map` | 16 | Alignment (BWA-MEM) |
| `dedup` | 16 | Duplicate marking |
| `gatk_haplotypecaller` | 1 | GATK HaplotypeCaller (whole-genome) |
| `gatk_haplotypecaller_interval` | 1 | GATK HaplotypeCaller (per-interval) |
| `joint_genomics_db_import` | 1 | GenomicsDB import (whole-genome) |
| `gatk_genomics_db_import` | 1 | GenomicsDB import (per-interval) |
| `compute_d4` | 6 | Mosdepth D4 depth |
| `clam_loci` | 6 | Callable loci computation |
| `sentieon_map` | 16 | Sentieon alignment |
| `sentieon_dedup` | 16 | Sentieon duplicate marking |
| `sentieon_haplotyper` | 32 | Sentieon HaplotypeCaller |
| `sentieon_combine_gvcf` | 32 | Sentieon gVCF merge |
| `bcftools_call` | 8 | bcftools mpileup/call |
| `deepvariant_call` | 8 | DeepVariant inference |
| `glnexus_joint` | 4 | GLnexus joint genotyping |
| `parabricks_haplotypecaller` | 16 | Parabricks HaplotypeCaller |

### Customizing per-rule resources

Per-rule resources can be overridden via `set-resources` in a profile.
See [Run on HPC](../how-to/run-hpc.md).

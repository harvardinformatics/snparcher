# Quickstart: run the test dataset

The snpArcher test dataset ships with five simulated samples and a small reference genome.
Running it end-to-end produces a joint-called, hard-filtered VCF and confirms that your installation is working before you commit time and compute to a real analysis.

**Time estimate:** ~5-10 minutes on a machine with 4 cores.

**Prerequisites:**

- snpArcher is cloned and the `snparcher` conda environment is installed.
  If not, follow the [installation guide](../how-to/install.md) first.

## Step 1: Activate the environment and navigate to snpArcher

Open a terminal and activate the snpArcher conda environment:

```bash
conda activate snparcher
```

Then navigate to the root of the snpArcher repository:

```bash
cd /path/to/snparcher  # <-- change this to your snpArcher clone location
```

The test dataset lives in the `example/` directory.
It contains a pre-configured sample sheet and config file that point to five simulated FASTQ pairs and a small reference genome.

## Step 2: Perform a dry run

A dry run verifies that the configuration files parse correctly and that Snakemake can build the full dependency graph, without executing any jobs.

```bash
snakemake --use-conda --dry-run --directory example/
```

You should see output listing all the rules that Snakemake would run, along with the total number of jobs.
Look for familiar rule names like `fastp`, `bwa_map`, `gatk_haplotypecaller`, and `genotype_gvcfs`.
If you see an error instead, double-check that you activated the `snparcher` conda environment and that you are in the snpArcher repository root.

!!! tip "What does `--dry-run` do?"
    The `--dry-run` (or `-n`) flag tells Snakemake to resolve the full dependency graph and report what it *would* do, without executing any jobs.
    This is a free sanity check. Use it before every real run.

## Step 3: Run the pipeline

Run the pipeline with `--cores 4`:

```bash
snakemake --use-conda --cores 4 --directory example/
```

!!! note "First run takes longer"
    On the first execution, Snakemake will create conda environments for each pipeline step (downloading and installing tools like BWA, GATK, samtools, etc.).
    This adds several minutes of one-time setup.
    Subsequent runs reuse the cached environments.

Snakemake will print progress to the terminal as each rule starts and finishes.
The entire run should complete in roughly five minutes once the conda environments are in place.

## Step 4: Check the outputs

A successful run produces output files under `example/results/`.
The most important output is the hard-filtered VCF.
Verify it exists and is not empty:

```bash
ls -lh example/results/vcfs/filtered.vcf.gz
```

You should see a file of at least a few kilobytes.
To confirm it contains variant records:

```bash
zgrep -c -v "^#" example/results/vcfs/filtered.vcf.gz
```

You should see a positive number of variant records.

!!! tip "Other output files"
    The test run also produces per-sample BAMs (`example/results/bams/`) and per-sample gVCFs (`example/results/gvcfs/`).
    For a complete list of output files and what they contain, see the [outputs reference](../reference/outputs.md).

## Step 5: Verify QC metrics

The core pipeline always produces a QC metrics summary, even when the interactive QC dashboard is not enabled.
Check that it was generated:

```bash
ls -lh example/results/qc_metrics/qc_report.tsv
```

This TSV file contains per-sample summary statistics including mapping rate and mean depth.

!!! note "QC dashboard"
    The interactive QC dashboard (`results/qc/qc_dashboard.html`) is produced only when `modules.qc.enabled` is set to `true` in the config.
    The test config has it disabled by default to keep the test fast.
    In a real project, we strongly recommend enabling it. See the [first project tutorial](first-project.md) and the [QC report how-to guide](../how-to/qc-report.md) for details.

## What just happened?

snpArcher executed the full variant calling pipeline:

1. **Trim** (fastp) and **align** (BWA-MEM) reads to the reference genome.
2. **Mark duplicates** in the aligned BAMs.
3. **Call variants** per sample (GATK HaplotypeCaller), splitting the genome at runs of Ns for parallelism.
4. **Joint genotype** across samples via GenomicsDB and GenotypeGVCFs.
5. **Hard filter** the joint VCF to produce the final call set.

For a deeper explanation of the pipeline architecture and the scatter-by-Ns parallelization strategy, see [Architecture](../explanation/architecture.md) and [Parallelization](../explanation/parallelization.md).

## Next steps

- **[Your first project](first-project.md)**: A real analysis from SRA data to a filtered VCF, including QC review and postprocessing.
- **[How to configure snpArcher](../how-to/configure.md)**: Customize variant calling settings for your organism.
- **[How to run on an HPC cluster](../how-to/run-hpc.md)**: Submit snpArcher jobs to SLURM or other schedulers.

# Your first project: from SRA to filtered VCF

We download 10 wild zebra finch (*Taeniopygia guttata*) whole-genome resequencing samples from the NCBI Sequence Read Archive, call variants against a reference genome, review the QC dashboard, and produce a filtered VCF ready for downstream population genomic analysis.

This dataset is from BioProject [PRJEB10586](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB10586) (Singhal et al. 2015, *Science*) and was used to benchmark snpArcher in [Mirchandani et al. 2024](https://doi.org/10.1093/molbev/msad270).

**Time estimate:** Several hours of compute time (varies with machine and network speed).
The hands-on configuration steps take about 20 minutes.

**Prerequisites:**

- snpArcher is cloned and the `snparcher` conda environment is installed.
  If not, follow the [installation guide](../how-to/install.md).
- You have successfully run the [Quickstart tutorial](quickstart.md) to confirm your installation works.
- At least 50 GB of free disk space for intermediate files and outputs.
- A machine with at least 4 CPU cores and 16 GB of RAM (more cores will speed up the run considerably).

!!! tip "HPC recommended for real projects"
    This tutorial uses local execution for simplicity.
    For datasets larger than this example (more samples, deeper coverage, or a larger genome), we strongly recommend running on an HPC cluster.
    We show both local and SLURM commands at the execution step. See [How to run on an HPC cluster](../how-to/run-hpc.md) for full details.

## Step 1: Set up the project directory

We recommend keeping each snpArcher analysis in its own project directory, separate from the snpArcher repository itself.
A single clone of snpArcher can serve any number of projects.

Assuming your snpArcher clone is at `~/snparcher` and you want to organize projects under a `~/projects` directory:

```bash
mkdir -p ~/projects/zebrafinch_reseq/config
```

```text
~/
├── snparcher/                  # snpArcher repository clone
└── projects/
    └── zebrafinch_reseq/
        └── config/             # will hold config.yaml and samples.csv
```

## Step 2: Create the sample sheet

The sample sheet tells snpArcher what samples to process and where their data are.
Since we are downloading data from the SRA, each row needs a `sample_id`, `input_type` set to `srr`, and an SRA run accession in the `input` column.

Create the file `~/projects/zebrafinch_reseq/config/samples.csv` with the following content:

```csv
sample_id,input_type,input
zf_26462,srr,ERR1013161
zf_28339,srr,ERR1013162
zf_28353,srr,ERR1013163
zf_26721,srr,ERR1013164
zf_28456,srr,ERR1013165
zf_28402,srr,ERR1013166
zf_26516,srr,ERR1013167
zf_28404,srr,ERR1013168
zf_26820,srr,ERR1013169
zf_26733,srr,ERR1013170
```

Column definitions:

- **`sample_id`**: A unique name for each sample. This will appear as the sample name in the final VCF. Must contain only letters, numbers, periods, underscores, and hyphens.
- **`input_type`**: Set to `srr` because we are downloading from the SRA. Other options are `fastq` (local FASTQ files), `bam` (aligned BAMs), and `gvcf` (per-sample gVCFs).
- **`input`**: The SRA run accession. snpArcher will download the reads automatically.

!!! note "Using local FASTQ files instead?"
    If your reads are already on disk, set `input_type` to `fastq` and provide the forward and reverse read paths joined with a semicolon:
    ```csv
    sample_id,input_type,input
    zf_26462,fastq,/storage/data/zf_26462_R1.fq.gz;/storage/data/zf_26462_R2.fq.gz
    ```
    Use absolute paths so the files can be found regardless of working directory.

For a full description of all sample sheet columns, see the [sample sheet reference](../reference/sample-sheet-schema.md).

## Step 3: Write the config file

Copy the default config from the snpArcher repository to your project:

```bash
cp ~/snparcher/config/config.yaml ~/projects/zebrafinch_reseq/config/config.yaml
```

Now open `~/projects/zebrafinch_reseq/config/config.yaml` in your editor and modify it to match the following:

```yaml
samples: "config/samples.csv"

reference:
  name: "zebrafinch"  # <-- short name for output filenames
  source: "GCF_048771995.1"  # <-- RefSeq accession; snpArcher downloads it automatically

variant_calling:
  expected_coverage: "low"
  tool: "gatk"
  ploidy: 2
  gatk:
    het_prior: 0.005  # raised from GATK default of 0.001; zebra finch has higher diversity than humans

intervals:
  enabled: true

callable_sites:
  generate_bed_file: true
  coverage:
    enabled: true
  mappability:
    enabled: true

modules:
  qc:
    enabled: true  # <-- we want the QC dashboard
    clusters: 3
  postprocess:
    enabled: false  # we will enable this in Step 9 after reviewing QC
```

Key choices in this config:

- **`reference.source`**: A RefSeq accession. snpArcher will download and index the genome automatically. You can also provide a local file path or a URL here.
- **`variant_calling.tool: "gatk"`**: GATK HaplotypeCaller is the default and best-tested variant caller in snpArcher.
- **`variant_calling.gatk.het_prior: 0.005`**: The GATK default (0.001) is calibrated for humans. Most non-model organisms have higher nucleotide diversity, so we raise this to 0.005 to improve variant detection sensitivity. See [Variant calling explained](../explanation/variant-calling.md) for more on this choice.
- **`variant_calling.ploidy: 2`**: Zebra finches are diploid. Set this correctly for your organism.
- **`modules.qc.enabled: true`**: Generates the interactive QC dashboard, which we will examine in Step 8.

!!! note "Default settings"
    The remaining options (`intervals`, `callable_sites`, and several `variant_calling` sub-blocks) use sensible defaults for most datasets.
    For a complete list of every config option, see the [config reference](../reference/config-schema.md).

## Step 4: Set up the execution profile

The execution profile controls computational resources (threads, memory) and Snakemake options.
Copy the default profile from the snpArcher repository to your project:

```bash
cp -r ~/snparcher/workflow-profiles ~/projects/zebrafinch_reseq/
```

For local execution, the defaults in `workflow-profiles/default/config.yaml` are a reasonable starting point.
The profile sets thread counts for each pipeline step (e.g., 16 threads for BWA alignment, 6 threads for read downloading) and memory allocation.

!!! tip "Adjusting for your machine"
    If your workstation has fewer than 16 cores, reduce the `bwa_map` thread count in `workflow-profiles/default/config.yaml`:
    ```yaml
    set-threads:
      bwa_map: 8  # <-- reduce if you have fewer cores
    ```
    The `--cores` flag you pass to Snakemake (in Step 6) sets the total pool of available cores.
    Snakemake will not start a 16-thread job if only 8 cores are available, but reducing per-rule threads lets more jobs run in parallel.

## Step 5: Dry run

Before committing compute time, verify that snpArcher has parsed the configuration correctly:

```bash
conda activate snparcher
snakemake \
  --snakefile ~/snparcher/workflow/Snakefile \
  --use-conda \
  --dry-run \
  --directory ~/projects/zebrafinch_reseq
```

Check the output for the following:

- **Sample IDs match**: You should see `zf_26462` through `zf_26733` appear in the planned job names.
- **Expected rules are queued**: Look for `get_fastq_pe` (SRA download), `fastp` (trimming), `bwa_map` (alignment), `gatk_haplotypecaller` (variant calling), and QC rules.
- **No errors**: If you see a schema validation error, double-check your sample sheet columns and config syntax.

!!! tip "Reading the dry-run output"
    The dry-run prints a summary at the end showing the total number of jobs and a count per rule.
    With 10 samples and interval-split calling, expect hundreds of jobs. This is normal and is how snpArcher achieves parallelism.

## Step 6: Run the pipeline

Since this will take a while, run inside a terminal multiplexer so the process survives if your terminal disconnects:

```bash
tmux new -s snparcher
```

Then start the pipeline:

```bash
conda activate snparcher
snakemake \
  --snakefile ~/snparcher/workflow/Snakefile \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --cores 8 \
  --directory ~/projects/zebrafinch_reseq \
  --workflow-profile ~/projects/zebrafinch_reseq/workflow-profiles/default
```

Flag breakdown:

- **`--snakefile`**: Points to the snpArcher Snakefile in your clone.
- **`--use-conda`**: Tells Snakemake to create and manage isolated conda environments for each pipeline step.
- **`--conda-prefix ~/snparcher_envs`**: Caches conda environments in a central location so they can be reused across projects. Without this, each project reinstalls all tools independently.
- **`--cores 8`**: Allow up to 8 concurrent cores. Adjust to match your machine.
- **`--directory`**: Points to your project directory where the config and sample sheet live.
- **`--workflow-profile`**: Points to the resource profile controlling threads and memory.

To detach from tmux without stopping the run, press `Ctrl-b` then `d`.
Reattach later with `tmux attach -t snparcher`.

??? example "SLURM cluster execution"
    If you are on an HPC cluster with SLURM, first install the executor plugin:

    ```bash
    conda activate snparcher
    pip install snakemake-executor-plugin-slurm
    ```

    Edit the workflow profile to uncomment and set `slurm_partition` and (if required) `slurm_account` under `default-resources`.

    Then submit the Snakemake controller as a SLURM job:

    ```bash
    #!/bin/bash
    #SBATCH --job-name=snparcher
    #SBATCH --time=7-00:00:00
    #SBATCH --mem=4G
    #SBATCH --cpus-per-task=1
    #SBATCH --output=snparcher_controller_%j.log

    conda activate snparcher
    snakemake \
      --executor slurm \
      --snakefile ~/snparcher/workflow/Snakefile \
      --use-conda \
      --conda-prefix ~/snparcher_envs \
      --directory ~/projects/zebrafinch_reseq \
      --workflow-profile ~/projects/zebrafinch_reseq/workflow-profiles/default
    ```

    The controller job needs minimal resources. It only orchestrates other jobs that run on compute nodes.
    See [How to run on an HPC cluster](../how-to/run-hpc.md) for full details.

## Step 7: Monitor progress

While the pipeline runs, Snakemake prints progress to the terminal.
You will see messages like:

```text
rule fastp:
    input: data/reads/zf_26462_R1.fq.gz, data/reads/zf_26462_R2.fq.gz
    output: results/reads/trimmed/zf_26462_R1.fq.gz, results/reads/trimmed/zf_26462_R2.fq.gz
    jobid: 42
```

If a job fails, Snakemake reports which rule and sample encountered the error.
To diagnose failures, check the log files in this order:

1. **Snakemake log** (`.snakemake/logs/`). Identifies which rule failed.
2. **Individual job log** (`logs/<rule_name>/`). Contains the actual error message from the underlying tool (e.g., BWA, GATK). This is usually the most informative.

!!! tip "Resuming after failure"
    If the pipeline stops due to a transient error (e.g., network timeout during SRA download), simply rerun the same command.
    Snakemake automatically detects which steps completed and picks up where it left off.

## Step 8: Examine the QC report

Once the pipeline completes, the QC dashboard is at:

```text
~/projects/zebrafinch_reseq/results/qc/qc_dashboard.html
```

Open this file in a web browser.
The dashboard is interactive. You can zoom, pan, and hover over points to see sample IDs.

Walk through each panel:

1. **Header summary**: Check the number of samples processed, total SNPs detected, estimated nucleotide diversity (Watterson's theta), and average genome-wide depth. For 10 zebra finch samples, expect Watterson's theta on the order of 0.001-0.01 and millions of SNPs given the ~1.1 Gb genome.

2. **PCA (Figure 1)**: Samples are colored by k-means cluster assignment (k=3 by default). Look for any unexpected groupings. With 10 wild zebra finch samples, you might expect a relatively tight cluster. Distinct groups could indicate population structure or a technical batch effect worth investigating.

3. **Depth vs. PC plots (Figure 2)**: Check whether sequencing depth correlates with any of the top principal components. A strong correlation indicates that depth variation is driving more of the observed genetic structure than real biological differences, and that depth normalization or stricter sample filtering is warranted.

4. **Depth vs. missingness (Figure 3)**: Samples with low depth tend to have high missingness. Identify outliers with unusually high missingness for exclusion in postprocessing.

5. **Mapping rate vs. depth (Figure 4)**: Mapping rates below ~80% may indicate contamination or that a sample is not the expected species. All zebra finch samples should have consistently high mapping rates against the zebra finch reference.

6. **Inbreeding coefficient (Figure 5)**: Look for outliers with strongly negative F values, which can indicate cross-contamination between samples. Positive values indicate excess homozygosity.

7. **Neighbor-joining tree (Figure 6)**: Provides a quick phylogenetic view. Branches are colored by PCA cluster, so you can assess whether PCA groupings and phylogenetic relationships agree.

8. **Relatedness heatmap (Figure 7)**: The KING kinship coefficient matrix flags closely related individuals. Values near 0.5 indicate duplicate or identical samples; values near 0.25 indicate first-degree relatives. Remove one member of each highly related pair to avoid biasing population-level estimates.

9. **ADMIXTURE (Figure 10)**: Shows model-based ancestry proportions for k=2 and k=3. Useful for a quick check of population structure, but formal ADMIXTURE analyses with cross-validation should be done separately.

!!! note "QC is for guidance, not final analysis"
    The QC report uses a random LD-pruned subset of ~100,000 SNPs after light filtering.
    It surfaces problems and guides decisions but is not publication-ready. See [How to interpret the QC report](../how-to/qc-report.md) for filter details.

## Step 9: Postprocess the VCF

After reviewing the QC dashboard, you may have identified samples to exclude (low depth, contaminated, duplicated, or highly related).
The postprocessing module re-filters the VCF with those exclusions applied and also performs site-level filtering.

### Mark samples for exclusion

Create a sample metadata file at `~/projects/zebrafinch_reseq/config/sample_metadata.csv`.
This file links back to sample IDs from your sample sheet and provides module-specific metadata.
For postprocessing, the key column is `exclude`:

```csv
sample_id,exclude
zf_26462,false
zf_28339,false
zf_28353,true
zf_26721,false
zf_28456,false
zf_28402,false
zf_26516,false
zf_28404,false
zf_26820,false
zf_26733,false
```

In this example, we are excluding `zf_28353` (perhaps it had low depth or high missingness in the QC report).
Set `exclude` to `true` for any sample you want to remove.

### Enable postprocessing in the config

Edit `~/projects/zebrafinch_reseq/config/config.yaml` and add the `sample_metadata` path and enable the postprocess module:

```yaml
sample_metadata: "config/sample_metadata.csv"  # <-- add this line

modules:
  qc:
    enabled: true
    clusters: 3
  postprocess:
    enabled: true  # <-- change from false to true
    filtering:
      contig_size: 10000     # exclude scaffolds smaller than 10 kb
      maf: 0.01              # minor allele frequency threshold
      missingness: 0.75      # keep sites genotyped in at least 75% of samples
      exclude_scaffolds: ""  # <-- comma-separated scaffold names to exclude, e.g. "mtDNA,W"
```

### Run postprocessing

Rerun snpArcher with the same command as Step 6.
Snakemake will detect that the core pipeline already completed and will only run the newly enabled postprocessing rules:

```bash
snakemake \
  --snakefile ~/snparcher/workflow/Snakefile \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --cores 8 \
  --directory ~/projects/zebrafinch_reseq \
  --workflow-profile ~/projects/zebrafinch_reseq/workflow-profiles/default
```

## Step 10: Verify the final output

Postprocessing produces several output files under `results/postprocess/`:

```bash
ls -lh ~/projects/zebrafinch_reseq/results/postprocess/
```

You should see:

- **`filtered.vcf.gz`**: The VCF after removing excluded samples and sites failing the hard filters or with reference N alleles or zero allele frequency.
- **`clean_snps.vcf.gz`**: Biallelic SNPs from the filtered VCF, after applying your missingness, MAF, callable-site, and scaffold-exclusion filters. This is typically the file you want for downstream analyses.
- **`clean_indels.vcf.gz`**: Indels with the same filters applied.

Count variants in the clean SNPs file to confirm it is populated:

```bash
bcftools stats ~/projects/zebrafinch_reseq/results/postprocess/clean_snps.vcf.gz | grep "^SN"
```

This prints summary numbers including total records, SNPs, and samples.
Confirm the sample count matches your expectation (original 10 minus any excluded).

You can also verify that excluded samples were actually removed:

```bash
bcftools query -l ~/projects/zebrafinch_reseq/results/postprocess/clean_snps.vcf.gz
```

This lists the sample names in the VCF.
`zf_28353` (or whichever samples you excluded) should not appear.

## Summary

You now have a complete snpArcher project: SRA accessions in, QC-reviewed and postprocessed VCF out, with a reusable directory structure and config for future analyses.

## Next steps

How-to guides for specific tasks:

- **[How to create a sample sheet](../how-to/sample-sheet.md)**: Handling multi-lane samples, mixed input types, and local FASTQ files.
- **[How to configure snpArcher](../how-to/configure.md)**: Deep dive into config options for different organisms and study designs.
- **[How to run on an HPC cluster](../how-to/run-hpc.md)**: SLURM profiles, resource tuning, and submitting the controller job.
- **[How to interpret the QC report](../how-to/qc-report.md)**: Detailed guidance on each QC panel and what to look for.
- **[How to postprocess and filter](../how-to/postprocess.md)**: Advanced filtering strategies including callable sites and SFS-guided filtering.

For deeper understanding of the design decisions behind snpArcher:

- **[Architecture](../explanation/architecture.md)**: How the pipeline is structured and why.
- **[Variant calling](../explanation/variant-calling.md)**: GATK vs. bcftools vs. DeepVariant: when and why to choose each.
- **[Parallelization](../explanation/parallelization.md)**: The scatter-by-Ns strategy explained.
- **[Filtering philosophy](../explanation/filtering.md)**: SFS-guided filtering with a worked burrowing owl example.

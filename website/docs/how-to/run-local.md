# How to run locally

Run snpArcher on a single machine: a laptop, workstation, or interactive cluster session.

## Prerequisites

You need:

1. [Installed snpArcher](install.md) and activated the `snparcher` conda environment
2. [Created a sample sheet](sample-sheet.md)
3. [Configured your run](configure.md) with a `config.yaml`

## Set up your project directory

Organize each analysis as a self-contained project directory:

```text
my_project/
├── config/
│   ├── config.yaml
│   └── samples.csv
```

Copy the default config from the snpArcher repository as a starting point:

```bash
mkdir -p my_project/config
cp /path/to/snpArcher/config/config.yaml my_project/config/config.yaml
```

Edit `my_project/config/config.yaml` to point to your sample sheet and reference genome.

!!! tip
    A single snpArcher clone can serve many projects.
    You do not need a copy of the repository in each project directory.

## Dry run first

Always perform a dry run before committing compute resources:

```bash
conda activate snparcher
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --dry-run
```

Check the output for:

- **Sample IDs match your sample sheet.** Sample identifiers in the job list should correspond to your `sample_id` values.
- **Expected rules are listed.** If you enabled the QC module, you should see QC rules; if postprocessing is disabled, those rules should be absent.
- **No errors about missing inputs or invalid config**.

## Run the pipeline

Once the dry run looks correct, launch the full run:

```bash
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --cores 8  # <-- change to the number of cores available
```

### Key flags

| Flag | Purpose |
|---|---|
| `--snakefile` / `-s` | Path to snpArcher's `workflow/Snakefile`. Not needed if you run from inside the snpArcher directory. |
| `--directory` / `-d` | Path to your project directory (where `config/config.yaml` lives). |
| `--use-conda` | Required. Tells Snakemake to create isolated environments for each pipeline step. |
| `--cores` | Maximum number of CPU cores to use concurrently. Set this to the number of available cores. |
| `--workflow-profile` | Path to a workflow profile directory for resource settings. Defaults to `workflow-profiles/default` in the snpArcher repo. |

### Share conda environments across projects

By default, Snakemake creates per-step conda environments inside the project directory under `.snakemake/conda/`, so each project gets its own copy of every tool.

To share environments across projects, set a central conda prefix:

```bash
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --cores 8
```

Any subsequent run using the same `--conda-prefix` will reuse existing environments instead of rebuilding them.

## Use a custom workflow profile

The workflow profile controls per-rule resource allocation (threads, memory).
snpArcher ships with a default profile at `workflow-profiles/default/config.yaml`.

To use a custom profile, copy it to your project and pass the directory path:

```bash
cp -r /path/to/snpArcher/workflow-profiles/default my_project/workflow-profile
# Edit my_project/workflow-profile/config.yaml as needed
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --workflow-profile my_project/workflow-profile \
  --cores 8
```

## Prevent disconnection from killing your run

Local runs can take hours.
Use a terminal multiplexer to keep the process alive after disconnecting:

```bash
tmux new -s snparcher
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --cores 8
# Detach: Ctrl-b, then d
# Reattach later: tmux attach -t snparcher
```

## Resume after interruption

Snakemake tracks which steps have completed.
If the pipeline is interrupted, re-run the same command and it will pick up where it left off.

```bash
# Same command as before; Snakemake resumes automatically
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --cores 8
```

!!! tip
    If a job failed partway through and left behind incomplete output files, add `--rerun-incomplete` to force those jobs to re-run:

    ```bash
    snakemake \
      --snakefile /path/to/snpArcher/workflow/Snakefile \
      --directory /path/to/my_project \
      --use-conda \
      --cores 8 \
      --rerun-incomplete
    ```

## Where to find outputs

All outputs are written to the `results/` directory inside your project.
Key output files include:

| Output | Path |
|---|---|
| Hard-filtered VCF | `results/vcfs/filtered.vcf.gz` |
| Raw (unfiltered) VCF | `results/vcfs/raw.vcf.gz` |
| Per-sample BAMs | `results/bams/markdup/{sample}.bam` |
| Per-sample gVCFs | `results/gvcfs/{sample}.g.vcf.gz` |
| QC dashboard | `results/qc/qc_dashboard.html` |
| Callable sites BED | `results/callable_sites/callable_sites.bed` |

For a complete listing, see the [outputs reference](../reference/outputs.md).

## Next steps

- [Interpret QC reports](qc-report.md) to check for problems
- [Filter and postprocess](postprocess.md) to prepare a clean VCF
- For larger datasets, consider [running on HPC](run-hpc.md)

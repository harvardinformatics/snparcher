# How to run on HPC

SLURM is the primary example here because it is the most common scheduler and the one snpArcher has been most thoroughly tested on.

## Install the SLURM executor plugin

In the conda environment where you installed Snakemake, install the SLURM plugin:

```bash
conda activate snparcher
pip install snakemake-executor-plugin-slurm
```

## Set up the workflow profile

The workflow profile is a YAML file specifying how much memory, time, and how many threads to allocate for each pipeline step.
snpArcher ships with a default profile at `workflow-profiles/default/config.yaml`.

Copy it into your project directory:

```bash
mkdir -p /path/to/my_project/workflow-profile
cp /path/to/snpArcher/workflow-profiles/default/config.yaml \
   /path/to/my_project/workflow-profile/config.yaml
```

Edit the copied file to configure SLURM-specific settings.
At minimum, uncomment and set the partition and account fields under `default-resources`:

```yaml
default-resources:
  mem_mb: attempt * 16000
  mem_mb_reduced: (attempt * 16000) * 0.9
  tmpdir: system_tmpdir
  slurm_partition: "short"  # <-- change to your cluster's partition
  slurm_account: "mylab"  # <-- change this, or remove if not required
  runtime: 720  # in minutes
```

!!! note
    Not all clusters require `slurm_account`.
    If yours does not, you can leave it blank or remove the line.

## Customize per-rule resources

The `set-threads` and `set-resources` sections let you override the defaults for specific rules.

### Threads

The profile ships with reasonable defaults for thread counts.
Adjust them if your cluster nodes have more or fewer cores:

```yaml
set-threads:
  bwa_map: 16  # <-- alignment benefits from many threads
  dedup: 16
  fastp: 6
  bcftools_call: 8
```

### Memory and partitions

Use `set-resources` to customize memory, runtime, or partition for individual rules.
This is useful when a specific step needs more resources than the default:

```yaml
set-resources:
  gatk_haplotypecaller:
    mem_mb: attempt * 8000
    mem_mb_reduced: (attempt * 8000) * 0.9
    slurm_partition: "long"  # <-- move long-running jobs to a different partition
    runtime: 600
```

!!! tip "Dynamic memory with retries"
    The `attempt * N` syntax increases memory on each retry.
    For example, `attempt * 8000` allocates 8 GB on the first attempt, 16 GB on the second, and 24 GB on the third.
    This requires setting retries in your Snakemake command (see below).

## Dry run

Before submitting jobs, validate your configuration with a dry run:

```bash
conda activate snparcher
snakemake \
  --executor slurm \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --workflow-profile /path/to/my_project/workflow-profile \
  --use-conda \
  --dry-run
```

## Submit the workflow

Two approaches:

### Option A: interactive session with `srun`

Request an interactive session on a compute node and run Snakemake from there:

```bash
srun --time=7-00:00:00 --mem=4G --cpus-per-task=1 --pty bash
conda activate snparcher
snakemake \
  --executor slurm \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --workflow-profile /path/to/my_project/workflow-profile \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --jobs 100 \
  --retries 3 \
  --latency-wait 20
```

### Option B: batch job with `sbatch`

Create a submission script (e.g., `run_snparcher.sh`):

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
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --workflow-profile /path/to/my_project/workflow-profile \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --jobs 100 \
  --retries 3 \
  --latency-wait 20
```

Submit it:

```bash
sbatch run_snparcher.sh
```

!!! tip "Controller job resources"
    The controller job (the Snakemake process itself) only orchestrates other jobs.
    It needs minimal resources: 4 GB of memory and 1 CPU is typically sufficient.
    Set the time limit generously to keep it alive for the entire run.

### Key flags for cluster execution

| Flag | Purpose |
|---|---|
| `--executor slurm` | Submit individual steps as SLURM jobs instead of running locally. |
| `--jobs` | Maximum number of jobs submitted to the queue at any time. |
| `--retries` | Number of times to retry a failed job. Works with `attempt * N` memory scaling. |
| `--latency-wait` | Seconds to wait for output files after a job finishes. Accounts for filesystem latency. |
| `--conda-prefix` | Central directory for conda environments. Avoids rebuilding envs for each project. |

## Where to find logs

See [Where to find logs](troubleshoot.md#where-to-find-logs) for debugging failed jobs.

## Other schedulers

Snakemake supports many schedulers through its [plugin catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html).
The general approach is the same:

1. Install the executor plugin (e.g., `pip install snakemake-executor-plugin-lsf`).
2. Set up a workflow profile with the appropriate resource fields for your scheduler.
3. Run with `--executor <scheduler_name>`.

Consult the plugin documentation for scheduler-specific resource field names and configuration.

## Next steps

- [Interpret QC reports](qc-report.md) once the pipeline completes
- [Troubleshoot errors](troubleshoot.md) if jobs fail

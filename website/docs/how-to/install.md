# How to install snpArcher

## Prerequisites

- **Operating system:** Linux or macOS (Windows users should use WSL)
- **Git:** for cloning the repository
- **Conda or Mamba:** for managing environments

If you are on an institutional cluster, conda may already be available.
Check by running:

```bash
conda --version
```

If you see a version number, you are set.
If you see "command not found", install conda through [Miniforge](https://github.com/conda-forge/miniforge):

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

## Install Snakemake

Create a dedicated conda environment with Snakemake:

```bash
conda create -c conda-forge -c bioconda -n snparcher "snakemake>=9"
```

Activate the environment:

```bash
conda activate snparcher
```

!!! note
    This is the only software you need to install manually.
    snpArcher uses Snakemake's `--use-conda` flag to create isolated environments for each pipeline step automatically.

## Clone snpArcher

Clone the repository:

```bash
git clone https://github.com/harvardinformatics/snpArcher.git
```

You can place this clone anywhere on your filesystem.
A single clone can serve multiple independent projects.

!!! tip "Pinning a release"
    To use a specific release, visit the [Releases page](https://github.com/harvardinformatics/snpArcher/releases) and download the version you want, or check out a tag after cloning:

    ```bash
    cd snpArcher
    git checkout v2.0.0  # <-- change this to the desired version
    ```

## Verify the installation

Run the bundled example dataset to confirm everything is working:

```bash
conda activate snparcher
snakemake --use-conda --cores 4 --directory snpArcher/example/
```

This processes five simulated samples against a small reference genome.
It should complete in about five minutes on a machine with four cores.
A successful run produces output files under `example/results/` with no error messages.

To do a faster check without actually running any jobs, use the dry-run flag:

```bash
snakemake --use-conda --dry-run --directory snpArcher/example/
```

This resolves the full dependency graph and reports which jobs would be run, without executing them.

## Optional: install executor plugins for HPC

If you plan to run snpArcher on a cluster, you need an executor plugin.
For SLURM (the most common scheduler):

```bash
conda activate snparcher
pip install snakemake-executor-plugin-slurm
```

For other schedulers, find the appropriate plugin in the [Snakemake Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html) and follow its installation instructions.

See [How to run on HPC](run-hpc.md) for full cluster setup details.

## Next steps

- [Create a sample sheet](sample-sheet.md) describing your samples
- [Configure your run](configure.md) with a `config.yaml` file
- [Run locally](run-local.md) or [run on HPC](run-hpc.md)

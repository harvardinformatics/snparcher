# How to troubleshoot errors

## Where to find logs

Check logs in this order:

1. **Snakemake log** at `.snakemake/logs/<timestamp>.snakemake.log` in your project directory. This identifies which rule failed and (on a cluster) the job ID.
2. **Cluster log** at `.snakemake/slurm_logs/<rule>/<wildcards>/<job_id>.log` (SLURM only). Shows resource-related failures like out-of-memory or time-limit exceeded.
3. **Tool log** under `logs/<rule_name>/<sample_id>/` in your project directory. Contains the actual standard output and standard error from the bioinformatics tool. This is usually the most informative.

## Snakemake will not start

### "Missing input files" error

**Symptom:** Snakemake reports missing input files and refuses to build the DAG.

**Causes and fixes:**

- **Sample sheet paths are wrong.** Check that the `input` column in your sample sheet contains correct, absolute paths. Relative paths are resolved from the project directory (`--directory`).
- **Reference genome not found.** Verify that `reference.source` in `config.yaml` points to a valid file, URL, or accession.
- **Sample sheet not found.** Confirm that `samples` in `config.yaml` points to your CSV file, relative to the project directory.

### Config validation errors

**Symptom:** Snakemake reports a schema validation error.

**Common causes:**

- **Wrong column names in sample sheet.** The v2 sample sheet uses `sample_id`, `input_type`, and `input`. If you are using v1 column names (`BioSample`, `Run`, `fq1`, `fq2`, etc.), the validation will fail.
- **Invalid characters in `sample_id`.** Only alphanumeric characters, periods, underscores, and hyphens are allowed (`^[A-Za-z0-9._-]+$`).
- **`input_type` not recognized.** Must be one of `srr`, `fastq`, `bam`, or `gvcf`.
- **FASTQ input missing semicolon.** For `input_type: fastq`, the `input` column must contain forward and reverse paths separated by a semicolon (e.g., `/data/R1.fq.gz;/data/R2.fq.gz`).
- **gVCF with incompatible caller.** gVCF input is only supported with `gatk` and `sentieon`. If your sample sheet includes `gvcf` rows and `variant_calling.tool` is set to `bcftools`, `deepvariant`, or `parabricks`, validation will fail.

### "No rule to produce" error

**Symptom:** `MissingRuleException: No rule to produce <target>`.

**Causes and fixes:**

- **Module not enabled.** If you expect postprocess or QC outputs, make sure the corresponding module is enabled in `config.yaml`.
- **Wrong Snakefile path.** If using `--snakefile`, verify it points to `workflow/Snakefile` in the snpArcher repository.

## Jobs fail on the cluster

### Out of memory (OOM)

**Symptom:** SLURM log shows `CANCELLED` with reason `OUT_OF_MEMORY`, or the tool log shows a Java `OutOfMemoryError`.

**Fixes:**

- **Increase memory for the failing rule** in your workflow profile under `set-resources`:

    ```yaml
    set-resources:
      gatk_haplotypecaller:
        mem_mb: attempt * 16000
        mem_mb_reduced: (attempt * 16000) * 0.9
    ```

- **Enable retries** with `--retries 3` in your Snakemake command. Combined with `attempt * N` memory scaling, this automatically increases memory on each retry.
- **Reduce concurrency.** On a local machine, lowering `--cores` reduces the number of memory-intensive jobs running simultaneously.

### Time limit exceeded

**Symptom:** SLURM log shows `TIMEOUT` or `CANCELLED` with reason `TIME_LIMIT`.

**Fixes:**

- **Increase runtime** for the failing rule in your workflow profile:

    ```yaml
    set-resources:
      gatk_haplotypecaller:
        runtime: 1440  # 24 hours, in minutes
    ```

- **Move long-running rules to a partition with longer time limits:**

    ```yaml
    set-resources:
      gatk_haplotypecaller:
        slurm_partition: "long"
    ```

### Conda environment creation fails

**Symptom:** errors during conda environment creation, often related to package resolution or network issues.

**Fixes:**

- **Check network access.** On some clusters, compute nodes cannot access the internet. Create conda environments from a node with network access first by running a dry run or a short test.
- **Use `--conda-prefix`** to share environments across projects and avoid repeated creation:

    ```bash
    snakemake --use-conda --conda-prefix ~/snparcher_envs ...
    ```

- **Clear corrupted environments.** If an environment was partially created, remove it and let Snakemake rebuild:

    ```bash
    rm -rf .snakemake/conda/<hash>
    ```

## Pipeline completes but output is unexpected

### Empty or very small VCF

**Possible causes:**

- **Reference genome mismatch.** If the reference genome is very different from your samples, most reads will not align, producing very few variants. Check mapping rates in the QC report.
- **Wrong ploidy.** If `variant_calling.ploidy` is set incorrectly, the caller may miss variants. Verify ploidy for your organism.
- **Heterozygosity prior too low.** The snpArcher default (0.005) may still be too low for highly diverse non-model organisms. Try raising `variant_calling.gatk.het_prior` to 0.01 or higher.

### QC report shows unexpected patterns

- **All samples cluster together in PCA but you expected structure:** this is normal for populations with low differentiation. Consider whether your k-means cluster count (`modules.qc.clusters`) is appropriate.
- **Depth drives the top PCs:** common in low-coverage datasets. Consider whether depth normalization or stricter depth filtering is warranted. See [Interpret QC reports](qc-report.md) for guidance.
- **Samples with very negative inbreeding coefficients:** likely cross-contamination. Consider excluding these samples.

### Postprocessed VCF is missing

- **Module not enabled.** Verify `modules.postprocess.enabled: true` in config.
- **No callable sites.** If `callable_sites.generate_bed_file: true` but both coverage and mappability are disabled, snpArcher warns and disables postprocessing. Enable at least one callable sites source.

## Resume after failures

Snakemake tracks completed steps automatically.
After fixing the issue, re-run the same command and Snakemake will pick up where it left off.

If a job failed partway through and left incomplete output files, add `--rerun-incomplete`:

```bash
snakemake \
  --snakefile /path/to/snpArcher/workflow/Snakefile \
  --directory /path/to/my_project \
  --use-conda \
  --cores 8 \
  --rerun-incomplete
```

## Getting help

If you cannot resolve an issue:

1. Search existing [GitHub Issues](https://github.com/harvardinformatics/snpArcher/issues). Your problem may already have a solution.
2. Open a new issue with:
    - The Snakemake command you ran
    - The relevant log output (Snakemake log, cluster log, and tool log)
    - Your `config.yaml` (with any sensitive paths redacted)
    - Your Snakemake and snpArcher versions (`snakemake --version` and `git log -1` in the snpArcher directory)

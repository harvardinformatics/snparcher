# How to add a custom module

snpArcher supports custom Snakemake modules that build on its outputs.

## Module requirements

To be included in snpArcher, a module must meet four criteria:

1. **Conda-defined environments.** The module's Snakefile must define all software dependencies using conda environment files, so that users do not need to install tools manually.
2. **Open source and documented.** The module must be freely distributed (e.g., on GitHub) with documentation sufficient for users to adapt it to their needs.
3. **Unit tested.** The module must include a test that runs on either the existing snpArcher test dataset or a module-specific minimal test dataset.
4. **Registered.** The module should be registered within the main snpArcher project so users can find it.

## Create the module

### Directory structure

A module lives in its own directory under `workflow/modules/`:

```text
workflow/modules/my_module/
├── Snakefile          # Snakemake rules
└── envs/
    └── my_tool.yaml   # Conda environment definition(s)
```

### Snakefile conventions

Your module's Snakefile should:

- Accept inputs from snpArcher's standard output paths (e.g., `results/vcfs/filtered.vcf.gz`, `results/bams/markdup/{sample}.bam`).
- Write outputs under a dedicated subdirectory (e.g., `results/my_module/`).
- Use conda environment files for all tool dependencies.
- Define a top-level target rule (e.g., `rule all`) that lists the module's final outputs.

Here is a minimal example:

```python
rule my_analysis:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
    output:
        result="results/my_module/output.tsv",
    conda:
        "envs/my_tool.yaml"
    shell:
        "my_tool --input {input.vcf} --output {output.result}"
```

### Conda environment file

Define the tools your module needs in a YAML file:

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - my_tool=1.2.3
```

## Register the module

To integrate your module with snpArcher's configuration system:

1. Add a toggle and any module-specific settings under `modules` in the config schema (`workflow/schemas/config.schema.yaml`).
2. Add the module import logic to the main `workflow/Snakefile` so it is conditionally loaded when enabled.
3. Add default values for module settings to `config/config.yaml`.

!!! note
    If you are developing a module for your own use rather than contributing it upstream, you can skip formal registration and simply invoke the module's Snakefile directly:

    ```bash
    snakemake \
      --snakefile workflow/modules/my_module/Snakefile \
      --use-conda \
      --configfile my_module_config.yaml \
      --cores 4
    ```

## Test the module

Write a test that runs the module on a small dataset:

1. Use the bundled snpArcher test data (in the `example/` directory) or provide a minimal test dataset with your module.
2. Verify that the module produces the expected output files and that they are non-empty.
3. Include the test in your module's repository so that others can validate it works.

## Contributing upstream

If you would like your module included in the main snpArcher distribution:

1. Open an issue on the [snpArcher GitHub repository](https://github.com/harvardinformatics/snpArcher/issues) describing your module.
2. Submit a pull request with the module code, tests, and documentation.
3. The snpArcher maintainers will review for compliance with the four criteria above.

For questions about module development, reach out via GitHub Issues. The maintainers are happy to discuss design decisions early in the process.

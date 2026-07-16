# How to filter and postprocess

The postprocess module re-filters the VCF after excluding flagged samples and applying site-level quality filters.

For the rationale behind each filtering step, see [Filtering philosophy](../explanation/filtering.md).

## Enable postprocessing

In your `config.yaml`, set the postprocess module to enabled:

```yaml
modules:
  postprocess:
    enabled: true
```

Then re-run snpArcher.
Snakemake detects that the postprocess rules are now required and runs only the new steps.

## Exclude samples

Sample exclusions are managed through a **sample metadata file**, separate from the main sample sheet.

### 1. Create the metadata file

Create a CSV file with at least `sample_id` and `exclude` columns:

```csv
sample_id,exclude
bird_1,false
bird_2,false
bird_3,true
bird_4,false
bird_5,true
```

Set `exclude` to `true` for any sample you want removed from the postprocessed VCF.
The `sample_id` values must match entries in your main sample sheet.

### 2. Point config to the metadata file

Add the `sample_metadata` path to your `config.yaml`:

```yaml
sample_metadata: "config/sample_metadata.csv"  # <-- change this
```

!!! tip
    The metadata file can also include optional columns used by other modules:

    - `lat`, `long` (numeric): geographic coordinates for the QC map panels.

## Configure site-level filters

The `modules.postprocess.filtering` block controls which sites are retained in the clean VCF:

```yaml
modules:
  postprocess:
    enabled: true
    filtering:
      contig_size: 10000  # <-- exclude contigs smaller than this (bp)
      maf: 0.01  # <-- minimum minor allele frequency
      missingness: 0.75  # <-- maximum fraction of missing genotypes per site
      exclude_scaffolds: "mtDNA,Y"  # <-- comma-separated scaffold names to remove
```

### Filter descriptions

| Filter | Default | Effect |
|---|---|---|
| `contig_size` | `10000` | SNPs on contigs/scaffolds this size or smaller are excluded. |
| `maf` | `0.01` | SNPs with minor allele frequency below this threshold are excluded. |
| `missingness` | `0.75` | SNPs where more than this fraction of genotypes are missing are excluded. |
| `exclude_scaffolds` | `"mtDNA,Y"` | SNPs on these scaffolds are excluded entirely. Comma-separated, no spaces. |

!!! tip
    For many organisms, you will want to exclude sex chromosomes and mitochondrial scaffolds.
    For example, for a ZW bird: `exclude_scaffolds: "Z,W,mtDNA"`.

## Intersect with callable sites

If callable sites are enabled (the default), the postprocess module automatically intersects the VCF with the callable sites BED file, restricting the clean VCF to regions that pass both coverage and mappability filters.

If you disabled callable sites, the postprocess module will warn and skip this intersection step.

To configure callable sites thresholds, see [How to configure your run](configure.md#configure-callable-sites).

## Output files

After postprocessing, the following files are produced in `results/postprocess/`:

| File | Description |
|---|---|
| `filtered.vcf.gz` | VCF after removing excluded samples and sites failing hard filters or with ref=N or AF=0. |
| `clean_snps.vcf.gz` | Biallelic SNPs from the filtered VCF, after applying missingness, MAF, callable-site, and scaffold-exclusion filters. |
| `clean_indels.vcf.gz` | Indels from the filtered VCF, after applying the same filters. |

The `clean_snps.vcf.gz` file is typically what you want for downstream population genomic analyses.

## Worked example

Suppose you ran snpArcher on 50 bird samples, reviewed the [QC report](qc-report.md), and decided to:

- Exclude 3 samples with very low depth (< 5x)
- Remove the Z chromosome and mitochondrial scaffolds
- Set a stricter MAF filter

1. Create `config/sample_metadata.csv`:

    ```csv
    sample_id,exclude
    bird_12,true
    bird_27,true
    bird_41,true
    ```

    (You only need to list the samples you want to exclude. You can optionally list all samples, setting `exclude` to `false` for those you want to keep.)

2. Update `config.yaml`:

    ```yaml
    sample_metadata: "config/sample_metadata.csv"

    modules:
      postprocess:
        enabled: true
        filtering:
          contig_size: 10000
          maf: 0.05  # <-- stricter MAF filter
          missingness: 0.75
          exclude_scaffolds: "Z,mtDNA"
    ```

3. Re-run snpArcher:

    ```bash
    snakemake \
      --snakefile /path/to/snpArcher/workflow/Snakefile \
      --directory /path/to/my_project \
      --use-conda \
      --cores 8
    ```

    Snakemake will only run the postprocessing steps, since the core pipeline outputs already exist.

## Next steps

- [Filtering philosophy](../explanation/filtering.md) for background on why these filters matter
- [Outputs reference](../reference/outputs.md) for a complete listing of all output files

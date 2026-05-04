# mkado Module Design Note

## Goal

Design a `snpArcher` module that uses [`mkado`](https://github.com/kr-colab/mkado) as the backend for McDonald-Kreitman analysis.

This note is for planning only. It records the current proposed shape so implementation can resume later without re-doing the design pass.

## Summary

`mkado` looks like a good fit for `snpArcher` because it can run directly from:

- an ingroup VCF
- an outgroup VCF
- a shared reference FASTA
- a GFF3 annotation with CDS features

The main integration problem is not the MK test itself. The main work is defining a stable contract between `snpArcher` outputs and `mkado` inputs.

The biggest constraints are:

- `mkado vcf` wants one multi-sample ingroup VCF plus one single-sample outgroup VCF
- `mkado` expects a GFF3, not a GTF
- `mkado` skips indels and multi-allelic sites unless we normalize upstream
- some existing `postprocess` filters in `snpArcher` are not appropriate as defaults for MK analysis because they alter the site-frequency spectrum

## Current `snpArcher` State

Relevant current workflow behavior:

- The default `all` target includes `results/vcfs/raw.vcf.gz`, callable-site outputs, optional postprocess outputs, and optional QC outputs.
- `call_variants` runs through `results/vcfs/filtered.vcf.gz`.
- The postprocess module currently starts from `results/vcfs/raw.vcf.gz`, not `results/vcfs/filtered.vcf.gz`.
- Sample metadata already supports `exclude` and `outgroup`.
- The reference workflow already normalizes the reference FASTA into `results/reference/{name}.fa.gz` and builds the `.fai`.
- There is currently no v2 module for MK analysis.
- There was an older removed `modules.mk` implementation based on `degenotate`, but it should be treated only as prior art, not as a contract.

## Recommended Module Name

Use `modules.mkado`, not `modules.mk`.

Reason:

- `modules.mk` was explicitly removed in v2
- `mkado` is a distinct backend with different assumptions from the old `degenotate` module
- a new name avoids accidental reuse of stale config expectations

## High-Level Module Contract

### Inputs from Main Workflow

The proposed main-workflow-to-module contract is:

- `samples`: main samples CSV
- `sample_metadata`: optional metadata CSV
- `ref`: normalized reference FASTA from `results/reference/{name}.fa.gz`
- `ref_fai`: matching FASTA index
- `vcf`: analysis source VCF, preferably `results/vcfs/filtered.vcf.gz`
- `callable_sites_bed`: optional callable BED from `results/callable_sites/callable_sites.bed`
- `gff`: resolved GFF3 annotation path

### Outputs from Module

Suggested public outputs:

- `results/mkado/ingroup.vcf.gz`
- `results/mkado/ingroup.vcf.gz.tbi`
- `results/mkado/outgroup.vcf.gz`
- `results/mkado/outgroup.vcf.gz.tbi`
- `results/mkado/results.tsv`
- `results/mkado/results.json`

Optional outputs, depending on mode:

- `results/mkado/volcano.png`
- `results/mkado/asymptotic_alpha.png`
- `results/mkado/genes_used.txt`
- `results/mkado/skipped_genes.tsv`

## Recommended Source VCF

Use `results/vcfs/filtered.vcf.gz` as the default source VCF for MKado, not `results/postprocess/clean_snps.vcf.gz`.

Reason:

- `postprocess` applies MAF and missingness filters intended for general downstream population-genetic use
- asymptotic and imputed MK analyses rely on allele-frequency information that should not be pre-distorted more than necessary
- `mkado` already has its own analysis-level frequency filtering options

`postprocess` is still useful as a source of implementation patterns:

- sample exclusion handling
- callable-region restriction
- contig/scaffold filtering
- SNP-only extraction
- removing SNPs that overlap indels

## Proposed Translation Pipeline

This is the proposed conversion path from `snpArcher` outputs to `mkado` inputs.

### 1. Resolve Sample Groups

From the main samples CSV and optional metadata:

- excluded samples = `exclude == true`
- outgroup samples = `outgroup == true`
- ingroup samples = all samples not excluded and not outgroup

### 2. Validate Sample Group Constraints

Recommended v1 rules:

- require at least one ingroup sample
- require exactly one outgroup sample
- fail if an outgroup sample is also excluded
- warn if metadata is absent and MKado is enabled

Why require exactly one outgroup in v1:

- `mkado vcf` expects a single-sample outgroup VCF
- allowing multiple outgroups would force `snpArcher` to invent a consensus or arbitration policy
- that policy should be explicit, tested, and probably deferred to a later version

### 3. Resolve Annotation

MKado needs a GFF3 with CDS features.

Preferred resolution order:

1. `reference.annotation.source` if provided
2. `modules.mkado.annotation` if a module-local override is allowed
3. auto-fetch from NCBI when `reference.source` is an accession and no annotation path is supplied
4. otherwise fail with a clear error

Important constraints:

- must be GFF3, not GTF
- CDS features must be present
- genes whose CDS length is not divisible by 3 will be skipped by MKado

## Proposed VCF Preparation Steps

Prepare dedicated MKado input VCFs in module-local rules.

### Ingroup/Outgroup Preparation

1. Start from `results/vcfs/filtered.vcf.gz`
2. Subset ingroup samples into `results/mkado/ingroup.raw.vcf.gz`
3. Subset the chosen outgroup sample into `results/mkado/outgroup.raw.vcf.gz`

### Normalize Variants

Normalize before running MKado:

- split multi-allelic sites with `bcftools norm -m -any`
- keep only biallelic SNPs
- drop indels
- drop SNPs overlapping indels, using the same basic pattern as `postprocess`

### Optional Region Restriction

If enabled and available:

- restrict both ingroup and outgroup VCFs to callable regions

This should probably default to `true`, but it should not be hard-required in v1.

### Optional Contig Filtering

Optionally:

- exclude named scaffolds/contigs
- exclude contigs below a minimum size threshold

These are reasonable defaults for cleaning up annotation and variant edge cases, but should remain configurable.

## Proposed Rule Flow

The module could be organized roughly like this:

1. `mkado_validate_inputs`
2. `mkado_resolve_annotation`
3. `mkado_write_sample_lists`
4. `mkado_subset_ingroup_vcf`
5. `mkado_subset_outgroup_vcf`
6. `mkado_prepare_callable_bed` (optional)
7. `mkado_normalize_ingroup_vcf`
8. `mkado_normalize_outgroup_vcf`
9. `mkado_filter_ingroup_vcf`
10. `mkado_filter_outgroup_vcf`
11. `mkado_run`
12. `mkado_collect_outputs`

This keeps the data-preparation logic explicit and makes it easier to test independently from the `mkado` invocation itself.

## Proposed Config Shape

This is the current recommended config shape.

```yaml
reference:
  name: "my_ref"
  source: "GCF_..."
  annotation:
    source: ""

modules:
  mkado:
    enabled: false
    outgroup_sample: ""
    source_vcf: "filtered"
    use_callable_sites: true
    filter:
      exclude_scaffolds: ""
      min_contig_size: 0
      max_missingness: null
    analysis:
      mode: "standard"
      per_gene: true
      gene_list: ""
      code_table: 1
      min_freq: null
      no_singletons: false
      bins: 10
      bootstrap: 100
      freq_cutoffs: [0.1, 0.9]
      workers: 0
      volcano: false
      plot_asymptotic: false
```

## Config Notes

### `reference.annotation.source`

This is the preferred place for the GFF3 path.

Reason:

- annotation is tied to the reference, not to the MK method
- other future annotation-dependent modules may need the same resource
- it avoids burying core reference metadata inside one downstream module

### `modules.mkado.outgroup_sample`

This should probably be optional, but used as an override.

Proposed behavior:

- if exactly one `outgroup=true` sample exists in metadata, use it
- if zero or multiple outgroup samples exist, require `outgroup_sample`
- fail if the named sample is not in the main sample sheet

### `modules.mkado.source_vcf`

Keep this small and explicit.

Suggested allowed values:

- `filtered`
- `raw`

Default should be `filtered`.

### `modules.mkado.analysis.mode`

Suggested allowed values:

- `standard`
- `asymptotic`
- `imputed`
- `alpha_tg`

Start by supporting one mode per run.

Running multiple modes from the same prepared inputs can be added later.

## Suggested Snakemake Import Pattern

If enabled, the main workflow would build a flat config dict and import the module similarly to QC and postprocess.

Example shape:

```python
_mkado_config = {
    "samples": config["samples"],
    "sample_metadata": config.get("sample_metadata", ""),
    "vcf": "results/vcfs/filtered.vcf.gz",
    "ref": REF_FILES["ref"],
    "ref_fai": REF_FILES["ref_fai"],
    "callable_sites_bed": "results/callable_sites/callable_sites.bed",
    "annotation": ...,
    **config["modules"]["mkado"],
}
```

## Environment and Dependency Notes

The module should use its own environment, not reuse the old removed `mk` env.

Likely dependencies:

- `mkado`
- `bcftools`
- `bedtools`
- `samtools`

Important note:

- MKado currently documents a Python 3.12 requirement
- the module env should pin a known working MKado release instead of following `main`

## What We Can Reuse from Existing `snpArcher`

Likely reusable logic:

- metadata parsing for `exclude` and `outgroup`
- callable BED generation from the main workflow
- reference normalization and indexing
- postprocess-style VCF restriction and SNP cleanup logic

Likely not reusable as-is:

- postprocess MAF filtering
- postprocess missingness filtering as a default
- the old `degenotate` module assumptions

## Open Questions

### 1. Should callable BED be mandatory?

Current recommendation:

- no, but default to using it if available

Reason:

- callable BED is useful for avoiding biased site comparisons
- making it mandatory would block MKado for gVCF-only or reduced workflows
- the right strictness may depend on how strongly we want MKado runs to reflect accessible coding sequence only

### 2. Should multiple outgroups be supported?

Current recommendation:

- no, not in v1

Reason:

- `mkado vcf` expects one outgroup VCF
- any multi-outgroup strategy would need explicit consensus rules

### 3. Should the module use `filtered` or `raw` VCF by default?

Current recommendation:

- use `filtered`

Reason:

- this keeps upstream caller hard filters
- it avoids importing postprocess population-genetic filters by default

### 4. Should module-local missingness filtering exist?

Current recommendation:

- allow it, but disable it by default

Reason:

- users may want a safeguard for pathological sites
- aggressive missingness filtering should not silently become part of the MK test definition

### 5. Where should annotation live in config?

Current recommendation:

- `reference.annotation.source`

Reason:

- it is reference metadata
- it will likely be useful beyond this one module

## Minimum Viable Implementation

A reasonable first implementation would do only this:

1. support exactly one outgroup sample
2. require a user-supplied GFF3 path
3. default to `results/vcfs/filtered.vcf.gz`
4. optionally restrict to callable regions
5. normalize to biallelic SNPs only
6. run `mkado vcf` in standard mode
7. emit TSV and JSON outputs

This would get a working backend in place without over-designing annotation retrieval or advanced mode orchestration.

## Good Follow-Up Work After MVP

- add asymptotic mode support
- add imputed mode support
- add volcano and asymptotic plot outputs
- allow auto-fetch of GFF3 for NCBI accession references
- add a dedicated test fixture with one ingroup population plus one outgroup sample
- decide whether skipped-gene reporting should be promoted to a public output
- decide whether future modules should share a common annotation-resolution helper

## Bottom Line

The cleanest plan is:

- create a new `modules.mkado` module
- treat `results/vcfs/filtered.vcf.gz` as the default source VCF
- generate dedicated ingroup/outgroup VCFs inside the module
- require or resolve a valid GFF3 annotation
- keep sample-group validation strict, especially around outgroup selection
- avoid inheriting postprocess MAF and missingness filtering by default

That gives `snpArcher` a clear MK backend contract without entangling the module with the old removed `mk` design.

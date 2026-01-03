# Initial V2 Brainstorming

Major goals of version 2 redesign are to:
- simplify and reorganize sample sheet and config
- normalize output paths to best practices
- remove unnecessary rule dependencies
- faciliate future moves to a cli and to add more options for each step

## Sample Sheet

The sample sheet currently controls too many things, including workflow behavior. 

The redesigned sample sheet should meet these requirements:

(1) A pure data manifest that captures only the minimum information needed to process each sample, specifically sample_id, input_type, and input fields, plus an optional library_id for cases where the same sample_id appears in more than one row. 

sample_id: defines the sample id in the final VCF and other locations; should have some formatting requirements (no spaces or other special characters, maximum length)
input_type: defines fastq, bam, gvcf, srr, though initially only fastq and srr are functional
library_id: if a sample_id appears on more than one row, the library design needs to be specified for dedup steps
input: path to fastq, srr accession, etc

(2) Other information that currently lives on the sample sheet should be removed. 

sample metaata (such as lat, lon, bioproject, organism, sex, population, etc) can be stored however the end user desires, although some aspects of this metadata may be used for downstream QC if available. The ultimate redesign of the QC steps needs to account for this.

reference genome is more appropriately stored in the config file; this means that multi-reference runs are no longer supported, but generally these were never a good idea to begin with. 

postprocessing information (exclusion/inclusion, mk output/ingroup) needs to be stored somewhere, perhaps in a postprocess-specific sample sheet

### Ref Genome - Considerations

Optionally, could add a ref_cache_dir and a ref_url to allow caching and allow alternate sources besides local fasta and NCBI. Have to be careful not to build a full-fledged reference manager, though (unless we intend to). 

Add a tiny registry abstraction: map ref_name -> {path, url, checksum}. Resolution order: explicit ref_path in manifest/config → ref_cache_dir/ref_name.{fna|fa|fasta} → download from url into cache → error. Always verify checksum and emit a single canonical location.


## Config

The config file currently has too many variables, organized by comments, and require commenting/uncommenting in order to use

A redesigned config should expose only parameters that are usefully changeable by the end user, and should simplify the process of updating the config for the end user. 

**Decision point**: There are some advantages to specifying resources in the config, instead of using set-resources, which can be confusing for users to edit when needed. Should we put the resources in the config? See for example [this](https://github.com/harvardinformatics/cactus-snakemake/blob/main/config-templates/config-template.yaml).

**Decision point**: Should "advanced options" be included at the end of the config?

**Decision point**: Flat with comments is often easier for end users to read; hierarchical is more logical if you are used to reading json or yaml but is more prone to formatting errors. I lean flat-with-comments. 

Minimally, we need to reorganize the config so that options don't require lots of commenting/uncommenting and other confusing issues.

## Sample Sheet and Config - Ideas

- Add explicit schema validation (JSON Schema/pydantic) for both sample sheet and config, with a snakemake --lint-manifest/--lint-config helper to catch missing fields, invalid characters in sample_id, and duplicate sample_id+library_id combos.

- Write a converter/backfill script that reads the current BioSample/Run/LibraryName/... sheet and emits the proposed minimal manifest plus a separate metadata/QC sheet, so existing users can adopt v2 without hand edits.

## Path Outputs

The output file names have evolved over time into a mishmash of different approaches. A clear, simple design goal for output names will help dramatically:

(1) Descriptive names for directories in each step
(2) Minimizing the use of wildcards, especially beyond sample_id, though some exceptions will be needed (e.g. for library_id in dedup logic)
(3) Better organization to facilitate the recommended directory setup in the docs

### Ideas

- Avoid deep nesting
- Define rule naming convention alongside paths (e.g., stage_action: q30_fastp, align_bwa, call_gatk_joint, qc_plink)
- Log/benchmark mirrors rule path: logs/{ref}/{rule}/..., benchmarks/{ref}/{rule}/
- Document which outputs are “public API” (stable) vs internal temps (important for module development)

## Rule Dependencies

It is currently difficult to support multiple entry points, and difficultt to support optional / variable paths through the DAG. 

Major design goal is to streamline rule dependencies as much as possible, and use target rules (inputs only) to help modularize the main workflow. 

## Design Architecture

Ultimately, we want to move towards a CLI interface that runs a variety of processing options, including allowing alternate aligners / snp callers and allowing easy run of partial steps (fastq -> bam; bam -> vcf; vcf -> qc; vcf -> postprocess, etc). Right now, the architectural design is brittle and does not easily facilitate this kind of modularity.  

Document pipeline “modes” that map to future CLI options (e.g., ingest-fastq, ingest-bam, call-intervals, qc-only, postprocess-only) so config keys and target rules align with modular entry points.

Add light import rules (copy/link/validate) that normalize into the same path layout, so downstream rules stay unchanged.

## Decisions Made

The following decisions have been finalized for v2:

### Sample Sheet Format

- **Format**: CSV (familiar to bio users, Excel-editable)
- **Columns**: `sample_id`, `input_type`, `input`, `library_id` (optional)
- **Paired FASTQ**: Semicolon-separated in `input` column: `r1.fq.gz;r2.fq.gz`
- **Input types**: `fastq`, `bam`, `gvcf`, `srr` (fastq/srr implemented first)
- **Sample ID validation**: `^[A-Za-z0-9_-]+$`
- **Library ID default**: Uses `sample_id` when not specified
- **Metadata**: Extra columns allowed (pass-through for QC module)

### Config Structure

- **Reference genome**: Nested under `reference:` with `name`, `path`, `accession`
- **Single reference per run**: Multi-reference runs no longer supported
- **Flat structure**: Keep flat-with-comments style for user readability
- **Schema validation**: JSON Schema via Snakemake's built-in `validate()`

### Wildcards

- **Simplified**: Only `{sample}` and `{library}` wildcards
- **No `{refGenome}`**: Removed from paths since single-reference per run

## Deferred Decisions

The following items are deferred for later discussion/implementation:

### Reference Genome Caching

The design doc proposes a `ref_cache_dir` and registry abstraction for reference management. This is deferred to avoid scope creep. Current implementation: simple `path` OR `accession` in config.

Future considerations:
- Cache directory for downloaded references
- Checksum verification
- Resolution order: config path → cache → download → error

### SRA RunSelector Integration

Need a helper script or documentation for converting SRA RunSelector output to v2 sample sheet format. This would ease migration for users reanalyzing existing SRA data.

### Postprocessing Sample Sheet

The postprocessing module needs sample inclusion/exclusion and MK ingroup/outgroup designation. Options:
- Separate `postprocess_samples.csv`
- Keep in main sample sheet with additional columns
- Defer until postprocessing module refactor

### Metadata Pass-through to QC

Extra columns in sample sheet (lat, long, population) should be available to QC module. Implementation details deferred until QC module refactor.

### CLI Interface

Future CLI will provide entry points for partial runs:
- `snparcher ingest-fastq`
- `snparcher call`
- `snparcher qc`
- `snparcher postprocess`

Design deferred until core workflow refactor is complete.

### Conda Environment Version Pinning

**Problem**: Conda env files (`workflow/envs/*.yml`) currently use unpinned versions
(`>=` instead of `==`) to support ARM (osx-arm64) testing. This is a tradeoff:
- Pinned versions: Better reproducibility, but many older bioinformatics packages
  lack ARM builds
- Unpinned versions: ARM-compatible, but builds may differ across platforms/time

**Recommended solution: Snakemake Wrappers**

[Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/) bundle tool + conda env
together, with versions managed by wrapper maintainers. Benefits:
- Wrapper maintainers handle cross-platform version compatibility
- Wrappers are automatically tested across platforms
- Reduces maintenance burden for snpArcher
- Provides standardized, well-tested implementations

**Migration plan**: Replace custom shell commands with wrappers where available:
- `bio/bwa/mem` for alignment
- `bio/fastp` for read filtering
- `bio/gatk/haplotypecaller` for variant calling
- `bio/samtools/faidx` for reference indexing
- etc.

This is tracked in `TODO.md` under "Conda Environment Versions".

# snpArcher v2 TODO

## Completed

- [x] Set up development environment with pixi (`pyproject.toml`)
- [x] Create AGENTS.md with development guidelines
- [x] Set up Material for MkDocs with Diataxis structure
- [x] Create JSON schemas for sample sheet and config validation
- [x] Document sample sheet specification (`docs/reference/sample-sheet.md`)
- [x] Document config options (`docs/reference/config.md`)
- [x] Define v2 sample sheet format (CSV with `sample_id`, `input_type`, `input`, `library_id`)
- [x] Define v2 config structure (`reference:` nested config)
- [x] Define v2 wildcard structure (`{sample}`, `{library}` only)
- [x] Implement new sample sheet parsing with validation in `rules/common.smk`
- [x] Reorganize rules into stage directories (`align/`, `call/`, `filter/`, `stats/`)
- [x] Refactor rules to use new wildcard structure (`{sample}`, `{unit}`)
- [x] Update main Snakefile for new structure
- [x] Create v2 test configs (`.test/fastq-to-vcf/`)
- [x] Run full integration test locally (ARM) - verified working
- [x] Set up pixi-based CI with platform-specific tasks
- [x] Update GitHub Actions workflow to use pixi

## High Priority

- [ ] Remove old v1 rule files after confirming v2 CI passes:
  - `workflow/rules/fastq2bam.smk`
  - `workflow/rules/bam2vcf_gatk.smk`
  - `workflow/rules/bam2vcf_gatk_intervals.smk`
  - `workflow/rules/intervals.smk`
  - `workflow/rules/sentieon.smk`
  - `workflow/rules/sumstats.smk`
  - `workflow/rules/cov_filter.smk`
  - `workflow/rules/mappability.smk`
- [ ] Remove `workflow/snparcher_utils/` after refactor complete
- [ ] Replace custom helper functions with Snakemake semantic helpers where possible
- [ ] Add Snakemake wrappers where available (see Conda Environment Versions below)

## Conda Environment Versions

**Current state:** Conda env files (`workflow/envs/*.yml`) use unpinned versions (`>=`)
to support ARM (osx-arm64) testing. This is not ideal for reproducibility.

**Why this matters:**
- Pinned versions ensure reproducible builds across platforms and time
- Unpinned versions may resolve differently on different platforms/dates
- ARM builds often lack older package versions available on x86_64

**Solution: Migrate to Snakemake Wrappers**
[Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/) bundle tool + conda env
together, with versions managed upstream. Benefits:
- Wrapper maintainers handle version compatibility
- Automatically tested across platforms
- Reduces maintenance burden for snpArcher

**Rules to migrate:**
- [ ] `align_bwa_map` -> `bio/bwa/mem` wrapper
- [ ] `fastq_fastp` -> `bio/fastp` wrapper
- [ ] `call_haplotypecaller_*` -> `bio/gatk/haplotypecaller` wrapper
- [ ] `reference_index` (bwa index) -> `bio/bwa/index` wrapper
- [ ] `reference_index` (samtools faidx) -> `bio/samtools/faidx` wrapper

## Design Issues to Address

- [ ] **Sumstats coupling** - stats rules force BAM generation even when starting from GVCF
  - Consider: separate target rule, conditional on input_type, or config flag to skip
  - Current design prevents efficient "GVCF-only" entry point
- [ ] **Module compatibility** - qc/postprocess/mk/trackhub modules need updates for v2
  - Currently disabled in main Snakefile
  - Need to update module common.smk files for v2 sample sheet format

## Deferred / To Discuss

- [ ] Reference genome caching system (`ref_cache_dir`, registry abstraction)
- [ ] SRA RunSelector integration / helper script for sample sheet generation
- [ ] Postprocessing sample sheet design (include/exclude, MK ingroup/outgroup)
- [ ] Metadata pass-through to QC module (lat, long, population)
- [ ] CLI interface design
- [ ] Migration script for v1 -> v2 sample sheet conversion

## Documentation

- [ ] Tutorial: First variant calling run
- [ ] How-to: Create sample sheet
- [ ] How-to: Run on SLURM cluster
- [ ] Explanation: Workflow architecture
- [ ] Migrate remaining content from `old-docs/`

## New v2 Rule Files

Created during refactor:
- `workflow/rules/common.smk` - Rewritten with v2 sample sheet parsing
- `workflow/rules/fastq.smk` - Updated for v2 wildcards
- `workflow/rules/reference.smk` - Updated for single reference (symlink for bgzip)
- `workflow/rules/align/bwa.smk` - BWA alignment rules
- `workflow/rules/align/sentieon.smk` - Sentieon alignment rules
- `workflow/rules/call/gatk.smk` - GATK non-interval calling
- `workflow/rules/call/gatk_intervals.smk` - GATK interval-based calling
- `workflow/rules/call/intervals.smk` - Interval generation checkpoints
- `workflow/rules/call/sentieon.smk` - Sentieon calling rules
- `workflow/rules/filter/coverage.smk` - Coverage filtering
- `workflow/rules/filter/mappability.smk` - Mappability filtering
- `workflow/rules/stats/sumstats.smk` - Summary statistics
- `workflow/Snakefile` - Main workflow entry point

## Test Configs

- `.test/fastq-to-vcf/config/config.yaml` - v2 format config
- `.test/fastq-to-vcf/config/samples.csv` - v2 format sample sheet

## CI Configuration

All platforms use `--use-conda` with unpinned conda environment versions.
GitHub Actions workflow (`.github/workflows/test.yaml`) uses `prefix-dev/setup-pixi`.

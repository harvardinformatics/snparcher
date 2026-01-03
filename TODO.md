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

## High Priority

- [ ] Implement new sample sheet parsing with validation in `rules/common.smk`
- [ ] Refactor rules to use new wildcard structure (`{sample}`, `{library}`)
- [ ] Move reference genome config from sample sheet to config.yaml
- [ ] Remove `workflow/snparcher_utils/` after refactor complete

## Medium Priority

- [ ] Reorganize rules into stage directories (`align/`, `call/`, `qc/`, `postprocess/`)
- [ ] Replace custom helper functions with Snakemake semantic helpers
- [ ] Add Snakemake wrappers where available
- [ ] Migrate docs content from `old-docs/` to Diataxis structure

## Deferred / To Discuss

- [ ] Reference genome caching system (`ref_cache_dir`, registry abstraction)
- [ ] SRA RunSelector integration / helper script for sample sheet generation
- [ ] Postprocessing sample sheet design (include/exclude, MK ingroup/outgroup)
- [ ] Metadata pass-through to QC module (lat, long, population)
- [ ] CLI interface design

## Documentation

- [ ] Tutorial: First variant calling run
- [ ] How-to: Create sample sheet
- [ ] How-to: Run on SLURM cluster
- [ ] Explanation: Workflow architecture

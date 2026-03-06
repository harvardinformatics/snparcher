# v2 Modularity and Maintainability Plan

This plan is aligned to the goals in `version-2-design.md`:
- simplify sample/config interfaces
- normalize outputs and rule naming
- reduce brittle dependencies
- support future CLI "modes" and partial pipeline execution

## High-priority findings from this review

1. P1: v1/v2 drift is high across config, docs, and modules.
   - `docs/setup.md` and `docs/examples.md` document legacy fields (for example `BioSample`, `Run`, `final_prefix`) that do not match the current v2 sample/config contract.
   - `workflow/modules/*` still depend on `snparcher_utils`, which is not present in this repo, and are not wired into the v2 main workflow.

2. P1: Configuration contract is fragmented.
   - Defaults are split between `workflow/rules/common.smk` and `workflow/schemas/config.schema.yaml`.
   - Keys currently used in workflow logic are not fully represented in schema defaults/structure (for example `callable_sites.coverage.enabled`), increasing risk of runtime key errors and silent drift.

3. P1: Resource profile thread keys are mostly stale.
   - `workflow-profiles/default/config.yaml` contains many `set-threads` keys that do not match current rule names, so user overrides often do not apply.

4. P2: Variant-calling logic is duplicated.
   - Filtering expressions and mapfile/import/genotype structure are repeated across `gatk.smk`, `gatk_intervals.smk`, and `sentieon.smk`, making behavior changes error-prone.

5. P2: Rule files mix orchestration and implementation details.
   - Long `run:` blocks and helper functions embedded in rule files reduce readability and testability.
   - Snakemake lint flags this pattern repeatedly.

## Recommended phased execution

## Phase 0: Define and lock contracts (first)
- [ ] Choose one canonical v2 data contract:
  - sample manifest columns
  - config keys and nesting
  - stable output paths ("public API")
- [ ] Make schema strict (`additionalProperties: false` where appropriate).
- [ ] Ensure every key used in workflow code has schema coverage and defaults.
- [ ] Add one migration document: v1 -> v2 sample/config mapping.

Acceptance criteria:
- `snakemake --lint` clean for config/sample contract issues.
- Running with missing optional keys does not raise key errors.

## Phase 1: Extract workflow core library
- [ ] Create `workflow/lib/` for reusable Python helpers:
  - `manifest.py` (sample parsing, library grouping, validation)
  - `config.py` (single source of defaults + typed accessors)
  - `paths.py` (central path builders)
  - `intervals.py` (checkpoint interval list handling)
- [ ] Refactor `workflow/rules/common.smk` to consume `workflow/lib/*` instead of owning most logic.
- [ ] Move long `run:` blocks into `workflow/scripts/` and unit-test them.

Acceptance criteria:
- No duplicated manifest/config parsing logic in rule files.
- Rule files are mostly declarative (inputs/outputs/resources/shell).

## Phase 2: Unify DAG paths and add explicit modes
- [ ] Define target rules for v2 modes (`setup`, `ingest`, `call`, `qc`, `callable-sites`, `full`).
- [ ] Collapse duplicated variant-calling steps into shared helpers/rules:
  - one source for filter expressions
  - one source for db mapfile generation
  - one source for interval list expansion logic
- [ ] Keep tool-specific implementations (GATK/Sentieon) behind clear interfaces.

Acceptance criteria:
- Mode targets can run independently with predictable outputs.
- GATK interval/non-interval behavior differences are data-driven, not copy-pasted files.

## Phase 3: Resolve legacy module strategy
- [ ] Decide per module (`qc`, `postprocess`, `mk`, `trackhub`):
  - port to v2 contracts now, or
  - archive as legacy with clear documentation
- [ ] If ported, introduce explicit module API:
  - required upstream artifacts
  - module-specific config namespace
  - tests per module

Acceptance criteria:
- No module references missing code (`snparcher_utils`).
- Module docs and tests reflect actual runnable behavior.

## Phase 4: Align docs, profiles, and CI guardrails
- [ ] Rewrite setup/examples docs for v2 contract only.
- [ ] Update `workflow-profiles/default/config.yaml` with current rule names.
- [ ] Add CI checks:
  - dry-run tests for mode targets
  - schema validation tests
  - lint check
  - docs config examples validated against schema

Acceptance criteria:
- New users can run from docs without guessing.
- Profile overrides affect real rules.

## Suggested PR sequence

1. Contract and schema lock (Phase 0)
2. `workflow/lib` extraction + script moves (Phase 1)
3. Variant path unification + mode targets (Phase 2)
4. Module port/archive decisions (Phase 3)
5. Docs/profile/CI alignment (Phase 4)

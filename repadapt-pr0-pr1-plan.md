# Implementation plan: PR 0 (fixture + golden output) and PR 1 (`tool: repadapt`)

Scope: the first two deliverables in `agents-plan.md` section 3. The design in `agents-plan.md` section 2 is unchanged. This plan fills in the file-level details, fixes places where the handoff doc and the code disagree, and orders the work.

**Order:** step 0 (dev setup) → PR 0 → PR 1.
- PR 1's code and dry-run tests don't depend on PR 0 and can start in parallel.
- Only PR 1's golden test needs PR 0's outputs.

## Where the handoff doc and the code disagree

| Handoff doc | What the code or this box actually has | What this plan does |
|---|---|---|
| Refers to `APPLY_GATK_HARD_FILTERS` | `common.smk:640` and `Snakefile:33` use `APPLY_HARD_FILTERS`. Separately, `common.smk:628` forces `GENERATE_FILTERED_VCF = False`, with a warning, for every caller outside the GATK family. So `rule all` would never request `filtered.vcf.gz` for `repadapt`. | Rename to `APPLY_GATK_HARD_FILTERS` (2 references), and exempt `repadapt` from the forced-off block. |
| Per-region commands can't be checked in a dry run; use `workflow_source` | `test_bcftools_long_contig_dry_run_uses_csi_indexes` pre-seeds `results/vcfs/regions/regions.tsv` and dry-runs one region target, and the full command then appears. | Assert the exact `mpileup \| call` command from dry-run output. |
| Edit `workflow-profiles/slurm/config.yaml` | Only `workflow-profiles/default/` exists. | Edit the default profile only. |
| Golden run needs Apptainer | This box (holybioinf) has **SingularityCE 4.5.1** and no `apptainer` binary. Nextflow and pixi aren't installed, and the system Java is 1.8. | `golden.config` uses whichever engine is installed (`singularity {}` here). Step 0 installs Nextflow and pixi. |
| `-C golden.config` replaces RepAdapt's config | Checked at `2077f6f`: `nextflow.config` holds only container paths and `apptainer { enabled = true }`. CPUs are set in the process files. `-C` is a top-level option: `nextflow run main.nf -C ...` fails with "Unknown option: -C". | `-C` loses nothing. Keep it, as `nextflow -C golden.config run main.nf`. `agents-plan.md` is fixed. |
| (not covered) | `bcftools_concat_regions` writes `temp(RAW_VCF)`. That's safe for bcftools only because `RAW_VCF` is its final VCF. For `repadapt`, `call_variants` → `filtered.vcf.gz` would delete the raw VCF after filtering. | `repadapt`'s raw VCF is **not** temp. |
| (not covered) | Newer htslib warns `MQ should be declared as Type=Float` when it reads 1.16's Integer MQ (seen with bcftools 1.24). Postprocess and QC use bcftools ≥ 1.19. | The warning is harmless. Document it; don't "fix" it. |
| (not covered) | Fixture BAMs in `tests/data/fixtures/results/bams/markdup/` have only `.bai`. Linking the directory makes `index_bam_csi` write `.csi` files into the fixtures. | The unit test links the BAM files one at a time. |

Already checked with bcftools 1.24 on a small hand-made VCF:
- The two-pass soft filter gives `PASS`, `AllHomAlt`, `LowMQ` and `AllHomAlt;LowMQ`, with no `PASS;LowMQ`.
- Its PASS set equals RepAdapt's one-pass `-e 'AC=AN || MQ < 30'`, including records with MQ missing.

It gets re-checked on 1.16 (step 1.1 and the golden test).

## Step 0: dev environment (one-time, this box)

1. **pixi** (user-level install into `~/.pixi`), then `pixi install -e dev`. All test commands below use `pixi run -e dev ...`.
2. **Nextflow 25.10.2** in its own conda env outside the repo: `mamba create -p $GOLDEN_WORK/nf-env -c conda-forge -c bioconda nextflow=25.10.2`. This brings a Java 17+ JDK; don't use the system Java 1.8.
3. **`GOLDEN_WORK`:** a directory outside the repo for the SIF images (about 2 GB), the Nextflow work dirs and the run outputs. Proposed: `/n/holylfs05/LABS/informatics/Lab/projects/tsackton/snparcher-dev/repadapt-golden/`.
4. Start each branch from an updated `feat/repadapt` (`git fetch && git checkout feat/repadapt && git merge --ff-only origin/feat/repadapt`).

## PR 0: fixture and golden output (`feat/repadapt-fixture`)

### 0.1 `tests/data/repadapt/make_fixture.py`

Pure Python and seeded (`random.Random(SEED)`). Writes gzip with `mtime=0` so reruns are byte-identical. It writes `reference.fasta`, `genes.gff` and `fastq/S{1..4}_{1,2}.fastq.gz` under `tests/data/repadapt/`.

Reference, about 25 kb:
- `ctgA`, about 15 kb, with an exact ~1 kb repeat. Reads there get MAPQ 0, are removed by `view -q 10`, and leave orphaned mates behind.
- `ctgB`, about 8 kb.
- `ctgC`, about 1.5 kb.
- **Added:** a near-identical 2 kb copy of a `ctgA` segment, placed in `ctgB`, with one difference every 500 bp. The exact repeat alone won't give LowMQ sites, because its MAPQ-0 reads are filtered out before calling.
  - A 2–4% diverged copy (the first idea) doesn't work either: bwa's paired-end scoring raises those reads to MAPQ ≈ 57.
  - A pair that spans exactly one difference gets MAPQ ≈ 27 (`raw_mapq(5) − 3`). That survives `-q 10` but keeps INFO/MQ under 30. Measured: MAPQ 27 or 0 across the copy.

Variants per sample (four diploid samples):
- 34 ordinary SNPs, het and hom-alt.
- 16 indels of 1–30 bp. 6 of them are one repeat unit inside a planted short tandem repeat; each indel has a SNP 8–20 bp downstream.
- **Added, so every filter outcome is exercised:**
  - 3 SNPs hom-alt in all four samples (→ `AllHomAlt`)
  - 3 SNPs inside the near-identical copy (→ `LowMQ`)
  - 1 SNP that is both (→ `AllHomAlt;LowMQ`)
- **Measured realignment effect:** GATK3 realigns only 1–6 reads per sample and changes one record (`BQBZ`). That's typical of bwa-mem input: gaps are already placed well, and soft-clipped ends aren't realigned. Getting more would take more reads than the 1,000-pair fastp limit allows. **PR 3 should therefore check directly that realignment ran (OC tags in the final BAMs)**, not rely on VCF differences.

Reads:
- 2×150 bp, insert 350 ± 50, about 0.5% errors, a realistic quality profile, about 5% duplicate pairs.
- **≤ 1,000 pairs per sample** (about 12× depth). This keeps fastp deterministic.
- **Added:** about 3% short-insert pairs (100–140 bp) with adapter read-through. This does nothing in PR 1, but PR 3's golden test then catches a stray `--detect_adapter_for_pe`.
- Names like `@SIM.<n>`, identical for both mates. That avoids fastp's polyG triggers and Picard's tile/x/y parsing.

Other outputs:
- `genes.gff`: a few `gene` lines. RepAdapt requires it; we don't use it.

### 0.2 `tests/repadapt/normalize_vcf.py`

One pure-Python normalizer, used both to produce the golden files and inside the tests, so both sides normalize the same way.
- Input: a VCF (plain or bgzipped; Python's `gzip` reads BGZF), the reference FASTA for contig order, and an optional `--pass-only`.
- It reorders the sample columns into sorted order and sorts records by (contig order from the FASTA, POS, REF, ALT).
- It prints CHROM through the last sample column, tab-separated, with no header.
- Callable as a function (tests load it with `importlib`) and as a CLI (`make_golden.sh`).

### 0.3 `tests/repadapt/make_golden.sh WORKDIR`

1. Check for Linux, `apptainer` or `singularity`, and Nextflow ≤ 25.10.2 on Java 17+.
2. Clone RepAdapt and check out `2077f6f311650c24ccc0df8c4088e616e624f5d8`.
3. `curl` the seven depot SIFs (`agents-plan.md` 5.1) into `WORKDIR/images/` and record their sha256.
4. Write `golden.config`:
   - a `withName:` → absolute-path SIF map for every process in the table in `agents-plan.md` section 3, PR 0
   - `singularity { enabled = true; autoMounts = true }`, or `apptainer {…}` when that binary is present
5. Run it twice, each run with its own `-work-dir` and `--outdir`:
   `LC_ALL=C nextflow run main.nf -C golden.config --ref_genome $FIX/reference.fasta --gff_file $FIX/genes.gff --reads "$FIX/fastq/*_{1,2}.fastq.gz" --outdir ...`
6. **Added: collect RepAdapt's raw calls.** Find the pre-filter `variants_chr_*.vcf.gz` files in the work dir and concatenate them in FASTA order with the bcftools SIF. This lets PR 1 test our `raw.vcf.gz` against RepAdapt's raw calls, not only the PASS set.
7. Check determinism across the two runs:
   - The normalized final and raw records are identical.
   - `samtools view` records (no header) of each realigned BAM are identical, using the samtools SIF.
8. **Check the fixture exercises the filter.** Golden raw minus golden final must include at least one `AC=AN` record and at least one `MQ<30` record, and final must have a reasonable number of PASS SNPs. If not, adjust the fixture parameters and rerun.
9. Copy into `tests/data/repadapt/golden/`:
   - `final_variants.vcf.gz`
   - `final_records.tsv`
   - `raw_records.tsv`
   - `S{1..4}_sorted_RG_dedup_realigned.bam`
   - `PROVENANCE.md`: RepAdapt commit, image URLs and sha256, Nextflow and Java versions, container engine and version, the exact command, date, and fixture seed

Keep the total under about 2 MB. At about 12× depth the four BAMs should come to roughly 0.5 MB.

### 0.4 Sample sheets and tests

- `tests/sample_sheets/repadapt_fastqs.csv`: S1–S4, fastq, one row each, `mark_duplicates: true`.
- `tests/sample_sheets/repadapt_golden_bams.csv`: S1–S4 with `input_type: bam`, pointing at the golden BAMs. Sample IDs must equal the BAMs' `SM` tags (RepAdapt sets `SM` to the file stem, so `S1`…`S4`).
- `tests/configs/repadapt.yaml` goes in PR 1, where it's first used, because `tool: repadapt` fails schema validation until then.
- Unmarked test `test_repadapt_fixture_is_reproducible` (runs in the dry-run CI group): regenerate the fixture into a temp dir and compare sha256 values with the committed files.

**PR 0 is done when:** the fixture is reproducible, the golden output is identical across two RepAdapt runs, the filter-coverage check passes, and the files fit in the size budget.

## PR 1: calling model (`feat/repadapt-calling`)

### 1.1 Pinned env (do first; `agents-plan.md` 7.1 and 7.7)

- `workflow/envs/repadapt/bcftools.yaml`: channels conda-forge and bioconda; `bcftools=1.16`.
- `workflow/envs/repadapt/bcftools.linux-64.pin.txt`: the Appendix A bcftools list, verbatim.
- Check:
  - `conda create -p <tmp> --file <pin>` succeeds
  - `bcftools --version` reports 1.16 with htslib 1.16
  - the soft-filter check gives the same result on 1.16
  - a Snakemake run logs that it used the pin file

  Each of these also becomes a test assertion (1.6) so it doesn't depend on one-off probing.

### 1.2 Shared region checkpoint

New `workflow/rules/variant_calling/contig_regions.smk`. Move these verbatim out of `bcftools.smk`:
- `bcftools_call_input`
- `checkpoint bcftools_regions` (name kept, because tests assert it)
- `_read_bcftools_regions`, `_get_bcftools_regions_file`, `get_bcftools_region_ids`, `get_bcftools_region_name`

Only one message changes: "for bcftools caller" becomes caller-neutral. `bcftools.smk` keeps `get_bcftools_region_vcfs` and `_tbis` and its two rules. The Snakefile includes `contig_regions.smk` in both the `bcftools` and `repadapt` branches.

### 1.3 `workflow/rules/variant_calling/repadapt.smk`

All three rules use `conda: "../../envs/repadapt/bcftools.yaml"`, plus a log and a benchmark like the bcftools rules.

**`repadapt_call`**
- Input: `unpack(bcftools_call_input)` and `regions_tsv`.
- Output: `temp("results/vcfs/regions/repadapt/{region_id}.vcf.gz")`, with `wildcard_constraints: region_id="L\d{6}"`.
- `threads: 1`.
- No per-region index; plain concat doesn't need one.
- Shell:
  ```
  bcftools mpileup -Ou -f {input.ref} -r {params.contig} -q 10 -I -a FMT/AD,FMT/DP {input.bams} 2> {log} \
    | bcftools call -G - -f GQ -mv --ploidy {params.ploidy} -Oz -o {output.vcf} - 2>> {log}
  ```
  No `-Q`, `-d` or `--threads`.

**`repadapt_concat_regions`**
- Input: the region VCFs for `get_bcftools_region_ids` (sorted IDs, which is `.fai` order).
- Output: `RAW_VCF` and `RAW_VCF_INDEX`, not temp.
- A `run:` block writes the list file in Python, which avoids argv limits on fragmented references, then runs:
  - `bcftools concat -f list -Oz -o RAW`
  - `bcftools index BCFTOOLS_INDEX_ARGS RAW`

**`repadapt_filter`**
- Defined only `if APPLY_REPADAPT_FILTER`, mirroring how `hard_filters.smk` is gated.
- Output: `FILTERED_VCF` and `FILTERED_VCF_INDEX`.
- Shell:
  ```
  bcftools filter -s AllHomAlt -e 'AC=AN' -Ou {input.vcf} 2> {log} \
    | bcftools filter -m + -s LowMQ -e 'MQ<30' -Oz -o {output.vcf} - 2>> {log}
  bcftools index {params.index_args} {output.vcf} 2>> {log}
  ```

### 1.4 `common.smk` and `Snakefile`

```python
REPADAPT_CALLER = VARIANT_TOOL == "repadapt"
FILTERING_CALLER = GATK_LINEAGE_CALLER or REPADAPT_CALLER
if GENERATE_FILTERED_VCF and not FILTERING_CALLER:   # was: not GATK_LINEAGE_CALLER
    logger.warning(...)                              # text: bcftools/deepvariant only
    GENERATE_FILTERED_VCF = False
APPLY_GATK_HARD_FILTERS = GATK_LINEAGE_CALLER and GENERATE_FILTERED_VCF
APPLY_REPADAPT_FILTER = REPADAPT_CALLER and GENERATE_FILTERED_VCF
FINAL_VCF = FILTERED_VCF if (APPLY_GATK_HARD_FILTERS or APPLY_REPADAPT_FILTER) else RAW_VCF
```

Other `common.smk` changes:
- Add `"repadapt"` to the gVCF-rejection set (line 472).
- New check: `tool: repadapt` requires ploidy 1 or 2, because of bcftools 1.16's `--ploidy` aliases. It sits next to the gVCF check.
- Add `repadapt` to the list of supported backends in the long-contig error message.
- New warning, placed after `POSTPROCESS_ENABLED` is resolved (about line 716): if `REPADAPT_CALLER and POSTPROCESS_ENABLED and POSTPROCESS_SPLIT_BY_TYPE`, warn that `clean_indels.vcf.gz` will be empty because `-I` skips indels.

`Snakefile` changes:
- An `elif VARIANT_TOOL == "repadapt":` branch.
- `if APPLY_GATK_HARD_FILTERS:`.
- The `rule all` comment and the `call_variants` docstring now mention repadapt.

Config, schema and profile:
- **Schema:** add `repadapt` to the tool enum. Rewrite the `generate_filtered_vcf` description. Add a description on the `bcftools` block saying it has no effect under `repadapt`.
- **`config/config.yaml`:** update the tool comment and the `generate_filtered_vcf` comment.
- **Default profile:** add `repadapt_call: 1` to `set-threads`, to state the single-threaded choice explicitly.

### 1.5 CI and pixi

- **`pyproject.toml`:** add a second invocation to `setup-test-envs` with `--config samples=config/samples.csv 'variant_calling={tool: repadapt}'`. The repadapt env is then built even though the per-region jobs sit behind the checkpoint, because `repadapt_concat_regions`/`repadapt_filter` are in the DAG and share the env.
- **`.github/workflows/test.yaml`:** add `'workflow/**/*.pin.txt'` to the three `hashFiles(...)` cache keys (lines 68, 94, 137).

### 1.6 Tests

Test names avoid the substrings `qc`, `postprocess` and `metadata`, so all of these land in the core CI groups.

**Dry-run tests (`tests/tests.py`)**

1. `test_repadapt_dry_run`:
   - `call_variants` schedules `bcftools_regions`, `repadapt_concat_regions` and `repadapt_filter`, and `filtered.vcf.gz` is the target.
   - It doesn't schedule `variant_filtration` or `bcftools_concat_regions`.
   - There's no "generate_filtered_vcf … disabling" warning.
2. `test_repadapt_region_command_flags`:
   - Pre-seed `regions.tsv`, link fixtures, and dry-run `results/vcfs/regions/repadapt/L000000.vcf.gz`.
   - Assert the exact `mpileup -Ou … -q 10 -I -a FMT/AD,FMT/DP` and `call -G - -f GQ -mv --ploidy 2` strings.
   - Assert `-Q`, `-d` and `--threads` are absent.
3. `test_repadapt_filter_commands`: the exact two-pass filter string and `bcftools index -f -t results/vcfs/filtered.vcf.gz`.
4. `test_repadapt_long_contig_uses_csi`: `write_long_contig_config(..., "repadapt")` gives `.csi` outputs and `-f -c`.
5. `test_repadapt_generate_filtered_vcf_false_uses_raw`: no `repadapt_filter`, and `call_variants` resolves to `raw.vcf.gz`.
6. `test_repadapt_rejects_unsupported_ploidy`: ploidy 3 fails with a clear message.
7. `test_repadapt_split_by_type_warns_empty_indels`: the warning appears with the postprocess config, and not with `split_by_type: false`.
8. Add `"repadapt"` to the parametrize list of `test_gvcf_input_rejected_for_new_callers`.
9. bcftools regression, in `test_bcftools_dry_run`:
   - `repadapt_call` isn't scheduled.
   - `bcftools.smk` doesn't contain `-G -`.

   The existing hard-filter tests must pass unchanged:
   - `test_generate_filtered_vcf_false_skips_hard_filtering`
   - `test_long_contig_hard_filters_*`
   - `test_short_contig_hard_filters_*`
   - `test_generate_filtered_vcf_auto_disabled_for_non_gatk`

**Unit test (`tests/unit_tests.py`)**

10. `test_repadapt_call_and_filter`:
    - Link `results/reference` and each fixture BAM file individually, not the directory. Use a repadapt config copy, and target `filtered.vcf.gz`.
    - The header has `##bcftools_callVersion=1.16+htslib-1.16`, which also catches a silent fallback from the pin file.
    - FORMAT AD, DP and GQ are declared.
    - There are no indel records.
    - FILTER values are all in {PASS, AllHomAlt, LowMQ, AllHomAlt;LowMQ}.
    - PASS records have AC<AN and MQ≥30 (or MQ missing).
    - `raw.vcf.gz` still exists afterwards.

**Full-run tests (`tests/tests.py`, need PR 0)**

11. `test_repadapt_golden_bams`: uses `tests/configs/repadapt.yaml` with `repadapt_golden_bams.csv`.
    - Normalized `-f PASS` records equal `final_records.tsv`.
    - Normalized raw records equal `raw_records.tsv`.
    - Each non-PASS label matches its own evaluation of `AC=AN` / `MQ<30`.
    - Both labels occur.
12. `test_full_pipeline_repadapt`: the default mapping pipeline on `get_samples_file()`, targets `all` and `call_variants`. Raw, filtered, `qc_report.tsv` and the callable-sites outputs all exist. This mirrors `test_full_pipeline_bcftools`.

### 1.7 Docs

- `reference/config-schema.md`: the tool enum, the `generate_filtered_vcf` text, and a note on the bcftools block.
- `explanation/variant-calling.md`, new RepAdapt section:
  - the command, and how it differs from RepAdapt's literal command (`-q 10` for `-q 5`, `--ploidy`)
  - why 1.16 is pinned
  - `-G -` semantics
  - no indels
  - Linux-only
- `explanation/filtering.md`: the soft filters and the PASS equivalence.
- `reference/outputs.md`.
- `how-to/configure.md`: the tool table.
- `reference/changelog.md`.

### 1.8 Before opening PR 1

```bash
pixi run -e dev pytest -v tests/tests.py --dry-run-only
pixi run -e dev setup-test-envs
pixi run -e dev pytest -v tests/unit_tests.py --conda-prefix $(pwd)/.snakemake/conda -k repadapt
pixi run -e dev pytest -v tests/tests.py -m full_run -k "repadapt or bcftools or hard_filters" --conda-prefix $(pwd)/.snakemake/conda
pixi run -e dev lint
git checkout tests/data/fixtures && git status   # no fixture changes committed
```

**PR 1 is done when:**
- All dry-run tests pass, and the GATK and bcftools filter tests are unchanged.
- The unit test and both full runs pass on this box.
- Lint is clean.
- The golden test passes (after PR 0 merges into `feat/repadapt`).

## Decisions (settled 2026-09-30)

1. `GOLDEN_WORK` is `.../snparcher-dev/repadapt-golden/`. It holds the RepAdapt clone, the SIFs, `nf-env` and the runs. pixi is installed in `~/.pixi/bin`, not on PATH.
2. All three additions are kept: the fixture features, the golden raw calls, and the dry-run per-region command check.
3. For cluster-sized runs, use the Cannon profile at `~/snpArcher/workflow-profiles/default` (by path; don't copy it into the repo). Tests run locally with `--cores 1`.

## Status

- **Step 0:** done.
- **PR 0:** merged (#347).
  - Golden run: two RepAdapt runs gave identical records. Raw outcomes are PASS 50, AllHomAlt 3, LowMQ 3, AllHomAlt;LowMQ 1.
  - `tests/data/repadapt/` totals 1.1 MB.
  - `test_repadapt_fixture_is_reproducible` and `test_repadapt_golden_records_are_consistent` pass.
- **PR 1:** implemented on `feat/repadapt-calling` (#348).
  - Verified on Linux (`agents-plan.md` 7.1 and 7.7): the pin file installs, Snakemake logs `Using pinnings from ...bcftools.linux-64.pin.txt`, and bcftools 1.16 gives the expected soft-filter labels.
  - `test_repadapt_golden_bams`: calling on RepAdapt's BAMs reproduces its raw and PASS records exactly.
  - **Real-data check** (`../repadapt-calbicans/`): 10 *C. albicans* isolates through RepAdapt's real pipeline, then `tool: repadapt` on its BAMs.
    - With RepAdapt's BAM order, all 296,311 raw and 293,429 PASS records are identical.
    - With the sample-sheet order, 4 records in the collapsed rDNA differ, because of BAM order at extreme depth (`agents-plan.md` gotcha 23).
- To regenerate: `G=.../repadapt-golden; NEXTFLOW=$G/nf-env/bin/nextflow JAVA_HOME=$G/nf-env/lib/jvm tests/repadapt/make_golden.sh $G`

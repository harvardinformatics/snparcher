# Implementation plan: PR 3 (`mapping.pipeline: repadapt`)

This is the last of the nested PRs into `feat/repadapt` (branch `feat/repadapt-mapping`). It adds the RepAdapt mapping pipeline: fastq → RepAdapt-equivalent final BAM, inside snpArcher's structure. The agreed design is `agents-plan.md` 2.4 and 2.5. This plan turns it into file-level steps, records what was checked on Linux while planning, and lists the few design choices made here.

Done means: `mapping.pipeline: repadapt` plus `tool: repadapt` on the fixture reads reproduces the golden RepAdapt records exactly, raw and PASS, and nothing else changes.

## Checked on Linux while planning (handoff section 7)

These were run with RepAdapt's own images, on snpArcher's reference files from the *C. albicans* run.

| Item | Result | Consequence |
|---|---|---|
| 7.4 GATK 3.8 and a bgzipped reference | Fails: "The GATK cannot process compressed (.gz) reference sequences" | Add a rule for an uncompressed `results/reference/{name}.fa` plus `.fai`. |
| 7.4 GATK 3.8 and snpArcher's `samtools dict` | Works with the uncompressed `.fa` (293 targets on a test region) | Reuse `results/reference/{name}.dict`; no Picard dictionary needed. |
| 7.5 GATK 3.8 and a CSI-only BAM | Fails: "not indexed" | Write a temp BAI for each realignment input. |
| 7.6 bwa 0.7.17 and snpArcher's index (built by bwa 0.7.19 from the `.fa.gz`) | Works | Reuse `REF_BWA_IDX`. |
| 7.8 fastp 0.20.1 JSON and threads | Has `summary.before_filtering` / `after_filtering` `total_reads`; caps at 16 threads with a warning | `collect_fastp_stats` works unchanged; the profile's 6 threads are fine. |
| 7.9 samtools 1.16.1 `flagstat -O tsv` | Same line layout as `parse_bam_stats` uses (0 total, 4 duplicates, 6 mapped, 10 paired in sequencing, 13 properly paired). "Properly paired %" is relative to paired in sequencing (72,081 / 72,726 = 99.11%). | The pre-filter QC sums lines 0, 6, 10 and 13 and recomputes the percentages. |
| 7.2 bwa 0.7.17 + samtools 1.16.1 in one env | Solvable. **Unpinned, the solver pairs samtools 1.16.1 with htslib 1.21**; RepAdapt's image has htslib 1.16. With `htslib=1.16=h6bc39ce_0` it solves to 27 packages. | Pin htslib 1.16 in the mapping env. |

Still to check in step 1:
- 7.1: the new pin files install, and Snakemake logs that it used them.
- 7.3: `gatk3 -T RealignerTargetCreator --help` runs non-interactively from the conda env.
- The Java versions: GATK on OpenJDK 8, Picard on 11.

## Design choices made here

These are open to review; none of them reopens `agents-plan.md` 2.7.

1. **Duplicate removal reads a sample's row BAMs directly.**
   - Picard runs with one `-INPUT` per row, and merges them internally, instead of a separate `samtools merge` step.
   - A single-row sample is then exactly RepAdapt's step. Multi-row samples, which RepAdapt can't take at all, are still deduplicated together, under one library tag.
   - This saves one BAM rewrite per sample, and it avoids a temp/final ambiguity that a separate merge output would create.
   - Samples with `mark_duplicates: false` are merged with `samtools merge`.
2. **No generic validation hook.**
   - The only startup check needing information the registry doesn't have is realignment versus long-contig mode. That becomes `resolve_repadapt_indel_realignment()`, placed in `common.smk` after `LONG_CONTIG_MODE`.
   - Caller compatibility is already enforced in `resolve_mapping_pipeline`: `tool: sentieon` with `repadapt` is an error.
3. **Pin htslib 1.16 in the mapping env**, per the 7.2 finding above.
4. **Pre-filter flagstat through a FIFO and `wait`, not `tee >(...)`.**
   - With process substitution, a failing flagstat doesn't fail the job.
   - The process can also still be writing when the rule ends.
   - With a FIFO, the job waits for flagstat and fails if it fails.

## 1. Envs (`workflow/envs/repadapt/`)

Each env is a `.yaml` with versions only, plus a `.linux-64.pin.txt` with exact builds. These go next to PR 1's `bcftools` env.

| Env | yaml | Pin source |
|---|---|---|
| `fastp` | `fastp=0.20.1` | Appendix A |
| `mapping` | `bwa=0.7.17`, `samtools=1.16.1`, `htslib=1.16` | Solve on Linux with `bwa=0.7.17=h5bf99c6_8 samtools=1.16.1=h6899075_0 htslib=1.16=h6bc39ce_0`, then `conda list --explicit --md5` |
| `picard` | `picard=2.26.3`, `openjdk=11.0.9.1` | Appendix A (106 packages) |
| `gatk3` | `gatk=3.8` | Appendix A (153 packages; build 9 bundles the 3.8-1 jar) |

To check each env, create it from its pin, then run:
- `fastp --version`
- `bwa 2>&1 | grep Version`
- `samtools --version`, which must report htslib 1.16
- `picard MarkDuplicates --version`
- `gatk3 -T RealignerTargetCreator --help`
- `java -version` in the picard and gatk3 envs

## 2. Config

- **Schema:**
  - `mapping.pipeline` enum `["default", "sentieon", "repadapt"]`.
  - New `mapping.repadapt` object (`additionalProperties: false`, `default: {}`) with `indel_realignment: oneOf [boolean, enum ["auto"]]`, default `"auto"`.
  - Description text for both.
- **`DEFAULTS` in `common.smk`:** `"mapping": {"pipeline": "default", "repadapt": {"indel_realignment": "auto"}}`. Both places must agree.
- **`common.smk`, after `LONG_CONTIG_MODE`:**
  ```python
  def resolve_repadapt_indel_realignment():
      setting = config["mapping"]["repadapt"]["indel_realignment"]
      if MAPPING_PIPELINE != "repadapt":
          return False
      if LONG_CONTIG_MODE:
          if setting is True:
              raise ValueError(...)  # GATK3 needs BAI, which can't index contigs > 2^29
          if setting == "auto":
              logger.warning(...)    # realignment skipped for long contigs
          return False
      return setting in (True, "auto")

  REPADAPT_INDEL_REALIGNMENT = resolve_repadapt_indel_realignment()
  ```
- **Registry entry:**
  ```python
  "repadapt": {
      "rules": "rules/mapping/repadapt.smk",
      "final_bam": _repadapt_final_bam,  # realigned/ if realigning, else dedup/ or merged/
      "qc_json": lambda sample: f"results/qc_metrics/repadapt/{sample}.json",
      "extra_qc": lambda: {},
  },
  ```
- **`config/config.yaml`:** the new pipeline and `mapping.repadapt.indel_realignment`, with comments.

## 3. Rules (`workflow/rules/mapping/repadapt.smk`)

All rules use the pinned envs above. Paths are under `results/bams/repadapt/` unless noted.

| Rule | Input → output | Command (result-affecting parts) |
|---|---|---|
| `fastp` (same name and outputs as default) | staged reads → `results/filtered_fastqs/...`, `results/fastp/...json` | `fastp -i R1 -I R2 -o O1 -O O2 -w {threads} -j {json} -h /dev/null`, without `--detect_adapter_for_pe` |
| `bwa_mem` (same name and row BAM path as default) | filtered reads → temp `results/bams/raw/{s}/{lib}/{unit}.bam`, plus temp `results/qc_metrics/repadapt/flagstat/{s}/{lib}/{unit}.tsv` | see below |
| `repadapt_remove_duplicates` (samples with `mark_duplicates: true`) | the sample's row BAMs → `dedup/{s}.bam` (temp when realigning), plus `results/qc_metrics/repadapt/{s}_duplicates.txt` | `picard -Xmx{mem_mb_reduced}m MarkDuplicates -INPUT row1 [-INPUT row2 ...] -OUTPUT OUT -METRICS_FILE M -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --TMP_DIR {resources.tmpdir}` |
| `repadapt_merge_sample` (`mark_duplicates: false`) | the sample's row BAMs → `merged/{s}.bam` (temp when realigning) | `samtools merge -o OUT rows...` |
| `repadapt_reference_fasta` (when realigning) | `{name}.fa.gz` → temp `results/reference/{name}.fa` and `.fa.fai` | `gunzip -c`, then `samtools faidx` |
| `repadapt_index_bai` (when realigning) | `dedup/` or `merged/` BAM → temp `.bam.bai` | `samtools index -b` |
| `repadapt_realigner_targets` (when realigning) | BAM, BAI, `.fa`, `.fa.fai`, `.dict` → temp `realign/{s}.intervals` | `gatk3 -Xmx{m}m -Djava.io.tmpdir={tmp} -T RealignerTargetCreator -R REF.fa -I IN -o T.intervals` |
| `repadapt_indel_realigner` (when realigning) | the same, plus the intervals → `realigned/{s}.bam`, plus temp `realigned/{s}.bai` (written by GATK) | `gatk3 ... -T IndelRealigner -R REF.fa -I IN -targetIntervals T --consensusDeterminationModel USE_READS -o OUT`; single thread, no `-nt` |
| `repadapt_mapping_qc` | pre-filter flagstats, Picard metrics, `bam_stats` coverage → `results/qc_metrics/repadapt/{s}.json` | see QC below |

**The `bwa_mem` command.** This keeps RepAdapt's order: map, MAPQ filter, name sort, fixmate, coordinate sort. It streams, with no SAM file.
```
fifo=$(mktemp -u {resources.tmpdir}/flagstat.XXXXXX); mkfifo "$fifo"
samtools flagstat -O tsv "$fifo" > {output.flagstat} &
flagstat_pid=$!
bwa mem -K 40000000 -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} \
  | tee "$fifo" \
  | samtools view -b -q 10 - \
  | samtools sort -n -@ {threads} -T {tmp}/n - \
  | samtools fixmate -m - - \
  | samtools sort -@ {threads} -T {tmp}/c -o {output.bam} - 2>> {log}
wait "$flagstat_pid"; rm -f "$fifo"
```

**Read group.** `@RG\tID:{library}.{input_unit}\tSM:{sample}\tLB:{sample}_LB\tPL:ILLUMINA`. The row-level ID is snpArcher's. The single per-sample LB is RepAdapt's, and it's what makes Picard remove duplicates per sample.

**After the final BAM, everything is shared:**
- `index_bam_csi`, which writes the CSI index.
- `bam_stats`, whose coverage the QC JSON uses.
- Callable sites and every caller except Sentieon's.

BAM inputs are used as they are and get the shared QC.

**QC JSON for samples this pipeline maps.** It has the same keys as `parse_bam_stats`, named `{sample}.json`.
- **From the per-row pre-filter flagstats, summed line by line:**
  - `total_reads` and `num_mapped` (lines 0 and 6)
  - `percent_mapped` = mapped / total × 100
  - `percent_properly_paired` = line 13 / line 10 × 100
- **From Picard's metrics:**
  - `num_duplicates` = `UNPAIRED_READ_DUPLICATES + 2 × READ_PAIR_DUPLICATES`
  - `percent_duplicates` = `PERCENT_DUPLICATION × 100`
  - both are 0 when duplicate removal is skipped
- **From `bam_stats` coverage of the final BAM:** `mean_depth` and `covered_bases`.

To compute the coverage values, move `parse_bam_stats`' coverage parsing into a module-level helper that both rules use. That's no behavior change, and `test_qc_metrics` covers it.

**Profiles.** In `workflow-profiles/default/config.yaml`:
- Threads: the existing `fastp: 6` and `bwa_mem: 16` apply. Add `repadapt_remove_duplicates`, `repadapt_realigner_targets` and `repadapt_indel_realigner` at 1.
- Resources:
  - Picard: `mem_mb` and `mem_mb_reduced`
  - GATK3: 16 GB and `runtime` 2880, matching RepAdapt's 48 h
- Heap size doesn't change results: Picard only spills to disk, and IndelRealigner's limits count reads, not bytes.

## 4. Tests

Test names avoid `qc`, `postprocess` and `metadata`, so they stay in the core CI groups.

**Dry-run**
- `test_repadapt_mapping_commands`: the exact commands.
  - fastp without `--detect_adapter_for_pe`
  - `bwa mem -K 40000000` without `-M`
  - `view -b -q 10 | sort -n | fixmate -m | sort`
  - `LB:{sample}_LB` with ID `{library}.{input_unit}`
  - Picard `-REMOVE_DUPLICATES true` with `-Xmx` from the profile
  - both GATK3 commands, on the uncompressed reference
- `test_repadapt_mapping_multirow_sample`: a multi-row, multi-library sample is mapped per row, then deduplicated in one Picard job with one `-INPUT` per row and one LB.
- `test_repadapt_mapping_no_dedup`: `mark_duplicates: false` gives no Picard job; the sample is merged, then realigned.
- `test_repadapt_mapping_realignment_settings`:
  - `false` gives no GATK jobs, and the final BAM is `dedup/`.
  - long contigs with `auto` skip realignment, with a warning.
  - long contigs with `true` are an error.
- `test_repadapt_mapping_bam_inputs`: a BAM-input sample is staged as is, with the shared QC JSON; fastq samples get the repadapt QC JSON in `combine_qc_metrics`' input.
- `test_repadapt_mapping_callers`:
  - `tool: sentieon` is rejected.
  - `gatk` and `bcftools` dry-run with `repadapt` mapping.
- Regression: the default and sentieon dry runs are unchanged. Rerun PR 2's 14-scenario snapshot script before and after.

**Full runs (Linux)**
- `test_repadapt_mapping_golden`:
  - Uses `mapping.pipeline: repadapt` and `tool: repadapt` on `repadapt_fastqs.csv`.
  - Raw records must equal `raw_records.tsv`, and PASS records `final_records.tsv`.
  - The fixture is at most 1,000 pairs, so fastp is deterministic at any thread count. `-K` fixes bwa's batching, and Picard and single-threaded GATK3 are deterministic.
- In the same test:
  - The final BAMs contain OC tags, so realignment ran. PR 0 showed the VCF barely changes. The check counts `OC:Z` tags in the decompressed BAM, as `b"OCZ"` followed by a digit, since an OC value is a CIGAR string; that's very unlikely to occur by chance in BAM's binary fields.
  - The QC report has sensible values. `percent_mapped` is below 100 because the fixture's junk mates are unmapped, which shows the values are pre-filter. `num_duplicates` is about 2 × the fixture's 50 duplicate pairs.
  - Callable sites are produced.
- Development-time check, not a test: compare our final BAM records with the golden BAMs, ignoring the `RG` and `PG` tags. The 11 mandatory fields should match exactly. If they don't, find the step that diverges before running the VCF golden test.

**CI**
- Add a third `setup-test-envs` pass with `'mapping={pipeline: repadapt}' 'variant_calling={tool: repadapt}'`.
- The picard and gatk3 envs are large (r-base and Java), so the first build after the cache key changes will be slow.

## 5. Docs

- **`explanation/` mapping section** (in `variant-calling.md`, or a new page):
  - what the pipeline does, step by step against RepAdapt
  - what's kept from snpArcher
  - one LB per sample, so duplicates are removed per sample
  - `mark_duplicates: false` skips duplicate removal, unlike RepAdapt
  - QC metric definitions, and how they differ from the default pipeline
  - limitations: Linux only; GATK3 is single-threaded and slow; it fails on heavily fragmented references (per RepAdapt's README); realignment is skipped for long contigs
- **`reference/config-schema.md`:** the enum and the `mapping.repadapt` block.
- **`how-to/configure.md`:** a row in the pipeline table, plus the RepAdapt-equivalent example from `repadapt-snparcher-plan.md`.
- **`reference/outputs.md`:** `results/bams/repadapt/...` and `results/qc_metrics/repadapt/...`.
- **`reference/changelog.md`.**

## 6. Real-data check (after the tests pass)

Reuse `../repadapt-calbicans/` with a new stage: snpArcher from reads, with `mapping.pipeline: repadapt` and `tool: repadapt`, compared with RepAdapt's VCF.
- **Not identical:** these samples have more than 1,000 pairs, so fastp's multi-threaded output order differs from run to run in both pipelines, and bwa's batches follow it. Report concordance instead:
  - shared PASS sites
  - genotype agreement
  - the distribution of differences
- **What a good result looks like:** near-total agreement, with differences scattered and not systematic. Systematic differences would point at a step.
- **Where to run it:** locally, as before, since Snakemake's per-job startup on the compute nodes was the bottleneck. Or on Slurm, if that turns out to have been transient.

## Order of work

1. Create the branch and the envs, and check them (section 1).
2. Config and registry. Take the snapshot baseline first.
3. Rules, profile entries and the QC rule.
4. Dry-run tests, then the development-time BAM comparison against the golden BAMs, then the full and golden tests.
5. Docs, CI, the snapshot regression, and the full suite: dry-run, unit, and the full runs for repadapt, gatk and bcftools.
6. The real-data check, then the PR.

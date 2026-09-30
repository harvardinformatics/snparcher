# Agent plan: RepAdapt support in snpArcher

A handoff document for the agent that implements this work. It records the agreed design, the research behind it, the codebase map, gotchas and the checks that still need doing. `repadapt-snparcher-plan.md` has the human-readable summary.

## 0. Status (2026-09-30)

- **Status:** planning is complete and no code has been written. The design was iterated with the user several times; section 2 is the final, agreed version. **Don't reopen the decisions in section 2.7 unless the user raises them.**
- **Branch:** `feat/repadapt` is the integration branch. It was cut from `main` at `09a58b7` (#345, the hard-filter fix) and so far contains only these two docs.
- **Delivery:** nested PRs, each targeting `feat/repadapt` (section 3). At the end, `feat/repadapt` merges to `main`.
- **Environment:** the user switched to an HPC dev box, which has Linux and presumably Apptainer and a scheduler. The work was planned on a Mac, where the pinned tools can't be installed.
- **Related:**
  - Issue harvardinformatics/snpArcher#346 revisits `bwa mem -M`/`-Y`/`-K` for the default pipeline. It's separate from this work.
  - A separate task suggestion was raised: DeepVariant's rule uses `container:` only, but nothing documents or checks that Apptainer is required. Not part of this work.
- **RepAdapt pipeline:** https://github.com/RepAdapt/nextflow_snp_calling_linux, pinned at commit `2077f6f311650c24ccc0df8c4088e616e624f5d8` (2026-07-08).

### User preferences that apply here
- Plan before implementing large features, and keep each PR focused.
- Validate changes with the pytest harness (`tests/tests.py`, `tests/unit_tests.py`), not ad-hoc probing. Ad-hoc cluster probes caused wrong conclusions before.
- **Never name an internal option or pipeline "snparcher".** The whole tool is snpArcher; use `default` or a descriptive name.
- The goal is equivalent **VCFs**. Byte-identical BAMs and RepAdapt's coverage tables are explicitly not goals.
- Keep snpArcher's structure and behavior (per-row mapping, QC report, callable sites, modules) wherever it doesn't change VCF results.
- Unit tests modify files under `tests/data/fixtures`; restore them with `git checkout tests/data/fixtures` before committing.

## 1. Goal

1. **A RepAdapt calling model** (`variant_calling.tool: repadapt`): an alternative BAM → VCF path that uses RepAdapt's bcftools settings. Raw calls go to `results/vcfs/raw.vcf.gz`. A separate RepAdapt filter step writes `results/vcfs/filtered.vcf.gz`.
2. **Pluggable mapping pipelines** (`mapping.pipeline: default | sentieon | repadapt`). The `repadapt` pipeline reproduces RepAdapt's fastq → BAM processing using RepAdapt's pinned tool versions and result-affecting arguments, inside snpArcher's structure. Its final BAMs feed the calling model (or any caller except Sentieon's), QC and callable sites.

Together, `mapping.pipeline: repadapt` plus `tool: repadapt` should produce VCF records equivalent to what RepAdapt's pipeline would produce from the same reads.

**Non-goals:**
- a literal Nextflow replica
- Apptainer containers at runtime
- byte-identical BAMs
- RepAdapt's per-gene, per-window and whole-genome depth tables (which need a GFF)
- RepAdapt's `errorStrategy 'ignore'`, which silently drops failed samples
- a separate BAM-input path for (b); BAM inputs simply go to the calling model

## 2. Agreed design (authoritative)

### 2.1 Config

```yaml
mapping:                        # NEW top-level block
  pipeline: default             # default | sentieon | repadapt
  repadapt:
    indel_realignment: auto     # auto | true | false
variant_calling:
  tool: repadapt                # NEW enum value
```

- **Schema:** edit `workflow/schemas/config.schema.yaml`. The top level has `additionalProperties: false`, so the new `mapping` object must be declared with `default: {}` and nested defaults. Model `indel_realignment` on `long_contig_mode`: `oneOf: [boolean, enum ["auto"]]`, default `"auto"`.
- **Resolving the pipeline:**
  - `tool: sentieon` with `pipeline: default` resolves to the `sentieon` pipeline. This keeps existing configs working, since Sentieon calling has always used Sentieon mapping.
  - `tool: sentieon` with `pipeline: repadapt` is an error.
  - `pipeline: sentieon` with any tool other than sentieon is an error in v1.
  - The value `default` can't be told apart from an explicit user setting, because schema defaults are filled in first. That's fine: `default` with `tool: sentieon` has never been a valid combination.
- **Ignored settings:** under `tool: repadapt`, `variant_calling.bcftools.*` (`min_mapq`, `min_baseq`, `max_depth`) has no effect. Document this rather than warning, because we can't tell user-set values from schema defaults.
- **Examples:** update `config/config.yaml` with examples and comments.

### 2.2 Calling model (`tool: repadapt`), PR 1

Per contig, on the final BAMs of every sample that has one (`SAMPLES_WITH_BAM` via `get_final_bam`):

```
bcftools mpileup -Ou -f REF -r CONTIG BAM... -q 10 -I -a FMT/AD,FMT/DP \
  | bcftools call -G - -f GQ -mv --ploidy PLOIDY -Oz -o OUT
```

- **Differences from RepAdapt's literal command** (`-q 5`, no `--ploidy`):
  - `-q 10` approximates RepAdapt's BAM-level `samtools view -q 10`. On `repadapt`-pipeline BAMs every read already has MAPQ ≥ 10, so `-q 10` and `-q 5` select exactly the same reads.
  - `--ploidy` behaves identically to omitting it for diploids. bcftools 1.16 accepts only its predefined aliases, which include `1` and `2`, so validate ploidy ∈ {1, 2} at startup.
- **Don't add `-Q`, `-d` or `--threads`.** RepAdapt uses the bcftools 1.16 defaults: **min BQ 1** (not 13), max depth 250, min MQ 0.
- **Pin bcftools to 1.16** (`workflow/envs/repadapt/bcftools.yaml` plus a linux-64 pin file; see appendix A). The reason: 1.16 writes INFO/MQ as an Integer, while newer versions write a Float, which can move sites across the `MQ < 30` filter.
- **Contigs and concat:** use one region per contig, reusing the bcftools region checkpoint logic. Concatenate the per-contig VCFs in `.fai` order with a plain `bcftools concat -f list` into `RAW_VCF`, then index with `BCFTOOLS_INDEX_ARGS` (TBI, or CSI in long-contig mode).
- **Filter** (`repadapt_filter`) writes `FILTERED_VCF` using two soft filters:
  ```
  bcftools filter -s AllHomAlt -e 'AC=AN' -Ou RAW \
    | bcftools filter -m + -s LowMQ -e 'MQ<30' -Oz -o FILTERED
  bcftools index BCFTOOLS_INDEX_ARGS FILTERED
  ```
  - The PASS set equals the records RepAdapt's `bcftools filter -e 'AC=AN || MQ < 30'` keeps.
  - Without `-s`, bcftools sets PASS on kept records whose FILTER was `.`, so RepAdapt's output has FILTER=PASS and our PASS records match it exactly.
  - Verify how `-m +` handles PASS; see section 7.
- **`generate_filtered_vcf`:** applies to `repadapt` too. When it's false, the raw VCF is the final call set.
- **Other wiring:**
  - reject gVCF inputs
  - allow long-contig mode (bcftools handles CSI)
  - warn that `postprocess.filtering.split_by_type` will produce an empty `clean_indels.vcf.gz`, because `-I` means no indels are called
  - QC and postprocess modules work unchanged; QC benefits from FORMAT/DP

### 2.3 Mapping pipelines (contract), PR 2

A mapping pipeline turns a sample's staged per-row reads into one coordinate-sorted final BAM, plus QC inputs. The shared code handles SRA downloads and fastq staging, external BAM inputs, CSI indexing, `bam_stats`, callable sites and calling.

**Contract.** Each pipeline provides:
1. **The final BAM** for each sample it maps. The shared `get_final_bam(sample)` keeps returning `results/bams/input/{sample}.bam` for BAM inputs and otherwise delegates to the pipeline.
2. **The per-sample mapping QC JSON path.** It must use the keys `sample, total_reads, num_mapped, percent_mapped, num_duplicates, percent_duplicates, percent_properly_paired, mean_depth, covered_bases`, and the basename must be `{sample}.json` (because `combine_qc_metrics` takes the sample name from the basename). Samples supplied as BAMs always use the shared `results/qc_metrics/bam/{sample}.json` from `parse_bam_stats`.
3. **Extra QC inputs**, optional. This replaces `if USE_SENTIEON` in `combine_qc_input`; Sentieon adds its insert-size JSONs this way.
4. **A validation hook**, run at startup, for things like long-contig support and caller compatibility.

**Registry.** Put a `MAPPING_PIPELINES` table in `common.smk` mapping each name to its rule file and hooks, plus `MAPPING_PIPELINE = resolve_mapping_pipeline(...)`. The Snakefile includes only the selected pipeline's file.

**File layout** (target):
- `workflow/rules/mapping/common.smk`: `stage_external_bam`, `index_bam_csi`, the shared merge rules, `bam_stats`, contract helpers.
- `workflow/rules/mapping/fastp.smk`: today's `fastp` rule, moved from `fastq.smk`. Included by `default` and `sentieon`.
- `workflow/rules/mapping/default.smk`: today's `bwa_mem` (`-M`) and `markdup_library` (sambamba), plus merges, moved verbatim.
- `workflow/rules/mapping/sentieon.smk`: today's `sentieon_map`, `sentieon_dedup_library` and `sentieon_bam_stats`, moved verbatim.
- `workflow/rules/mapping/repadapt.smk`: new, in PR 3.
- `workflow/rules/fastq.smk` keeps `download_sra`, `fastp_input` and `collect_fastp_stats`.

**Constraints on the refactor:**
- **Keep every existing rule name** so profiles and tests don't change: `fastp, bwa_mem, markdup_library, merge_library_bams, merge_dedup_libraries, merge_library_level_bams, bam_stats, sentieon_map, sentieon_dedup_library, sentieon_bam_stats, stage_external_bam, index_bam_csi`.
- `USE_SENTIEON` currently serves both mapping and calling. After the refactor, mapping code keys off `MAPPING_PIPELINE`, and calling code keeps using `VARIANT_TOOL`.
- **Extending later:** a variant that only changes flags (for example, preserving today's `-M` behavior as `legacy` if #346 changes the default) is a new table entry that reuses an existing rule file with a different settings dict. A pipeline with different steps is a new rule file.

### 2.4 The `repadapt` mapping pipeline, PR 3

Per row, then per sample:

| # | Step | Tool (pinned) | Command / behavior | Kept from snpArcher |
|---|---|---|---|---|
| 1 | Trim | fastp 0.20.1 | `fastp -i R1 -I R2 -o O1 -O O2 -w {threads} -j {json} -h /dev/null`. No `--detect_adapter_for_pe`. | Per row, same rule name `fastp`, same outputs, snpArcher threads, JSON report. |
| 2 | Map + filter + fixmate | bwa 0.7.17 + samtools 1.16.1, one env | `bwa mem -K 40000000 -t {threads} -R '{rg}' REF R1 R2` → tee into `samtools flagstat -O tsv` (pre-filter QC) → `samtools view -b -q 10` → `samtools sort -n` → `samtools fixmate -m - -` → `samtools sort -o OUT` | Per row, streamed with no SAM file, read group set at mapping, snpArcher threads. |
| 3 | Merge | samtools | Merge **all rows of the sample** (not per library). | |
| 4 | Remove duplicates | Picard 2.26.3 | `picard -Xmx{mem}m MarkDuplicates -INPUT IN -OUTPUT OUT -METRICS_FILE M -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT` (a `--TMP_DIR` may be added). Once per sample. Skipped if the sheet says `mark_duplicates: false`. | |
| 5 | BAI index | samtools | `samtools index` (BAI, temp). GATK3 can't use CSI. | |
| 6 | Realign | GATK 3.8-1 | `gatk3 -Xmx{mem}m -T RealignerTargetCreator -R REF.fa -I IN -o T.intervals`, then `gatk3 -Xmx{mem}m -T IndelRealigner -R REF.fa -I IN -targetIntervals T.intervals --consensusDeterminationModel USE_READS -o OUT`. Single thread (no `-nt`). | Skipped automatically in long-contig mode (`auto`). |
| 7 | Final BAM | shared | `index_bam_csi`, `bam_stats` (coverage), mosdepth/callable sites, calling | all |

**Read group** for every row: `@RG\tID:{library}.{input_unit}\tSM:{sample}\tLB:{sample}_LB\tPL:ILLUMINA`.
- The library tag is **one per sample**, matching RepAdapt's `-RGLB ${id}_LB`. Picard groups duplicates by LB, so this makes duplicate removal per sample, the way RepAdapt does it.
- The row-level ID stays, as in snpArcher.

**Why each step's arguments are what they are:**
- **`-K 40000000`:** RepAdapt runs bwa with `-t 4`. bwa's batch size is 10,000,000 × threads unless `-K` is given, and per-batch insert-size estimation changes results. `-K 40000000` reproduces that batching at any thread count, which matters because tests run with `--cores 1`.
- **No `-M`:** RepAdapt doesn't use it. Supplementary alignments stay supplementary, and mpileup counts them.
- **`view -q 10`, then `sort -n`, then `fixmate -m`:** keep RepAdapt's exact order.
  - fixmate clears PAIRED and PROPER_PAIR on a read whose mate was filtered out, so mpileup then counts it as single-end.
  - `sort -n` also sets the tie order of the later coordinate sort, which affects Picard's tie-breaking and mpileup's depth cap.
  - Don't drop `sort -n`, even though bwa output is already grouped by name.
- **Picard `REMOVE_DUPLICATES=true`:** this is RepAdapt's flag. Picard also adds `PG:Z` tags to reads by default; that's fine.
- **GATK3 inputs:**
  - It needs an **uncompressed** reference. Add a rule that decompresses `results/reference/{name}.fa.gz` to `results/reference/{name}.fa` and runs `samtools faidx`. GATK 3.8's htsjdk almost certainly can't read bgzipped FASTA; verify.
  - It looks for the dictionary at `results/reference/{name}.dict`, the same path as snpArcher's existing `samtools dict` output. Verify GATK 3.8 accepts that file. If it doesn't, generate one with Picard 2.26.3 under a separate basename.
  - It needs a BAI for its input BAM.
  - It writes `{out}.bai` itself; declare that file as temp or delete it.
- **Dropped RepAdapt steps**, none of which change the VCF:
  - AddOrReplaceReadGroups (the read group is set at mapping)
  - indexing the sorted BAM
  - the SAM and temp files and `rm` lines
  - the depth tables
  - Nextflow file naming
  - per-sample fastq handling (we map per row)

**QC for samples the pipeline maps** (rule e.g. `repadapt_mapping_qc`, output `results/qc_metrics/repadapt/{sample}.json`):
- `total_reads`, `num_mapped`, `percent_mapped` and `percent_properly_paired`: sum the per-row pre-filter flagstat TSVs line by line and recompute the percentages:
  - mapped % = mapped / total
  - properly paired % = properly paired / paired in sequencing (check against samtools)
- `num_duplicates` = `UNPAIRED_READ_DUPLICATES + 2 × READ_PAIR_DUPLICATES`, and `percent_duplicates` = `PERCENT_DUPLICATION × 100`, both from Picard's metrics. Both are 0 if duplicate removal was skipped.
- `mean_depth` and `covered_bases`: from `samtools coverage` of the final BAM (`bam_stats` already produces it), parsed the same way as `parse_bam_stats`. These reflect the reads actually used for calling; document that.
- **Metric definitions differ** from the default pipeline (duplicates among MAPQ ≥ 10 primary reads, and so on). Document this.
- fastp 0.20.1's JSON contains `summary.before_filtering.total_reads` and `summary.after_filtering.total_reads`, the fields `collect_fastp_stats` reads. Verify.

**Validation and behavior:**
- `indel_realignment`:
  - `auto`: on, but off with a warning when `LONG_CONTIG_MODE`
  - `true` with long contigs: an error
  - `false`: off
  - When realignment is off, the final BAM is the deduplicated BAM (or the merged BAM if duplicates were skipped).
- `tool: sentieon` is an error with this pipeline. All other callers are allowed.
- BAM inputs are used as they are (no realignment) and get the shared QC. gVCF inputs follow the caller's rules.
- Add resources and threads for the new rules to `workflow-profiles/default/config.yaml` and `workflow-profiles/slurm/config.yaml`. Java rules use `mem_mb_reduced` for `-Xmx`.
- **Known limitations:** GATK3 is slow (single-threaded per sample; RepAdapt allows 48 h), and RepAdapt's README says GATK3 fails on heavily fragmented references. Document both. Scattering the realignment is a possible later optimization, but it must be shown to give identical output first.

### 2.5 Envs (conda only, no containers)

Put all RepAdapt envs under `workflow/envs/repadapt/`, each as a `.yaml` plus a `.linux-64.pin.txt`. Snakemake uses the pin file automatically.

| Env | yaml spec | Pin source |
|---|---|---|
| `bcftools` | `bcftools=1.16` | appendix A: the exact package list from the RepAdapt image (31 packages, conda-forge and bioconda only) |
| `fastp` | `fastp=0.20.1` | appendix A |
| `mapping` | `bwa=0.7.17=h5bf99c6_8`, `samtools=1.16.1=h6899075_0` | **Solve on Linux**, then export `conda list --explicit --md5`. It can't be copied from an image, because the two tools come from separate images; appendix A lists both images' packages for reference. |
| `picard` | `picard=2.26.3`, `openjdk=11.0.9.1` | appendix A (the image used openjdk 11.0.9.1 and r-base 4.1.1) |
| `gatk3` | `gatk=3.8=9` | appendix A (openjdk 8.0.265, python 3.9.0, r-base 3.6.3). The package **bundles the GATK 3.8-1 jar**, so no `gatk3-register` step is needed. |

- In the yaml, pin versions and leave builds unpinned, so macOS can at least attempt to solve. Exact builds live in the pin files.
- None of the exact builds exist for osx-arm64, and most don't exist for osx-64 either, so full runs are Linux-only. Use `skip_if_arm64_packages_unavailable` in tests.

### 2.6 Validation summary

**Errors:**
- `tool: sentieon` with `pipeline: repadapt`
- `pipeline: sentieon` without `tool: sentieon`
- `indel_realignment: true` with long contigs
- `tool: repadapt` with ploidy other than 1 or 2
- `tool: repadapt` with gVCF inputs

**Warnings:**
- realignment automatically skipped because of long contigs
- `split_by_type` with `tool: repadapt`

### 2.7 Decision log (settled with the user; don't reopen)

| Decision | Choice | Reason |
|---|---|---|
| Keep the BAM-level MAPQ ≥ 10 filter in the repadapt pipeline | **Yes** | It changes duplicate removal, realignment input and orphan handling, not just calling. Mapping-rate QC comes from the pre-filter flagstat. In general (default pipeline, `tool: repadapt`), the user considers BAM filtering wrong; there the approximation is mpileup `-q 10`. |
| Pin bcftools 1.16 for the calling model | **Yes** | INFO/MQ Integer vs Float changes MQ<30 filtering. Accept Linux-only. |
| fastp threads | **snpArcher's thread count** | Accept nondeterministic read order, which RepAdapt also has with `-w 4`. |
| Duplicate removal unit | **Per sample** (one LB per sample) | Matches RepAdapt, as a RepAdapt user would get by concatenating a multi-library sample. It only matters for samples with more than one `library_id`: roughly 0.01% of read pairs plus a few percent of orphaned reads at 30×, clustered near repeats. |
| Config shape | `mapping.pipeline: default\|sentieon\|repadapt` + `variant_calling.tool: repadapt` | Mapping paths should be extensible. |
| Name of the standard pipeline | `default` | Not "snparcher". |
| Sentieon | Moves into the mapping-pipeline structure (PR 2) | |
| PR structure | Nested PRs into `feat/repadapt` | |
| Rejected: a literal Nextflow replica with Apptainer, RepAdapt output naming, depth tables, byte-identical BAMs, GFF input | | The user judged these over-indexed on RepAdapt's processing quirks. |

## 3. PR breakdown (each targets `feat/repadapt`)

| PR | Branch | Depends on | Summary |
|---|---|---|---|
| 0 | `feat/repadapt-fixture` | none | Simulated fixture, a RepAdapt golden-output generator script, the golden outputs |
| 1 | `feat/repadapt-calling` | PR 0 (golden test only) | `tool: repadapt`, bcftools 1.16 env, filter, filtering wiring, tests, docs |
| 2 | `feat/mapping-pipelines` | none | Contract and registry, move `default` and `sentieon` with no behavior change, `mapping` config block |
| 3 | `feat/repadapt-mapping` | PRs 1 and 2 | `repadapt` mapping pipeline, envs, QC wiring, golden tests, docs |
| final | `feat/repadapt` → `main` | all | |

PRs 1 and 2 can proceed in parallel. When starting a sub-branch, update `feat/repadapt` from `origin` first.

### PR 0: fixture and golden output
- **`tests/data/repadapt/make_fixture.py`:** seeded and pure Python, with no numpy dependency. It writes:
  - **`reference.fasta`**, about 25 kb, with contigs:
    - `ctgA`, about 15 kb, containing an embedded ~1 kb segment present twice (a repeat that produces low-MAPQ and orphaned reads)
    - `ctgB`, about 8 kb
    - `ctgC`, about 1.5 kb
  - **Four diploid samples** `S1`–`S4`, each a single fastq pair:
    - genotypes at about 40 SNPs (a mix of het and hom-alt) and about 6 short indels (1–10 bp) to exercise realignment
    - 2×150 bp reads with a 350 ± 50 bp insert and about 0.5% sequencing error, with qualities
    - about 5% PCR-duplicate pairs
    - **at most 1,000 pairs per sample** (see the fastp note below)
  - **Read names:** SRA-like, e.g. `@SIM.123`. Avoid names starting `@NS`, `@NB` or `@A0`, which switch on fastp's polyG trimming. Tile/x/y-style names switch on Picard's optical-duplicate logic; both are best avoided.
  - **`genes.gff`**, a minimal GFF with a few `gene` lines, needed **only** because RepAdapt's pipeline requires `--gff_file`.
- **Sample sheets:** `tests/sample_sheets/repadapt_fastqs.csv` (fastq rows), and later `repadapt_golden_bams.csv` (`input_type: bam` pointing at the golden BAMs). **Config:** `tests/configs/repadapt.yaml`.
- **`tests/repadapt/make_golden.sh`** runs the real RepAdapt pipeline on the fixture:
  - **Requirements:** Linux, Apptainer, Nextflow **≤ 25.10.2** (per RepAdapt's README) and Java 17+.
  - **Steps:**
    1. Clone RepAdapt at `2077f6f`.
    2. Pull the seven depot images (section 5.1) and record their sha256.
    3. Write `golden.config` mapping each process to its image with `process { withName:trimSequences { container = '/abs/fastp.sif' } ... }` and `apptainer { enabled = true; autoMounts = true }` (or a `singularity {}` scope where only SingularityCE is installed, as on the holybioinf dev box). The table below lists the process names.
    4. Run with **`-C golden.config`** (capital C, which uses only this config). `-C` is a top-level option, so it goes **before** `run`. RepAdapt's own `nextflow.config` hard-codes the author's image paths.
    5. `LC_ALL=C nextflow -C golden.config run main.nf --ref_genome $FIX/reference.fasta --gff_file $FIX/genes.gff --reads "$FIX/fastq/*_{1,2}.fastq.gz" --outdir OUT`
  - **The reads pattern matters.** RepAdapt's default pattern is `./*{1,2}.fastq.gz`, which names `S1_1.fastq.gz` as sample `S1_`, with a trailing underscore. Use `*_{1,2}`.
  - **File suffixes:** the reference must end in `.fasta` and the GFF in `.gff`.
  - **Run twice** and confirm the VCF records and BAM records match between runs. This checks that RepAdapt is deterministic on the fixture.
- **Commit** to `tests/data/repadapt/golden/`:
  - the normalized `final_variants` records (text)
  - `final_variants.vcf.gz`
  - the realigned BAMs `{S}_sorted_RG_dedup_realigned.bam`
  - `PROVENANCE.md` (RepAdapt commit, image URLs and sha256, Nextflow version, command, date)

  Keep the total under about 2 MB.
- **Where to run it:** directly on the HPC dev box is simplest. A CI job is optional; if you add one, `workflow_dispatch` only works once the workflow file is on the default branch, so trigger it on `push` to the fixture branch.

  | Nextflow process | Image |
  |---|---|
  | `trimSequences` | fastp |
  | `fastaIndex`, `samtoolsSort`, `samtoolsRealignedIndex`, `samtoolsDedupIndex`, `calculateDepth` | samtools |
  | `gatkIndex`, `addRG`, `dupRemoval` | picard |
  | `bwaIndex`, `bwaMap` | bwa |
  | `realignIndel` | gatk |
  | `calculateGenesDepth`, `calculateWindowsDepth`, `calculateWgDepth` | bedtools |
  | `snpCalling`, `concatVCFs` | bcftools |

  `prepareDepth` and `joinDepth` have no container and run on the host.
- **Normalizing for comparisons:**
  - Put samples in a fixed sorted order (`bcftools view -s`) and contigs in `.fai` order, sort, then compare `bcftools view -H` records for CHROM through the samples.
  - Take our PASS records with `-f PASS`. The golden records are already all PASS.
  - Ignore the header: bcftools adds `; Date=` stamps, and RepAdapt's sample and contig order depends on Nextflow's collect order, which is arbitrary.

### PR 1: calling model
- **Files:** `workflow/rules/variant_calling/repadapt.smk` (new) with `repadapt_call` (per region), `repadapt_concat_regions` and `repadapt_filter`; env `workflow/envs/repadapt/bcftools.yaml` plus pin file.
- **Region checkpoint:**
  - Reuse the `bcftools_regions` checkpoint by moving it and its helpers (`_read_bcftools_regions`, `_get_bcftools_regions_file`, `get_bcftools_region_ids`, `get_bcftools_region_name`) into a shared file included by both callers.
  - Keep the checkpoint name `bcftools_regions`, because tests assert it.
  - Region IDs are `L{idx:06d}` in `.fai` order, so sorted IDs give `.fai` order.
  - Write repadapt per-region outputs to `results/vcfs/regions/repadapt/{region_id}.vcf.gz`, separate from bcftools'.
- **`common.smk` filtering** (about lines 619–649; #345 just made `generate_filtered_vcf` authoritative):
  ```
  REPADAPT_CALLER = VARIANT_TOOL == "repadapt"
  APPLY_GATK_HARD_FILTERS = GATK_LINEAGE_CALLER and GENERATE_FILTERED_VCF
  APPLY_REPADAPT_FILTER = REPADAPT_CALLER and GENERATE_FILTERED_VCF
  FINAL_VCF = FILTERED_VCF if (APPLY_GATK_HARD_FILTERS or APPLY_REPADAPT_FILTER) else RAW_VCF
  ```
  - Keep the "generate_filtered_vcf ignored" warning for bcftools and deepvariant only.
  - The Snakefile includes `hard_filters.smk` on `APPLY_GATK_HARD_FILTERS`.
  - **GATK behavior must not change.** Guarded by `test_generate_filtered_vcf_false_skips_hard_filtering`, `test_long_contig_hard_filters_*` and `test_short_contig_hard_filters_*`.
- **Other `common.smk` changes:** add `"repadapt"` to the gVCF-rejection set (about line 472), the ploidy check, and the `split_by_type` warning. The Snakefile tool dispatch is about lines 14–31.
- **Tests:**
  - dry-run: rules scheduled; `call_variants` → `filtered.vcf.gz`; no `variant_filtration`; `generate_filtered_vcf: false` gives raw; gVCF rejection (add `repadapt` to the parametrize of `test_gvcf_input_rejected_for_new_callers`); long-contig CSI; the `split_by_type` warning
  - bcftools regression: unchanged, and its source doesn't contain `-G -`
  - Per-region commands **don't appear in dry-run output** until the checkpoint has run. Assert flags via `workflow_source("rules","variant_calling","repadapt.smk")`, as `test_bcftools_long_contig_dry_run` does.
  - full run (Linux): header has FORMAT AD, DP and GQ; no indel records; FILTER values are only PASS, AllHomAlt, LowMQ or AllHomAlt;LowMQ; PASS records satisfy AC<AN and MQ≥30
  - golden: calling on the golden realigned BAMs gives PASS records equal to the golden records
- **CI:**
  - `pyproject.toml`'s `setup-test-envs` only builds the envs used by `config/config.yaml`'s DAG (`tool: gatk`). Add a second invocation with a repadapt config, and later one with `mapping.pipeline: repadapt`, so those envs are built and cached.
  - In `.github/workflows/test.yaml`, the conda cache key only hashes `workflow/**/*.yml` and `*.yaml`. **Add `workflow/**/*.pin.txt`**, or changes to pin files won't invalidate the cache.
- **Docs:**
  - `website/docs/reference/config-schema.md`: tool enum, `generate_filtered_vcf` text, a note on the bcftools block
  - `website/docs/explanation/variant-calling.md`: new RepAdapt section
  - `website/docs/explanation/filtering.md`, `reference/outputs.md`, `how-to/configure.md` (tool table), `reference/changelog.md`
  - `config/config.yaml` comments

### PR 2: mapping-pipeline refactor (no behavior change)
- Implement section 2.3.
- **Guard with snapshots.** Before touching code, capture dry-run output (`snakemake -n -p`: rule names and shell commands) for:
  - `tests/configs/local_genome.yaml` with `local_fastqs.csv`
  - the multi-row multi-library sheet
  - the no-dedup sheet
  - a mixed SRA + fastq sheet
  - `tool: sentieon`
  - a long-contig config

  After the refactor, assert the output is identical, after normalizing job IDs and ordering. A one-off script is fine; the existing tests stay the permanent guard.
- Profiles and existing tests must pass unchanged.
- Schema: the `mapping` block, defaulting to `default`. Update docs (`config-schema.md`, `configure.md`).

### PR 3: `repadapt` mapping pipeline
- Implement sections 2.4 and 2.5.
- **Tests:**
  - dry-run exact commands: no `-M`, `-K 40000000`, `view -b -q 10`, `sort -n`, `fixmate -m`, `LB:{sample}_LB`, Picard `REMOVE_DUPLICATES true`, and the realigner commands
  - realignment skipped with a warning in long-contig mode; `true` with long contigs errors; `tool: sentieon` errors
  - full run on the fixture: QC report present with sensible numbers (mapped % taken before filtering), callable sites produced
  - golden: `mapping.pipeline: repadapt` + `tool: repadapt` gives PASS records equal to the golden records
- Add the envs to the CI env build. Update profiles and docs; the docs should explain the QC metric definitions, per-sample duplicate removal, the one library tag per sample, and the limitations.

## 4. Codebase map (at `09a58b7`; line numbers approximate)

- **`workflow/Snakefile`:**
  - includes: common, reference, intervals, fastq, mapping, qc_metrics, callable_sites
  - tool dispatch (lines 14–31); hard filters included if `APPLY_HARD_FILTERS` (line 33)
  - `rule all` (lines 38–59)
  - module imports: postprocess and qc get `FINAL_VCF`, `REF_FILES`, and so on
  - targets `setup`, `map_samples`, `call_variants`, `qc_report`, `callable_sites`, `gvcfs`
- **`workflow/rules/common.smk`:**
  - `VARIANT_TOOL`/`USE_SENTIEON` (334–335)
  - `GATK_HARD_FILTERS` (339)
  - `mark_duplicates` parsing (373–445)
  - `MIXED_READ_INPUT_TYPES`/`SINGLE_ROW_INPUT_TYPES` (446–447)
  - gVCF rejection for bcftools/deepvariant/parabricks (472)
  - `REF_FILES`: `results/reference/{REF_NAME}.fa.gz` (bgzipped), `.fa.gz.fai`, `{REF_NAME}.dict`; `REF_BWA_IDX` = multiext of `.fa.gz` (502–512)
  - long-contig logic `_resolve_long_contig_mode`, `LONG_CONTIG_MODE`, `BCFTOOLS_INDEX_ARGS` (`-f -c` or `-f -t`), `get_compressed_vcf_index` (516–600)
  - `RAW_VCF`, `FILTERED_VCF` and work paths (604–614)
  - filter logic `GATK_LINEAGE_TOOLS`, `GENERATE_FILTERED_VCF`, `APPLY_HARD_FILTERS`, `FINAL_VCF` (619–649)
  - sample lists `SAMPLES_WITH_BAM`, `SAMPLES_WITH_FASTQ`, `SAMPLES_NEED_ALIGNMENT`, `SAMPLES_WITH_GVCF` (654–667)
  - callable-sites flags (669–735)
  - helpers `get_sample_libraries`, `get_sample_rows`, `get_sample_mark_duplicates`, `get_sample_inputs`, `BAM_INDEX_SUFFIX=".csi"`, `get_bam_index`, `get_external_bam` (739–812)
  - `get_final_bam` (814): input BAM → `results/bams/input/{s}.bam`; mark_dups → `results/bams/markdup/{s}.bam`; otherwise `results/bams/merged/{s}.bam`
- **`workflow/rules/fastq.smk`:** `fastp_input`, `download_sra` (sam-dump extraction, #343), `fastp` (about line 284; `--detect_adapter_for_pe`, `-j json -h /dev/null`, output `results/filtered_fastqs/{sample}/{library}/{input_unit}_{1,2}.fastq.gz` and `results/fastp/{sample}/{library}/{input_unit}.json`), `collect_fastp_stats` (323; sums before/after `total_reads`).
- **`workflow/rules/mapping.smk`:**
  - `bwa_mem_input`, `get_read_group` (`ID:{library}.{input_unit}`, `SM:{sample}`, `LB:{library}`, `PL:ILLUMINA`), merge-input helpers
  - `stage_external_bam`, `index_bam_csi`
  - `if USE_SENTIEON:` `sentieon_map`, `sentieon_dedup_library`, `sentieon_bam_stats`; `else:` `bwa_mem` (`bwa mem -M -t -R | samtools sort`, output `results/bams/raw/{sample}/{library}/{input_unit}.bam`), `markdup_library` (sambamba markdup to `results/bams/library_markdup/{s}/{lib}.bam`)
  - `merge_library_bams` (`results/bams/library/...`), `merge_dedup_libraries` (`results/bams/markdup/{s}.bam`), `merge_library_level_bams` (`results/bams/merged/{s}.bam`)
  - `bam_stats`: `samtools coverage` and `samtools flagstat -O tsv` on `get_final_bam`, to `results/qc_metrics/bam/{s}_coverage.txt` and `_flagstat.txt`
- **`workflow/rules/qc_metrics.smk`:**
  - `parse_bam_stats` reads flagstat TSV by **line index**: 0 total, 4 duplicates, 6 mapped, 7 mapped %, 14 properly paired %. From coverage it reads field 4 (covbases) and field 6 (meandepth), weighted by contig length.
  - `parse_sentieon_stats`; `combine_qc_input` (uses `USE_SENTIEON`); `combine_qc_metrics` (sample name from the JSON basename)
- **`workflow/rules/reference.smk`:** `prepare_reference` (local, URL or accession, bgzipped `.fa.gz`); `index_reference` (`samtools faidx`, `samtools dict`, `bwa index` on the `.fa.gz`; env `reference.yaml`).
- **`workflow/rules/variant_calling/bcftools.smk`:** `checkpoint bcftools_regions` (`results/vcfs/regions/regions.tsv`, `L%06d` IDs), `bcftools_call` (`-q/-Q/-d` from config, `-r contig`, `call -m --ploidy -v`), `bcftools_concat_regions` (`concat -D -a | sort`, output `temp(RAW_VCF)`; `RAW_VCF` is also a `rule all` target, so it's kept).
- **`workflow/rules/variant_calling/hard_filters.smk`:** GATK `variant_filtration`. Under long-contig mode it filters the uncompressed work VCF, then `compress_filtered_vcf`.
- **`workflow/rules/callable_sites.smk`:** `mosdepth --d4` on the final BAM plus CSI; clam; genmap.
- **`workflow/schemas/config.schema.yaml`:** `additionalProperties: false` at the top level and in nested blocks; `variant_calling.tool` enum; `long_contig_mode` oneOf bool/"auto".
- **Profiles:** `workflow-profiles/default/config.yaml` (`use-conda: True`; `set-threads` e.g. `fastp: 6`, `bwa_mem: 16`, `markdup_library: 16`, `bcftools_call: 8`; `set-resources`; `mem_mb_reduced` for Java) and `workflow-profiles/slurm/`.
- **Tests:**
  - `tests/conftest.py`:
    - `SnakemakeRunner(workdir, use_conda, conda_prefix)` with `.run(target, configfile, samples, extra_args, config_overrides)` and `.dry_run(...)`
    - runs with **`--cores 1`** and `--workflow-profile workflow-profiles/default`
    - markers `dry_run`, `full_run`, `unit`; `--dry-run-only`
  - `tests/tests.py` helpers: `get_config_file()` and `get_samples_file()` (arm64 uses the no-dedup sheet and the no-clam config), `write_config_for_tool` (text-replaces `tool: "gatk"`), `write_long_contig_config`, `scheduled_rule_present`, `workflow_source`, `iter_vcf_records`, `get_vcf_samples`, `skip_if_arm64_packages_unavailable`.
  - Existing bcftools tests are about lines 1567, 1695 and 2054; gVCF rejection is about line 1928.
- **CI** (`.github/workflows/test.yaml`):
  - dry-run matrix groups by `-k`: core = `not metadata and not postprocess and not qc`
  - `build-conda-envs` (pixi `setup-test-envs`), unit tests (4 splits), full runs (core in 3 splits, plus qc)
  - **Test names containing `qc`, `postprocess` or `metadata` land in those groups.**
- **pixi** (`pyproject.toml`): `pixi run -e dev pytest ...`; tasks `setup-test-envs`, `test-dry-run`, `lint` (ruff) and `format` (ruff + snakefmt). Snakemake 9.14.5 (`pixi.lock`), which supports pin files.

## 5. Research findings (verified from source or registries unless marked)

### 5.1 RepAdapt's pinned images
- RepAdapt's docs (RepAdapt/singularity `RepAdaptSingularity.imagelocations.md`, repo HEAD `ec3792b`) pull images from `https://depot.galaxyproject.org/singularity/<tag>`. Each depot SIF was built from `quay.io/biocontainers/<tag>`, and all the URLs return 200.
- `nextflow.config` points at the author's local SIF paths.
- All images use a **BusyBox** userland. BusyBox `sort` is always byte order, whatever the locale. This only matters for RepAdapt's depth tables, which we dropped.

| Tag | quay.io manifest digest | Bioconda build | Platforms |
|---|---|---|---|
| `fastp:0.20.1--h8b12597_0` | `sha256:56ca79fc827c1e9f48120cfa5adb654c029904d8e0b75d01d5f86fdd9b567bc5` | `fastp=0.20.1=h8b12597_0` | linux-64 only |
| `samtools:1.16.1--h6899075_0` | `sha256:3ee6abdc8f4d842ad82c4e81498609e6b3fbc5d0be51e7b459eabcbbe5d36a78` | `samtools=1.16.1=h6899075_0` | linux-64 only |
| `picard:2.26.3--hdfd78af_0` | `sha256:1fa3a74c7c9ce445c5c7c75135f62dea3e52029fbb499ca963626db19531ee46` | `picard=2.26.3=hdfd78af_0` | noarch |
| `bwa:0.7.17--h5bf99c6_8` | `sha256:f8494324de6da332792dc8e4acc2549152375e1966c96163087d6ff6d42ff48c` | `bwa=0.7.17=h5bf99c6_8` | linux-64 only |
| `gatk:3.8--9` | `sha256:e07c301b41224bd79f114438945677e3c62339d84e96659be29315c2b6d6c5db` | `gatk=3.8=9` | noarch |
| `bedtools:2.27.1--0` | `sha256:47a727b5eb7bf5ced0cf660973fe839671224c368d141a8dca44d7d229e7380a` | (not used) | |
| `bcftools:1.16--hfe4b78e_1` | `sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa` | `bcftools=1.16=hfe4b78e_1` | linux-64 only |

- **GATK image:** it contains GATK **3.8-1** (`opt/gatk-3.8/GenomeAnalysisTK.jar`, reporting `3.8-1-0-gf15c1c3ef`) on openjdk 8.0.265.
  - The wrapper is `gatk3`, a Python script. Its default heap is `-Xms512m -Xmx1g`; any `-Xm*` argument replaces the default, and `-D`/`-XX` arguments pass to the JVM.
  - **If `_JAVA_OPTIONS` is set at all, the defaults are skipped.**
  - It uses `$JAVA_HOME/bin/java` if that exists.
  - So RepAdapt effectively runs GATK3 with a 1 GB heap.
- **Picard wrapper:** default heap `-Xms512m -Xmx2g`. `-Xm*` arguments anywhere on the command line go to the JVM. It exports `LC_ALL=en_US.UTF-8` and prefers `$JAVA_HOME/bin/java`. The image used openjdk 11.0.9.1.
- **bedtools image:** its runtime libraries come from Anaconda's `defaults` channel, which has licensing terms. It isn't needed now.
- **Snakemake pin files:** `envs/X.linux-64.pin.txt` next to `envs/X.yaml`, in conda explicit format. Supported since 7.8. If creating from the pin fails, Snakemake warns and falls back to the yaml. On platforms without a pin file, it solves the yaml.
- **Combining `conda:` and `container:` on one rule** under `--sdm conda apptainer` builds the conda env inside the container, which fails for biocontainers because they have no conda. That doesn't matter now, since we use conda only.

### 5.2 Tool behavior that affects results
- **bwa 0.7.17** (`fastmap.c:359`, `bwamem.c:73`, `bwa.c:77`):
  - A batch is `-K` bases if given, otherwise 10,000,000 × threads, and it closes at the first even read count once that size is reached.
  - Insert-size statistics are estimated per batch (`mem_pestat`, `bwamem.c:1226`).
  - Random tie-breaks hash the global read index (`bwamem.c:527`), so **output depends on input read order**.
  - Output order is deterministic at any thread count, and `@PG CL` records the full argv.
  - `-t 8 -K 40000000` gives the same alignments as `-t 4`.
- **fastp 0.20.1:**
  - With two or more workers (the default is 2), it writes packs of 1,000 pairs in the **order they finish**. That is nondeterministic, which is why RepAdapt itself isn't bit-reproducible on inputs over 1,000 pairs. `-w 1` keeps input order.
  - It reads gzip with zlib `gzread`, which handles concatenated multi-member gzip files; gzip is detected only by the `.gz` suffix.
  - In PE mode with no options: adapter trimming by overlap only (sequence detection needs `--detect_adapter_for_pe`); quality filter `-q 15 -u 40 -n 5`; length filter `-l 15`; a pair is kept only if both reads pass.
  - PolyG trimming is automatic only when the first R1 header starts with `@NS`, `@NB` or `@A0`.
  - It always writes `fastp.json` and `fastp.html` to the current directory unless `-j`/`-h` are given, so always pass them.
- **samtools 1.16.1:**
  - `sort` gives the same order at any thread count or memory setting: key (tid, pos, strand), with ties stable in input order.
  - `sort -n` uses natural name order, then READ1/READ2.
  - **`fixmate` on a read whose mate is missing** clears 0x1, 0x2 and 0x20, sets mate tid/pos to -1 and TLEN to 0, and does **not** set 0x8 (`bam_mate.c:365-366`).
  - `depth` excludes UNMAP, SECONDARY, QCFAIL and DUP, applies no MAPQ filter, and counts overlapping mates twice.
- **bcftools 1.16:**
  - **mpileup defaults:** skips UNMAP, SECONDARY, QCFAIL and DUP (0x704). It skips "anomalous" reads, meaning PAIRED set and PROPER_PAIR not set, unless `-A`. That's why fixmate orphans, with PAIRED cleared, **are counted**.
  - More mpileup defaults: **min BQ 1**, min MQ 0, max depth 250 per file (order-dependent), partial BAQ, smart overlaps. mpileup writes `##bcftoolsCommand` with no date and `##reference=file://<path as given>`.
  - **`filter -e` without `-s`** drops failing records and sets PASS on passing records whose FILTER was `.` (`vcffilter.c:694-699`). **With `-s NAME`**, passing records get PASS and failing ones get NAME; `-m +` appends.
  - `##bcftools_*Command` header lines always get `; Date=` appended, unless `--no-version`.
  - **`call -m -G -`:** each sample is its own group. Each sample's allele-frequency estimate is its own normalized FORMAT/QS, or FORMAT/AD if QS is absent, and it's an **error without either**. The site's alleles are the union across samples, QUAL is the maximum per-sample QUAL, and the genotype prior comes from the sample's own read fraction. Behavior matches current master except INFO/MQ, which is **Integer in 1.16 and Float in master**.
  - The mutation prior θ (default 1.1e-3) is scaled by a Watterson factor over the whole cohort's chromosome count, even under `-G -`.
  - `--ploidy` accepts only predefined aliases (`GRCh37`, `GRCh38`, `X`, `Y`, `1`, `2`).
- **GATK 3.8:**
  - RealignerTargetCreator downsamples to 1000× per sample by default. IndelRealigner randomly subsamples reads for consensus building above 120 reads. The seed is fixed (`47382911`) unless `-ndrs`, so **single-threaded runs are deterministic** given the same input.
  - It adds OC/OP tags to realigned reads. The `@PG CL` holds the arguments as typed, with no date. Java 8 is required; Java 9+ likely fails.
- **Picard 2.26.3 MarkDuplicates:**
  - Deterministic: the comparator is a total order, and ties go to the lowest tile/x/y and then file order. Heap size only affects spilling to disk.
  - It groups duplicates **by library (LB)**. Fragments, including fixmate orphans, are marked as duplicates when they coincide with a paired read end in the same library.
  - It adds `@PG` and, by default, `PG:Z:MarkDuplicates` on every read. The metrics file has a date.

### 5.3 Background discussed with the user (context only)
- **`-G -` versus GATK:**
  - bcftools' default mode pools samples: it detects sites across the cohort and applies an HWE prior from the cohort frequency.
  - `-G -` does per-sample site detection (the site is kept if any sample calls it, which gives more false positives at low coverage) with a per-sample prior. Heterozygous calls need reads from both alleles, so genotypes lean homozygous at low depth.
  - GATK GenotypeGVCFs pools site detection (an allele-frequency model with `--heterozygosity` as the prior) but assigns genotypes by maximum likelihood with a flat prior (`USE_PLS_TO_ASSIGN`). snpArcher's `het_prior` therefore only affects QUAL and which sites are called.
- **bwa `-M`/`-Y`/`-K`:** Broad's WARP pipelines use `bwa mem -K 100000000 -p -v 3 -t 16 -Y` with no `-M`, and the CCDG standard says not to use `-M`. snpArcher's default `bwa_mem` has no `-K`, so its results depend on thread count; the Sentieon path uses `-K 10000000`. All of this is issue #346, and it doesn't block this work.

## 6. Gotchas

1. **The snpArcher reference is bgzipped** (`.fa.gz` + `.gzi` + `.fai`). GATK 3.8 needs an uncompressed FASTA. bwa and bcftools read the bgzipped one fine.
2. **BAM indexes are CSI throughout snpArcher** (#326). GATK3 needs BAI. Long contigs (>2^29) can't have a BAI, so realignment is skipped in long-contig mode.
3. `long_contig_mode: auto` inspects an existing `.fai` at DAG-build time. VCF index type follows `BCFTOOLS_INDEX_ARGS` and `get_compressed_vcf_index`.
4. **#345 made `generate_filtered_vcf` authoritative.** Don't change GATK filtering behavior when generalizing `APPLY_HARD_FILTERS`.
5. **Checkpoint jobs don't show their commands in dry runs.** Test per-region flags through `workflow_source`.
6. **Tests run with `--cores 1`**, so `bwa -t 1`. `-K 40000000` keeps batching the same.
7. **`parse_bam_stats` depends on the flagstat TSV line layout.** Any synthetic or summed flagstat must keep it, or the repadapt QC must write the JSON itself (recommended).
8. **`combine_qc_metrics` takes the sample name from the JSON basename**, so name the repadapt QC JSONs `{sample}.json`.
9. **fastp writes `fastp.json`/`fastp.html` to the working directory by default**, which would clobber the project root. Always pass `-j`/`-h`.
10. **Keep `samtools sort -n` before `fixmate`**, even though bwa output is already grouped by name. It sets the tie order used by the later steps.
11. **Picard's wrapper passes any `-Xm*` argument to the JVM.** Use `resources.mem_mb_reduced`. If `_JAVA_OPTIONS` is set, the wrapper's defaults are skipped. The GATK3 wrapper works the same way.
12. **GATK3 writes its own `{out}.bai`.** Declare it temp or remove it; snpArcher's final index is CSI via `index_bam_csi`.
13. **The pinned tools have no macOS builds**, so local Mac runs can only be dry runs.
14. **CI's conda cache key doesn't include `*.pin.txt`**; add it. `setup-test-envs` must also build the repadapt envs.
15. **CI groups tests by `-k` name substrings** (`qc`, `postprocess`, `metadata`). Name tests deliberately.
16. **Unit tests modify `tests/data/fixtures`.** Run `git checkout tests/data/fixtures` before committing.
17. **RepAdapt's default `--reads` glob adds a trailing underscore** to sample names (`S1_`). Use `*_{1,2}.fastq.gz` for golden runs.
18. **RepAdapt is not run-to-run deterministic above 1,000 read pairs per sample** because of fastp's output order. Keep fixture samples at or under 1,000 pairs, and compare against real data only statistically.
19. **`bcftools call -G -` needs FORMAT/AD or QS** from mpileup. Our `-a FMT/AD,FMT/DP` provides AD.
20. **Sample-sheet semantics in the repadapt pipeline:**
    - `library_id` doesn't affect duplicate removal (one library tag per sample).
    - `mark_duplicates: false` still skips duplicate removal. That's snpArcher behavior and a documented deviation from RepAdapt, which always removes duplicates.
    - Multi-row samples are mapped per row and then merged. RepAdapt can't take them at all, so equivalence is defined for single-row samples.
21. **Rule-name reuse.** The repadapt pipeline defines `fastp` and `bwa_mem` with the same names and outputs as `default`, which works because only one pipeline file is included. That keeps `collect_fastp_stats`, the profiles and `merge_*` wiring working.
22. **`include:` with a computed path** (e.g. `include: MAPPING_PIPELINES[...]["rules"]`) should work in Snakemake, because the path is a Python expression. Verify it, including with `snakefmt`/lint.

## 7. Still to verify on Linux (early in each PR)

1. The pin files install (`conda create --file X.linux-64.pin.txt`) and Snakemake picks them up.
2. bwa 0.7.17 (h5bf99c6_8) and samtools 1.16.1 (h6899075_0) install into one env. If they don't, split the mapping step at a pipe boundary.
3. `gatk=3.8=9` runs non-interactively: `gatk3 -T RealignerTargetCreator --help`.
4. GATK 3.8 can't read the bgzipped `.fa.gz` (so the uncompressed copy is needed), and it accepts snpArcher's `samtools dict` output.
5. GATK 3.8 needs a BAI, not CSI.
6. bwa 0.7.17 reads snpArcher's bwa index, which was built by a newer bwa from the `.fa.gz`. The index format is believed stable and the BWT is identical whichever version builds it.
7. The two-step soft filter gives: PASS when neither filter fails; `AllHomAlt` only; `LowMQ` only (**check that `-m +` removes PASS rather than leaving `PASS;LowMQ`**); `AllHomAlt;LowMQ` when both fail.
8. fastp 0.20.1's JSON contains `summary.before_filtering.total_reads` and `summary.after_filtering.total_reads`. Also check fastp 0.20.1's maximum `-w`; the profile uses 6.
9. The flagstat TSV layout in samtools 1.16.1 matches the `parse_bam_stats` indices, and the denominator for properly-paired % is "paired in sequencing".
10. The PR 2 dry-run snapshots are identical before and after the refactor.
11. The golden run is deterministic across two RepAdapt runs, and both golden tests pass.

## 8. Handy commands

```bash
pixi run -e dev pytest -v tests/tests.py --dry-run-only -k repadapt
pixi run -e dev setup-test-envs
pixi run -e dev pytest -v tests/tests.py -m full_run -k repadapt --conda-prefix $(pwd)/.snakemake/conda
pixi run -e dev lint
```

Commit and PR attribution: follow the current session's instructions.

## 9. RepAdapt reference (commit `2077f6f`, verbatim script blocks)

Process CPUs: fastp 4, bwa 4, samtoolsSort 4, everything else 1. The memory settings don't reach the Java heap; see section 5.1.

```bash
# trimSequences (fastp_trimming.nf)
fastp -w $task.cpus -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_1_trimmed.fastq.gz -O ${sample_id}_2_trimmed.fastq.gz
# bwaMap (bwa_mapping.nf)
bwa mem -t $task.cpus $reference ${trimmed_reads[0]} ${trimmed_reads[1]} > ${sample_id}.sam
# samtoolsSort (samtools_sort.nf)
samtools view -Sb -q 10 $sample_sam > temp666.bam
samtools sort -n -o temp777.bam temp666.bam
samtools fixmate -m temp777.bam temp888.bam
samtools sort --threads $task.cpus temp888.bam > ${sample_sam.baseName}_sorted.bam
samtools index ${sample_sam.baseName}_sorted.bam
# addRG (picard_add_read_groups.nf); id = sorted_bam baseName minus _sorted
picard AddOrReplaceReadGroups -INPUT ${sorted_bam[0]} -OUTPUT ${sorted_bam[0].baseName}_RG.bam -RGID ${id} -RGLB ${id}_LB -RGPL ILLUMINA -RGPU unit1 -RGSM ${id} --VALIDATION_STRINGENCY SILENT
# dupRemoval (picard_duplicates_removal.nf)
picard MarkDuplicates -INPUT $rg_bam -OUTPUT ${rg_bam[0].baseName}_dedup.bam -METRICS_FILE ${rg_bam[0].baseName}_DUP_metrics.txt -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT
# samtoolsDedupIndex
samtools index $dedup_bam
# realignIndel (gatk3_indel_realignment.nf)
gatk3 -T RealignerTargetCreator -R $reference -I $dedup_bam -o ${dedup_bam[0].baseName}_intervals.intervals
gatk3 -T IndelRealigner -R $reference -I $dedup_bam -targetIntervals ${dedup_bam[0].baseName}_intervals.intervals --consensusDeterminationModel USE_READS  -o ${dedup_bam[0].baseName}_realigned.bam
# samtoolsRealignedIndex
samtools index $realigned_bam
# snpCalling (bcftools_snp_calling.nf), one job per chromosome from the .fai
bcftools mpileup -Ou -f ${reference} -r $chr ${bam_files} -q 5 -I -a FMT/AD,FMT/DP | \
bcftools call -G - -f GQ -mv -Oz > variants_chr_${chr}.vcf.gz
bcftools filter -e 'AC=AN || MQ < 30' variants_chr_${chr}.vcf.gz -Oz > final_variants_chr_${chr}.vcf.gz
# concatVCFs
bcftools concat $vcfs -Oz > final_variants.vcf.gz
tabix -p vcf final_variants.vcf.gz
# reference prep: fastaIndex / gatkIndex / bwaIndex
samtools faidx $reference
picard CreateSequenceDictionary R=$reference O=${reference.baseName}.dict
bwa index -a bwtsw $reference
```

The published outputs are the realigned BAMs, `final_variants.vcf.gz` with its `.tbi`, and `combined_{windows,genes,wg}.tsv` (not replicated). Nextflow's `collect()` order is arbitrary, so RepAdapt's VCF sample order and contig order are arbitrary.

## Appendix A: conda explicit package lists from RepAdapt's images

These were reconstructed from each image's `usr/local/conda-meta` (in dependency order). Each tool's own package URL was rewritten to anaconda.org; its md5 was verified to match. All URLs returned HTTP 200 when checked. None of these lists has been test-installed yet (section 7).

To regenerate on a machine with Apptainer: `apptainer exec docker://quay.io/biocontainers/<tag> ls /usr/local/conda-meta`, then read the `url` and `md5` from each JSON. Without Apptainer, the quay.io registry API works anonymously: get a token from `https://quay.io/v2/auth?service=quay.io&scope=repository:biocontainers/<n>:pull`, fetch the manifest and the last layer blob, and read `usr/local/conda-meta/*.json` from the tarball.

### bcftools (`bcftools:1.16--hfe4b78e_1`, 31 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-12.1.0-h8d9b700_16.tar.bz2#f013cf7749536ce43d82afbffdf499ab
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-2_gnu.tar.bz2#73aaf86a425cc6e73fcf236a5a46396d
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-12.1.0-h8d9b700_16.tar.bz2#4f05bc9844f7c101e6e147dab3c88d5c
https://conda.anaconda.org/conda-forge/linux-64/libgfortran5-12.1.0-hdcd56e2_16.tar.bz2#b02605b875559ff99f04351fd5040760
https://conda.anaconda.org/conda-forge/linux-64/libgfortran-ng-12.1.0-h69a702a_16.tar.bz2#6bf15e29a20f614b18ae89368260d0a2
https://conda.anaconda.org/conda-forge/linux-64/libopenblas-0.3.21-pthreads_h78a6416_3.tar.bz2#8c5963a49b6035c40646a763293fbb35
https://conda.anaconda.org/conda-forge/linux-64/libblas-3.9.0-16_linux64_openblas.tar.bz2#d9b7a8639171f6c6fa0a983edabcfe2b
https://conda.anaconda.org/conda-forge/linux-64/libcblas-3.9.0-16_linux64_openblas.tar.bz2#20bae26d0a1db73f758fc3754cab4719
https://conda.anaconda.org/conda-forge/linux-64/gsl-2.7-he838d99_0.tar.bz2#fec079ba39c9cca093bf4c00001825de
https://conda.anaconda.org/conda-forge/linux-64/bzip2-1.0.8-h7f98852_4.tar.bz2#a1fd65c7ccbf10880423d82bca54eb54
https://conda.anaconda.org/conda-forge/linux-64/keyutils-1.6.1-h166bdaf_0.tar.bz2#30186d27e2c9fa62b45fb1476b7200e3
https://conda.anaconda.org/conda-forge/linux-64/ncurses-6.3-h27087fc_1.tar.bz2#4acfc691e64342b9dae57cf2adc63238
https://conda.anaconda.org/conda-forge/linux-64/libedit-3.1.20191231-he28a2e2_2.tar.bz2#4d331e44109e3f0e19b4cb8f9b82f3e1
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-12.1.0-ha89aaad_16.tar.bz2#6f5ba041a41eb102a1027d9e68731be7
https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2022.9.24-ha878542_0.tar.bz2#41e4e87062433e283696cf384f952ef6
https://conda.anaconda.org/conda-forge/linux-64/openssl-1.1.1q-h166bdaf_0.tar.bz2#07acc367c7fc8b716770cd5b36d31717
https://conda.anaconda.org/conda-forge/linux-64/krb5-1.19.3-h3790be6_0.tar.bz2#7d862b05445123144bec92cb1acc8ef8
https://conda.anaconda.org/conda-forge/linux-64/c-ares-1.18.1-h7f98852_0.tar.bz2#f26ef8098fab1f719c91eb760d63381a
https://conda.anaconda.org/conda-forge/linux-64/libev-4.33-h516909a_1.tar.bz2#6f8720dff19e17ce5d48cfe7f3d2f0a3
https://conda.anaconda.org/conda-forge/linux-64/libzlib-1.2.12-h166bdaf_4.tar.bz2#6a2e5b333ba57ce7eec61e90260cbb79
https://conda.anaconda.org/conda-forge/linux-64/libnghttp2-1.47.0-hdcd2b5c_1.tar.bz2#6fe9e31c2b8d0b022626ccac13e6ca3c
https://conda.anaconda.org/conda-forge/linux-64/libssh2-1.10.0-haa6b8db_3.tar.bz2#89acee135f0809a18a1f4537390aa2dd
https://conda.anaconda.org/conda-forge/linux-64/libcurl-7.85.0-h7bff187_0.tar.bz2#054fb5981fdbe031caeb612b71d85f84
https://conda.anaconda.org/conda-forge/linux-64/libdeflate-1.13-h166bdaf_0.tar.bz2#4b5bee2e957570197327d0b20a718891
https://conda.anaconda.org/conda-forge/linux-64/xz-5.2.6-h166bdaf_0.tar.bz2#2161070d867d1b1204ea749c8eec4ef0
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.12-h166bdaf_4.tar.bz2#995cc7813221edbc25a3db15357599a0
https://conda.anaconda.org/bioconda/linux-64/htslib-1.16-h6bc39ce_0.tar.bz2#503d5675802731748f348164d978c6de
https://conda.anaconda.org/conda-forge/linux-64/libnsl-2.0.0-h7f98852_0.tar.bz2#39b1328babf85c7c3a61636d9cd50206
https://conda.anaconda.org/conda-forge/linux-64/perl-5.32.1-2_h7f98852_perl5.tar.bz2#09ba115862623f00962e9809ea248f1a
https://conda.anaconda.org/bioconda/linux-64/bcftools-1.16-hfe4b78e_1.tar.bz2#be03e1f9478c474f6645b3841a3c23da
```

### fastp (`fastp:0.20.1--h8b12597_0`, 7 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-9.2.0-h24d8f2e_2.tar.bz2#7c6154111414b737f850a5d02b9c9380
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-0_gnu.tar.bz2#3053d142226093f417abe884f2ab6620
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-9.2.0-h24d8f2e_2.tar.bz2#51817a9a1e6064d0a065a93c2f4986f5
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-9.2.0-hdf63c60_2.tar.bz2#ca372846bf35fc54a853f4175ae5c776
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.11-h516909a_1006.tar.bz2#eb59ca0a6123517e35bae003be4b4bfe
https://conda.anaconda.org/bioconda/linux-64/fastp-0.20.1-h8b12597_0.tar.bz2#a30b35cd2f8693a10a56b40ee3ac8496
```

### bwa (`bwa:0.7.17--h5bf99c6_8`, 7 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-9.3.0-h2828fa1_18.tar.bz2#fc7a2a7e6a741c8afdd764715ac7039d
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-1_gnu.tar.bz2#561e277319a41d4f24f5c05a9ef63c04
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-9.3.0-h2828fa1_18.tar.bz2#5a9490c49a3505a6d19bda012cde6ad3
https://conda.anaconda.org/conda-forge/linux-64/perl-5.32.0-h36c2ea0_0.tar.bz2#79bf76d543579f7644ff59e01ed482d4
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.11-h516909a_1010.tar.bz2#339cc5584e6d26bc73a875ba900028c3
https://conda.anaconda.org/bioconda/linux-64/bwa-0.7.17-h5bf99c6_8.tar.bz2#296699e49a5c7d784cbe9dd788848fb5
```

### samtools (`samtools:1.16.1--h6899075_0`, 23 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-12.1.0-h8d9b700_16.tar.bz2#f013cf7749536ce43d82afbffdf499ab
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-2_gnu.tar.bz2#73aaf86a425cc6e73fcf236a5a46396d
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-12.1.0-h8d9b700_16.tar.bz2#4f05bc9844f7c101e6e147dab3c88d5c
https://conda.anaconda.org/conda-forge/linux-64/bzip2-1.0.8-h7f98852_4.tar.bz2#a1fd65c7ccbf10880423d82bca54eb54
https://conda.anaconda.org/conda-forge/linux-64/c-ares-1.18.1-h7f98852_0.tar.bz2#f26ef8098fab1f719c91eb760d63381a
https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2022.9.24-ha878542_0.tar.bz2#41e4e87062433e283696cf384f952ef6
https://conda.anaconda.org/conda-forge/linux-64/keyutils-1.6.1-h166bdaf_0.tar.bz2#30186d27e2c9fa62b45fb1476b7200e3
https://conda.anaconda.org/conda-forge/linux-64/ncurses-6.3-h27087fc_1.tar.bz2#4acfc691e64342b9dae57cf2adc63238
https://conda.anaconda.org/conda-forge/linux-64/libedit-3.1.20191231-he28a2e2_2.tar.bz2#4d331e44109e3f0e19b4cb8f9b82f3e1
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-12.1.0-ha89aaad_16.tar.bz2#6f5ba041a41eb102a1027d9e68731be7
https://conda.anaconda.org/conda-forge/linux-64/openssl-1.1.1q-h166bdaf_0.tar.bz2#07acc367c7fc8b716770cd5b36d31717
https://conda.anaconda.org/conda-forge/linux-64/krb5-1.19.3-h3790be6_0.tar.bz2#7d862b05445123144bec92cb1acc8ef8
https://conda.anaconda.org/conda-forge/linux-64/libev-4.33-h516909a_1.tar.bz2#6f8720dff19e17ce5d48cfe7f3d2f0a3
https://conda.anaconda.org/conda-forge/linux-64/libzlib-1.2.12-h166bdaf_3.tar.bz2#29b2d63b0e21b765da0418bc452538c9
https://conda.anaconda.org/conda-forge/linux-64/libnghttp2-1.47.0-hdcd2b5c_1.tar.bz2#6fe9e31c2b8d0b022626ccac13e6ca3c
https://conda.anaconda.org/conda-forge/linux-64/libssh2-1.10.0-haa6b8db_3.tar.bz2#89acee135f0809a18a1f4537390aa2dd
https://conda.anaconda.org/conda-forge/linux-64/libcurl-7.83.1-h7bff187_0.tar.bz2#d0c278476dba3b29ee13203784672ab1
https://conda.anaconda.org/conda-forge/linux-64/libdeflate-1.13-h166bdaf_0.tar.bz2#4b5bee2e957570197327d0b20a718891
https://conda.anaconda.org/conda-forge/linux-64/xz-5.2.6-h166bdaf_0.tar.bz2#2161070d867d1b1204ea749c8eec4ef0
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.12-h166bdaf_3.tar.bz2#76c717057865201aa2d24b79315645bb
https://conda.anaconda.org/bioconda/linux-64/htslib-1.16-h6bc39ce_0.tar.bz2#503d5675802731748f348164d978c6de
https://conda.anaconda.org/bioconda/linux-64/samtools-1.16.1-h6899075_0.tar.bz2#89ddc96b0bd29a5e34e683326c2201f2
```

### picard (`picard:2.26.3--hdfd78af_0`, 106 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-11.2.0-h1d223b6_10.tar.bz2#2cb1281aba496f83cd7ff83e1a455bbf
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-1_gnu.tar.bz2#561e277319a41d4f24f5c05a9ef63c04
https://conda.anaconda.org/conda-forge/noarch/_r-mutex-1.0.1-anacondar_1.tar.bz2#19f9db5f4f1b7f5ef5f6d67207f25f38
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-11.2.0-h1d223b6_10.tar.bz2#d363bc6e9ed11c7e0e0440c9229fffa4
https://conda.anaconda.org/conda-forge/linux-64/alsa-lib-1.2.3-h516909a_0.tar.bz2#1378b88874f42ac31b2f8e4f6975cb7b
https://conda.anaconda.org/conda-forge/linux-64/ld_impl_linux-64-2.36.1-hea4e1c9_2.tar.bz2#bd4f2e711b39af170e7ff15163fe87ee
https://conda.anaconda.org/conda-forge/noarch/kernel-headers_linux-64-2.6.32-he073ed8_14.tar.bz2#40c41dffc04c17136f02498538db1d2b
https://conda.anaconda.org/conda-forge/noarch/sysroot_linux-64-2.12-he073ed8_14.tar.bz2#78c8c32c25226732442d101d4fe1d785
https://conda.anaconda.org/conda-forge/linux-64/binutils_impl_linux-64-2.36.1-h193b22a_2.tar.bz2#32aae4265554a47ea77f7c09f86aeb3b
https://conda.anaconda.org/conda-forge/linux-64/binutils_linux-64-2.36-hf3e587d_1.tar.bz2#7d750bafc3cc7121170387583ecf0ff1
https://conda.anaconda.org/conda-forge/linux-64/libzlib-1.2.11-h36c2ea0_1013.tar.bz2#dcddf696ff5dfcab567100d691678e18
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.11-h36c2ea0_1013.tar.bz2#cf7190238072a41e9579e4476a6a60b8
https://conda.anaconda.org/conda-forge/linux-64/tk-8.6.11-h27826a3_1.tar.bz2#84e76fb280e735fec1efd2d21fd9cb27
https://conda.anaconda.org/conda-forge/linux-64/bwidget-1.9.14-ha770c72_0.tar.bz2#8e5a537ea36e426631695472f949e529
https://conda.anaconda.org/conda-forge/linux-64/bzip2-1.0.8-h7f98852_4.tar.bz2#a1fd65c7ccbf10880423d82bca54eb54
https://conda.anaconda.org/conda-forge/linux-64/c-ares-1.17.2-h7f98852_0.tar.bz2#a25871010e5104556045aa01850fbddf
https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2021.10.8-ha878542_0.tar.bz2#575611b8a84f45960e87722eeb51fa26
https://conda.anaconda.org/conda-forge/linux-64/libpng-1.6.37-h21135ba_2.tar.bz2#b6acf807307d033d4b7e758b4f44b036
https://conda.anaconda.org/conda-forge/linux-64/freetype-2.10.4-h0708190_1.tar.bz2#4a06f2ac2e5bfae7b6b245171c3f07aa
https://conda.anaconda.org/conda-forge/linux-64/libuuid-2.32.1-h7f98852_1000.tar.bz2#772d69f030955d9646d3d0eaf21d859d
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-11.2.0-he4da1e4_10.tar.bz2#4ab58cef1aea5e9e436f3e59f9982ba1
https://conda.anaconda.org/conda-forge/linux-64/icu-68.1-h58526e2_0.tar.bz2#fc7a4271dc2a7f4fd78cd63695baf7c3
https://conda.anaconda.org/conda-forge/linux-64/libiconv-1.16-h516909a_0.tar.bz2#5c0f338a513a2943c659ae619fca9211
https://conda.anaconda.org/conda-forge/linux-64/xz-5.2.5-h516909a_1.tar.bz2#33f601066901f3e1a85af3522a8113f9
https://conda.anaconda.org/conda-forge/linux-64/libxml2-2.9.12-h72842e0_0.tar.bz2#bd14fdf5b9ee5568056a40a6a2f41866
https://conda.anaconda.org/conda-forge/linux-64/fontconfig-2.13.1-hba837de_1005.tar.bz2#fd3611672eb91bc9d24fd6fb970037eb
https://conda.anaconda.org/conda-forge/linux-64/libffi-3.4.2-h9c3ff4c_4.tar.bz2#dea515312db9d4788e52c2edfc657635
https://conda.anaconda.org/conda-forge/linux-64/gettext-0.19.8.1-h73d1719_1008.tar.bz2#af49250eca8e139378f8ff0ae9e57251
https://conda.anaconda.org/conda-forge/linux-64/pcre-8.45-h9c3ff4c_0.tar.bz2#c05d1820a6d34ff07aaaab7a9b7eddaa
https://conda.anaconda.org/conda-forge/linux-64/libglib-2.68.4-h174f98d_1.tar.bz2#bbdd1d97559e052c98edd08a408218b4
https://conda.anaconda.org/conda-forge/linux-64/pthread-stubs-0.4-h36c2ea0_1001.tar.bz2#22dad4df6e8630e8dff2428f6f6a7036
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxau-1.0.9-h7f98852_0.tar.bz2#bf6f803a544f26ebbdc3bfff272eb179
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxdmcp-1.1.3-h7f98852_0.tar.bz2#be93aabceefa2fac576e971aef407908
https://conda.anaconda.org/conda-forge/linux-64/libxcb-1.13-h7f98852_1003.tar.bz2#a9371e9e40aded194dcba1447606c9a1
https://conda.anaconda.org/conda-forge/linux-64/pixman-0.40.0-h36c2ea0_0.tar.bz2#660e72c82f2e75a6b3fe6a6e75c79f19
https://conda.anaconda.org/conda-forge/linux-64/xorg-libice-1.0.10-h7f98852_0.tar.bz2#d6b0b50b49eccfe0be0373be628be0f3
https://conda.anaconda.org/conda-forge/linux-64/xorg-libsm-1.2.3-hd9c2040_1000.tar.bz2#9e856f78d5c80d5a78f61e72d1d473a3
https://conda.anaconda.org/conda-forge/linux-64/xorg-kbproto-1.0.7-h7f98852_1002.tar.bz2#4b230e8381279d76131116660f5a241a
https://conda.anaconda.org/conda-forge/linux-64/xorg-xproto-7.0.31-h7f98852_1007.tar.bz2#b4a4381d54784606820704f7b5f05a15
https://conda.anaconda.org/conda-forge/linux-64/xorg-libx11-1.7.2-h7f98852_0.tar.bz2#12a61e640b8894504326aadafccbb790
https://conda.anaconda.org/conda-forge/linux-64/xorg-xextproto-7.3.0-h7f98852_1002.tar.bz2#1e15f6ad85a7d743a2ac68dae6c82b98
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxext-1.3.4-h7f98852_1.tar.bz2#536cc5db4d0a3ba0630541aec064b5e4
https://conda.anaconda.org/conda-forge/linux-64/xorg-renderproto-0.11.1-h7f98852_1002.tar.bz2#06feff3d2634e3097ce2fe681474b534
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxrender-0.9.10-h7f98852_1003.tar.bz2#f59c1242cc1dd93e72c2ee2b360979eb
https://conda.anaconda.org/conda-forge/linux-64/cairo-1.16.0-h6cf1ce9_1008.tar.bz2#a43fb47d15e116f8be4be7e6b17ab59f
https://conda.anaconda.org/conda-forge/linux-64/ncurses-6.2-h58526e2_4.tar.bz2#509f2a21c4a09214cd737a480dfd80c9
https://conda.anaconda.org/conda-forge/linux-64/libedit-3.1.20191231-he28a2e2_2.tar.bz2#4d331e44109e3f0e19b4cb8f9b82f3e1
https://conda.anaconda.org/conda-forge/linux-64/openssl-3.0.0-h7f98852_1.tar.bz2#784c93a4856e2693d2a3ed337e04f2f0
https://conda.anaconda.org/conda-forge/linux-64/krb5-1.19.2-h48eae69_2.tar.bz2#2be27c612724e0f84c61565d8520e603
https://conda.anaconda.org/conda-forge/linux-64/libev-4.33-h516909a_1.tar.bz2#6f8720dff19e17ce5d48cfe7f3d2f0a3
https://conda.anaconda.org/conda-forge/linux-64/libnghttp2-1.43.0-ha19adfc_1.tar.bz2#f012f7ad2b1a6607dc2b863304387787
https://conda.anaconda.org/conda-forge/linux-64/libssh2-1.10.0-ha35d2d1_2.tar.bz2#e0adb0915b3e6971a84e347f998ca837
https://conda.anaconda.org/conda-forge/linux-64/libcurl-7.79.1-h494985f_1.tar.bz2#82c9d3315f9f356188e2e57aa810ee58
https://conda.anaconda.org/conda-forge/linux-64/curl-7.79.1-h494985f_1.tar.bz2#721a0c2ff787c2b683732f4082ef47d2
https://conda.anaconda.org/conda-forge/noarch/font-ttf-dejavu-sans-mono-2.37-hab24e00_0.tar.bz2#0c96522c6bdaed4b1566d11387caaf45
https://conda.anaconda.org/conda-forge/noarch/font-ttf-inconsolata-3.000-h77eed37_0.tar.bz2#34893075a5c9e55cdafac56607368fc6
https://conda.anaconda.org/conda-forge/noarch/font-ttf-source-code-pro-2.038-h77eed37_0.tar.bz2#4d59c254e01d9cde7957100457e2d5fb
https://conda.anaconda.org/conda-forge/noarch/font-ttf-ubuntu-0.83-hab24e00_0.tar.bz2#19410c3df09dfb12d1206132a1d357c5
https://conda.anaconda.org/conda-forge/noarch/fonts-conda-forge-1-0.tar.bz2#f766549260d6815b0c52253f1fb1bb29
https://conda.anaconda.org/conda-forge/noarch/fonts-conda-ecosystem-1-0.tar.bz2#fee5683a3f04bd15cbd8318b096a27ab
https://conda.anaconda.org/conda-forge/linux-64/fribidi-1.0.10-h36c2ea0_0.tar.bz2#ac7bc6a654f8f41b352b38f4051135f8
https://conda.anaconda.org/conda-forge/linux-64/libgcc-devel_linux-64-9.4.0-hd854feb_10.tar.bz2#273198c3a2d4fcbcc4223ae77980d5de
https://conda.anaconda.org/conda-forge/linux-64/libsanitizer-9.4.0-h79bfe98_10.tar.bz2#3d123a3e8b3a64f9ad3c467630aefebd
https://conda.anaconda.org/conda-forge/linux-64/gcc_impl_linux-64-9.4.0-h03d3576_10.tar.bz2#6e528e317f6bc53c84fcee1421476bb3
https://conda.anaconda.org/conda-forge/linux-64/gcc_linux-64-9.4.0-h391b98a_1.tar.bz2#1dc00eab8ce4b3341aae9f54c19b810a
https://conda.anaconda.org/conda-forge/linux-64/libgfortran5-11.2.0-h5c6108e_10.tar.bz2#ade52c46464c0632a92c73e66fb0aeed
https://conda.anaconda.org/conda-forge/linux-64/gfortran_impl_linux-64-9.4.0-h0003116_10.tar.bz2#75c922b7db6569e3f465fdaf49130d56
https://conda.anaconda.org/conda-forge/linux-64/gfortran_linux-64-9.4.0-hf0ab688_1.tar.bz2#a78f3dc7f1318710899565c50734c0f0
https://conda.anaconda.org/conda-forge/linux-64/giflib-5.2.1-h36c2ea0_2.tar.bz2#626e68ae9cc5912d6adb79d318cf962d
https://conda.anaconda.org/conda-forge/linux-64/graphite2-1.3.13-h58526e2_1001.tar.bz2#8c54672728e8ec6aa6db90cf2806d220
https://conda.anaconda.org/conda-forge/linux-64/libgfortran-ng-11.2.0-h69a702a_10.tar.bz2#eae1ff8390aef33d6341274a48228213
https://conda.anaconda.org/conda-forge/linux-64/libopenblas-0.3.17-pthreads_h8fe5266_1.tar.bz2#7f96c04618e952e0f9d94d5e07545a71
https://conda.anaconda.org/conda-forge/linux-64/libblas-3.9.0-11_linux64_openblas.tar.bz2#b8a498e2cac5746b808d5961cb584a13
https://conda.anaconda.org/conda-forge/linux-64/libcblas-3.9.0-11_linux64_openblas.tar.bz2#59bf439337c9ec59297f701e4ee97e09
https://conda.anaconda.org/conda-forge/linux-64/gsl-2.7-he838d99_0.tar.bz2#fec079ba39c9cca093bf4c00001825de
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-devel_linux-64-9.4.0-hd854feb_10.tar.bz2#02d38e05c1c5e1b925e923dea4495366
https://conda.anaconda.org/conda-forge/linux-64/gxx_impl_linux-64-9.4.0-h03d3576_10.tar.bz2#d8d6c89c5429eda816db0d2a6c31fb4c
https://conda.anaconda.org/conda-forge/linux-64/gxx_linux-64-9.4.0-h0316aca_1.tar.bz2#8c020fd66512f59c968546a840af6112
https://conda.anaconda.org/conda-forge/linux-64/harfbuzz-2.9.1-h83ec7ef_1.tar.bz2#9a9e823b2e31e84e5ce06f54ffce9d70
https://conda.anaconda.org/conda-forge/linux-64/jbig-2.1-h7f98852_2003.tar.bz2#1aa0cee79792fa97b7ff4545110b60bf
https://conda.anaconda.org/conda-forge/linux-64/jpeg-9d-h36c2ea0_0.tar.bz2#ea02ce6037dbe81803ae6123e5ba1568
https://conda.anaconda.org/conda-forge/linux-64/lerc-2.2.1-h9c3ff4c_0.tar.bz2#ea833dcaeb9e7ac4fac521f1a7abec82
https://conda.anaconda.org/conda-forge/linux-64/libdeflate-1.7-h7f98852_5.tar.bz2#10e242842cd30c59c12d79371dc0f583
https://conda.anaconda.org/conda-forge/linux-64/libwebp-base-1.2.1-h7f98852_0.tar.bz2#90607c4c0247f04ec98b48997de71c1a
https://conda.anaconda.org/conda-forge/linux-64/lz4-c-1.9.3-h9c3ff4c_1.tar.bz2#fbe97e8fa6f275d7c76a09e795adc3e6
https://conda.anaconda.org/conda-forge/linux-64/zstd-1.5.0-ha95c52a_0.tar.bz2#b56f94865e2de36abf054e7bfa499034
https://conda.anaconda.org/conda-forge/linux-64/libtiff-4.3.0-hf544144_1.tar.bz2#a65a4158716bd7d95bfa69bcfd83081c
https://conda.anaconda.org/conda-forge/linux-64/lcms2-2.12-hddcbb42_0.tar.bz2#797117394a4aa588de6d741b06fad80f
https://conda.anaconda.org/conda-forge/linux-64/liblapack-3.9.0-11_linux64_openblas.tar.bz2#00d3680586af1f0689398b080e273cbb
https://conda.anaconda.org/conda-forge/linux-64/make-4.3-hd18ef5c_1.tar.bz2#4049ebfd3190b580dffe76daed26155a
https://conda.anaconda.org/conda-forge/linux-64/xorg-inputproto-2.3.2-h7f98852_1002.tar.bz2#bcd1b3396ec6960cbc1d2855a9e60b2b
https://conda.anaconda.org/conda-forge/linux-64/xorg-fixesproto-5.0-h7f98852_1002.tar.bz2#65ad6e1eb4aed2b0611855aff05e04f6
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxfixes-5.0.3-h7f98852_1004.tar.bz2#e9a21aa4d5e3e5f1aed71e8cefd46b6a
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxi-1.7.10-h7f98852_0.tar.bz2#e77615e5141cad5a2acaa043d1cf0ca5
https://conda.anaconda.org/conda-forge/linux-64/xorg-recordproto-1.14.2-h7f98852_1002.tar.bz2#2f835e6c386e73c6faaddfe9eda67e98
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxtst-1.2.3-h7f98852_1002.tar.bz2#a220b1a513e19d5cb56c1311d44f12e6
https://conda.anaconda.org/conda-forge/linux-64/openjdk-11.0.9.1-h5cc2fde_1.tar.bz2#173813c547d6b3652fbf1db24cdc0c1b
https://conda.anaconda.org/conda-forge/linux-64/pango-1.48.10-hb8ff022_1.tar.bz2#f67c24bfd760cd50c285556ee7507853
https://conda.anaconda.org/conda-forge/linux-64/pcre2-10.37-h032f7d1_0.tar.bz2#6469e4602e914febe6f057ad2271a54e
https://conda.anaconda.org/conda-forge/linux-64/readline-8.1-h46c0cb4_0.tar.bz2#5788de3c8d7a7d64ac56c784c4ef48e6
https://conda.anaconda.org/conda-forge/linux-64/sed-4.8-he412f7d_0.tar.bz2#7362f0042e95681f5d371c46c83ebd08
https://conda.anaconda.org/conda-forge/linux-64/tktable-2.10-hb7b940f_3.tar.bz2#ea4d0879e40211fa26f38d8986db1bbe
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxt-1.2.1-h7f98852_2.tar.bz2#60d6eec5273f1c9af096c10c268912e3
https://conda.anaconda.org/conda-forge/linux-64/r-base-4.1.1-hb93adac_1.tar.bz2#393858be0a91b3df02d3cba3f0ad4b60
https://conda.anaconda.org/bioconda/noarch/picard-2.26.3-hdfd78af_0.tar.bz2#571c10d9a3768ee5b94c3c1ac90a2807
```

### gatk (`gatk:3.8--9`, 153 packages)

```
# reconstructed from usr/local/conda-meta in container layer (dependency order)
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-9.3.0-h5dbcf3e_17.tar.bz2#8fd587013b9da8b52050268d50c12305
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-1_gnu.tar.bz2#561e277319a41d4f24f5c05a9ef63c04
https://conda.anaconda.org/conda-forge/noarch/_r-mutex-1.0.1-anacondar_1.tar.bz2#19f9db5f4f1b7f5ef5f6d67207f25f38
https://conda.anaconda.org/conda-forge/linux-64/ld_impl_linux-64-2.35-h769bd43_9.tar.bz2#e91fb361f3d158f06546dc87cbe55739
https://conda.anaconda.org/conda-forge/noarch/kernel-headers_linux-64-2.6.32-h77966d4_13.tar.bz2#182b3bbe97ca334be3ccb50b80810bb1
https://conda.anaconda.org/conda-forge/noarch/sysroot_linux-64-2.12-h77966d4_13.tar.bz2#e411486a18c4f61c59083c5792c1ce3b
https://conda.anaconda.org/conda-forge/linux-64/binutils_impl_linux-64-2.35-h18a2f87_9.tar.bz2#f0f95ebd15e1eca2d43d5ab2167045e4
https://conda.anaconda.org/conda-forge/linux-64/binutils_linux-64-2.35-hc3fd857_29.tar.bz2#17d622904723a89e59ac9251aa432078
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-9.3.0-h5dbcf3e_17.tar.bz2#fc9f5adabc4d55cd4b491332adc413e0
https://conda.anaconda.org/conda-forge/linux-64/zlib-1.2.11-h516909a_1010.tar.bz2#339cc5584e6d26bc73a875ba900028c3
https://conda.anaconda.org/conda-forge/linux-64/tk-8.6.10-hed695b0_1.tar.bz2#7ef837cd455bd0f19f49b8b62d4cb568
https://conda.anaconda.org/conda-forge/linux-64/bwidget-1.9.14-ha770c72_0.tar.bz2#8e5a537ea36e426631695472f949e529
https://conda.anaconda.org/conda-forge/linux-64/bzip2-1.0.8-h516909a_3.tar.bz2#a05ea3fc1a51cf629bd49b481f729ebd
https://conda.anaconda.org/conda-forge/linux-64/c-ares-1.16.1-h516909a_3.tar.bz2#8d0f54b0a09bb496dea3f8dae0c551e4
https://conda.anaconda.org/conda-forge/linux-64/ca-certificates-2020.6.20-hecda079_0.tar.bz2#1b1cca86e95c416a8e7eb6062af6d503
https://conda.anaconda.org/conda-forge/linux-64/libpng-1.6.37-h21135ba_2.tar.bz2#b6acf807307d033d4b7e758b4f44b036
https://conda.anaconda.org/conda-forge/linux-64/freetype-2.10.4-h7ca028e_0.tar.bz2#1e6bf409916fdf455a5bb8627b626285
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-9.3.0-h2ae2ef3_17.tar.bz2#342f3c931d0a3a209ab09a522469d20c
https://conda.anaconda.org/conda-forge/linux-64/icu-67.1-he1b5a44_0.tar.bz2#7ced6a5e5c94726af797d2b5a2b09228
https://conda.anaconda.org/conda-forge/linux-64/libuuid-2.32.1-h14c3975_1000.tar.bz2#39c6326f6ee5297632c47db6520546fe
https://conda.anaconda.org/conda-forge/linux-64/libiconv-1.16-h516909a_0.tar.bz2#5c0f338a513a2943c659ae619fca9211
https://conda.anaconda.org/conda-forge/linux-64/xz-5.2.5-h516909a_1.tar.bz2#33f601066901f3e1a85af3522a8113f9
https://conda.anaconda.org/conda-forge/linux-64/libxml2-2.9.10-h68273f3_2.tar.bz2#0315cae0468a1e17f1e7fad5b13d53f8
https://conda.anaconda.org/conda-forge/linux-64/fontconfig-2.13.1-h7e3eb15_1002.tar.bz2#24b40f20bd46b6f75872dbff650fa129
https://conda.anaconda.org/conda-forge/linux-64/libffi-3.2.1-he1b5a44_1007.tar.bz2#11389072d7d6036fd811c3d9460475cd
https://conda.anaconda.org/conda-forge/linux-64/gettext-0.19.8.1-hf34092f_1004.tar.bz2#5582e1349bee4a25705adca745bf6845
https://conda.anaconda.org/conda-forge/linux-64/pcre-8.44-he1b5a44_0.tar.bz2#e647d89cd5cdf62760cf283a001841ff
https://conda.anaconda.org/conda-forge/linux-64/libglib-2.66.2-hbe7bbb4_0.tar.bz2#9dbfa33d92fd7e5f6c0e51802c47c1a3
https://conda.anaconda.org/conda-forge/linux-64/pthread-stubs-0.4-h14c3975_1001.tar.bz2#e0b9987b65bba8016d21f5bdbe0a7923
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxau-1.0.9-h14c3975_0.tar.bz2#ffa7c2b7a2c7dc779ed9e38b10a93c3c
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxdmcp-1.1.3-h516909a_0.tar.bz2#e95a160e60b2a327309a6d323a4d780e
https://conda.anaconda.org/conda-forge/linux-64/libxcb-1.13-h14c3975_1002.tar.bz2#ca4eb860b5528d5c6de8d97021d9ef78
https://conda.anaconda.org/conda-forge/linux-64/pixman-0.40.0-h36c2ea0_0.tar.bz2#660e72c82f2e75a6b3fe6a6e75c79f19
https://conda.anaconda.org/conda-forge/linux-64/xorg-libice-1.0.10-h516909a_0.tar.bz2#4dfda1ccfd0cc90c4fe6786acece9a30
https://conda.anaconda.org/conda-forge/linux-64/xorg-libsm-1.2.3-h84519dc_1000.tar.bz2#56cc238b81624db9e3e36b2b6a482cc1
https://conda.anaconda.org/conda-forge/linux-64/xorg-kbproto-1.0.7-h14c3975_1002.tar.bz2#6dfe5dbe10d55266e4a5e89287eed578
https://conda.anaconda.org/conda-forge/linux-64/xorg-xproto-7.0.31-h14c3975_1007.tar.bz2#a45d8cd411bdf8f08ced463f68986b62
https://conda.anaconda.org/conda-forge/linux-64/xorg-libx11-1.6.12-h516909a_0.tar.bz2#cdea71c3e27611e012419faea02c56bd
https://conda.anaconda.org/conda-forge/linux-64/xorg-xextproto-7.3.0-h14c3975_1002.tar.bz2#f08999859c405bad87c4bf9b6cdc7bbb
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxext-1.3.4-h516909a_0.tar.bz2#aba40ec20012e0a2641ced275a0abdca
https://conda.anaconda.org/conda-forge/linux-64/xorg-renderproto-0.11.1-h14c3975_1002.tar.bz2#fbcb7fa11dee1a5d3df4371cc55bb229
https://conda.anaconda.org/conda-forge/linux-64/xorg-libxrender-0.9.10-h516909a_1002.tar.bz2#bf3514891c6743ed808602fa9f4f1508
https://conda.anaconda.org/conda-forge/linux-64/cairo-1.16.0-h9f066cc_1006.tar.bz2#9f8ef44d205f35056dccc1a92fbe6592
https://conda.anaconda.org/conda-forge/linux-64/ncurses-6.2-he1b5a44_2.tar.bz2#48af29b07d7336037d00cdea3e11775d
https://conda.anaconda.org/conda-forge/linux-64/openssl-1.1.1h-h516909a_0.tar.bz2#3a99e0cb8f325dbf8f616da2d2fb6d4f
https://conda.anaconda.org/conda-forge/linux-64/readline-8.0-he28a2e2_2.tar.bz2#4d0ae8d473f863696088f76800ef9d38
https://conda.anaconda.org/conda-forge/linux-64/sqlite-3.33.0-h4cf870e_1.tar.bz2#1968ff6c4b8fbd2903f6672e292d932d
https://conda.anaconda.org/conda-forge/noarch/tzdata-2020d-h516909a_0.tar.bz2#153abe529c372caadf52bde5e07656e3
https://conda.anaconda.org/conda-forge/linux-64/python-3.9.0-h2a148a8_4_cpython.tar.bz2#cb273019874ded6e4dc57713f7c89bf6
https://conda.anaconda.org/conda-forge/linux-64/python_abi-3.9-1_cp39.tar.bz2#9ee3692d976902241a3392495768fe98
https://conda.anaconda.org/conda-forge/linux-64/certifi-2020.6.20-py39h079e4ff_2.tar.bz2#a5532022d85afa3f5cbc93dcd491495c
https://conda.anaconda.org/conda-forge/linux-64/libedit-3.1.20191231-he28a2e2_2.tar.bz2#4d331e44109e3f0e19b4cb8f9b82f3e1
https://conda.anaconda.org/conda-forge/linux-64/krb5-1.17.1-hfafb76e_3.tar.bz2#b9c0993124fbf5f4ccf37fd00d6a3705
https://conda.anaconda.org/conda-forge/linux-64/libev-4.33-h516909a_1.tar.bz2#6f8720dff19e17ce5d48cfe7f3d2f0a3
https://conda.anaconda.org/conda-forge/linux-64/libnghttp2-1.41.0-h8cfc5f6_2.tar.bz2#726ca0fed4bde95b056ef26df1efaf60
https://conda.anaconda.org/conda-forge/linux-64/libssh2-1.9.0-hab1572f_5.tar.bz2#18aaa1bd2238ae2b5e89591046973123
https://conda.anaconda.org/conda-forge/linux-64/libcurl-7.71.1-hcdd3856_8.tar.bz2#2c81fb9f0d82c04b08617fc73eb615af
https://conda.anaconda.org/conda-forge/linux-64/curl-7.71.1-he644dc0_8.tar.bz2#7fadfa44182e1cc0dbbe0c169c0b15cc
https://conda.anaconda.org/conda-forge/linux-64/fribidi-1.0.10-h36c2ea0_0.tar.bz2#ac7bc6a654f8f41b352b38f4051135f8
https://conda.anaconda.org/conda-forge/linux-64/openjdk-8.0.265-h516909a_0.tar.bz2#2c42874842fd46ba127b443119c8464f
https://conda.anaconda.org/conda-forge/linux-64/libgcc-devel_linux-64-9.3.0-hfd08b2a_17.tar.bz2#0134fffc3c28ccfef3d42ecd8671613a
https://conda.anaconda.org/conda-forge/linux-64/gcc_impl_linux-64-9.3.0-h28f5a38_17.tar.bz2#40235a65140e7f9a10716f555f7ce409
https://conda.anaconda.org/conda-forge/linux-64/gcc_linux-64-9.3.0-h7247604_29.tar.bz2#081695d5c895c0491ce4d8e8d4ade215
https://conda.anaconda.org/conda-forge/linux-64/libgfortran5-9.3.0-he4bcb1c_17.tar.bz2#0c15349375fc3d0cb2114fcabe2f0aba
https://conda.anaconda.org/conda-forge/linux-64/gfortran_impl_linux-64-9.3.0-h2bb4189_17.tar.bz2#33f393e20a07fb36469482b5fe2cbb3b
https://conda.anaconda.org/conda-forge/linux-64/gfortran_linux-64-9.3.0-ha1c937c_29.tar.bz2#87b2ba97e70fe7801c36b8472390df5b
https://conda.anaconda.org/conda-forge/linux-64/libgfortran-ng-9.3.0-he4bcb1c_17.tar.bz2#f92019e2b944dc3e5d33d0efba4a3461
https://conda.anaconda.org/conda-forge/linux-64/libopenblas-0.3.12-pthreads_h4812303_1.tar.bz2#2c7126a584f05e7bc8885d64dad3d21a
https://conda.anaconda.org/conda-forge/linux-64/libblas-3.9.0-2_openblas.tar.bz2#b0bc5ed57c1e8fbfaf25a2a011e89392
https://conda.anaconda.org/conda-forge/linux-64/libcblas-3.9.0-2_openblas.tar.bz2#e842cfe553dd36f59b2770be0b7d3cd2
https://conda.anaconda.org/conda-forge/linux-64/gsl-2.6-hf94e986_0.tar.bz2#21e191a5561d4d22602112e2cbbc5197
https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-devel_linux-64-9.3.0-h4084dd6_17.tar.bz2#88192ce1087609c6e27712ae81776d88
https://conda.anaconda.org/conda-forge/linux-64/gxx_impl_linux-64-9.3.0-h53cdd4c_17.tar.bz2#eb5cff7db946ed93fea593dc8eff1e32
https://conda.anaconda.org/conda-forge/linux-64/gxx_linux-64-9.3.0-h0d07fa4_29.tar.bz2#be1da0f94d85f1125363cc252c446542
https://conda.anaconda.org/conda-forge/linux-64/jpeg-9d-h36c2ea0_0.tar.bz2#ea02ce6037dbe81803ae6123e5ba1568
https://conda.anaconda.org/conda-forge/linux-64/liblapack-3.9.0-2_openblas.tar.bz2#51be6ef77752c43c1371008a1ed787f1
https://conda.anaconda.org/conda-forge/linux-64/libwebp-base-1.1.0-h36c2ea0_3.tar.bz2#5d2fa0f54b5cc9b500269437a62cba3b
https://conda.anaconda.org/conda-forge/linux-64/lz4-c-1.9.2-he1b5a44_3.tar.bz2#b2e54aad8640e7a877d2280d3ebfe85b
https://conda.anaconda.org/conda-forge/linux-64/zstd-1.4.5-h6597ccf_2.tar.bz2#d60d50f369d40f9787878cdd866fc9d3
https://conda.anaconda.org/conda-forge/linux-64/libtiff-4.1.0-h4f3a223_6.tar.bz2#26f3a89355b93d26418088a749fd9880
https://conda.anaconda.org/conda-forge/linux-64/make-4.3-hd18ef5c_1.tar.bz2#4049ebfd3190b580dffe76daed26155a
https://conda.anaconda.org/conda-forge/linux-64/graphite2-1.3.13-h58526e2_1001.tar.bz2#8c54672728e8ec6aa6db90cf2806d220
https://conda.anaconda.org/conda-forge/linux-64/harfbuzz-2.7.2-hb1ce69c_1.tar.bz2#6ad4d845564ca39dc15f0a5f2b676c2b
https://conda.anaconda.org/conda-forge/linux-64/pango-1.42.4-h80147aa_5.tar.bz2#8dc418e9db17340d77ca9dda93110013
https://conda.anaconda.org/conda-forge/linux-64/sed-4.8-hbfbb72e_0.tar.bz2#dd25e60fc346d08361d109d181d8b467
https://conda.anaconda.org/conda-forge/linux-64/tktable-2.10-hb7b940f_3.tar.bz2#ea4d0879e40211fa26f38d8986db1bbe
https://conda.anaconda.org/conda-forge/linux-64/r-base-3.6.3-hc603457_4.tar.bz2#8ea60b6caca09e4dc615c2ead4e8ac0e
https://conda.anaconda.org/conda-forge/linux-64/r-digest-0.6.27-r36h1b71b39_0.tar.bz2#64426e42f24294f6a6f20f2728b68c25
https://conda.anaconda.org/conda-forge/linux-64/r-glue-1.4.2-r36hcdcec82_0.tar.bz2#d8c13b22bbca3c3c2ec34956ca3fccdf
https://conda.anaconda.org/conda-forge/noarch/r-gtable-0.3.0-r36h6115d3f_3.tar.bz2#9888efe05ffb140c2256296cde37225f
https://conda.anaconda.org/conda-forge/linux-64/r-rcpp-1.0.5-r36h0357c0b_0.tar.bz2#1743ff0bd984ad2ec304fcfbcec6bb5c
https://conda.anaconda.org/conda-forge/linux-64/r-brio-1.1.0-r36h9e2df91_1.tar.bz2#6b24231ddaf10da54ef12d88bf6c512a
https://conda.anaconda.org/conda-forge/linux-64/r-ps-1.4.0-r36h0eb13af_0.tar.bz2#3c0eb479845fea39b293845747b587b8
https://conda.anaconda.org/conda-forge/noarch/r-r6-2.5.0-r36hc72bb7e_0.tar.bz2#f72d5aeb259db2c40f4305257202dfbc
https://conda.anaconda.org/conda-forge/linux-64/r-processx-3.4.4-r36hcdcec82_0.tar.bz2#0ea2b823427da665735b697d605ce407
https://conda.anaconda.org/conda-forge/noarch/r-callr-3.5.1-r36h142f84f_0.tar.bz2#25570fb40a0d29f0c8531ad9a54287a2
https://conda.anaconda.org/conda-forge/noarch/r-assertthat-0.2.1-r36h6115d3f_2.tar.bz2#89894bf1427eea74d7adb24bf93858dc
https://conda.anaconda.org/conda-forge/noarch/r-crayon-1.3.4-r36h6115d3f_1003.tar.bz2#f744faf7c73cc8abbe1b87a25a24878e
https://conda.anaconda.org/conda-forge/linux-64/r-fansi-0.4.1-r36hcdcec82_1.tar.bz2#22b36278eae6a6a917d9eba12325c4d3
https://conda.anaconda.org/conda-forge/noarch/r-cli-2.1.0-r36h142f84f_0.tar.bz2#723e86c45cb8aaf3fa7f3a274a7fb329
https://conda.anaconda.org/conda-forge/linux-64/r-backports-1.2.0-r36h9e2df91_0.tar.bz2#5257cc47bb62f7df7fd2b8ccb4ac9b9b
https://conda.anaconda.org/conda-forge/noarch/r-rprojroot-1.3_2-r36h6115d3f_1003.tar.bz2#20938c6c697d99831822775109a94306
https://conda.anaconda.org/conda-forge/noarch/r-desc-1.2.0-r36h6115d3f_1003.tar.bz2#ad3140f8a680f33dfe59fd4fec7201b2
https://conda.anaconda.org/conda-forge/linux-64/r-rlang-0.4.8-r36h9e2df91_0.tar.bz2#58d00701a8f4465fe87a8c029cab4e62
https://conda.anaconda.org/conda-forge/linux-64/r-ellipsis-0.3.1-r36hcdcec82_0.tar.bz2#041d3c7b6c127aaa50ce32e11e2f4b74
https://conda.anaconda.org/conda-forge/noarch/r-evaluate-0.14-r36h6115d3f_2.tar.bz2#de3bec436333836c3729452b91161d2b
https://conda.anaconda.org/conda-forge/linux-64/r-jsonlite-1.7.1-r36hcdcec82_0.tar.bz2#abc9f72c705463f1484f53b4c52f137c
https://conda.anaconda.org/conda-forge/noarch/r-lifecycle-0.2.0-r36h6115d3f_1.tar.bz2#ff33d80f55b80fda34f0f8980af58b7b
https://conda.anaconda.org/conda-forge/noarch/r-magrittr-1.5-r36h6115d3f_1003.tar.bz2#89a866d76b41e98cd5bfe458007877e8
https://conda.anaconda.org/conda-forge/noarch/r-prettyunits-1.1.1-r36h6115d3f_1.tar.bz2#c6ff88b560efa942ba8d0ae60b2727af
https://conda.anaconda.org/conda-forge/noarch/r-withr-2.3.0-r36h6115d3f_0.tar.bz2#5aa918d15567c119298541b7b2afa38a
https://conda.anaconda.org/conda-forge/noarch/r-pkgbuild-1.1.0-r36h6115d3f_0.tar.bz2#bbc33a8ed43a342d24bc2346d2b25b31
https://conda.anaconda.org/conda-forge/noarch/r-rstudioapi-0.11-r36h6115d3f_1.tar.bz2#8ed5ad583302bed8cf4b1cfe8858ef48
https://conda.anaconda.org/conda-forge/linux-64/r-pkgload-1.1.0-r36h0357c0b_0.tar.bz2#2b8830644030cdd9dba29758dc0b51d0
https://conda.anaconda.org/conda-forge/noarch/r-praise-1.0.0-r36h6115d3f_1004.tar.bz2#c1c81352a7102756ef1fa1838d25ad95
https://conda.anaconda.org/conda-forge/linux-64/r-diffobj-0.3.2-r36h9e2df91_1.tar.bz2#e5ea8e90fb3b699b38bb72ebf4cfe718
https://conda.anaconda.org/conda-forge/linux-64/r-utf8-1.1.4-r36hcdcec82_1003.tar.bz2#1d49c6ca95781f9a58907f5a46199ad1
https://conda.anaconda.org/conda-forge/noarch/r-zeallot-0.1.0-r36h6115d3f_1002.tar.bz2#12fd0aa909c5641bf8ae36ab462a8ccc
https://conda.anaconda.org/conda-forge/linux-64/r-vctrs-0.3.4-r36hcdcec82_0.tar.bz2#3519560dd4589460711b3a9bca0c2b66
https://conda.anaconda.org/conda-forge/noarch/r-pillar-1.4.6-r36h6115d3f_0.tar.bz2#0562c96c3feab8bc4212374cdad86071
https://conda.anaconda.org/conda-forge/noarch/r-pkgconfig-2.0.3-r36h6115d3f_1.tar.bz2#d8c0bd8e2b4de78157be2866e6cc320b
https://conda.anaconda.org/conda-forge/linux-64/r-tibble-3.0.4-r36h0eb13af_0.tar.bz2#93c5e53bd0debc53798c82b006e67ba5
https://conda.anaconda.org/conda-forge/noarch/r-rematch2-2.1.2-r36h6115d3f_1.tar.bz2#ce9c3aa34838ec179fd0703123b02392
https://conda.anaconda.org/conda-forge/noarch/r-waldo-0.2.2-r36hc72bb7e_0.tar.bz2#5df3512a74210ae30603035b2fc6b8d5
https://conda.anaconda.org/conda-forge/linux-64/r-testthat-3.0.0-r36he524a50_0.tar.bz2#551c84032f79428d878a341ff96c9601
https://conda.anaconda.org/conda-forge/linux-64/r-isoband-0.2.2-r36h0357c0b_0.tar.bz2#7a0b07a45ca6e00f779918bac81be137
https://conda.anaconda.org/conda-forge/linux-64/r-mass-7.3_53-r36hcdcec82_0.tar.bz2#c3b5fb00c4abd18691a929d1c5b33f33
https://conda.anaconda.org/conda-forge/linux-64/r-lattice-0.20_41-r36hcdcec82_2.tar.bz2#b8bdf1e73795dab6125e9c32e7ef749e
https://conda.anaconda.org/conda-forge/linux-64/r-matrix-1.2_18-r36h7fa42b6_3.tar.bz2#df7f26b2aa8df312ec0ad53deed7deed
https://conda.anaconda.org/conda-forge/linux-64/r-nlme-3.1_150-r36h580db52_0.tar.bz2#51978e7e545886fef674d9d7799ebfc1
https://conda.anaconda.org/conda-forge/linux-64/r-mgcv-1.8_33-r36h7fa42b6_0.tar.bz2#2986bd026e3eff8ef6091559e191ce56
https://conda.anaconda.org/conda-forge/linux-64/r-farver-2.0.3-r36h0357c0b_1.tar.bz2#e979b708e9dbef4cb0d5728aa385d4e2
https://conda.anaconda.org/conda-forge/noarch/r-labeling-0.4.2-r36h142f84f_0.tar.bz2#0d9680449e3264b41b924b8fa52fdf90
https://conda.anaconda.org/conda-forge/linux-64/r-colorspace-1.4_1-r36hcdcec82_2.tar.bz2#7af73a139316d9ce8e46fba398928bc1
https://conda.anaconda.org/conda-forge/noarch/r-munsell-0.5.0-r36h6115d3f_1003.tar.bz2#6ea109782455a5738b7093889138931d
https://conda.anaconda.org/conda-forge/noarch/r-rcolorbrewer-1.1_2-r36h6115d3f_1003.tar.bz2#2ba73e5efd87bf6c2e483f6cdedeae91
https://conda.anaconda.org/conda-forge/noarch/r-viridislite-0.3.0-r36h6115d3f_1003.tar.bz2#9bd35b3bef17761b2608772f9c85cff4
https://conda.anaconda.org/conda-forge/noarch/r-scales-1.1.1-r36h6115d3f_0.tar.bz2#92945b973d607aec324d9bb4cad2e772
https://conda.anaconda.org/conda-forge/noarch/r-ggplot2-3.3.2-r36h6115d3f_0.tar.bz2#831e6f24739fccf44f2750d538227b58
https://conda.anaconda.org/conda-forge/linux-64/r-bitops-1.0_6-r36hcdcec82_1004.tar.bz2#85c5a1f6b7be5a1f33132435459711d0
https://conda.anaconda.org/conda-forge/linux-64/r-catools-1.18.0-r36h0357c0b_1.tar.bz2#e8d8d549014c61b4a80a2b6f735429e5
https://conda.anaconda.org/conda-forge/linux-64/r-gtools-3.8.2-r36hcdcec82_1.tar.bz2#9e7577788c2fceafe315c2fb4e3b8c11
https://conda.anaconda.org/conda-forge/noarch/r-gdata-2.18.0-r36h6115d3f_1003.tar.bz2#fd5f680b2403a9e0d5f665fe8f48e40d
https://conda.anaconda.org/conda-forge/linux-64/r-kernsmooth-2.23_18-r36h742201e_0.tar.bz2#7d38f23a89e8a594b04b5261a1299e0a
https://conda.anaconda.org/conda-forge/noarch/r-gplots-3.1.0-r36h6115d3f_0.tar.bz2#0b26b785542578b85b6d0b78de889b63
https://conda.anaconda.org/conda-forge/noarch/r-gsalib-2.1-r36_1002.tar.bz2#74d56a3d14e51f5844f92ba6e1a74299
https://conda.anaconda.org/conda-forge/linux-64/r-plyr-1.8.6-r36h0357c0b_1.tar.bz2#666de5043337b67331958586b4063f32
https://conda.anaconda.org/conda-forge/linux-64/r-reshape-0.8.8-r36hcdcec82_2.tar.bz2#55d99d916bb652d66775454c2bf5e81f
https://conda.anaconda.org/bioconda/noarch/gatk-3.8-9.tar.bz2#abbf30321eca326fdbb8f417d17e8673
https://conda.anaconda.org/conda-forge/linux-64/setuptools-49.6.0-py39h079e4ff_2.tar.bz2#7cebf33355b537e5706bfe0a546bb193
https://conda.anaconda.org/conda-forge/noarch/wheel-0.35.1-pyh9f0ad1d_0.tar.bz2#126827869be32f21872a2b30ebe2b038
https://conda.anaconda.org/conda-forge/noarch/pip-20.2.4-py_0.tar.bz2#d2c0e7b7ca15440dc445e725f1e79ccf
```


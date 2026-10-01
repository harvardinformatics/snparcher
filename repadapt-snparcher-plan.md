# RepAdapt support in snpArcher: goals and approach

This is a human-readable summary. The agent-facing version, with research details, the codebase map, gotchas and pinned package lists, is `agents-plan.md`.

## Why

[RepAdapt](https://github.com/RepAdapt/nextflow_snp_calling_linux) calls variants with a Nextflow pipeline: fastp, bwa, a MAPQ filter, Picard duplicate removal, GATK3 indel realignment, then bcftools mpileup/call with a minimal filter. We want two things:

1. Run RepAdapt's **calling model** inside snpArcher, on any snpArcher BAMs.
2. Go from a **snpArcher sample sheet to VCFs equivalent to RepAdapt's**, while keeping what snpArcher does better: per-row mapping, handling several libraries per sample, SRA downloads, the QC report, callable sites, the postprocess and QC modules, and failing loudly instead of silently dropping samples.

**The goal is equivalent VCFs.** Byte-identical BAMs, RepAdapt's depth-of-coverage tables (which need a GFF, and many species don't have one) and a literal Nextflow replica are not goals.

## What we're building

### 1. A RepAdapt calling model: `variant_calling.tool: repadapt`

This is an alternative BAM → VCF path.

- **Per contig**, RepAdapt's bcftools command: `mpileup -I -a FMT/AD,FMT/DP`, then `call -G - -f GQ -mv`.
  - It uses `-q 10` in place of RepAdapt's BAM-level MAPQ ≥ 10 filter.
  - bcftools is pinned to **1.16**, RepAdapt's version.
- **Raw calls** go to `results/vcfs/raw.vcf.gz`.
- **A separate RepAdapt filter step** writes `results/vcfs/filtered.vcf.gz` with two soft filters: `AllHomAlt` (every called genotype is alt) and `LowMQ` (MQ < 30). The PASS records are exactly what RepAdapt keeps.
- **Everything downstream keeps working:** QC report, callable sites, and the postprocess and QC modules.

### 2. Pluggable mapping pipelines: `mapping.pipeline`

"fastq → final BAM" becomes a named, pluggable stage:

| Pipeline | What it is |
|---|---|
| `default` | Today's snpArcher mapping: fastp, `bwa mem -M`, sambamba duplicate marking per library |
| `sentieon` | Today's Sentieon mapping, moved into the new structure. Chosen automatically with `tool: sentieon`. |
| `repadapt` | RepAdapt-equivalent processing, below |

Every pipeline delivers a final BAM per sample, plus the per-sample mapping QC in the same format. Callable sites, QC and variant calling don't need to know which pipeline ran. New variants, for example keeping today's `-M` behavior available if issue #346 changes the default, can be added as new pipeline names.

**What the `repadapt` pipeline does.** These are RepAdapt's pinned tools and result-affecting arguments, run inside snpArcher's structure:

| Step | RepAdapt setting we reproduce | snpArcher structure we keep |
|---|---|---|
| Trim | fastp 0.20.1 with RepAdapt's defaults (no extra adapter detection) | per row, snpArcher threads, JSON report |
| Map | bwa 0.7.17, no `-M`, batch size matched to RepAdapt's (`-K 40000000`) | per row, read groups at mapping, snpArcher threads |
| Filter | keep MAPQ ≥ 10, then `fixmate`, exactly as RepAdapt | streamed, with no SAM files |
| Duplicates | Picard 2.26.3, duplicates **removed**, **once per sample** | the sheet's `mark_duplicates` is still honored |
| Realign | GATK 3.8-1 indel realignment | skipped automatically, with a warning, for very long contigs |

**QC stays correct.** Mapped and properly-paired percentages come from the alignments *before* the MAPQ filter, and the duplication rate comes from Picard's metrics. Coverage in the QC report and callable sites reflects the reads actually used for calling.

### Example configs

RepAdapt calling on snpArcher's usual BAMs:
```yaml
variant_calling:
  tool: repadapt
```

RepAdapt-equivalent from reads to VCF:
```yaml
mapping:
  pipeline: repadapt
  repadapt:
    indel_realignment: auto   # auto | true | false
variant_calling:
  tool: repadapt
```

## Key decisions

| Decision | Choice | Why |
|---|---|---|
| BAM-level MAPQ filter in the repadapt pipeline | **Keep** | It changes what duplicate removal and realignment see, not only calling. In the general calling model, BAM filtering is avoided; mpileup `-q 10` stands in. |
| bcftools version for the calling model | **Pin 1.16** | Versions differ in how MQ is stored (integer or float), which can move sites across the MQ < 30 cutoff. |
| fastp threads | **snpArcher's thread count** | Thread count only changes read order, and RepAdapt's own order varies between runs too. |
| Duplicate removal | **Once per sample** | This matches RepAdapt. It only matters for samples with several libraries, where the effect is small (roughly 0.1% of reads, near repeats). |
| Depth tables / GFF | **Dropped** | Many species have no GFF, and snpArcher's coverage outputs cover this need. |
| Software | **Pinned conda envs** with the exact builds from RepAdapt's images | No containers needed at runtime. The trade-off is that full runs work on Linux only. |

## What "equivalent" means, and the differences that remain

With `mapping.pipeline: repadapt` and `tool: repadapt`, on single-row, single-library samples, we expect the same VCF records as RepAdapt, apart from:

- **Order:** sample and contig order. RepAdapt's is arbitrary; ours follows the sample sheet and the reference index.
- **Headers:** bcftools writes dates into header lines, so headers always differ.
- **Run-to-run variation:** RepAdapt isn't exactly reproducible on large inputs, because fastp's multi-threaded output order changes between runs. We share that property.
- **Extreme-depth sites:** at sites with extreme depth (collapsed repeats far above mpileup's 250-read cap), bcftools' PL and GQ depend on the order of the BAMs on the command line. Ours is the sample-sheet order; RepAdapt's is Nextflow's arbitrary collect order. On 10 *C. albicans* isolates called from RepAdapt's BAMs, 4 of 296,311 records differed, all in the collapsed rDNA, and none differed once we used RepAdapt's order.
- **Situations RepAdapt can't handle,** where snpArcher's behavior applies:
  - several rows per sample (mapped per row, then merged)
  - BAM inputs (used as they are)
  - `mark_duplicates: false` (duplicate removal skipped)

The calling model alone, used with the default mapping pipeline, is "RepAdapt-adjacent". The remaining differences are all upstream of calling: `-M`, no BAM-level filter, sambamba instead of Picard, and no realignment.

## Limitations

- The pinned tools only exist for Linux, so full runs don't work on macOS; dry runs do.
- GATK3 can't handle contigs longer than about 537 Mb (2^29 bp), so realignment is skipped for them. It's also slow (single-threaded per sample) and, per RepAdapt's own README, struggles with very fragmented references.

## How we'll deliver it

The PRs are nested, each targeting the `feat/repadapt` integration branch, which merges to `main` at the end:

0. **Test fixture and golden output:** a small simulated dataset, plus the real RepAdapt pipeline run once on it to produce reference results.
1. **Calling model:** `tool: repadapt`.
2. **Mapping-pipeline refactor,** with no behavior change: `default` and `sentieon` moved into the new structure.
3. **The `repadapt` mapping pipeline.**

**Testing:** PASS records from our pipeline must match RepAdapt's output on the fixture. That's checked twice:
- calling on RepAdapt's own BAMs, which checks the calling model
- calling from reads, which checks the whole repadapt path

Dry-run tests pin the exact commands.

## Background from the planning discussion

- **What `bcftools call -G -` does:** it calls each sample independently. A site is kept if any sample supports it, and each sample's genotype prior comes from its own reads. This avoids assuming one randomly mating population, which suits RepAdapt's structured sample sets. The costs are more false-positive sites at low coverage and genotypes that lean homozygous at low depth.
- **How GATK compares:** GATK joint genotyping detects sites using the whole cohort but assigns each sample's genotype on its own likelihoods, with no population prior.
- **bwa `-M`:** it's a GATK3-era convention. Current Broad pipelines drop it, add `-Y`, and fix `-K` for reproducibility. snpArcher's default mapping currently uses `-M` without `-K`, which makes results depend on thread count. That's tracked separately in issue #346.

## Possible future work

- A `legacy` mapping pipeline, if #346 changes the default flags.
- Allowing the Sentieon mapping pipeline with other callers.
- Speeding up GATK3 realignment by splitting it into pieces, once that's proven to give identical output.
- Adding AD/DP/GQ annotations to the plain `bcftools` caller.
- Documenting and checking that DeepVariant needs Apptainer (tracked separately).

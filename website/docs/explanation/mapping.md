# Mapping pipelines

A mapping pipeline turns each sample's reads into one coordinate-sorted final BAM, plus per-sample mapping QC.
Everything around it is shared and works the same whichever pipeline runs: SRA downloads, read staging, samples supplied as BAMs, callable sites, variant calling and the QC report.
You choose the pipeline with `mapping.pipeline`.

| Pipeline | What it does | Use with |
|---|---|---|
| `default` | fastp, `bwa mem -M`, per-library duplicate marking with sambamba | Any caller except Sentieon's |
| `sentieon` | Sentieon's `bwa mem` and Dedup | `tool: sentieon`, which selects it automatically |
| `repadapt` | RepAdapt's processing with its pinned tools, below | Any caller except Sentieon's; with `tool: repadapt` for RepAdapt-equivalent calls |

Samples supplied as BAMs (`input_type: bam`) are used as they are, whichever pipeline is set.

## The `repadapt` pipeline

[RepAdapt](https://github.com/RepAdapt/nextflow_snp_calling_linux)'s SNP-calling pipeline processes reads with a fixed set of tool versions and arguments.
`mapping.pipeline: repadapt` reproduces the steps that affect the BAMs, with RepAdapt's pinned tools, inside snpArcher's structure: rows are mapped separately, read groups are set at mapping, and snpArcher's thread counts apply.

| Step | Tool | What runs |
|---|---|---|
| Trim | fastp 0.20.1 | RepAdapt's defaults, without `--detect_adapter_for_pe` |
| Map | bwa 0.7.17 | `bwa mem -K 40000000`, without `-M`. `-K` reproduces RepAdapt's batch size at any thread count; bwa estimates insert sizes per batch. |
| Filter | samtools 1.16.1 | Keep MAPQ >= 10, then name sort, `fixmate -m` and coordinate sort, in RepAdapt's order. Reads whose mate was filtered out become single-end reads. |
| Remove duplicates | Picard 2.26.3 | `MarkDuplicates -REMOVE_DUPLICATES true`, once per sample |
| Realign indels | GATK 3.8-1 | `RealignerTargetCreator`, then `IndelRealigner --consensusDeterminationModel USE_READS`, single-threaded |

Pair it with `variant_calling.tool: repadapt` to get calls equivalent to RepAdapt's:

```yaml
mapping:
  pipeline: repadapt
  repadapt:
    indel_realignment: auto # auto | true | false
variant_calling:
  tool: repadapt
```

On single-row samples, this reproduces RepAdapt's VCF records, apart from sample order and header lines.
At sites with extreme depth, such as collapsed repeats, calls also depend on the order of the BAMs; see [the RepAdapt calling model](variant-calling.md#repadapt-calling-model).
RepAdapt itself isn't exactly reproducible on large inputs, because fastp's multi-threaded output order changes between runs, and this pipeline shares that property.

**What snpArcher does differently:**

- **One library tag per sample.** Every row's reads get `LB:{sample}_LB`, as RepAdapt sets it. Picard groups duplicates by library, so duplicates are removed per sample, across all of a sample's rows and libraries.
- **`mark_duplicates: false` skips duplicate removal.** RepAdapt always removes duplicates.
- **Several rows per sample are mapped separately, then deduplicated together.** RepAdapt takes one pair of FASTQs per sample.
- **Realignment is skipped in long-contig mode.** With `indel_realignment: auto`, a warning says so; `true` is an error. GATK3 needs BAI indexes, which can't index contigs longer than 2^29 bp. `false` always skips realignment.
- **Depth tables aren't produced.** RepAdapt's per-gene and per-window depth tables need a GFF. snpArcher's [callable sites](architecture.md#callable-sites-a-parallel-track) cover coverage instead.

### Mapping QC

The QC report has the same columns for every pipeline, but in the `repadapt` pipeline some are measured differently:

| Column | Measured on |
|---|---|
| `total_reads`, `percent_mapped`, `percent_properly_paired` | bwa's output before the MAPQ filter, summed over rows. Properly paired is relative to reads paired in sequencing. |
| `num_duplicates`, `percent_duplicates` | Picard's duplication metrics. The final BAM has no duplicates left to count. |
| `mean_depth`, `covered_bases` | The final BAM, so they reflect the reads that are actually called. |

Because the MAPQ filter removes reads before duplicate removal, Picard's duplicate rate is computed among MAPQ >= 10 reads.

### Limitations

- **Linux only.** The pinned tool builds exist only for linux-64. Dry runs work on macOS.
- **GATK3 realignment is slow.** It is single-threaded per sample, and RepAdapt allows it 48 hours.
- **Fragmented references.** RepAdapt's README notes that GATK3 struggles with very fragmented references.

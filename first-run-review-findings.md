# snpArcher first-run review findings

Reviewed on November 5, 2026 against `main` at `5002bb4`, with `feat/repadapt` at `227b8e5` for comparison. All eight bugs below can affect a run started in a fresh working directory. The numbers match the original review. Findings 1 and 8 are grouped separately as BAM input validation.

## Confirmed bugs

### 2. The QC dashboard counts PASS variants as filtered

The dashboard counts only the VCF filter value `.` as retained and sums every other category as filtered. On a three-record fixture where every record passed, vcftools reported `PASS 3`. The dashboard therefore counts all three as filtered and has no retained count for its diversity estimate. This also affects RepAdapt, which explicitly labels passing calls `PASS`. [Dashboard calculation](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/modules/qc/scripts/qc_dashboard_interactive.Rmd#L155-L189)

Sum the counts for both `.` and `PASS`, consistent with downstream filtering. Test a PASS-only summary and a mixed summary containing both retained categories and a named filter.

### 3. Whole-genome GATK joint calling fails when intervals are disabled

With `intervals.enabled: false`, the workflow passes the reference FASTA index to `GenomicsDBImport -L`. GATK 4.6.2.0 rejected the fixture’s `.fai` as an unsupported interval format. The dry run succeeds because it only constructs the command. [Joint import rule](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/joint_gvcf.smk#L59-L76)

Generate a whole-reference BED or interval list from the FAI. Exercise the joint import on a small reference with intervals disabled.

### 4. GenMap can overwrite and delete a reference FASTA

GenMap decompresses its input to `results/reference/{name}.fa` and deletes that file on exit. On `main`, a successful GenMap run deleted a user-supplied FASTA at that path. On the RepAdapt branch, the realignment rule declares that exact path as its reference. With mappability and realignment enabled together, GenMap can truncate or remove the FASTA while GATK3 needs it. [GenMap scratch file](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/callable_sites.smk#L168-L206) · [RepAdapt realignment input](https://github.com/harvardinformatics/snpArcher/blob/227b8e5/workflow/rules/mapping/repadapt.smk#L164-L190)

Decompress into a job-private temporary directory. Verify that an existing FASTA at the shared path survives byte-for-byte, and exercise mappability and realignment together with multiple cores before integrating RepAdapt.

### 5. An unsupported postprocess setting overrides the selected VCF

Setting `modules.postprocess.filtering.vcf` to a different VCF produces a warning that the key is being ignored. The generated postprocess command nevertheless reads that file. The parent workflow sets `vcf` to `FINAL_VCF`, then expands the entire user-supplied filtering dictionary over it. A registered pipeline that fixes the call set would have the same bypass. [Postprocess configuration](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/Snakefile#L62-L75)

Forward only supported filtering keys to the module. Test that supplying an unsupported `vcf` key either fails validation or leaves `FINAL_VCF` as the actual rule input.

### 6. Ploidy validation allows values the selected caller cannot honor

The schema accepts any positive integer. With `ploidy: 4`, plain bcftools reaches the calling step and bcftools 1.23.1 rejects `--ploidy 4`. DeepVariant has a different failure mode: startup accepts a non-diploid setting, but neither DeepVariant command uses it. The RepAdapt caller already checks the numerical ploidies supported by its pinned bcftools version. [Config schema](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/schemas/config.schema.yaml) · [bcftools command](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/bcftools.smk#L99-L113) · [DeepVariant commands](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/deepvariant.smk) · [RepAdapt validation](https://github.com/harvardinformatics/snpArcher/blob/227b8e5/workflow/rules/common.smk#L595-L603)

Validate ploidy against the selected caller at startup. Apply the same check when a registered pipeline supplies the value. Test both an unsupported value and the command generated for each supported value.

### 7. Callable-site filtering can require one sample too many

The minimum count is calculated as `ceil(fraction * total_samples)`. With 100 samples and a configured fraction of `0.55`, binary floating-point arithmetic produces `55.00000000000001` and the ceiling becomes 56. A site callable in exactly 55 samples is incorrectly dropped. [Threshold calculation](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/scripts/callable_zarr_to_bed.py#L263-L277)

Calculate the ceiling from the decimal fraction. A boundary test should retain a site with 55 callable samples and reject one with 54.

### 9. The bcftools final VCF loses its index

On a fresh `call_variants` run, bcftools produced both `raw.vcf.gz` and `raw.vcf.gz.tbi`. Snakemake then removed the index because the gather rule marks it temporary and the target requests only the VCF. A region query requires the user to rebuild the index. [Gather outputs](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/bcftools.smk#L118-L137) · [Documented final outputs](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/website/docs/reference/outputs.md#L8-L15)

Keep the final VCF and its matching TBI or CSI index together. Verify this with an indexed region query after normal workflow cleanup.

### 10. Sentieon rules reference the wrong environment filename on main

The rules request `workflow/envs/sentieon.yaml`, while `main` ships `workflow/envs/sentieon.yml`. This breaks Conda-backed execution before the licensed tool runs. The RepAdapt branch already carries the rename. [Rule declaration](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/sentieon.smk#L24-L34) · [Environment on main](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/envs/sentieon.yml) · [Environment on RepAdapt](https://github.com/harvardinformatics/snpArcher/blob/227b8e5/workflow/envs/sentieon.yaml)

Bring the filename fix to `main`. Checking that each Sentieon rule resolves its environment file requires no license.

## BAM input validation (1 and 8)

Both cases involve inconsistent or unsupported external BAM inputs. The proposed checks would identify the input problem early enough for the user to correct it.

**1. Conflicting sample identities.** The sheet supplies a `sample_id`. In a BAM, the read-group `SM` field identifies the sample, and several read groups can legitimately share that value within one sample. The collision occurs when two **distinct sheet samples** have BAMs carrying the same `SM`. In the reproduced case, bcftools pooled both BAMs into one VCF sample column. [BAM staging](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/mapping.smk#L78-L104) · [bcftools calling inputs](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/variant_calling/bcftools.smk#L4-L10)

The proposed default is to stop before calling and report both sheet IDs, both BAM paths, and the shared `SM`. The user can then resolve the identity conflict. Whether to permit a one-to-one alias between a sheet ID and a BAM sample name is a separate, open policy question.

**8. BAMs with no reads flagged as paired.** snpArcher currently supports paired-end data. On an external BAM with 24 reads and zero flagged as paired, `samtools flagstat` returned `N/A` for percent properly paired. The QC parser checked total reads and attempted `float("N/A")`. [BAM QC parser](https://github.com/harvardinformatics/snpArcher/blob/5002bb4/workflow/rules/qc_metrics.smk#L30-L39)

For a nonempty BAM with zero reads flagged as paired, report the sample, BAM path, observed counts, and paired-end requirement. Base the check on the **paired-in-sequencing count**. A BAM with paired reads and **0% properly paired** has a valid QC measurement and should continue through reporting.

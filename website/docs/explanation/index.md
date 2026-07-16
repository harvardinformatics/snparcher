# Concepts

These pages cover the design decisions, trade-offs, and biology behind snpArcher.
None of this is required to use the pipeline, but it will help you make better choices for your project.

- [Pipeline architecture](architecture.md): How the pipeline is structured and why.
- [Variant calling](variant-calling.md): GATK vs. bcftools vs. DeepVariant, and when to choose each.
- [Parallelization](parallelization.md): The scatter-by-Ns strategy.
- [QC metrics](qc-metrics.md): What each QC figure means biologically.
- [Filtering philosophy](filtering.md): SFS-guided filtering with rationale and worked example.
- [Non-model organisms](non-model.md): Challenges specific to non-model organisms and how snpArcher addresses them.

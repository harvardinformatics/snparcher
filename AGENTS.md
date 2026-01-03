# AGENTS.md - snpArcher Development Guidelines

## Overview

snpArcher v2 is a Snakemake 9+ workflow for population genomic variant calling.
This document defines coding standards and practices for contributors.

**Key Resources:**

- [Snakemake 9 Documentation](https://snakemake.readthedocs.io/en/stable/)
- [Snakemake Best Practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
- [Snakemake Wrappers](https://snakemake-wrappers.readthedocs.io/)

## Development Environment

### Setup

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup
git clone https://github.com/harvardinformatics/snpArcher.git
cd snpArcher
pixi install -e dev  # Install dev environment
```

### Available Environments

| Environment | Purpose |
|-------------|---------|
| `default` | Runtime environment with Snakemake and dependencies |
| `dev` | Development tools (ruff, snakefmt) |
| `docs` | Documentation tools (mkdocs, mkdocs-material) |

### Available Tasks

```bash
# Default environment tasks
pixi run lint-snakemake    # Run Snakemake linter
pixi run test-dry-run      # Validate workflow with dry-run

# Dev environment tasks
pixi run -e dev lint       # Check code with ruff
pixi run -e dev format     # Format Python and Snakemake files
pixi run -e dev check      # Run all linters

# Docs environment tasks
pixi run -e docs docs-serve  # Serve documentation locally
pixi run -e docs docs-build  # Build documentation
```

## Python Code Standards

### Version & Tools

- Python 3.11+ required
- Use `ruff` for linting and formatting
- Configuration in `pyproject.toml`

### Type Hints

Use modern syntax (PEP 604, PEP 585):

```python
# Good
def process(items: list[str], config: dict[str, int] | None = None) -> bool: ...

# Avoid
from typing import List, Dict, Optional
def process(items: List[str], config: Optional[Dict[str, int]] = None) -> bool: ...
```

### Docstrings

Use Google style:

```python
def align_reads(sample_id: str, reference: Path) -> Path:
    """Align reads to reference genome.

    Args:
        sample_id: Unique sample identifier.
        reference: Path to reference FASTA.

    Returns:
        Path to output BAM file.

    Raises:
        FileNotFoundError: If reference doesn't exist.
    """
```

### No Custom Utility Modules

Do not create custom utility modules. Use:

- Snakemake built-in functions and helpers
- Standard library
- Well-maintained third-party packages

## Snakemake Standards

### Version

- Target Snakemake 9.x
- Use `min_version("9.0")` in Snakefile

### Formatting

- Use `snakefmt` for all `.smk` and `Snakefile` files
- Run `snakemake --lint` before committing

### Rule Naming Convention

Use `{stage}_{action}` format:

```
align_bwa_mem2
call_gatk_haplotypecaller
qc_fastqc
filter_bcftools
```

### Prefer Snakemake Wrappers

Use [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/) over custom
shell commands when available:

```python
rule align_bwa:
    input:
        reads=["reads/{sample}_1.fq.gz", "reads/{sample}_2.fq.gz"],
        idx=multiext("genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "aligned/{sample}.bam",
    wrapper:
        "v5.0.0/bio/bwa/mem"
```

### Semantic Helpers

Use built-in helpers instead of custom functions. See
[Snakemake semantic helpers documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-semantic-helpers).

```python
# Use collect() for aggregation
rule aggregate:
    input:
        collect("results/{sample}.txt", sample=SAMPLES)

# Use lookup() for sample sheet queries
rule process:
    input:
        lookup(query="sample_id == '{sample}'", within=samples, cols="fastq_path")

# Use branch() for conditional inputs
rule conditional:
    input:
        branch(
            lookup(dpath="use_filtered", within=config),
            then="filtered/{sample}.bam",
            otherwise="raw/{sample}.bam",
        )
```

### Config and Sample Sheet Validation

Use JSON Schema for validation:

```python
configfile: "config/config.yaml"
validate(config, "schemas/config.schema.json")

samples = pd.read_csv(config["samples"])
validate(samples, "schemas/samples.schema.json")
```

### Wildcard Constraints

Always define constraints:

```python
wildcard_constraints:
    sample="[A-Za-z0-9_-]+",
    refGenome="[A-Za-z0-9_.-]+",
```

### Logging

Always capture stdout/stderr:

```python
rule example:
    output:
        "out.txt",
    log:
        "logs/example.log",
    shell:
        "command {input} > {output} 2> {log}"
```

### Conda Environments

Every rule should specify a conda environment:

```python
rule align:
    input:
        ...
    output:
        ...
    conda:
        "../envs/align.yml"
    shell:
        ...
```

## Project Structure

```
snpArcher/
├── pyproject.toml          # Project config, pixi, ruff
├── AGENTS.md               # This file
├── mkdocs.yml              # Documentation config
├── docs/
│   ├── index.md
│   ├── tutorials/          # Learning-oriented
│   ├── how-to/             # Task-oriented
│   ├── reference/          # Information-oriented
│   └── explanation/        # Understanding-oriented
├── workflow/
│   ├── Snakefile           # Main entry, target rules
│   ├── schemas/
│   │   ├── config.schema.json
│   │   └── samples.schema.json
│   ├── rules/
│   │   ├── common.smk      # Shared config, helpers
│   │   ├── align/          # Alignment stage
│   │   │   └── bwa.smk
│   │   ├── call/           # Variant calling stage
│   │   │   ├── gatk.smk
│   │   │   └── sentieon.smk
│   │   ├── qc/             # Quality control
│   │   └── postprocess/    # Post-processing
│   ├── envs/               # Conda environment yamls
│   └── scripts/            # Python/R scripts
├── config/
│   └── config.yaml
└── tests/
    └── integration/
```

## Documentation Standards

### Framework

Follow the [Diataxis](https://diataxis.fr/) documentation framework:

| Type | Purpose | User Need |
|------|---------|-----------|
| **Tutorials** | Learning-oriented | "Help me get started" |
| **How-to Guides** | Task-oriented | "Help me solve a problem" |
| **Reference** | Information-oriented | "Help me find information" |
| **Explanation** | Understanding-oriented | "Help me understand" |

### Tools

- Use [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
- Docs source in `docs/`
- Serve locally: `pixi run -e docs docs-serve`

### Document as You Refactor

When modifying code:

1. Update relevant documentation immediately
2. Add docstrings to new functions
3. Include code references using `file_path:line_number` format

### Code References

When referencing code in documentation, use `file_path:line_number` format:

```markdown
The alignment is performed in `workflow/rules/align/bwa.smk:15`.
```

## Testing

### Integration Tests

- Test data lives in `.test/`
- Run dry-run validation: `pixi run test-dry-run`

### CI Checks

The following must pass before merging:

- `snakemake --lint`
- `snakemake -n` (dry-run)
- `ruff check`

## Commit Guidelines

Prefer conventional commits:

```
feat: add support for BWA-MEM2 aligner
fix: correct interval handling for small genomes
docs: add tutorial for custom reference genomes
refactor: reorganize rules by stage
```

Reference issues when applicable: `fix: handle empty VCF files (#123)`

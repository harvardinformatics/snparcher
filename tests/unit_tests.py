import tempfile
from pathlib import Path

import pytest

from conftest import SnakemakeRunner, FIXTURES_DIR

SAMPLES = FIXTURES_DIR / "config" / "samples.csv"


@pytest.mark.unit
def test_setup(request):
    """Setup from scratch: reference prep + indexing + intervals."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data")

        result = smk.run(target="setup", samples=SAMPLES)
        result.assert_success()
        result.assert_output_exists(
            "results/reference/my_organism.fa.gz",
            "results/reference/my_organism.fa.gz.fai",
            "results/reference/my_organism.dict",
            "results/intervals/gvcf/intervals.txt",
            "results/intervals/db/intervals.txt",
        )


@pytest.mark.unit
def test_fastp(request):
    """Just fastp for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data", "results/reference", "results/intervals")

        result = smk.run(
            target=[
                "results/filtered_fastqs/sample0/sample0/u1_1.fastq.gz",
                "results/filtered_fastqs/sample0/sample0/u1_2.fastq.gz",
            ],
            samples=SAMPLES,
        )
        result.print_log()
        result.assert_success()
        result.assert_output_exists(
            "results/filtered_fastqs/sample0/sample0/u1_1.fastq.gz",
        )


@pytest.mark.unit
def test_bwa_mem(request):
    """Just bwa mem for one sample — fastp output provided."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/filtered_fastqs/sample0",
        )

        result = smk.run(
            target=[
                "results/bams/raw/sample0/sample0/u1.bam",
                "results/bams/raw/sample0/sample0/u1.bam.bai",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/raw/sample0/sample0/u1.bam",
        )


@pytest.mark.unit
def test_markdup(request):
    """Just merge + markdup for one sample — raw BAMs provided."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/raw/sample0",
        )

        result = smk.run(
            target=[
                "results/bams/markdup/sample0.bam",
                "results/bams/markdup/sample0.bam.bai",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/markdup/sample0.bam",
        )


@pytest.mark.unit
def test_gvcfs(request):
    """Per-sample gVCF generation for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target=[
                "results/gvcfs/sample0.g.vcf.gz",
                "results/gvcfs/sample0.g.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/gvcfs/sample0.g.vcf.gz",
            "results/gvcfs/sample0.g.vcf.gz.tbi",
        )


@pytest.mark.unit
def test_qc_metrics(request):
    """QC metrics for one sample + combine."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
            "results/fastp",
        )

        result = smk.run(
            target="results/qc_metrics/bam/sample0.json",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc_metrics/bam/sample0.json",
        )


@pytest.mark.unit
def test_callable_sites_mappability(request):
    """Mappability only — no per-sample work needed."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
        )

        result = smk.run(
            target="results/callable_sites/mappability.bed",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/mappability.bed",
        )


@pytest.mark.unit
def test_callable_sites_coverage(request):
    """Coverage for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target="results/callable_sites/depths/sample0.per-base.d4",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/depths/sample0.per-base.d4",
        )


@pytest.mark.unit
def test_haplotypecaller(request):
    """HaplotypeCaller for one sample, one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target=[
                "results/interval_gvcfs/sample0/0000.g.vcf.gz",
                "results/interval_gvcfs/sample0/0000.g.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/interval_gvcfs/sample0/0000.g.vcf.gz",
        )


@pytest.mark.unit
def test_genomics_db_import(request):
    """GenomicsDB import for one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
        )

        result = smk.run(
            target="results/gatk_genomics_db/L0000.tar",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/gatk_genomics_db/L0000.tar",
        )


@pytest.mark.unit
def test_genotype_gvcfs(request):
    """GenotypeGVCFs for one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
            "results/gatk_genomics_db",
        )

        result = smk.run(
            target=[
                "results/vcfs/intervals/L0000.vcf.gz",
                "results/vcfs/intervals/L0000.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/intervals/L0000.vcf.gz",
        )


@pytest.mark.unit
def test_collect_fastp_stats(request):
    """Collect fastp stats for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/fastp",
        )

        result = smk.run(
            target="results/qc_metrics/fastp/sample0.json",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc_metrics/fastp/sample0.json",
        )


@pytest.mark.unit
def test_mosdepth(request):
    """Mosdepth for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target="results/callable_sites/depths/sample0.mosdepth.summary.txt",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/depths/sample0.mosdepth.summary.txt",
        )
@pytest.mark.unit
def test_coverage_bed(request):
    """Coverage BED generation from callable loci zarr."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/callable_sites/callable_loci.zarr",
        )

        result = smk.run(
            target="results/callable_sites/coverage.bed",
            samples=SAMPLES,
        )
        result.print_log()
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/coverage.bed",
        )

        # Verify BED file is non-empty and well-formed
        bed_path = Path(tmpdir) / "results" / "callable_sites" / "coverage.bed"
        lines = bed_path.read_text().strip().splitlines()
        assert len(lines) > 0, "coverage.bed should not be empty"
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 3, f"BED line should have 3 fields: {line}"
            assert fields[0] == "chr2l", f"Expected contig chr2l: {line}"
            assert int(fields[1]) < int(fields[2]), f"Start should be < end: {line}"

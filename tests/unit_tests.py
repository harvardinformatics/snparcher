import re
import subprocess
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
                "results/bams/markdup/sample0.bam.csi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/markdup/sample0.bam",
            "results/bams/markdup/sample0.bam.csi",
        )


@pytest.mark.unit
def test_external_bam_staging_and_indexing(request):
    """External BAM inputs are staged and indexed with workflow-managed CSI."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        samples = tmpdir / "external_bam_samples.csv"
        external_bam = FIXTURES_DIR / "results/bams/markdup/sample0.bam"
        samples.write_text(
            "sample_id,input_type,input\n"
            f"sample_ext,bam,{external_bam}\n"
        )

        smk = SnakemakeRunner(tmpdir, use_conda=not no_conda, conda_prefix=conda_prefix)

        result = smk.run(
            target="results/bams/input/sample_ext.bam.csi",
            samples=samples,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/input/sample_ext.bam",
            "results/bams/input/sample_ext.bam.csi",
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


@pytest.mark.unit
def test_callable_sites_bed_is_intersection(request):
    """The final callable BED keeps only sites passing coverage AND mappability.

    The sources are written directly, and --allowed-rules stops Snakemake from
    rebuilding them, so only callable_sites_bed runs. A union of the two
    sources would yield chr2l:0-300 rather than the two overlapping pieces.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data", "results/reference")

        callable_dir = Path(tmpdir) / "results" / "callable_sites"
        callable_dir.mkdir(parents=True)
        # Coverage BED is written in reference order, not sorted.
        (callable_dir / "coverage.bed").write_text("chr2l\t200\t300\nchr2l\t0\t100\n")
        (callable_dir / "mappability.bed").write_text("chr2l\t50\t250\n")

        result = smk.run(
            target="results/callable_sites/callable_sites.bed",
            samples=SAMPLES,
            extra_args=["--allowed-rules", "callable_sites_bed"],
        )
        result.print_log()
        result.assert_success()

        bed = (callable_dir / "callable_sites.bed").read_text()
        assert bed == "chr2l\t50\t100\nchr2l\t200\t250\n", bed


@pytest.mark.unit
def test_download_sra_completeness(request):
    """download_sra emits every read the archive holds.

    Exercises the completeness assertion in the rule: the extractor's output is
    counted and compared against the archive's own spot count, rather than the
    rule merely checking that the files exist. A truncated extractor, or a pigz
    that compresses nothing and exits 0, fails here.

    Requires network access: prefetch pulls SRR000001 (~312 MB) from NCBI.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    samples = Path(__file__).parent / "sample_sheets" / "local_fastqs_and_srr_same_library.csv"

    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data")

        r1 = "results/fastqs/sample_mixed_reads/libA/u2/SRR000001_1.fastq.gz"
        r2 = "results/fastqs/sample_mixed_reads/libA/u2/SRR000001_2.fastq.gz"

        result = smk.run(target=[r1, r2], samples=samples)
        result.print_log()
        result.assert_success()
        result.assert_output_exists(r1, r2)

        # The rule logs what it counted; assert it actually ran the check
        # rather than passing because the block was skipped.
        log = (
            Path(tmpdir)
            / "logs"
            / "download_sra"
            / "sample_mixed_reads"
            / "libA"
            / "u2"
            / "SRR000001.txt"
        ).read_text()
        assert "completeness:" in log, f"completeness check did not run:\n{log}"

        emitted, expected = re.search(
            r"completeness: emitted=(\d+) expected=(\d+)", log
        ).groups()
        assert emitted == expected, f"emitted {emitted} != expected {expected}"


@pytest.mark.unit
def test_csra_expected_reads():
    """The cSRA census predicts what sam-dump | collate | fastq emits.

    Each spot is one vdb-dump tab line: READ_LEN, PRIMARY_ALIGNMENT_ID and
    ALIGNMENT_COUNT, one comma-separated value per read. The layouts are the
    ones found in ficAlb's ERR2167036 and ERR2167049, where the prediction
    matched the extracted counts exactly.
    """
    script = Path(__file__).parents[1] / "workflow" / "scripts" / "csra_expected_reads.awk"
    spots = [
        "151, 151\t12, 9\t1, 3",  # both primary: a pair
        "151, 151\t5, 0\t1, 0",  # half-aligned: the unaligned mate is emitted
        "151, 151\t0, 0\t0, 0",  # fully unaligned: a pair
        "151, 151\t0, 0\t0, 2",  # secondary-only beside unaligned: still visited
        "151, 151\t7, 0\t1, 2",  # secondary-only mate of a primary: a singleton
        "151, 151\t0, 0\t2, 1",  # both secondary-only: nothing emitted
        "151, 0\t8, 0\t1, 0",  # mate-less BAM record: unpaired
        "0, 151\t0, 0\t0, 0",  # mate-less and unaligned: unpaired
        "151\t3\t1",  # one-read spot: not modelled
    ]
    result = subprocess.run(
        ["awk", "-f", str(script)],
        input="\n".join(spots) + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    census = dict(line.split("=") for line in result.stdout.split())
    assert census == {
        "spots": "9",
        "spots_unmodelled": "1",
        "lines_malformed": "0",
        "reads_primary": "5",
        "reads_zero_length": "2",
        "reads_lost_secondary_only": "3",
        "pairs": "4",
        "singletons": "1",
        "unpaired": "2",
    }

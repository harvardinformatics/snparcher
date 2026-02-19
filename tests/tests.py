import platform
from pathlib import Path
import tempfile

import pytest

from conftest import SnakemakeRunner

TEST_DIR = Path(__file__).parent
CONFIGS_DIR = TEST_DIR / "configs"
SAMPLES_DIR = TEST_DIR / "sample_sheets"


def get_samples_file():
    """Get appropriate sample sheet based on platform."""
    if platform.machine() == "arm64":
        return SAMPLES_DIR / "local_fastqs_no_dedup.csv"
    return SAMPLES_DIR / "local_fastqs.csv"

def get_config_file():
    """Get appropriate sample sheet based on platform."""
    if platform.machine() == "arm64":
        return CONFIGS_DIR / "local_genome_no_clam.yaml"
    return CONFIGS_DIR / "local_genome.yaml"

@pytest.mark.dry_run
def test_setup_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        expected_rules = [
            "prepare_reference",
            "index_reference",
            "picard_intervals",
            "create_gvcf_intervals",
            "create_db_intervals",
        ]
        output = result.stdout + result.stderr
        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found"


@pytest.mark.dry_run
def test_full_pipeline_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.dry_run
def test_multirow_same_library_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_multirow_same_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "InputFunctionException" not in output
        assert "merge_library_bams" in output
        assert "markdup_library" in output
        assert "merge_dedup_libraries" in output
        assert "input_unit=u1" in output
        assert "input_unit=u2" in output


@pytest.mark.dry_run
def test_multirow_multi_library_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_multirow_multi_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "InputFunctionException" not in output
        assert "merge_library_bams" in output
        assert "markdup_library" in output
        assert "merge_dedup_libraries" in output
        assert "library=libA" in output
        assert "library=libB" in output


@pytest.mark.full_run
def test_setup(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()
        result.assert_output_exists(
            "results/reference/test_genome.fa.gz",
            "results/reference/test_genome.fa.gz.fai",
            "results/reference/test_genome.dict",
            "results/intervals/gvcf/intervals.txt",
            "results/intervals/db/intervals.txt",
        )


@pytest.mark.full_run
def test_full_pipeline(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/raw.vcf.gz",
            "results/qc_metrics/qc_report.tsv",
            "results/callable_sites/mappability.bed",
        )


@pytest.mark.full_run
def test_incremental(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        # Run setup
        result = smk.run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        # Dry run all - setup rules should not appear
        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "prepare_reference" not in output
        assert "index_reference" not in output

        # Run all
        result = smk.run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

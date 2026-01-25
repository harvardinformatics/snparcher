from pathlib import Path
import tempfile

import pytest

TEST_DIR = Path(__file__).parent
CONFIGS_DIR = TEST_DIR / "configs"
SAMPLES_DIR = TEST_DIR / "sample_sheets"


@pytest.fixture(scope="module")
def setup_workdir():
    """Shared workdir for setup tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture(scope="module")
def setup_runner(setup_workdir):
    """Shared runner for setup tests."""
    from conftest import SnakemakeRunner

    return SnakemakeRunner(setup_workdir)


@pytest.mark.dependency()
def test_setup_dry_run(setup_runner):
    result = setup_runner.dry_run(
        target="setup",
        configfile=CONFIGS_DIR / "local_genome.yaml",
        samples=SAMPLES_DIR / "local_fastqs.csv",
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
        assert rule in output, f"Expected rule '{rule}' not found in dry run output"


@pytest.mark.dependency(depends=["test_setup_dry_run"])
def test_setup_run(setup_runner):
    result = setup_runner.run(
        target="setup",
        configfile=CONFIGS_DIR / "local_genome.yaml",
        samples=SAMPLES_DIR / "local_fastqs.csv",
    )

    result.assert_success()
    result.assert_output_exists(
        "results/reference/test_genome.fa.gz",
        "results/reference/test_genome.fa.gz.fai",
        "results/reference/test_genome.dict",
        "results/intervals/gvcf/intervals.txt",
        "results/intervals/db/intervals.txt",
    )

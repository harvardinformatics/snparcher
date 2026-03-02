import platform
import re
import shutil
from pathlib import Path
import tempfile

import pytest

from conftest import SnakemakeRunner, WORKFLOW_DIR, TEST_DATA_DIR

TEST_DIR = Path(__file__).parent
CONFIGS_DIR = TEST_DIR / "configs"
SAMPLES_DIR = TEST_DIR / "sample_sheets"
METADATA_DIR = TEST_DIR / "sample_metadata"


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


def get_multistage_config_file():
    """Get multistage concat config file based on platform."""
    if platform.machine() == "arm64":
        return CONFIGS_DIR / "local_genome_no_clam_multistage.yaml"
    return CONFIGS_DIR / "local_genome_multistage.yaml"

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


@pytest.mark.full_run
def test_multistage_interval_concat(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        configfile = get_multistage_config_file()
        samples = get_samples_file()

        result = smk.run(
            target="setup",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()

        result = smk.dry_run(
            target="results/vcfs/raw.vcf.gz",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "concat_interval_gvcfs_stage" in output
        assert "concat_interval_vcfs_stage" in output
        gvcf_stage_match = re.search(r"concat_interval_gvcfs_stage\s+(\d+)", output)
        assert gvcf_stage_match is not None
        assert int(gvcf_stage_match.group(1)) >= 1

        vcf_stage_match = re.search(r"concat_interval_vcfs_stage\s+(\d+)", output)
        assert vcf_stage_match is not None
        assert int(vcf_stage_match.group(1)) > 1

        result = smk.run(
            target="results/vcfs/raw.vcf.gz",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/raw.vcf.gz",
            "results/vcfs/raw.vcf.gz.tbi",
            "logs/concat_interval_gvcfs/staged/sample1/r1/c0.txt",
            "logs/concat_interval_gvcfs/staged/sample2/r1/c0.txt",
            "logs/concat_interval_vcfs/staged/r2/c0.txt",
        )


# --- Metadata tests ---

@pytest.mark.dry_run
def test_no_metadata_dry_run(request):
    """Pipeline works without sample_metadata configured."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.dry_run
def test_valid_metadata_dry_run(request):
    """Pipeline works with valid sample_metadata."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "valid_metadata.csv")},
        )
        result.assert_success()


@pytest.mark.dry_run
def test_minimal_metadata_dry_run(request):
    """Pipeline works with metadata that only has sample_id and exclude."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "minimal_metadata.csv")},
        )
        result.assert_success()


@pytest.mark.dry_run
def test_partial_metadata_warns(request):
    """Pipeline warns when some samples lack metadata rows."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "partial_metadata.csv")},
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "sample2" in output, "Expected warning about sample2 missing from metadata"


@pytest.mark.dry_run
def test_unknown_sample_in_metadata_fails(request):
    """Pipeline fails when metadata contains sample_ids not in the sample sheet."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "unknown_sample.csv")},
        )
        assert not result.succeeded, "Expected failure for unknown sample in metadata"
        output = result.stdout + result.stderr
        assert "unknown_sample" in output


# --- Postprocess module tests ---

def get_postprocess_config():
    """Config with postprocess module enabled."""
    return CONFIGS_DIR / "local_genome_postprocess.yaml"


@pytest.mark.dry_run
def test_postprocess_dry_run(request):
    """Dry run postprocess basic_filter target (no callable_sites.bed dependency)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        # Target filtered.vcf.gz directly — this exercises filter_individuals
        # and basic_filter without needing callable_sites.bed (which requires
        # the zarr→BED conversion rule not yet added).
        result = smk.dry_run(
            target="results/postprocess/filtered.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "postprocess_filter_individuals" in output, \
            "Expected postprocess_filter_individuals rule in DAG"
        assert "postprocess_basic_filter" in output, \
            "Expected postprocess_basic_filter rule in DAG"


@pytest.mark.dry_run
def test_postprocess_disabled_no_rules(request):
    """When postprocess is disabled, no postprocess rules appear."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "postprocess_" not in output, \
            "Postprocess rules should not appear when module is disabled"


@pytest.mark.dry_run
def test_postprocess_with_metadata_dry_run(request):
    """Postprocess works with sample metadata (exclude column)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/postprocess/filtered.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "exclude_and_outgroup.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "postprocess_filter_individuals" in output


# --- QC module tests ---

def get_qc_config():
    """Config with QC module enabled."""
    return CONFIGS_DIR / "local_genome_qc.yaml"


@pytest.mark.dry_run
def test_qc_dry_run(request):
    """Dry run QC dashboard target — exercises all QC rules."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_check_fai" in output, \
            "Expected qc_check_fai rule in DAG"
        assert "qc_vcftools_individuals" in output, \
            "Expected qc_vcftools_individuals rule in DAG"
        assert "qc_subsample_snps" in output, \
            "Expected qc_subsample_snps rule in DAG"
        assert "qc_plink" in output, \
            "Expected qc_plink rule in DAG"
        assert "qc_admixture" in output, \
            "Expected qc_admixture rule in DAG"
        assert "qc_qc_dashboard" in output, \
            "Expected qc_qc_dashboard rule in DAG"
        # copy_qc_report should appear since main workflow provides qc_report
        assert "qc_copy_qc_report" in output, \
            "Expected qc_copy_qc_report rule in DAG"


@pytest.mark.dry_run
def test_qc_disabled_no_rules(request):
    """When QC is disabled, no QC rules appear."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_check_fai" not in output, \
            "QC rules should not appear when module is disabled"
        assert "qc_plink" not in output, \
            "QC rules should not appear when module is disabled"


@pytest.mark.dry_run
def test_qc_no_coords_without_metadata(request):
    """QC works without metadata — generate_coords_file should not run."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        # No metadata → no coords → no generate_coords_file rule
        assert "qc_generate_coords_file" not in output, \
            "generate_coords_file should not run without metadata"


@pytest.mark.dry_run
def test_qc_with_coords_metadata(request):
    """QC generates coords when metadata has lat/long columns."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "valid_metadata.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_generate_coords_file" in output, \
            "generate_coords_file should run when metadata has lat/long"


@pytest.mark.full_run
def test_qc_standalone_full_run(request):
    """Full execution of QC module as standalone workflow against test fixtures."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
                "vcf": str(TEST_DATA_DIR / "qc" / "raw.vcf.gz"),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
                "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
            },
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc/qc_dashboard.html",
            "results/qc/plink.eigenvec",
            "results/qc/plink.3.Q",
            "results/qc/coords.txt",
            "results/qc/qc_report.tsv",
        )

        # Copy dashboard HTML out before tmpdir cleanup so CI can upload it
        artifacts_dir = Path("test-artifacts")
        artifacts_dir.mkdir(exist_ok=True)
        dashboard = Path(tmpdir) / "results" / "qc" / "qc_dashboard.html"
        if dashboard.exists():
            shutil.copy2(dashboard, artifacts_dir / "qc_dashboard.html")

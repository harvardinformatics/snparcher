import functools
import gzip
import http.server
import platform
import re
import shutil
import socket
import tempfile
import threading
from contextlib import contextmanager
from pathlib import Path

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


def write_config_for_tool(base_config, out_dir, tool, parabricks_image=None):
    """Write a config copy with variant_calling.tool overridden."""
    text = Path(base_config).read_text()
    text = text.replace('tool: "gatk"', f'tool: "{tool}"', 1)
    if parabricks_image is not None:
        text = text.replace(
            'container_image: "/tmp/parabricks.sif"',
            f'container_image: "{parabricks_image}"',
            1,
        )

    out_path = Path(out_dir) / f"config_{tool}.yaml"
    out_path.write_text(text)
    return out_path

def write_callable_sites_config(base_config, out_dir, *, coverage_enabled, mappability_enabled):
    """Write a config copy with callable-sites toggles overridden."""
    text = Path(base_config).read_text()
    new_block = f"""callable_sites:
  generate_bed_file: false
  coverage:
    enabled: {"true" if coverage_enabled else "false"}
    stdev: 2
    merge_distance: 100
  mappability:
    enabled: {"true" if mappability_enabled else "false"}
    min_score: 1
    kmer: 150
    merge_distance: 100
"""
    pattern = re.compile(
        r"callable_sites:\n"
        r"  generate_bed_file: false\n"
        r"  coverage:\n"
        r"    enabled: (?:true|false)\n"
        r"    stdev: 2\n"
        r"    merge_distance: 100\n"
        r"  mappability:\n"
        r"    enabled: (?:true|false)\n"
        r"    min_score: 1\n"
        r"    kmer: 150\n"
        r"    merge_distance: 100\n"
    )
    if not pattern.search(text):
        raise AssertionError("Expected callable_sites block not found in config")

    out_path = Path(out_dir) / (
        f"config_callable_sites_{int(coverage_enabled)}_{int(mappability_enabled)}.yaml"
    )
    out_path.write_text(pattern.sub(new_block, text, count=1))
    return out_path

def write_intervals_config(base_config, out_dir, *, enabled):
    """Write a config copy with intervals.enabled overridden."""
    text = Path(base_config).read_text()
    pattern = re.compile(r"intervals:\n  enabled: (?:true|false)\n")
    if not pattern.search(text):
        raise AssertionError("Expected intervals block not found in config")

    out_path = Path(out_dir) / f"config_intervals_{int(enabled)}.yaml"
    out_path.write_text(
        pattern.sub(
            f"intervals:\n  enabled: {'true' if enabled else 'false'}\n",
            text,
            count=1,
        )
    )
    return out_path


def write_reference_source_config(base_config, out_dir, *, source):
    """Write a config copy with reference.source overridden."""
    text = Path(base_config).read_text()
    pattern = re.compile(
        r'(^reference:\n  name: ".*"\n  source: ").*(".*$)',
        re.MULTILINE,
    )
    if not pattern.search(text):
        raise AssertionError("Expected reference block not found in config")

    out_path = Path(out_dir) / "config_reference_source.yaml"
    out_path.write_text(pattern.sub(rf'\1{source}\2', text, count=1))
    return out_path


def write_gvcf_sample_sheet(out_dir, *, sample_id, gvcf_path):
    """Write a sample sheet for one external gVCF input."""
    out_path = Path(out_dir) / "external_gvcf_samples.csv"
    out_path.write_text(
        "sample_id,input_type,input\n"
        f"{sample_id},gvcf,{gvcf_path}\n"
    )
    return out_path


@contextmanager
def serve_directory(directory):
    """Serve a directory over HTTP for reference URL tests."""
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        host, port = sock.getsockname()

    server = http.server.ThreadingHTTPServer((host, port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://{host}:{port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


@pytest.mark.dry_run
def test_v1_style_config_exits_with_migration_message(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=CONFIGS_DIR / "legacy_v1_style.yaml",
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "Detected a v1-style snpArcher config" in output
        assert "docs/v2-migration.md" in output
        assert "sentieon" in output


@pytest.mark.dry_run
def test_extra_config_keys_warn_but_do_not_fail(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=CONFIGS_DIR / "local_genome_extra_keys.yaml",
            samples=get_samples_file(),
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert "Ignoring unsupported config key(s):" in output
        assert "reads.sort" in output
        assert "remote_reads" in output
        assert "variant_calling.gatk.min_pruning" in output
        assert "reference.path" in output


@pytest.mark.dry_run
def test_global_mark_duplicates_config_default_is_respected(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=CONFIGS_DIR / "local_genome_global_markdup.yaml",
            samples=SAMPLES_DIR / "local_fastqs_global_markdup.csv",
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert "merge_library_level_bams" in output
        assert "markdup_library" not in output
        assert "merge_dedup_libraries" not in output

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
def test_gvcfs_target_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="gvcfs",
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


@pytest.mark.dry_run
@pytest.mark.parametrize(
    ("coverage_enabled", "mappability_enabled"),
    [
        (True, True),
        (False, True),
        (True, False),
    ],
)
def test_callable_sites_target_dry_run(request, coverage_enabled, mappability_enabled):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=coverage_enabled,
            mappability_enabled=mappability_enabled,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert ("clam_loci" in output) == coverage_enabled
        assert ("mappability_bed" in output) == mappability_enabled


@pytest.mark.dry_run
def test_gatk_without_intervals_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_intervals_config(get_config_file(), tmpdir, enabled=False)

        result = smk.dry_run(
            target="call_variants",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "gatk_haplotypecaller" in output
        assert "joint_genomics_db_import" in output


@pytest.mark.full_run
@pytest.mark.parametrize("compressed", [False, True])
def test_reference_url_sources(request, compressed):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        source_dir = tmp_path / "reference_server"
        source_dir.mkdir()

        gz_ref = TEST_DATA_DIR / "reference" / "test_genome.fa.gz"
        plain_ref = source_dir / "test_genome.fa"
        plain_ref.write_text(gzip.decompress(gz_ref.read_bytes()).decode())
        shutil.copy2(gz_ref, source_dir / "test_genome.fa.gz")

        filename = "test_genome.fa.gz" if compressed else "test_genome.fa"

        with serve_directory(source_dir) as base_url:
            smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
            cfg = write_reference_source_config(
                get_config_file(),
                tmpdir,
                source=f"{base_url}/{filename}",
            )

            result = smk.run(
                target="setup",
                configfile=cfg,
                samples=get_samples_file(),
            )
            result.assert_success()
            result.assert_output_exists(
                "results/reference/test_genome.fa.gz",
                "results/reference/test_genome.fa.gz.fai",
                "results/reference/test_genome.dict",
            )


@pytest.mark.full_run
@pytest.mark.parametrize("intervals_enabled", [False, True])
def test_create_db_mapfile_preserves_external_gvcf_sample_id(request, intervals_enabled):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gvcf_path = tmp_path / "external_name.g.vcf.gz"
        gvcf_path.write_text("")

        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=False,
            mappability_enabled=True,
        )
        cfg = write_intervals_config(cfg, tmpdir, enabled=intervals_enabled)
        samples = write_gvcf_sample_sheet(
            tmpdir,
            sample_id="sample_gvcf",
            gvcf_path=gvcf_path,
        )

        result = smk.run(
            target="create_db_mapfile",
            configfile=cfg,
            samples=samples,
        )
        result.assert_success()

        mapfile = (tmp_path / "results/genomics_db/mapfile.txt").read_text().strip()
        assert mapfile == f"sample_gvcf\t{gvcf_path}"


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


@pytest.mark.dry_run
def test_bcftools_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "bcftools")

        result = smk.dry_run(
            target="call_variants",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "bcftools_regions" in output
        assert "bcftools_concat_regions" in output
        assert "variant_filtration" in output


@pytest.mark.dry_run
def test_deepvariant_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "deepvariant")

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "deepvariant_call" in output
        assert "glnexus_joint" in output


@pytest.mark.dry_run
def test_parabricks_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "parabricks")

        result = smk.dry_run(
            target=["call_variants", "results/gvcfs/sample1.g.vcf.gz"],
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "parabricks_haplotypecaller" in output
        assert "concat_interval_vcfs" in output


@pytest.mark.dry_run
@pytest.mark.parametrize("tool", ["bcftools", "deepvariant", "parabricks"])
def test_gvcf_input_rejected_for_new_callers(request, tool):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, tool)

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=SAMPLES_DIR / "local_gvcf.csv",
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "does not support samples with input_type='gvcf'" in output


@pytest.mark.dry_run
def test_callable_sites_coverage_rejected_for_gvcf_only_inputs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=True,
            mappability_enabled=True,
        )

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=SAMPLES_DIR / "local_gvcf.csv",
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "callable_sites.coverage.enabled requires at least one BAM-backed sample" in output


@pytest.mark.dry_run
def test_parabricks_requires_container_image(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(
            get_config_file(), tmpdir, "parabricks", parabricks_image=""
        )

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "parabricks.container_image is required" in output


@pytest.mark.full_run
def test_full_pipeline_bcftools(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "bcftools")

        result = smk.run(
            target=["all", "call_variants"],
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/raw.vcf.gz",
            "results/vcfs/filtered.vcf.gz",
            "results/qc_metrics/qc_report.tsv",
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


# --- MK module tests ---

def get_mk_config():
    """Config with MK module enabled."""
    return CONFIGS_DIR / "local_genome_mk.yaml"


@pytest.mark.dry_run
def test_mk_dry_run(request):
    """Dry run MK module — exercises all MK rules."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/mk/mk_table.tsv",
            configfile=get_mk_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "mk_decompress_ref" in output, \
            "Expected mk_decompress_ref rule in DAG"
        assert "mk_split_samples" in output, \
            "Expected mk_split_samples rule in DAG"
        assert "mk_degenotate" in output, \
            "Expected mk_degenotate rule in DAG"
        assert "mk_copy_gff" in output, \
            "Expected mk_copy_gff rule in DAG"


@pytest.mark.dry_run
def test_mk_disabled_no_rules(request):
    """When MK is disabled, no MK rules appear."""
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
        assert "mk_degenotate" not in output, \
            "MK rules should not appear when module is disabled"
        assert "mk_split_samples" not in output, \
            "MK rules should not appear when module is disabled"


@pytest.mark.dry_run
def test_mk_with_metadata_dry_run(request):
    """MK works with sample metadata (outgroup/exclude columns)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/mk/mk_table.tsv",
            configfile=get_mk_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "exclude_and_outgroup.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "mk_split_samples" in output

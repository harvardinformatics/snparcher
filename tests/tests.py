import functools
import gzip
import http.server
import platform
import re
import shutil
import socket
import subprocess
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

def write_callable_sites_config(
    base_config,
    out_dir,
    *,
    generate_bed_file=True,
    coverage_enabled,
    mappability_enabled,
    fraction="1.0",
    min_coverage="auto",
    max_coverage="auto",
):
    """Write a config copy with callable-sites toggles overridden."""
    text = Path(base_config).read_text()
    new_block = f"""callable_sites:
  generate_bed_file: {"true" if generate_bed_file else "false"}
  coverage:
    enabled: {"true" if coverage_enabled else "false"}
    fraction: {fraction}
    min_coverage: {min_coverage}
    max_coverage: {max_coverage}
    merge_distance: 100
  mappability:
    enabled: {"true" if mappability_enabled else "false"}
    min_score: 1
    kmer: 150
    merge_distance: 100
"""
    pattern = re.compile(
        r"callable_sites:\n"
        r"  generate_bed_file: (?:true|false)\n"
        r"  coverage:\n"
        r"    enabled: (?:true|false)\n"
        r"    fraction: .*\n"
        r"    min_coverage: .*\n"
        r"    max_coverage: .*\n"
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


def write_qc_config_without_exclude_scaffolds(base_config, out_dir):
    """Write a QC config copy without modules.qc.exclude_scaffolds."""
    text = Path(base_config).read_text()
    pattern = re.compile(r"\n    exclude_scaffolds: \".*\"\n")
    if not pattern.search(text):
        raise AssertionError("Expected modules.qc.exclude_scaffolds line not found in config")

    out_path = Path(out_dir) / "config_qc_no_exclude_scaffolds.yaml"
    out_path.write_text(pattern.sub("\n", text, count=1))
    return out_path


def iter_vcf_records(path):
    """Yield parsed VCF records from a plain-text or gzipped VCF."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            yield line.rstrip().split("\t")


def get_vcf_samples(path):
    """Return sample names from the VCF header."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]
    raise AssertionError(f"Missing #CHROM header in {path}")


def skip_if_arm64_packages_unavailable(result, *package_markers):
    """Skip full-run tests on osx-arm64 when conda cannot provide required packages."""
    if platform.machine() != "arm64":
        return

    stderr = result.stderr
    if "PackagesNotFoundError" not in stderr or "Platform: osx-arm64" not in stderr:
        return
    if package_markers and not any(marker in stderr for marker in package_markers):
        return

    missing = []
    collect_missing = False
    for line in stderr.splitlines():
        if line.strip().startswith("PackagesNotFoundError:"):
            collect_missing = True
            continue
        if not collect_missing:
            continue

        stripped = line.strip()
        if stripped.startswith("- "):
            missing.append(stripped[2:])
            continue
        if missing and not stripped:
            break

    details = ", ".join(missing) if missing else "required conda packages"
    pytest.skip(f"Conda package(s) unavailable on osx-arm64: {details}")


def write_numeric_qc_inputs(out_dir):
    """Write numeric-only QC fixtures derived from the existing standalone QC data."""
    source_vcf = TEST_DATA_DIR / "qc" / "raw.vcf.gz"
    source_fai = TEST_DATA_DIR / "qc" / "ref.fai"
    out_dir = Path(out_dir)
    out_vcf = out_dir / "numeric_raw.vcf.gz"
    out_fai = out_dir / "numeric_ref.fai"

    contigs = []
    with open(source_fai) as handle:
        for line in handle:
            fields = line.strip().split()
            if fields:
                contigs.append(fields[0])
    if len(contigs) < 2:
        raise AssertionError("Expected at least two contigs in QC FAI fixture")

    contig_map = {contig: str(i) for i, contig in enumerate(contigs, start=1)}

    with open(source_fai) as src, open(out_fai, "w") as dst:
        for line in src:
            fields = line.strip().split()
            if not fields:
                continue
            fields[0] = contig_map[fields[0]]
            dst.write("\t".join(fields) + "\n")

    sample_names = None
    contig_header_pattern = re.compile(r"^(##contig=<ID=)([^,>]+)(.*)$")
    with gzip.open(source_vcf, "rt") as src, gzip.open(out_vcf, "wt") as dst:
        for line in src:
            if line.startswith("##contig=<ID="):
                match = contig_header_pattern.match(line.rstrip("\n"))
                if match and match.group(2) in contig_map:
                    line = f"{match.group(1)}{contig_map[match.group(2)]}{match.group(3)}\n"
            elif line.startswith("#CHROM"):
                sample_names = line.rstrip().split("\t")[9:]
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if fields[0] in contig_map:
                    fields[0] = contig_map[fields[0]]
                line = "\t".join(fields) + "\n"
            dst.write(line)

        if sample_names is None:
            raise AssertionError("Missing #CHROM header in QC VCF fixture")

        extra_variant = [
            contig_map[contigs[1]],
            "1",
            ".",
            "A",
            "G",
            "60",
            "PASS",
            "AF=0.0104;ReadPosRankSum=0;FS=0;SOR=0;MQ=60;MQRankSum=0",
            "GT",
            *(["0/1", "0/1"] + ["0/0"] * (len(sample_names) - 2)),
        ]
        dst.write("\t".join(extra_variant) + "\n")

    return out_vcf, out_fai


def get_vcf_contig_headers(path):
    """Return contig IDs declared in the VCF header."""
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = []
    pattern = re.compile(r"^##contig=<ID=([^,>]+)")
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                break
            match = pattern.match(line.rstrip())
            if match:
                contigs.append(match.group(1))
    return contigs


def extract_r_function_source(path, function_name):
    """Extract an R function body from an Rmd file by balanced braces."""
    text = Path(path).read_text()
    marker = f"{function_name} <- function"
    start = text.find(marker)
    if start == -1:
        raise AssertionError(f"Function '{function_name}' not found in {path}")

    brace_start = text.find("{", start)
    if brace_start == -1:
        raise AssertionError(f"Could not find opening brace for '{function_name}' in {path}")

    depth = 0
    end = None
    for idx in range(brace_start, len(text)):
        char = text[idx]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                end = idx + 1
                break

    if end is None:
        raise AssertionError(f"Could not find closing brace for '{function_name}' in {path}")

    return text[start:end]


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


def write_fastq_sample_sheet(out_dir, *, sample_id, library_id):
    """Write a one-row FASTQ sample sheet with explicit sample and library IDs."""
    out_path = Path(out_dir) / "numeric_id_fastqs.csv"
    out_path.write_text(
        "sample_id,input_type,input,library_id,mark_duplicates\n"
        f"{sample_id},fastq,tests/data/fastq/sample1_1.fastq.gz;tests/data/fastq/sample1_2.fastq.gz,{library_id},true\n"
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
def test_mixed_srr_and_fastq_same_sample_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_and_srr_same_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "unsupported mixed input_type values" not in output
        assert "download_sra" in output
        assert "while IFS= read -r url; do" in output
        assert 'curl -fSL "$url"' in output
        assert "merge_library_bams" in output
        assert "library=libA" in output
        assert "input_unit=u1" in output
        assert "input_unit=u2" in output


@pytest.mark.dry_run
@pytest.mark.parametrize(
    ("generate_bed_file", "coverage_enabled", "mappability_enabled"),
    [
        (True, True, True),
        (True, False, True),
        (True, True, False),
        (False, True, True),
    ],
)
def test_callable_sites_target_dry_run(
    request,
    generate_bed_file,
    coverage_enabled,
    mappability_enabled,
):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=generate_bed_file,
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
        assert ("callable_coverage_thresholds" in output) == coverage_enabled
        assert ("clam_loci" in output) == coverage_enabled
        assert ("coverage_bed" in output) == (generate_bed_file and coverage_enabled)
        assert ("mappability_bed" in output) == mappability_enabled
        assert ("callable_sites_bed" in output) == (
            generate_bed_file and (coverage_enabled or mappability_enabled)
        )


@pytest.mark.dry_run
def test_callable_sites_disabled_sources_warn_and_skip_final_bed(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=False,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "Skipping results/callable_sites/callable_sites.bed generation" in output
        assert "callable_sites_bed" not in output
        assert "clam_loci" not in output
        assert "mappability_bed" not in output


@pytest.mark.dry_run
def test_callable_sites_numeric_thresholds_validate(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
            fraction="0.75",
            min_coverage="5",
            max_coverage="40",
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "coverage_bed" in output


@pytest.mark.dry_run
def test_callable_sites_invalid_fraction_fails_validation(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
            fraction="1.5",
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "callable_sites.coverage.fraction" in output


@pytest.mark.dry_run
def test_clam_loci_uses_per_sample_and_thresholds(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
            min_coverage="5",
            max_coverage="40",
        )

        result = smk.dry_run(
            target="results/callable_sites/callable_loci.zarr",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "--per-sample" in output
        assert '--min-coverage 5' in output or "min_coverage\t5" in output
        assert '--max-coverage 40' in output or "max_coverage\t40" in output


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
@pytest.mark.parametrize("sample_id", ["sample_gvcf", "00123"])
def test_create_db_mapfile_preserves_external_gvcf_sample_id(request, intervals_enabled, sample_id):
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
            sample_id=sample_id,
            gvcf_path=gvcf_path,
        )

        result = smk.run(
            target="create_db_mapfile",
            configfile=cfg,
            samples=samples,
        )
        result.assert_success()

        mapfile = (tmp_path / "results/genomics_db/mapfile.txt").read_text().strip()
        assert mapfile == f"{sample_id}\t{gvcf_path}"


@pytest.mark.dry_run
def test_fastq_dry_run_accepts_numeric_like_sample_and_library_ids(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        samples = write_fastq_sample_sheet(tmpdir, sample_id="00123", library_id="456")

        result = smk.dry_run(
            target=[
                "results/filtered_fastqs/00123/456/u1_1.fastq.gz",
                "results/filtered_fastqs/00123/456/u1_2.fastq.gz",
            ],
            configfile=get_config_file(),
            samples=samples,
        )

        result.assert_success()


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
            "results/callable_sites/callable_sites.bed",
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
    """Dry run postprocess target with callable-sites BED generation enabled."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/postprocess/clean_snps.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "callable_sites_bed" in output
        assert "postprocess_filter_individuals" in output, \
            "Expected postprocess_filter_individuals rule in DAG"
        assert "postprocess_basic_filter" in output, \
            "Expected postprocess_basic_filter rule in DAG"
        assert "postprocess_update_bed" in output, \
            "Expected postprocess_update_bed rule in DAG"


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
            target="results/postprocess/clean_snps.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "exclude_and_outgroup.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "callable_sites_bed" in output
        assert "postprocess_filter_individuals" in output


@pytest.mark.dry_run
def test_postprocess_warns_and_disables_without_callable_bed(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_postprocess_config(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=False,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "Disabling postprocess" in output
        assert "postprocess_" not in output


@pytest.mark.full_run
def test_postprocess_standalone_full_run(request):
    """Full execution of postprocess module against dedicated standalone fixtures."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        fixture_dir = TEST_DATA_DIR / "postprocess"
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                "samples": str(fixture_dir / "samples.csv"),
                "sample_metadata": str(fixture_dir / "sample_metadata.csv"),
                "vcf": str(fixture_dir / "raw.vcf"),
                "ref_fai": str(fixture_dir / "ref.fai"),
                "callable_sites_bed": str(fixture_dir / "callable_sites.bed"),
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        result.assert_success()
        result.assert_output_exists(
            "results/postprocess/filtered.vcf.gz",
            "results/postprocess/clean_snps.vcf.gz",
            "results/postprocess/clean_indels.vcf.gz",
        )

        filtered_vcf = Path(tmpdir) / "results" / "postprocess" / "filtered.vcf.gz"
        clean_snps_vcf = Path(tmpdir) / "results" / "postprocess" / "clean_snps.vcf.gz"
        clean_indels_vcf = Path(tmpdir) / "results" / "postprocess" / "clean_indels.vcf.gz"

        assert get_vcf_samples(filtered_vcf) == ["sample1"]
        assert get_vcf_samples(clean_snps_vcf) == ["sample1"]
        assert get_vcf_samples(clean_indels_vcf) == ["sample1"]

        snp_records = list(iter_vcf_records(clean_snps_vcf))
        indel_records = list(iter_vcf_records(clean_indels_vcf))
        assert snp_records, "Expected at least one SNP record in clean_snps.vcf.gz"
        assert indel_records, "Expected at least one indel record in clean_indels.vcf.gz"
        assert all(len(record[3]) == 1 and all(len(alt) == 1 for alt in record[4].split(",")) for record in snp_records)
        assert all(
            len(record[3]) != 1 or any(len(alt) != 1 for alt in record[4].split(","))
            for record in indel_records
        )


# --- QC module tests ---

def get_qc_config():
    """Config with QC module enabled."""
    return CONFIGS_DIR / "local_genome_qc.yaml"


def get_qc_underscore_ref_config():
    """Config with QC enabled and a local reference name containing an underscore."""
    return CONFIGS_DIR / "local_genome_ref_name_underscore.yaml"


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
        assert "qc_contig_map" in output, \
            "Expected qc_contig_map rule in DAG"
        assert "qc_vcftools_individuals" in output, \
            "Expected qc_vcftools_individuals rule in DAG"
        assert "qc_subsample_snps" in output, \
            "Expected qc_subsample_snps rule in DAG"
        assert "qc_prepare_plink_inputs" in output, \
            "Expected qc_prepare_plink_inputs rule in DAG"
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
def test_default_target_with_qc_enabled_runs_full_pipeline(request):
    """No explicit target should still resolve to the main workflow's all rule."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target=[],
            configfile=get_qc_underscore_ref_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "rule all:" in output, "Expected the main all rule to be the default target"
        assert "bwa_mem" in output, "Expected core mapping/variant DAG rules in the default target"
        assert "qc_qc_dashboard" in output, "Expected the QC dashboard DAG in the default target"


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
        assert "qc_contig_map" not in output, \
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


@pytest.mark.dry_run
def test_qc_defaults_exclude_scaffolds_when_omitted(request):
    """QC dry-run should succeed when modules.qc.exclude_scaffolds is omitted."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_qc_config_without_exclude_scaffolds(get_qc_config(), tmpdir)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()


def test_write_numeric_qc_inputs_rewrites_contigs():
    """The synthetic numeric QC fixture should rewrite both header and record contigs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        numeric_vcf, numeric_fai = write_numeric_qc_inputs(tmpdir)

        headers = get_vcf_contig_headers(numeric_vcf)
        assert "1" in headers
        assert "2" in headers
        assert "JAKDEW010000001.1" not in headers
        assert "JAKDEW010000002.1" not in headers

        records = list(iter_vcf_records(numeric_vcf))
        chroms = {record[0] for record in records}
        assert "1" in chroms
        assert "2" in chroms

        fai_lines = Path(numeric_fai).read_text().strip().splitlines()
        assert fai_lines[:2] == ["1\t1\t3", "2\t1\t3"]


def test_qc_dashboard_helper_preserves_numeric_like_ids():
    """QC dashboard helper should preserve leading zeros for headered and headerless tables."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not available")

    helper_source = extract_r_function_source(
        WORKFLOW_DIR / "modules" / "qc" / "scripts" / "qc_dashboard_interactive.Rmd",
        "read_table_preserve_ids",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        headerless = tmp_path / "dist.id"
        headerless.write_text("0 000123\n0 000456\n")
        headered = tmp_path / "depth.tsv"
        headered.write_text("INDV\tMEAN_DEPTH\n000123\t5.2\n")

        script = tmp_path / "validate_read_table_preserve_ids.R"
        script.write_text(
            "\n".join(
                [
                    "args <- commandArgs(trailingOnly = TRUE)",
                    "headerless <- args[[1]]",
                    "headered <- args[[2]]",
                    helper_source,
                    "headerless_df <- read_table_preserve_ids(headerless, id_cols = c('V1', 'V2'))",
                    "if (!is.character(headerless_df$V1) || !is.character(headerless_df$V2)) stop('Headerless ID columns were not read as character')",
                    "if (!identical(headerless_df$V2[[1]], '000123')) stop('Headerless leading zeros were not preserved')",
                    "headered_df <- read_table_preserve_ids(headered, header = TRUE, sep = '\\t', id_cols = c('INDV'))",
                    "if (!is.character(headered_df$INDV)) stop('Headered ID column was not read as character')",
                    "if (!identical(headered_df$INDV[[1]], '000123')) stop('Headered leading zeros were not preserved')",
                ]
            )
        )

        result = subprocess.run(
            ["Rscript", str(script), str(headerless), str(headered)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, (result.stdout + result.stderr).strip()


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
        skip_if_arm64_packages_unavailable(
            result,
            "admixture",
            "plink2",
            "plink",
            "r-ggmap",
            "r-ape",
            "bioconductor-ggtree",
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc/qc_dashboard.html",
            "results/qc/contig_map.tsv",
            "results/qc/plink.eigenvec",
            "results/qc/plink.bim_fixed",
            "results/qc/plink.3.Q",
            "results/qc/coords.txt",
            "results/qc/qc_report.tsv",
        )

        contig_map = Path(tmpdir) / "results" / "qc" / "contig_map.tsv"
        contig_lines = contig_map.read_text().strip().splitlines()
        assert contig_lines[0].split("\t") == ["original_contig", "plink_contig", "admixture_id"]
        assert [line.split("\t")[2] for line in contig_lines[1:]] == ["1", "2"]

        # Copy dashboard HTML out before tmpdir cleanup so CI can upload it
        artifacts_dir = Path("test-artifacts")
        artifacts_dir.mkdir(exist_ok=True)
        dashboard = Path(tmpdir) / "results" / "qc" / "qc_dashboard.html"
        if dashboard.exists():
            shutil.copy2(dashboard, artifacts_dir / "qc_dashboard.html")


@pytest.mark.full_run
def test_qc_numeric_contigs_full_run(request):
    """QC should normalize numeric-only contig labels and still run end-to-end."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        numeric_vcf, numeric_fai = write_numeric_qc_inputs(tmpdir)
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
                "vcf": str(numeric_vcf),
                "fai": str(numeric_fai),
                "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
            },
        )
        skip_if_arm64_packages_unavailable(
            result,
            "admixture",
            "plink2",
            "plink",
            "r-ggmap",
            "r-ape",
            "bioconductor-ggtree",
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc/contig_map.tsv",
            "results/qc/plink_input.vcf.gz",
            "results/qc/plink.bim_fixed",
            "results/qc/qc_dashboard.html",
        )

        contig_map = Path(tmpdir) / "results" / "qc" / "contig_map.tsv"
        rows = [line.split("\t") for line in contig_map.read_text().strip().splitlines()[1:]]
        assert [row[1] for row in rows] == ["qcctg1", "qcctg2"]
        assert [row[2] for row in rows] == ["1", "2"]

        bim_fixed = Path(tmpdir) / "results" / "qc" / "plink.bim_fixed"
        chroms = []
        for line in bim_fixed.read_text().strip().splitlines():
            chrom = line.split("\t", 1)[0]
            if chrom not in chroms:
                chroms.append(chrom)
        assert chroms[:2] == ["1", "2"]

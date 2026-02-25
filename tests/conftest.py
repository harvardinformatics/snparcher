# tests/conftest.py

import os
import subprocess
from pathlib import Path

import pytest

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"
CONDA_PREFIX = Path(__file__).parent / ".conda-envs"
TEST_DATA_DIR = Path(__file__).parent / "data"


def pytest_addoption(parser):
    parser.addoption(
        "--no-conda",
        action="store_true",
        default=False,
        help="Run without --use-conda (use system/pixi environment)",
    )
    parser.addoption(
        "--dry-run-only",
        action="store_true",
        default=False,
        help="Only run dry-run tests",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "dry_run: mark test as dry-run only")
    config.addinivalue_line("markers", "full_run: mark test as requiring full execution")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--dry-run-only"):
        skip_full = pytest.mark.skip(reason="--dry-run-only specified")
        for item in items:
            if "full_run" in [m.name for m in item.iter_markers()]:
                item.add_marker(skip_full)


class SnakemakeRunner:
    """Simple runner for testing snakemake workflows."""

    def __init__(self, workdir, use_conda=True):
        self.snakefile = WORKFLOW_DIR / "Snakefile"
        self.workdir = Path(workdir)
        self.use_conda = use_conda
        self.conda_prefix = CONDA_PREFIX

        # Symlink test data into workdir
        test_data_link = self.workdir / "tests" / "data"
        test_data_link.parent.mkdir(parents=True, exist_ok=True)
        if not test_data_link.exists():
            test_data_link.symlink_to(TEST_DATA_DIR)

    def run(self, target, configfile, samples=None, extra_args=None, config_overrides=None):
        if isinstance(target, str):
            targets = [target]
        else:
            targets = list(target)

        runtime_cache = self.workdir / ".snakemake-runtime-cache"
        runtime_cache.mkdir(parents=True, exist_ok=True)

        cmd = [
            "snakemake",
            "--show-failed-logs",
            "-s",
            str(self.snakefile),
            "--configfile",
            str(configfile),
            "--cores",
            "1",
            "--directory",
            str(self.workdir),
            "--runtime-source-cache-path",
            str(runtime_cache),
            "-p",
            *targets,
        ]

        # Build a single --config with all overrides
        config_pairs = []
        if samples:
            config_pairs.append(f"samples={samples}")
        if config_overrides:
            for key, value in config_overrides.items():
                config_pairs.append(f"{key}={value}")
        if config_pairs:
            cmd.extend(["--config"] + config_pairs)

        if self.use_conda:
            cmd.extend(["--use-conda", "--conda-prefix", str(self.conda_prefix)])

        # Optional shell override for test runs (e.g. /bin/zsh).
        shell_exec = os.environ.get("SNPARCHER_TEST_SHELL_EXECUTABLE")
        if shell_exec:
            cmd.extend(["--default-resources", f"shell_exec={shell_exec}"])

        if extra_args:
            cmd.extend(extra_args)

        env = os.environ.copy()
        env.setdefault("XDG_CACHE_HOME", str(self.workdir / ".cache"))
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        return SnakemakeResult(result, self.workdir)

    def dry_run(self, target, configfile, samples=None, extra_args=None, config_overrides=None):
        all_args = ["--dry-run"]
        if extra_args:
            all_args.extend(extra_args)
        return self.run(target, configfile, samples, extra_args=all_args, config_overrides=config_overrides)


class SnakemakeResult:
    def __init__(self, result, workdir):
        self.workdir = Path(workdir)
        self.returncode = result.returncode
        self.stdout = result.stdout
        self.stderr = result.stderr

    @property
    def succeeded(self):
        return self.returncode == 0

    def assert_success(self):
        assert self.succeeded, f"Snakemake failed:\nstdout: {self.stdout}\nstderr: {self.stderr}"

    def assert_output_exists(self, *paths):
        for path in paths:
            full_path = self.workdir / path
            assert full_path.exists(), f"Expected output not found: {path}"


@pytest.fixture
def temp_workdir(tmp_path):
    return tmp_path


@pytest.fixture
def runner(temp_workdir, request):
    no_conda = request.config.getoption("--no-conda")
    return SnakemakeRunner(temp_workdir, use_conda=not no_conda)

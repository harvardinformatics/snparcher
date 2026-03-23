import os
import subprocess
from pathlib import Path

import pytest

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"
TEST_DATA_DIR = Path(__file__).parent / "data"
FIXTURES_DIR = (TEST_DATA_DIR / "fixtures").resolve()
# Default conda prefix — can be overridden via --conda-prefix
CONDA_PREFIX = (Path(__file__).parent.parent / ".snakemake" / "conda").resolve()


def pytest_addoption(parser):
    parser.addoption("--no-conda", action="store_true", default=False)
    parser.addoption("--dry-run-only", action="store_true", default=False)
    parser.addoption("--conda-prefix", action="store", default=None)


def pytest_configure(config):
    config.addinivalue_line("markers", "dry_run: mark test as dry-run only")
    config.addinivalue_line("markers", "full_run: mark test as requiring full execution")
    config.addinivalue_line("markers", "unit: mark test as a unit test using fixtures")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--dry-run-only"):
        skip = pytest.mark.skip(reason="--dry-run-only specified")
        for item in items:
            markers = [m.name for m in item.iter_markers()]
            if "full_run" in markers or "unit" in markers:
                item.add_marker(skip)


class SnakemakeRunner:
    def __init__(self, workdir, use_conda=True, snakefile=None, conda_prefix=None):
        self.snakefile = Path(snakefile) if snakefile else WORKFLOW_DIR / "Snakefile"
        self.workdir = Path(workdir)
        self.use_conda = use_conda
        if conda_prefix:
            self.conda_prefix = Path(conda_prefix).resolve()
        else:
            self.conda_prefix = CONDA_PREFIX

    def link_fixtures(self, *paths):
        """Symlink specific fixture paths into the workdir.

        Each path is relative to FIXTURES_DIR. Examples:
            "config"                -> workdir/config
            "data"                  -> workdir/data
            "results/reference"     -> workdir/results/reference
            "results/bams/markdup"  -> workdir/results/bams/markdup
        """
        for rel_path in paths:
            src = FIXTURES_DIR / rel_path
            if not src.exists():
                raise FileNotFoundError(f"Fixture not found: {src}")

            dest = self.workdir / rel_path
            dest.parent.mkdir(parents=True, exist_ok=True)

            if dest.exists() or dest.is_symlink():
                continue

            dest.symlink_to(src.resolve())

    def run(self, target, configfile=None, samples=None, extra_args=None, config_overrides=None):
        if isinstance(target, str):
            targets = [target]
        else:
            targets = list(target)

        if configfile is None:
            configfile = FIXTURES_DIR / "config" / "config.yaml"

        runtime_cache = self.workdir / ".snakemake-runtime-cache"
        runtime_cache.mkdir(parents=True, exist_ok=True)

        cmd = [
            "snakemake",
            "--show-failed-logs",
            "-s",
            str(self.snakefile.resolve()),
            "--configfile",
            str(Path(configfile).resolve()),
            "--cores",
            "1",
            "--directory",
            str(self.workdir),
            "--runtime-source-cache-path",
            str(runtime_cache),
            "-p",
            *targets,
        ]

        config_pairs = []
        if samples:
            config_pairs.append(f"samples={Path(samples).resolve()}")
        if config_overrides:
            for key, value in config_overrides.items():
                config_pairs.append(f"{key}={value}")
        if config_pairs:
            cmd.extend(["--config"] + config_pairs)

        if self.use_conda:
            cmd.extend(["--use-conda", "--conda-prefix", str(self.conda_prefix)])

        if extra_args:
            cmd.extend(extra_args)

        env = os.environ.copy()
        env.setdefault("XDG_CACHE_HOME", str(self.workdir / ".cache"))
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        return SnakemakeResult(result, self.workdir)

    def dry_run(self, target, **kwargs):
        extra = kwargs.pop("extra_args", None) or []
        extra.append("--dry-run")
        return self.run(target, extra_args=extra, **kwargs)


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

    def print_log(self):
        if self.stdout:
            print(f"\n--- stdout ---\n{self.stdout}")
        if self.stderr:
            print(f"\n--- stderr ---\n{self.stderr}")
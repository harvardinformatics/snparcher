import subprocess
from pathlib import Path

import pytest

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"
CONDA_PREFIX = Path(__file__).parent / ".conda-envs"
TEST_DATA_DIR = Path(__file__).parent / "data"


class SnakemakeRunner:
    """Simple runner for testing snakemake workflows."""

    def __init__(self, workdir, use_conda=True):
        self.snakefile = WORKFLOW_DIR / "Snakefile"
        self.workdir = Path(workdir)
        self.use_conda = use_conda
        self.conda_prefix = CONDA_PREFIX

        test_data_link = self.workdir / "tests" / "data"
        test_data_link.parent.mkdir(parents=True, exist_ok=True)
        if not test_data_link.exists():
            test_data_link.symlink_to(TEST_DATA_DIR)

    def run(self, target, configfile, samples=None, extra_args=None):
        """
        Run snakemake with given target and config file.

        Args:
            target: Rule name or output file(s)
            configfile: Path to config YAML
            samples: Optional path to samples CSV (overrides config)
            extra_args: Additional snakemake arguments
        """
        if isinstance(target, str):
            targets = [target]
        else:
            targets = list(target)

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
            "-p",
            *targets,
        ]

        if samples:
            cmd.extend(["--config", f"samples={samples}"])

        if self.use_conda:
            cmd.extend(["--use-conda", "--conda-prefix", str(self.conda_prefix)])

        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(cmd, capture_output=True, text=True)
        return SnakemakeResult(result, self.workdir)

    def dry_run(self, target, configfile, samples=None):
        return self.run(target, configfile, samples, extra_args=["--dry-run"])


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
def runner(temp_workdir):
    return SnakemakeRunner(temp_workdir)

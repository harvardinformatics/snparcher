import importlib.util
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "callable_zarr_to_bed.py"
SPEC = importlib.util.spec_from_file_location("callable_zarr_to_bed", SCRIPT)
assert SPEC is not None
callable_zarr_to_bed = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(callable_zarr_to_bed)


class FakeArray:
    def __init__(self, passing: list[bool], samples: int = 2, row_chunk: int = 4):
        self._data = np.repeat(np.asarray(passing, dtype=bool)[:, None], samples, axis=1)
        self.shape = self._data.shape
        self.chunks = (row_chunk, samples)

    def __getitem__(self, key):
        return self._data[key]


def write_intervals(
    tmp_path: Path,
    passing_by_contig: dict[str, list[bool]],
    contig_order: list[str],
    merge_distance: int,
) -> list[str]:
    output_bed = write_intervals_to_path(
        tmp_path=tmp_path,
        output_name="callable.bed",
        passing_by_contig=passing_by_contig,
        contig_order=contig_order,
        merge_distance=merge_distance,
    )

    return output_bed.read_text(encoding="utf-8").splitlines()


def write_intervals_to_path(
    tmp_path: Path,
    output_name: str,
    passing_by_contig: dict[str, list[bool]],
    contig_order: list[str],
    merge_distance: int,
) -> Path:
    group = {name: FakeArray(passing) for name, passing in passing_by_contig.items()}
    contigs = [{"name": name, "length": len(passing_by_contig[name])} for name in contig_order]
    output_bed = tmp_path / output_name

    callable_zarr_to_bed.write_bed_intervals(
        group=group,
        contigs=contigs,
        total_samples=2,
        min_callable_samples=2,
        merge_distance=merge_distance,
        output_bed=output_bed,
    )

    return output_bed


def mask_from_runs(length: int, runs: list[tuple[int, int]]) -> list[bool]:
    mask = [False] * length
    for start, end in runs:
        mask[start:end] = [True] * (end - start)
    return mask


def read_nonempty_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def test_write_bed_intervals_merges_distance_across_chunk_boundaries(tmp_path):
    passing = [
        True,
        True,
        False,
        True,
        True,
        True,
        False,
        False,
        True,
        True,
        False,
        True,
        False,
        True,
        True,
    ]

    lines = write_intervals(
        tmp_path=tmp_path,
        passing_by_contig={"chr1": passing},
        contig_order=["chr1"],
        merge_distance=1,
    )

    assert lines == ["chr1\t0\t6", "chr1\t8\t15"]


def test_write_bed_intervals_keeps_bedtools_merge_zero_semantics(tmp_path):
    passing = [
        True,
        True,
        False,
        True,
        True,
        True,
        False,
        False,
        True,
        True,
        False,
        True,
        False,
        True,
        True,
    ]

    lines = write_intervals(
        tmp_path=tmp_path,
        passing_by_contig={"chr1": passing},
        contig_order=["chr1"],
        merge_distance=0,
    )

    assert lines == [
        "chr1\t0\t2",
        "chr1\t3\t6",
        "chr1\t8\t10",
        "chr1\t11\t12",
        "chr1\t13\t15",
    ]


def test_write_bed_intervals_preserves_contig_order_and_does_not_merge_between_contigs(
    tmp_path,
):
    lines = write_intervals(
        tmp_path=tmp_path,
        passing_by_contig={
            "chrB": [False, True, True],
            "chrA": [True, False, True],
        },
        contig_order=["chrB", "chrA"],
        merge_distance=100,
    )

    assert lines == ["chrB\t1\t3", "chrA\t0\t3"]


def test_direct_merge_matches_old_bedtools_pipeline_on_nonempty_records(tmp_path):
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        pytest.skip("bedtools is required for old-pipeline parity test")

    passing_by_contig = {
        "chr1": mask_from_runs(
            423,
            [
                (0, 5),
                (105, 110),  # gap 100: merge
                (211, 215),  # gap 101: do not merge
                (315, 320),  # gap 100: merge
                (421, 422),  # gap 101: do not merge
            ],
        ),
        "chr2": mask_from_runs(
            216,
            [
                (3, 4),  # single-base run
                (103, 108),  # gap 99: merge
                (208, 210),  # gap 100: merge
                (211, 215),  # gap 1: merge
            ],
        ),
        "chr3": [False] * 20,  # empty contig
        "chr4": [True] * 18,  # whole-contig interval
        "chr5": mask_from_runs(
            220,
            [
                (0, 10),
                (110, 115),  # gap 100: merge
                (216, 220),  # gap 101: do not merge
            ],
        ),
    }
    contig_order = ["chr1", "chr2", "chr3", "chr4", "chr5"]

    direct_bed = write_intervals_to_path(
        tmp_path=tmp_path,
        output_name="direct.bed",
        passing_by_contig=passing_by_contig,
        contig_order=contig_order,
        merge_distance=100,
    )
    premerge_bed = write_intervals_to_path(
        tmp_path=tmp_path,
        output_name="premerge.bed",
        passing_by_contig=passing_by_contig,
        contig_order=contig_order,
        merge_distance=0,
    )

    sorted_bed = tmp_path / "sorted.bed"
    with sorted_bed.open("w", encoding="utf-8") as stdout:
        subprocess.run([bedtools, "sort", "-i", str(premerge_bed)], check=True, stdout=stdout)

    old_pipeline_bed = tmp_path / "old-pipeline.bed"
    with old_pipeline_bed.open("w", encoding="utf-8") as stdout:
        subprocess.run(
            [bedtools, "merge", "-d", "100", "-i", str(sorted_bed)],
            check=True,
            stdout=stdout,
        )

    expected_records = [
        "chr1\t0\t110",
        "chr1\t211\t320",
        "chr1\t421\t422",
        "chr2\t3\t215",
        "chr4\t0\t18",
        "chr5\t0\t115",
        "chr5\t216\t220",
    ]
    direct_records = read_nonempty_records(direct_bed)

    assert direct_records == expected_records
    assert direct_records == read_nonempty_records(old_pipeline_bed)


def test_validate_merge_distance_rejects_negative_values():
    with pytest.raises(ValueError, match="non-negative"):
        callable_zarr_to_bed.validate_merge_distance(-1)

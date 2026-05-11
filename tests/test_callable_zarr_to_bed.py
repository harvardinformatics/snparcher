import importlib.util
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
    group = {name: FakeArray(passing) for name, passing in passing_by_contig.items()}
    contigs = [{"name": name, "length": len(passing_by_contig[name])} for name in contig_order]
    output_bed = tmp_path / "callable.bed"

    callable_zarr_to_bed.write_bed_intervals(
        group=group,
        contigs=contigs,
        total_samples=2,
        min_callable_samples=2,
        merge_distance=merge_distance,
        output_bed=output_bed,
    )

    return output_bed.read_text(encoding="utf-8").splitlines()


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


def test_validate_merge_distance_rejects_negative_values():
    with pytest.raises(ValueError, match="non-negative"):
        callable_zarr_to_bed.validate_merge_distance(-1)

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "lib"))

from staged_concat import final_stage_file, stage_chunk_counts, stage_inputs  # noqa: E402


def stage_path(round_idx, chunk_idx):
    return f"r{round_idx}/c{chunk_idx}"


FILES = [f"f{i}" for i in range(7)]


@pytest.mark.parametrize(
    "num_files, batch_size, expected",
    [
        (1, 3, []),
        (3, 3, [1]),
        (4, 3, [2, 1]),
        (7, 3, [3, 1]),
        (10, 3, [4, 2, 1]),
    ],
)
def test_stage_chunk_counts(num_files, batch_size, expected):
    assert stage_chunk_counts(num_files, batch_size, max_rounds=5) == expected


def test_stage_chunk_counts_rejects_too_many_rounds():
    with pytest.raises(ValueError, match="concat_max_rounds=2"):
        stage_chunk_counts(10, 3, max_rounds=2)


def test_stage_chunk_counts_rejects_no_files():
    with pytest.raises(ValueError, match="at least one input file"):
        stage_chunk_counts(0, 3, max_rounds=5)


def test_stage_inputs_first_round_batches_base_files():
    assert stage_inputs(FILES, 1, 0, 3, 5, stage_path) == ["f0", "f1", "f2"]
    assert stage_inputs(FILES, 1, 2, 3, 5, stage_path) == ["f6"]


def test_stage_inputs_later_rounds_gather_previous_chunks():
    # 7 files in batches of 3: round 1 has 3 chunks, round 2 gathers them.
    assert stage_inputs(FILES, 2, 0, 3, 5, stage_path) == ["r1/c0", "r1/c1", "r1/c2"]


@pytest.mark.parametrize(
    "round_idx, chunk_idx, message",
    [
        (0, 0, "Valid rounds are 1..2"),
        (3, 0, "Valid rounds are 1..2"),
        (1, 3, "Valid chunks are 0..2"),
        (2, 1, "Valid chunks are 0..0"),
    ],
)
def test_stage_inputs_rejects_out_of_range_chunks(round_idx, chunk_idx, message):
    with pytest.raises(ValueError, match=message):
        stage_inputs(FILES, round_idx, chunk_idx, 3, 5, stage_path)


def test_stage_inputs_rejects_a_single_file():
    with pytest.raises(ValueError, match="only one base file"):
        stage_inputs(["f0"], 1, 0, 3, 5, stage_path)


def test_final_stage_file():
    assert final_stage_file(["f0"], 3, 5, stage_path) == "f0"
    assert final_stage_file(FILES, 3, 5, stage_path) == "r2/c0"
    with pytest.raises(ValueError, match="No files available"):
        final_stage_file([], 3, 5, stage_path)

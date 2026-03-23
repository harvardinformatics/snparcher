#!/usr/bin/env python3
"""
Convert a clam per-sample callable Zarr into merged BED intervals.

Example:
    python scripts/callable_zarr_to_bed.py callable_masks.zarr callable.bed
    python scripts/callable_zarr_to_bed.py callable_masks.zarr callable.bed --fraction 0.8
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a clam --per-sample callable Zarr into BED intervals for sites "
            "where at least a given fraction of samples are callable."
        )
    )
    parser.add_argument(
        "callable_zarr",
        type=Path,
        help="Path to a clam callable Zarr produced with `clam loci --per-sample`.",
    )
    parser.add_argument(
        "output_bed",
        type=Path,
        help="Path to the output BED file.",
    )
    parser.add_argument(
        "--fraction",
        type=float,
        default=1.0,
        help="Minimum fraction of samples that must be callable at a site. Default: 1.0",
    )
    return parser.parse_args()


def die(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 1


def validate_fraction(fraction: float) -> None:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"--fraction must be between 0 and 1 inclusive, got {fraction}")


def import_zarr() -> Any:
    try:
        import zarr
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Python package `zarr` is required to read clam Zarr stores. "
            "Install it in the runtime environment and retry."
        ) from exc

    return zarr


def open_group(zarr_path: Path) -> Any:
    if not zarr_path.exists():
        raise FileNotFoundError(f"Zarr path does not exist: {zarr_path}")
    if not zarr_path.is_dir():
        raise ValueError(f"Zarr path is not a directory: {zarr_path}")

    zarr = import_zarr()

    try:
        return zarr.open_group(str(zarr_path), mode="r")
    except Exception as exc:
        raise ValueError(f"Failed to open Zarr group at {zarr_path}: {exc}") from exc


def get_clam_metadata(group: Any) -> dict[str, Any]:
    metadata = group.attrs.get("clam_metadata")
    if metadata is None:
        raise ValueError(
            "Missing `clam_metadata` in Zarr root. This does not look like a clam callable Zarr."
        )
    if not isinstance(metadata, dict):
        raise ValueError("`clam_metadata` is present but not a JSON object.")
    return metadata


def validate_callable_metadata(metadata: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    loci_type = metadata.get("callable_loci_type")
    if loci_type != "sample_masks":
        if loci_type == "population_counts":
            raise ValueError(
                "This callable Zarr stores population counts. Re-run `clam loci` with "
                "`--per-sample` and use that output instead."
            )
        if loci_type is None:
            raise ValueError(
                "Missing `callable_loci_type` in clam metadata. This may be a depth Zarr or "
                "an unsupported older layout."
            )
        raise ValueError(f"Unsupported callable_loci_type `{loci_type}`; expected `sample_masks`.")

    contigs = metadata.get("contigs")
    if not isinstance(contigs, list):
        raise ValueError("Missing or invalid `contigs` in clam metadata.")
    for contig in contigs:
        if not isinstance(contig, dict):
            raise ValueError("Invalid contig entry in clam metadata.")
        if not isinstance(contig.get("name"), str):
            raise ValueError("Contig metadata entry is missing a string `name`.")
        if not isinstance(contig.get("length"), int):
            raise ValueError("Contig metadata entry is missing an integer `length`.")

    column_names = metadata.get("column_names")
    if not isinstance(column_names, list) or not all(
        isinstance(column, str) for column in column_names
    ):
        raise ValueError("Missing or invalid `column_names` in clam metadata.")

    return contigs, column_names


def iter_true_runs(mask: np.ndarray) -> Iterable[tuple[int, int]]:
    if mask.ndim != 1:
        raise ValueError("Expected a 1D boolean mask.")
    if mask.size == 0:
        return

    padded = np.concatenate(([False], mask, [False]))
    transitions = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(transitions == 1)
    ends = np.flatnonzero(transitions == -1)

    for start, end in zip(starts, ends):
        yield int(start), int(end)


def write_bed_intervals(
    group: Any,
    contigs: list[dict[str, Any]],
    total_samples: int,
    min_callable_samples: int,
    output_bed: Path,
) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    with output_bed.open("w", encoding="utf-8") as handle:
        for contig in contigs:
            name = contig["name"]
            expected_length = contig["length"]

            try:
                array = group[name]
            except Exception as exc:
                raise ValueError(f"Contig array `{name}` missing from Zarr group.") from exc

            if len(array.shape) != 2:
                raise ValueError(
                    f"Contig array `{name}` must be 2D, found shape {array.shape}."
                )
            if array.shape[0] != expected_length:
                raise ValueError(
                    f"Contig `{name}` length mismatch: metadata says {expected_length}, "
                    f"array has {array.shape[0]} rows."
                )
            if array.shape[1] != total_samples:
                raise ValueError(
                    f"Contig `{name}` sample count mismatch: metadata says {total_samples}, "
                    f"array has {array.shape[1]} columns."
                )

            row_chunk = array.chunks[0] if getattr(array, "chunks", None) else array.shape[0]
            current_start: int | None = None
            current_end = 0

            for offset in range(0, expected_length, row_chunk):
                chunk = np.asarray(array[offset : offset + row_chunk, :], dtype=bool)
                callable_counts = np.count_nonzero(chunk, axis=1)
                passing = callable_counts >= min_callable_samples

                for rel_start, rel_end in iter_true_runs(passing):
                    abs_start = offset + rel_start
                    abs_end = offset + rel_end

                    if current_start is None:
                        current_start = abs_start
                        current_end = abs_end
                    elif abs_start == current_end:
                        current_end = abs_end
                    else:
                        handle.write(f"{name}\t{current_start}\t{current_end}\n")
                        current_start = abs_start
                        current_end = abs_end

            if current_start is not None:
                handle.write(f"{name}\t{current_start}\t{current_end}\n")


def main() -> int:
    args = parse_args()

    try:
        validate_fraction(args.fraction)
        group = open_group(args.callable_zarr)
        metadata = get_clam_metadata(group)
        contigs, column_names = validate_callable_metadata(metadata)
        total_samples = len(column_names)
        min_callable_samples = math.ceil(args.fraction * total_samples)

        write_bed_intervals(
            group=group,
            contigs=contigs,
            total_samples=total_samples,
            min_callable_samples=min_callable_samples,
            output_bed=args.output_bed,
        )
    except Exception as exc:
        return die(str(exc))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

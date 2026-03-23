#!/usr/bin/env python3
"""Compute cohort-wide callable coverage thresholds from mosdepth summaries."""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute cohort-wide min/max callable coverage thresholds from mosdepth "
            "summary files."
        )
    )
    parser.add_argument(
        "summaries",
        nargs="+",
        type=Path,
        help="Mosdepth summary files to aggregate.",
    )
    parser.add_argument(
        "--min-coverage",
        required=True,
        help="Minimum coverage threshold or 'auto'.",
    )
    parser.add_argument(
        "--max-coverage",
        required=True,
        help="Maximum coverage threshold or 'auto'.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path.",
    )
    return parser.parse_args()


def die(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 1


def parse_total_mean(summary_path: Path) -> float:
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary file does not exist: {summary_path}")

    total_mean = None
    fallback_mean = None
    with summary_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            if fields[0] == "chrom":
                continue
            try:
                mean = float(fields[3])
            except ValueError:
                continue
            fallback_mean = mean
            if fields[0].lower() == "total":
                total_mean = mean

    if total_mean is not None:
        return total_mean
    if fallback_mean is not None:
        return fallback_mean
    raise ValueError(f"Could not parse mean coverage from mosdepth summary: {summary_path}")


def parse_threshold(raw_value: str, *, auto_value: float, label: str) -> float:
    if raw_value == "auto":
        return auto_value

    try:
        value = float(raw_value)
    except ValueError as exc:
        raise ValueError(f"{label} must be numeric or 'auto', got {raw_value!r}") from exc

    if value < 0:
        raise ValueError(f"{label} must be >= 0, got {value}")
    return value


def format_threshold(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return str(value)


def main() -> int:
    args = parse_args()

    try:
        means = [parse_total_mean(path) for path in args.summaries]
        cohort_mean = statistics.fmean(means)
        auto_min = float(max(1, math.floor(cohort_mean / 2)))
        auto_max = float(math.ceil(cohort_mean * 2))
        min_coverage = parse_threshold(
            args.min_coverage,
            auto_value=auto_min,
            label="--min-coverage",
        )
        max_coverage = parse_threshold(
            args.max_coverage,
            auto_value=auto_max,
            label="--max-coverage",
        )
        if max_coverage < min_coverage:
            raise ValueError(
                f"Computed max coverage {max_coverage} is lower than min coverage {min_coverage}"
            )

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as handle:
            handle.write("cohort_mean_coverage\tmin_coverage\tmax_coverage\n")
            handle.write(
                f"{cohort_mean:.6f}\t{format_threshold(min_coverage)}\t"
                f"{format_threshold(max_coverage)}\n"
            )
    except Exception as exc:
        return die(str(exc))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

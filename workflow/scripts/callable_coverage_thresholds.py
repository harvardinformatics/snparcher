#!/usr/bin/env python3
"""Choose the clam loci thresholds for the callable-sites coverage mask.

A site is coverage-callable when
  * its depth summed over all samples is at most the upper cutoff, and
  * at least `fraction` of samples have at least `min_sample_depth` reads.

The upper cutoff is anchored on the genome-wide median of summed depth over
mappable sites with data. Most of a genome is single-copy, so that median is
the cohort's single-copy depth. With --max-depth auto the cutoff is an outlier
test: the median plus K robust standard deviations, measured on the
square-root scale where read-count noise is roughly constant. A number
instead sets the cutoff to that multiple of the median.

clam applies the thresholds. This script only picks them, and writes them as
a key/value TSV that the clam_loci rule reads and users can inspect.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

# 1.4826 * MAD estimates a normal standard deviation.
MAD_TO_SD = 1.4826
# sqrt(Poisson) has standard deviation of about 1/2 at any mean, so summed
# depth can never be less dispersed than this. It keeps the cutoff sensible
# when the observed spread is degenerate (e.g. near-constant depth).
MIN_SQRT_SD = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("depths_zarr", type=Path, help="Depth store from `clam collect`.")
    parser.add_argument("--mappability-bed", type=Path, help="Restrict the depth distribution to these regions.")
    parser.add_argument("--max-depth", required=True, help="'auto' or a multiple of the median summed depth.")
    parser.add_argument("--max-depth-k", type=float, required=True, help="Robust SDs above the median for 'auto'.")
    parser.add_argument("--min-sample-depth", type=int, required=True)
    parser.add_argument("--fraction", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_bed(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Sorted, merged BED intervals per contig."""
    intervals: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            contig, start, end = line.split("\t")[:3]
            intervals.setdefault(contig, []).append((int(start), int(end)))
    return {
        contig: (np.array([s for s, _ in sorted(ivs)]), np.array([e for _, e in sorted(ivs)]))
        for contig, ivs in intervals.items()
    }


def interval_mask(starts: np.ndarray, ends: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Boolean mask of [lo, hi) positions covered by the intervals."""
    mask = np.zeros(hi - lo, dtype=bool)
    first = np.searchsorted(ends, lo, side="right")
    last = np.searchsorted(starts, hi, side="left")
    for start, end in zip(starts[first:last], ends[first:last]):
        mask[max(start, lo) - lo : min(end, hi) - lo] = True
    return mask


def summed_depth_histogram(depths_zarr: Path, mappable) -> tuple[np.ndarray, int]:
    """Histogram of per-site summed depth over mappable sites with data."""
    import zarr

    group = zarr.open_group(str(depths_zarr), mode="r")
    metadata = group.attrs["clam_metadata"]
    n_samples = len(metadata["column_names"])
    hist = np.zeros(1, dtype=np.int64)

    for contig in metadata["contigs"]:
        name, length = contig["name"], contig["length"]
        array = group[name]
        rows = array.chunks[0]
        for lo in range(0, length, rows):
            hi = min(lo + rows, length)
            summed = array[lo:hi, :].sum(axis=1, dtype=np.uint64)
            keep = summed > 0
            if mappable is not None:
                if name not in mappable:
                    continue
                keep &= interval_mask(*mappable[name], lo, hi)
            counts = np.bincount(summed[keep].astype(np.int64))
            if len(counts) > len(hist):
                hist = np.pad(hist, (0, len(counts) - len(hist)))
            hist[: len(counts)] += counts

    return hist, n_samples


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, cumulative[-1] / 2)])


def main() -> int:
    args = parse_args()
    try:
        if not 0.0 <= args.fraction <= 1.0:
            raise ValueError(f"--fraction must be between 0 and 1, got {args.fraction}")

        mappable = read_bed(args.mappability_bed) if args.mappability_bed else None
        hist, n_samples = summed_depth_histogram(args.depths_zarr, mappable)
        sites = int(hist.sum())
        if sites == 0:
            raise ValueError("No mappable sites with read data; cannot estimate the depth distribution.")

        depth_values = np.flatnonzero(hist)
        weights = hist[depth_values]
        median = weighted_median(depth_values.astype(float), weights)
        sqrt_median = math.sqrt(median)
        sqrt_mad = weighted_median(np.abs(np.sqrt(depth_values) - sqrt_median), weights)
        sqrt_sd = max(MAD_TO_SD * sqrt_mad, MIN_SQRT_SD)

        # Fraction keeps the numeric bounds exact. In binary floating point
        # 2.3 * 50 is 114.99999999999999, which floors to 114, and 0.55 * 200
        # is 110.00000000000001, whose ceiling is 111. The median is a whole
        # number, since it is a sum of read counts.
        if args.max_depth == "auto":
            cutoff = (sqrt_median + args.max_depth_k * sqrt_sd) ** 2
        else:
            multiple = Fraction(args.max_depth)
            if multiple < 1:
                raise ValueError(f"--max-depth must be 'auto' or >= 1, got {args.max_depth}")
            cutoff = multiple * int(median)
        summed_cutoff = math.floor(cutoff)

        min_samples = math.ceil(Fraction(str(args.fraction)) * n_samples)

        # clam compares in floating point, and rounds the -d boundary the same
        # wrong way (it drops a site with exactly 110 of 200 samples at 0.55).
        # Both thresholds are integers, so pass each half a unit inside its
        # boundary, where floating-point error cannot land.
        thresholds = {
            "samples": n_samples,
            "sites_in_distribution": sites,
            "median_summed_depth": median,
            "sqrt_scale_sd": round(sqrt_sd, 6),
            "max_depth": args.max_depth,
            "max_depth_k": args.max_depth_k,
            "summed_depth_cutoff": summed_cutoff,
            "cutoff_over_median": round(summed_cutoff / median, 4),
            "min_sample_depth": args.min_sample_depth,
            "fraction": args.fraction,
            "min_samples_with_data": min_samples,
            "clam_max_mean_depth": (summed_cutoff + 0.5) / n_samples,
            "clam_min_proportion": max(0.0, (min_samples - 0.5) / n_samples),
        }

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as handle:
            for key, value in thresholds.items():
                handle.write(f"{key}\t{value}\n")
        print(json.dumps(thresholds, indent=2), file=sys.stderr)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

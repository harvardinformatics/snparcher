from re import I
from snakemake.exceptions import WorkflowError
import os

"""
Reads output file from ScatterIntervalsByNs and puts intervals into (relatively) specified num of equal groups.
Writes intervals groups to individual files for use by HaplotypeCaller
"""


def make_intervals(in_file: str, num_intervals: int, output_dir: str, int_output_file: str) -> None:
    intervals = []

    with open(in_file, "r") as f:
        for line in f:
            if not line.startswith("@"):
                line = line.strip().split()
                (
                    chrom,
                    start,
                    end,
                ) = (
                    line[0],
                    int(line[1]),
                    int(line[2]),
                )
                size = end - start
                intervals.append((chrom, start, end, size))

    if num_intervals > len(intervals):
        num_intervals = len(intervals)

    groups = [[] for i in range(num_intervals)]
    sums = {i: 0 for i in range(num_intervals)}
    c = 0
    for chrom, start, end, size in sorted(intervals, key=lambda x: x[3]):
        for i in sums:
            if c == sums[i]:
                groups[i].append((chrom, start, end))
                break
        sums[i] += size
        c = min(sums.values())

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(int_output_file, "w") as out:
        for i, group in enumerate(groups):
            file = os.path.join(output_dir, f"{i}.list")
            with open(file, "w") as f:
                for chrom, start, end in group:
                    print(f"{chrom}:{start}-{end}", file=f)
                    print(f"{chrom}:{start}-{end}", file=out)


def main():
    make_intervals(
        snakemake.input["in_file"],
        snakemake.params["max_intervals"],
        snakemake.output["out_dir"],
        snakemake.output["intervals"],
    )


if __name__ == "__main__":
    main()

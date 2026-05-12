#!/usr/bin/env python3
"""Utilities for Picard/GATK interval_list filtering and shard splitting."""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class IntervalRecord:
    raw: str
    contig: str


@dataclass(frozen=True)
class IntervalList:
    header: list[str]
    records: list[IntervalRecord]
    contig_lengths: dict[str, int]


def parse_interval_list(path: Path) -> IntervalList:
    header: list[str] = []
    records: list[IntervalRecord] = []
    contig_lengths: dict[str, int] = {}

    with path.open() as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped.startswith("@"):
                header.append(line)
                if stripped.startswith("@SQ"):
                    fields = stripped.split("\t")
                    contig = None
                    length = None
                    for field in fields[1:]:
                        if field.startswith("SN:"):
                            contig = field[3:]
                        elif field.startswith("LN:"):
                            try:
                                length = int(field[3:])
                            except ValueError as err:
                                raise ValueError(
                                    f"Invalid @SQ length in {path}: {field}"
                                ) from err
                    if contig is not None and length is not None:
                        contig_lengths[contig] = length
                continue

            if not stripped:
                continue

            fields = stripped.split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid interval record in {path}: {stripped}")
            records.append(IntervalRecord(raw=line, contig=fields[0]))

    return IntervalList(header=header, records=records, contig_lengths=contig_lengths)


def write_interval_list(path: Path, header: list[str], records: list[IntervalRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.writelines(header)
        for record in records:
            handle.write(record.raw)


def filter_intervals(args: argparse.Namespace) -> int:
    interval_list = parse_interval_list(args.input)

    if args.min_contig_length <= 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.input, args.output)
        retained = len(interval_list.records)
    else:
        retained_records: list[IntervalRecord] = []
        for record in interval_list.records:
            if record.contig not in interval_list.contig_lengths:
                raise ValueError(
                    f"Cannot filter {args.input}: contig {record.contig!r} has no @SQ LN entry"
                )
            if interval_list.contig_lengths[record.contig] >= args.min_contig_length:
                retained_records.append(record)

        write_interval_list(args.output, interval_list.header, retained_records)
        retained = len(retained_records)

    if retained == 0:
        raise ValueError(
            f"Interval filtering removed all intervals from {args.input}; "
            "lower intervals.min_contig_length or disable it with 0"
        )

    print(
        f"Retained {retained} of {len(interval_list.records)} interval(s) "
        f"with min contig length {args.min_contig_length}"
    )
    return 0


def would_exceed_limit(
    records: list[IntervalRecord],
    contigs: set[str],
    next_record: IntervalRecord,
    max_intervals: int,
    max_contigs: int,
) -> bool:
    if max_intervals > 0 and len(records) + 1 > max_intervals:
        return True

    return (
        max_contigs > 0
        and next_record.contig not in contigs
        and len(contigs) + 1 > max_contigs
    )


def split_records(
    records: list[IntervalRecord],
    max_intervals: int,
    max_contigs: int,
) -> list[list[IntervalRecord]]:
    if max_intervals <= 0 and max_contigs <= 0:
        return [records]

    chunks: list[list[IntervalRecord]] = []
    current: list[IntervalRecord] = []
    current_contigs: set[str] = set()

    for record in records:
        if current and would_exceed_limit(
            current, current_contigs, record, max_intervals, max_contigs
        ):
            chunks.append(current)
            current = []
            current_contigs = set()

        current.append(record)
        current_contigs.add(record.contig)

    if current:
        chunks.append(current)

    return chunks


def remove_existing_final_shards(output_dir: Path, fof: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in output_dir.glob("*-scattered.interval_list"):
        if path.is_file():
            path.unlink()
    if fof.exists():
        fof.unlink()


def split_db_intervals(args: argparse.Namespace) -> int:
    input_files = sorted(args.input_dir.glob("*-scattered.interval_list"))
    if not input_files:
        raise ValueError(f"No GATK scattered interval lists found in {args.input_dir}")

    output_chunks: list[tuple[list[str], list[IntervalRecord]]] = []
    total_records = 0
    for path in input_files:
        interval_list = parse_interval_list(path)
        total_records += len(interval_list.records)
        if not interval_list.records:
            continue
        for chunk in split_records(
            interval_list.records,
            args.max_intervals_per_shard,
            args.max_contigs_per_shard,
        ):
            output_chunks.append((interval_list.header, chunk))

    if not output_chunks:
        raise ValueError(f"No interval records found in {args.input_dir}")

    remove_existing_final_shards(args.output_dir, args.fof)

    width = max(4, len(str(len(output_chunks) - 1)))
    output_paths: list[Path] = []
    for index, (header, records) in enumerate(output_chunks):
        output_path = args.output_dir / f"{index:0{width}d}-scattered.interval_list"
        write_interval_list(output_path, header, records)
        output_paths.append(output_path)

    with args.fof.open("w") as handle:
        for path in output_paths:
            handle.write(f"{path}\n")

    print(
        f"Wrote {len(output_paths)} DB interval shard(s) from {len(input_files)} "
        f"GATK shard(s) and {total_records} interval record(s)"
    )
    return 0


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    filter_parser = subparsers.add_parser("filter")
    filter_parser.add_argument("--input", type=Path, required=True)
    filter_parser.add_argument("--output", type=Path, required=True)
    filter_parser.add_argument("--min-contig-length", type=nonnegative_int, required=True)
    filter_parser.set_defaults(func=filter_intervals)

    split_parser = subparsers.add_parser("split-db")
    split_parser.add_argument("--input-dir", type=Path, required=True)
    split_parser.add_argument("--output-dir", type=Path, required=True)
    split_parser.add_argument("--fof", type=Path, required=True)
    split_parser.add_argument(
        "--max-intervals-per-shard", type=nonnegative_int, required=True
    )
    split_parser.add_argument("--max-contigs-per-shard", type=nonnegative_int, required=True)
    split_parser.set_defaults(func=split_db_intervals)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

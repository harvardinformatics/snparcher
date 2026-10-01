#!/usr/bin/env python3
"""Normalize VCF records to compare snpArcher's RepAdapt calling model with RepAdapt.

Each record is printed as CHROM through the last sample column, tab-separated,
with the sample columns in sorted order. Records are sorted by (contig order in
the reference FASTA, POS, REF, ALT). Headers are dropped: bcftools stamps dates
into them, and RepAdapt's sample and contig order depends on Nextflow's
collect() order, which is arbitrary.

Several input VCFs may be given (for example RepAdapt's per-chromosome calls);
their records are pooled. Plain and bgzipped VCFs are both read.

The same code writes the golden record files (tests/repadapt/make_golden.sh)
and normalizes snpArcher's output in the tests, so both sides of a comparison
are normalized identically.
"""

import argparse
import gzip
import sys
from pathlib import Path


def contig_order(fasta):
    """Return {contig: index} in FASTA order."""
    order = {}
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                order[line[1:].split()[0]] = len(order)
    return order


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_records(vcf):
    """Yield records with sample columns in sorted order, as field lists."""
    samples = None
    sample_order = None
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                sample_order = sorted(range(len(samples)), key=lambda i: samples[i])
                continue
            if sample_order is None:
                raise ValueError(f"{vcf}: record before #CHROM header")
            yield fields[:9] + [fields[9 + i] for i in sample_order]


def normalized_records(vcfs, fasta, pass_only=False):
    """Return normalized record lines (no trailing newline) for the VCFs."""
    order = contig_order(fasta)
    records = []
    for vcf in vcfs:
        for fields in read_records(vcf):
            if pass_only and fields[6] != "PASS":
                continue
            if fields[0] not in order:
                raise ValueError(f"{vcf}: contig {fields[0]!r} not in {fasta}")
            records.append(fields)
    records.sort(key=lambda f: (order[f[0]], int(f[1]), f[3], f[4]))
    return ["\t".join(f) for f in records]


def parse_info(info):
    """Return {key: value} for an INFO column; flags map to True."""
    if info == ".":
        return {}
    out = {}
    for item in info.split(";"):
        key, sep, value = item.partition("=")
        out[key] = value if sep else True
    return out


def repadapt_failed_filters(fields):
    """Return the RepAdapt filters a record fails, as a sorted tuple.

    Mirrors bcftools' evaluation of RepAdapt's ``AC=AN || MQ < 30``:
    ``AC=AN`` is true if any per-ALT AC equals AN, and a missing MQ never
    satisfies ``MQ < 30``.
    """
    info = parse_info(fields[7])
    failed = []
    an = info.get("AN")
    ac = info.get("AC")
    if (
        an is not None
        and ac is not None
        and ac is not True
        and any(value == an for value in ac.split(","))
    ):
        failed.append("AllHomAlt")
    mq = info.get("MQ")
    if mq not in (None, True, ".") and float(mq) < 30:
        failed.append("LowMQ")
    return tuple(failed)


def expected_repadapt_filter(fields):
    """Return the FILTER value snpArcher's repadapt_filter should give a record."""
    return ";".join(repadapt_failed_filters(fields)) or "PASS"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("vcf", nargs="+", type=Path)
    parser.add_argument("--reference", required=True, type=Path, help="FASTA giving contig order")
    parser.add_argument("--pass-only", action="store_true", help="Keep only FILTER=PASS records")
    parser.add_argument("-o", "--output", type=Path, help="Output file (default: stdout)")
    args = parser.parse_args()

    lines = normalized_records(args.vcf, args.reference, pass_only=args.pass_only)
    text = "".join(f"{line}\n" for line in lines)
    if args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()

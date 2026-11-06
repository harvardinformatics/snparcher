import functools
import gzip
import http.server
import math
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
from contextlib import contextmanager
from pathlib import Path

import pytest
from conftest import (
    FIXTURES_DIR,
    TEST_DATA_DIR,
    WORKFLOW_DIR,
    WORKFLOW_PROFILE_DIR,
    SnakemakeRunner,
)

TEST_DIR = Path(__file__).parent
CONFIGS_DIR = TEST_DIR / "configs"
SAMPLES_DIR = TEST_DIR / "sample_sheets"
METADATA_DIR = TEST_DIR / "sample_metadata"
INTERVAL_LIST_TOOLS = WORKFLOW_DIR / "scripts" / "interval_list_tools.py"


def get_samples_file():
    """Get appropriate sample sheet based on platform."""
    if platform.machine() == "arm64":
        return SAMPLES_DIR / "local_fastqs_no_dedup.csv"
    return SAMPLES_DIR / "local_fastqs.csv"


def get_config_file():
    """Get appropriate sample sheet based on platform."""
    if platform.machine() == "arm64":
        return CONFIGS_DIR / "local_genome_no_clam.yaml"
    return CONFIGS_DIR / "local_genome.yaml"


def write_profile_with_genomicsdb_heap_override(out_dir, mem_mb_reduced):
    """Copy the default profile and override intervalized GenomicsDB Java heap."""
    profile_dir = Path(out_dir) / "workflow-profile"
    shutil.copytree(WORKFLOW_PROFILE_DIR, profile_dir)

    config_path = profile_dir / "config.yaml"
    old = (
        "  gatk_genomics_db_import:\n"
        "    mem_mb: attempt * 32000\n"
        "    mem_mb_reduced: attempt * 24000\n"
    )
    new = (
        "  gatk_genomics_db_import:\n"
        "    mem_mb: attempt * 32000\n"
        f"    mem_mb_reduced: {mem_mb_reduced}\n"
    )
    profile_config = config_path.read_text()
    if old not in profile_config:
        raise AssertionError("Expected gatk_genomics_db_import resource block not found")
    config_path.write_text(profile_config.replace(old, new))
    return profile_dir


def scheduled_rule_present(output, rule):
    return f"rule {rule}:" in output or f"checkpoint {rule}:" in output


def workflow_source(*parts):
    return (WORKFLOW_DIR.joinpath(*parts)).read_text()


def profile_resource_block(profile_dir, rule):
    config_path = Path(profile_dir) / "config.yaml"
    match = re.search(
        rf"(?ms)^  {re.escape(rule)}:\n(?P<block>(?:    .+\n)+)",
        config_path.read_text(),
    )
    if not match:
        raise AssertionError(f"Expected resource block for {rule} in {config_path}")
    return match.group("block")


def assert_profile_resource(profile_dir, rule, key, value):
    block = profile_resource_block(profile_dir, rule)
    assert f"    {key}: {value}\n" in block


def assert_profile_resource_scales_by_attempt(profile_dir, rule, key):
    block = profile_resource_block(profile_dir, rule)
    assert re.search(rf"^    {re.escape(key)}: attempt\s*\*\s*\d+\n", block, re.M)


def assert_gatk_rule_uses_profile_heap(source):
    assert "--java-options '-Xmx{resources.mem_mb_reduced}m'" in source
    assert "-Xmx3072m" not in source
    assert "mem_mb=4096" not in source


def assert_haplotypecaller_uses_explicit_bam_csi_index(source, expected_count):
    assert "**get_indexed_final_bam_input(wildcards.sample)" in source
    assert source.count("--read-index {input.bam_index}") == expected_count


def assert_glnexus_rule_uses_profile_memory(source):
    assert "glnexus_mem_gbytes=$(( {resources.mem_mb_reduced} / 1024 ))" in source
    assert '--mem-gbytes "$glnexus_mem_gbytes"' in source


def get_multistage_config_file():
    """Get multistage concat config file based on platform."""
    if platform.machine() == "arm64":
        return CONFIGS_DIR / "local_genome_no_clam_multistage.yaml"
    return CONFIGS_DIR / "local_genome_multistage.yaml"


def write_interval_list(path, contig_lengths, records):
    """Write a small Picard/GATK interval_list for script tests."""
    lines = ["@HD\tVN:1.6\n"]
    lines.extend(f"@SQ\tSN:{contig}\tLN:{length}\n" for contig, length in contig_lengths.items())
    lines.extend(f"{contig}\t{start}\t{end}\t+\tACGTmer\n" for contig, start, end in records)
    path.write_text("".join(lines))


def read_interval_records(path):
    return [
        line
        for line in path.read_text().splitlines()
        if line and not line.startswith("@")
    ]


def read_interval_header(path):
    return [line for line in path.read_text().splitlines() if line.startswith("@")]


def run_interval_list_tool(*args):
    result = subprocess.run(
        [sys.executable, str(INTERVAL_LIST_TOOLS), *map(str, args)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout
    return result


def write_interval_complexity_config(base_config, out_dir):
    """Write a config copy with the new interval complexity keys explicitly set."""
    text = Path(base_config).read_text()
    pattern = re.compile(
        r"(intervals:\n"
        r"  enabled: (?:true|false)\n"
        r"  min_nmer: \d+\n"
        r"  num_gvcf_intervals: \d+\n"
        r"  db_scatter_factor: [0-9.]+\n)"
    )
    if not pattern.search(text):
        raise AssertionError("Expected intervals block not found in config")

    out_path = Path(out_dir) / "config_interval_complexity.yaml"
    out_path.write_text(
        pattern.sub(
            r"\1"
            "  min_contig_length: 1000\n"
            "  db_max_intervals_per_shard: 200\n"
            "  db_max_contigs_per_shard: 200\n",
            text,
            count=1,
        )
    )
    return out_path


def test_db_interval_split_caps_fragmented_tail(tmp_path):
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    raw_dir.mkdir()
    contig_lengths = {f"ctg{i:04d}": 1000 for i in range(1005)}
    records = [(contig, 1, length) for contig, length in contig_lengths.items()]
    raw_path = raw_dir / "0591-scattered.interval_list"
    write_interval_list(raw_path, contig_lengths, records)

    run_interval_list_tool(
        "split-db",
        "--input-dir",
        raw_dir,
        "--output-dir",
        out_dir,
        "--fof",
        out_dir / "intervals.txt",
        "--max-intervals-per-shard",
        200,
        "--max-contigs-per-shard",
        200,
    )

    shard_paths = [Path(line) for line in (out_dir / "intervals.txt").read_text().splitlines()]
    assert len(shard_paths) == 6

    combined_records = []
    expected_header = read_interval_header(raw_path)
    for shard_path in shard_paths:
        shard_records = read_interval_records(shard_path)
        shard_contigs = {record.split("\t", 1)[0] for record in shard_records}
        assert len(shard_records) <= 200
        assert len(shard_contigs) <= 200
        assert read_interval_header(shard_path) == expected_header
        combined_records.extend(shard_records)

    assert combined_records == read_interval_records(raw_path)


def test_interval_filter_uses_contig_length_not_interval_length(tmp_path):
    input_path = tmp_path / "input.interval_list"
    output_path = tmp_path / "filtered.interval_list"
    contig_lengths = {"short_contig": 900, "long_contig": 4000}
    records = [
        ("short_contig", 1, 900),
        ("long_contig", 1, 50),
        ("long_contig", 1000, 2000),
    ]
    write_interval_list(input_path, contig_lengths, records)

    run_interval_list_tool(
        "filter",
        "--input",
        input_path,
        "--output",
        output_path,
        "--min-contig-length",
        1000,
    )

    assert read_interval_header(output_path) == read_interval_header(input_path)
    assert read_interval_records(output_path) == [
        "long_contig\t1\t50\t+\tACGTmer",
        "long_contig\t1000\t2000\t+\tACGTmer",
    ]


def test_db_interval_split_disabled_preserves_raw_shard_contents(tmp_path):
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    raw_dir.mkdir()

    first_path = raw_dir / "0010-scattered.interval_list"
    second_path = raw_dir / "0020-scattered.interval_list"
    write_interval_list(
        first_path,
        {"ctg1": 1000, "ctg2": 1000},
        [("ctg1", 1, 1000), ("ctg2", 1, 1000)],
    )
    write_interval_list(
        second_path,
        {"ctg3": 1000},
        [("ctg3", 1, 1000)],
    )

    run_interval_list_tool(
        "split-db",
        "--input-dir",
        raw_dir,
        "--output-dir",
        out_dir,
        "--fof",
        out_dir / "intervals.txt",
        "--max-intervals-per-shard",
        0,
        "--max-contigs-per-shard",
        0,
    )

    shard_paths = [Path(line) for line in (out_dir / "intervals.txt").read_text().splitlines()]
    assert [path.name for path in shard_paths] == [
        "0000-scattered.interval_list",
        "0001-scattered.interval_list",
    ]
    assert shard_paths[0].read_text() == first_path.read_text()
    assert shard_paths[1].read_text() == second_path.read_text()


def test_db_interval_split_rewrites_pathological_contig_shards(tmp_path):
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    raw_dir.mkdir()
    input_path = raw_dir / "0000-scattered.interval_list"
    write_interval_list(
        input_path,
        {"ctg1": 1000, "ctg2": 2000, "ctg3": 3000},
        [
            ("ctg1", 10, 100),
            ("ctg1", 500, 900),
            ("ctg2", 20, 200),
            ("ctg3", 30, 300),
        ],
    )

    run_interval_list_tool(
        "split-db",
        "--input-dir",
        raw_dir,
        "--output-dir",
        out_dir,
        "--fof",
        out_dir / "intervals.txt",
        "--max-intervals-per-shard",
        0,
        "--max-contigs-per-shard",
        3,
        "--merge-contigs-threshold",
        2,
    )

    shard_paths = [Path(line) for line in (out_dir / "intervals.txt").read_text().splitlines()]
    assert [path.name for path in shard_paths] == [
        "0000-scattered.interval_list",
    ]
    assert read_interval_records(shard_paths[0]) == [
        "ctg1\t1\t1000\t+\tctg1",
        "ctg2\t1\t2000\t+\tctg2",
        "ctg3\t1\t3000\t+\tctg3",
    ]
    assert read_interval_header(shard_paths[0]) == read_interval_header(input_path)


def split_db_shards(raw_dir, out_dir, max_intervals, max_contigs, merge_threshold):
    """Run split-db and return the output shard paths in order."""
    run_interval_list_tool(
        "split-db",
        "--input-dir",
        raw_dir,
        "--output-dir",
        out_dir,
        "--fof",
        out_dir / "intervals.txt",
        "--max-intervals-per-shard",
        max_intervals,
        "--max-contigs-per-shard",
        max_contigs,
        "--merge-contigs-threshold",
        merge_threshold,
    )
    return [Path(line) for line in (out_dir / "intervals.txt").read_text().splitlines()]


def assert_shards_disjoint(shard_paths):
    """Fail if two shard records overlap, since that sequence is genotyped twice."""
    spans = sorted(
        (contig, int(start), int(end), path.name)
        for path in shard_paths
        for contig, start, end, *_ in (
            record.split("\t") for record in read_interval_records(path)
        )
    )
    overlaps = [
        f"{contig} overlaps between {shard} and {next_shard}"
        for (contig, _, end, shard), (next_contig, next_start, _, next_shard) in zip(
            spans, spans[1:]
        )
        if contig == next_contig and next_start <= end
    ]
    assert not overlaps, "\n".join(overlaps)


def test_db_interval_split_rewrite_leaves_contigs_split_by_gatk(tmp_path):
    """Only contigs wholly inside a many-contig shard are rewritten.

    In a chromosome-level assembly GATK often puts the last chromosome's tail in
    the same shard as the unplaced scaffolds. Rewriting that tail to the whole
    chromosome would genotype the head of the chromosome in two shards.
    """
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    raw_dir.mkdir()
    scaffolds = [f"scaf{i:03d}" for i in range(100)]
    contig_lengths = {"chr10": 100_000, **dict.fromkeys(scaffolds, 2000), "big1": 80_000}
    write_interval_list(
        raw_dir / "0000-scattered.interval_list",
        contig_lengths,
        [("chr10", 1, 90_000)],
    )
    write_interval_list(
        raw_dir / "0001-scattered.interval_list",
        contig_lengths,
        [
            ("chr10", 90_001, 100_000),
            *[(scaffold, 1, 1500) for scaffold in scaffolds],
            ("big1", 1, 40_000),
        ],
    )
    write_interval_list(
        raw_dir / "0002-scattered.interval_list",
        contig_lengths,
        [("big1", 40_001, 80_000)],
    )

    shard_paths = split_db_shards(raw_dir, out_dir, 200, 200, 50)

    assert_shards_disjoint(shard_paths)
    assert [read_interval_records(path) for path in shard_paths] == [
        ["chr10\t1\t90000\t+\tACGTmer"],
        ["chr10\t90001\t100000\t+\tACGTmer"],
        [f"{scaffold}\t1\t2000\t+\t{scaffold}" for scaffold in scaffolds],
        ["big1\t1\t40000\t+\tACGTmer"],
        ["big1\t40001\t80000\t+\tACGTmer"],
    ]
    # The scaffold shard is still all whole contigs, so GenomicsDB merges them.
    result = run_interval_list_tool(
        "genomicsdb-merge-contigs-arg",
        "--input",
        shard_paths[2],
        "--threshold",
        50,
    )
    assert result.stdout.strip() == "--merge-contigs-into-num-partitions 50"


def test_db_interval_split_rewrite_leaves_contigs_split_by_cap(tmp_path):
    """The per-shard interval cap can split a gappy contig; it stays partial too."""
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "out"
    raw_dir.mkdir()
    scaffolds = [f"scaf{i:03d}" for i in range(60)]
    gappy = [("gappy", start, start + 49) for start in range(1, 25_000, 100)]
    write_interval_list(
        raw_dir / "0000-scattered.interval_list",
        {**dict.fromkeys(scaffolds, 2000), "gappy": 100_000},
        [*[(scaffold, 1, 1500) for scaffold in scaffolds], *gappy],
    )

    shard_paths = split_db_shards(raw_dir, out_dir, 200, 200, 50)

    assert_shards_disjoint(shard_paths)
    gappy_records = [f"gappy\t{start}\t{end}\t+\tACGTmer" for _, start, end in gappy]
    assert len(gappy_records) == 250
    assert [read_interval_records(path) for path in shard_paths] == [
        [f"{scaffold}\t1\t2000\t+\t{scaffold}" for scaffold in scaffolds],
        gappy_records[:140],
        gappy_records[140:],
    ]


def test_genomicsdb_merge_contigs_arg_only_for_whole_contig_pathology(tmp_path):
    whole_path = tmp_path / "whole.interval_list"
    partial_path = tmp_path / "partial.interval_list"
    fai_path = tmp_path / "ref.fa.fai"
    contig_lengths = {"ctg1": 1000, "ctg2": 2000, "ctg3": 3000}
    write_interval_list(
        whole_path,
        contig_lengths,
        [("ctg1", 1, 1000), ("ctg2", 1, 2000), ("ctg3", 1, 3000)],
    )
    write_interval_list(
        partial_path,
        contig_lengths,
        [("ctg1", 1, 1000), ("ctg2", 20, 200), ("ctg3", 1, 3000)],
    )

    result = run_interval_list_tool(
        "genomicsdb-merge-contigs-arg",
        "--input",
        whole_path,
        "--threshold",
        2,
    )
    assert result.stdout.strip() == "--merge-contigs-into-num-partitions 2"

    result = run_interval_list_tool(
        "genomicsdb-merge-contigs-arg",
        "--input",
        partial_path,
        "--threshold",
        2,
    )
    assert result.stdout.strip() == ""
    assert "not all records are whole-contig intervals" in result.stderr

    fai_path.write_text("ctg1\t1000\t0\t80\t81\nctg2\t2000\t1020\t80\t81\nctg3\t3000\t3040\t80\t81\n")
    result = run_interval_list_tool(
        "genomicsdb-merge-contigs-arg",
        "--input",
        fai_path,
        "--threshold",
        2,
    )
    assert result.stdout.strip() == "--merge-contigs-into-num-partitions 2"


def write_config_for_tool(base_config, out_dir, tool, parabricks_image=None):
    """Write a config copy with variant_calling.tool overridden."""
    text = Path(base_config).read_text()
    text = text.replace('tool: "gatk"', f'tool: "{tool}"', 1)
    if parabricks_image is not None:
        text = text.replace(
            'container_image: "/tmp/parabricks.sif"',
            f'container_image: "{parabricks_image}"',
            1,
        )

    out_path = Path(out_dir) / f"config_{tool}.yaml"
    out_path.write_text(text)
    return out_path


def write_long_contig_config(base_config, out_dir, tool="gatk", mode="true"):
    """Write a config copy with long-contig mode enabled or set to auto."""
    text = Path(base_config).read_text()
    text = text.replace('tool: "gatk"', f'tool: "{tool}"', 1)
    text = text.replace('tool: "gatk" #', f'tool: "{tool}" #', 1)
    pattern = re.compile(r'(variant_calling:\n(?:  .+\n)*?  tool: "[^"]+".*\n)')
    if not pattern.search(text):
        raise AssertionError("Expected variant_calling.tool entry not found")
    text = pattern.sub(rf"\1  long_contig_mode: {mode}\n", text, count=1)
    if tool == "parabricks":
        text = text.replace(
            'container_image: ""',
            'container_image: "/tmp/parabricks.sif"',
            1,
        )

    out_path = Path(out_dir) / f"config_{tool}_long_{mode}.yaml"
    out_path.write_text(text)
    return out_path


def write_auto_long_contig_reference_config(base_config, out_dir, contig_length, tool="sentieon"):
    """Write a config whose local reference has a pre-existing .fai."""
    out_dir = Path(out_dir)
    ref = out_dir / f"auto_{contig_length}.fa"
    ref.write_text(">ctg0\nA\n")
    Path(f"{ref}.fai").write_text(f"ctg0\t{contig_length}\t0\t80\t81\n")

    cfg = write_long_contig_config(base_config, out_dir, tool=tool, mode="auto")
    text = cfg.read_text()
    text = re.sub(
        r'  source: "[^"]+".*',
        f'  source: "{ref}"',
        text,
        count=1,
    )
    cfg.write_text(text)
    return cfg

def write_callable_sites_config(
    base_config,
    out_dir,
    *,
    generate_bed_file=True,
    coverage_enabled,
    mappability_enabled,
    fraction="1.0",
    min_sample_depth="1",
    max_depth="auto",
):
    """Write a config copy with callable-sites toggles overridden."""
    text = Path(base_config).read_text()
    new_block = f"""callable_sites:
  generate_bed_file: {"true" if generate_bed_file else "false"}
  coverage:
    enabled: {"true" if coverage_enabled else "false"}
    fraction: {fraction}
    min_sample_depth: {min_sample_depth}
    max_depth: {max_depth}
    merge_distance: 100
  mappability:
    enabled: {"true" if mappability_enabled else "false"}
    min_score: 1
    kmer: 150
    merge_distance: 100
"""
    pattern = re.compile(
        r"callable_sites:\n"
        r"  generate_bed_file: (?:true|false)\n"
        r"  coverage:\n"
        r"    enabled: (?:true|false)\n"
        r"    fraction: .*\n"
        r"    merge_distance: 100\n"
        r"  mappability:\n"
        r"    enabled: (?:true|false)\n"
        r"    min_score: 1\n"
        r"    kmer: 150\n"
        r"    merge_distance: 100\n"
    )
    if not pattern.search(text):
        raise AssertionError("Expected callable_sites block not found in config")

    out_path = Path(out_dir) / (
        f"config_callable_sites_{int(coverage_enabled)}_{int(mappability_enabled)}.yaml"
    )
    out_path.write_text(pattern.sub(new_block, text, count=1))
    return out_path

def write_intervals_config(base_config, out_dir, *, enabled):
    """Write a config copy with intervals.enabled overridden."""
    text = Path(base_config).read_text()
    pattern = re.compile(r"intervals:\n  enabled: (?:true|false)\n")
    if not pattern.search(text):
        raise AssertionError("Expected intervals block not found in config")

    out_path = Path(out_dir) / f"config_intervals_{int(enabled)}.yaml"
    out_path.write_text(
        pattern.sub(
            f"intervals:\n  enabled: {'true' if enabled else 'false'}\n",
            text,
            count=1,
        )
    )
    return out_path


def write_qc_config_without_exclude_scaffolds(base_config, out_dir):
    """Write a QC config copy without modules.qc.exclude_scaffolds."""
    text = Path(base_config).read_text()
    pattern = re.compile(r"\n    exclude_scaffolds: \".*\"\n")
    if not pattern.search(text):
        raise AssertionError("Expected modules.qc.exclude_scaffolds line not found in config")

    out_path = Path(out_dir) / "config_qc_no_exclude_scaffolds.yaml"
    out_path.write_text(pattern.sub("\n", text, count=1))
    return out_path


def iter_vcf_records(path):
    """Yield parsed VCF records from a plain-text or gzipped VCF."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            yield line.rstrip().split("\t")


def get_vcf_samples(path):
    """Return sample names from the VCF header."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip().split("\t")[9:]
    raise AssertionError(f"Missing #CHROM header in {path}")


def skip_if_arm64_packages_unavailable(result, *package_markers):
    """Skip full-run tests on osx-arm64 when conda cannot provide required packages."""
    if platform.machine() != "arm64":
        return

    stderr = result.stderr
    if "PackagesNotFoundError" not in stderr or "Platform: osx-arm64" not in stderr:
        return
    if package_markers and not any(marker in stderr for marker in package_markers):
        return

    missing = []
    collect_missing = False
    for line in stderr.splitlines():
        if line.strip().startswith("PackagesNotFoundError:"):
            collect_missing = True
            continue
        if not collect_missing:
            continue

        stripped = line.strip()
        if stripped.startswith("- "):
            missing.append(stripped[2:])
            continue
        if missing and not stripped:
            break

    details = ", ".join(missing) if missing else "required conda packages"
    pytest.skip(f"Conda package(s) unavailable on osx-arm64: {details}")


def write_numeric_qc_inputs(out_dir):
    """Write numeric-only QC fixtures derived from the existing standalone QC data."""
    source_vcf = TEST_DATA_DIR / "qc" / "raw.vcf.gz"
    source_fai = TEST_DATA_DIR / "qc" / "ref.fai"
    out_dir = Path(out_dir)
    out_vcf = out_dir / "numeric_raw.vcf.gz"
    out_fai = out_dir / "numeric_ref.fai"

    contigs = []
    with open(source_fai) as handle:
        for line in handle:
            fields = line.strip().split()
            if fields:
                contigs.append(fields[0])
    if len(contigs) < 2:
        raise AssertionError("Expected at least two contigs in QC FAI fixture")

    contig_map = {contig: str(i) for i, contig in enumerate(contigs, start=1)}

    with open(source_fai) as src, open(out_fai, "w") as dst:
        for line in src:
            fields = line.strip().split()
            if not fields:
                continue
            fields[0] = contig_map[fields[0]]
            dst.write("\t".join(fields) + "\n")

    sample_names = None
    contig_header_pattern = re.compile(r"^(##contig=<ID=)([^,>]+)(.*)$")
    with gzip.open(source_vcf, "rt") as src, gzip.open(out_vcf, "wt") as dst:
        for line in src:
            if line.startswith("##contig=<ID="):
                match = contig_header_pattern.match(line.rstrip("\n"))
                if match and match.group(2) in contig_map:
                    line = f"{match.group(1)}{contig_map[match.group(2)]}{match.group(3)}\n"
            elif line.startswith("#CHROM"):
                sample_names = line.rstrip().split("\t")[9:]
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if fields[0] in contig_map:
                    fields[0] = contig_map[fields[0]]
                line = "\t".join(fields) + "\n"
            dst.write(line)

        if sample_names is None:
            raise AssertionError("Missing #CHROM header in QC VCF fixture")

        extra_variant = [
            contig_map[contigs[1]],
            "1",
            ".",
            "A",
            "G",
            "60",
            "PASS",
            "AF=0.0133;ReadPosRankSum=0;FS=0;SOR=0;MQ=60;MQRankSum=0",
            "GT",
            # 4 hets in 150 samples (AF 0.0133), so it passes MAF >= 0.01:
            # QC recomputes AF from genotypes rather than trusting INFO/AF.
            *(["0/1"] * 4 + ["0/0"] * (len(sample_names) - 4)),
        ]
        dst.write("\t".join(extra_variant) + "\n")

    return out_vcf, out_fai


def write_qc_vcf_with_half_call(out_dir):
    """Write a QC VCF fixture with one GT half-call for PLINK import tests."""
    source_vcf = TEST_DATA_DIR / "qc" / "raw.vcf.gz"
    out_vcf = Path(out_dir) / "half_call_raw.vcf.gz"
    changed = False

    with gzip.open(source_vcf, "rt") as src, gzip.open(out_vcf, "wt") as dst:
        for line in src:
            if changed or line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            format_fields = fields[8].split(":")
            gt_i = format_fields.index("GT")
            sample_fields = fields[9].split(":")
            sample_fields[gt_i] = "0/."
            fields[9] = ":".join(sample_fields)
            dst.write("\t".join(fields) + "\n")
            changed = True

    if not changed:
        raise AssertionError("Expected at least one variant in QC VCF fixture")

    return out_vcf


def write_qc_vcf_without_gatk_annotations(out_dir):
    """Write a small QC VCF whose header lacks GATK-specific INFO annotations."""
    out_dir = Path(out_dir)
    out_vcf = out_dir / "missing_gatk_annotations_raw.vcf.gz"
    sample_names = [
        line.strip()
        for line in (TEST_DATA_DIR / "qc" / "samples.csv").read_text().splitlines()[1:]
        if line.strip()
    ]
    with open(TEST_DATA_DIR / "qc" / "ref.fai") as handle:
        contig = handle.readline().split()[0]

    header = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        f"##contig=<ID={contig}>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names),
    ]
    genotypes = [
        f"{'0/1' if sample_idx < 4 else '0/0'}:10"
        for sample_idx, _sample in enumerate(sample_names)
    ]
    record = [
        contig,
        "1",
        ".",
        "A",
        "G",
        "60",
        "PASS",
        "AF=0.0133;MQ=50",
        "GT:DP",
        *genotypes,
    ]

    with gzip.open(out_vcf, "wt") as handle:
        handle.write("\n".join(header))
        handle.write("\n")
        handle.write("\t".join(record))
        handle.write("\n")

    return out_vcf, sample_names


def write_qc_cohort(out_dir, *, contig_names, n_samples, filter_label=None):
    """Derive QC inputs from the standalone QC fixture.

    Renames the fixture's two contigs to `contig_names`, keeps the `n_samples`
    samples with the highest mean FORMAT/DP (well above min_depth), and can
    relabel unfiltered records (FILTER ".") as `filter_label`. Returns the VCF,
    FAI and sample sheet paths, and the number of records.
    """
    out_dir = Path(out_dir)
    source_vcf = TEST_DATA_DIR / "qc" / "raw.vcf.gz"
    with open(TEST_DATA_DIR / "qc" / "ref.fai") as handle:
        source_contigs = [line.split()[0] for line in handle if line.strip()]
    rename = dict(zip(source_contigs, contig_names, strict=True))

    with gzip.open(source_vcf, "rt") as handle:
        lines = handle.read().splitlines()
    header = [line for line in lines if line.startswith("#")]
    records = [line.split("\t") for line in lines if not line.startswith("#")]
    samples = header[-1].split("\t")[9:]

    dp_index = records[0][8].split(":").index("DP")
    depth = [0] * len(samples)
    for fields in records:
        for i, call in enumerate(fields[9:]):
            parts = call.split(":")
            if len(parts) > dp_index and parts[dp_index] not in ("", "."):
                depth[i] += int(parts[dp_index])
    keep = sorted(sorted(range(len(samples)), key=lambda i: -depth[i])[:n_samples])

    out_vcf = out_dir / "cohort_raw.vcf.gz"
    with gzip.open(out_vcf, "wt") as dst:
        for line in header[:-1]:
            for old, new in rename.items():
                line = line.replace(f"##contig=<ID={old},", f"##contig=<ID={new},")
            dst.write(line + "\n")
        dst.write("\t".join(header[-1].split("\t")[:9] + [samples[i] for i in keep]) + "\n")
        for fields in records:
            fields = fields[:9] + [fields[9 + i] for i in keep]
            fields[0] = rename[fields[0]]
            if filter_label and fields[6] == ".":
                fields[6] = filter_label
            dst.write("\t".join(fields) + "\n")

    out_fai = out_dir / "cohort_ref.fai"
    with open(TEST_DATA_DIR / "qc" / "ref.fai") as src, open(out_fai, "w") as dst:
        for line in src:
            fields = line.split()
            if fields:
                dst.write("\t".join([rename[fields[0]]] + fields[1:]) + "\n")

    out_samples = out_dir / "cohort_samples.csv"
    out_samples.write_text("sample_id\n" + "".join(f"{samples[i]}\n" for i in keep))
    return out_vcf, out_fai, out_samples, len(records)


def write_ad_only_qc_vcf(out_dir):
    """Write a QC VCF with FORMAT/AD but no per-genotype FORMAT/DP."""
    source_vcf = TEST_DATA_DIR / "qc" / "raw.vcf.gz"
    out_vcf = Path(out_dir) / "ad_only_raw.vcf.gz"

    with gzip.open(source_vcf, "rt") as src, gzip.open(out_vcf, "wt") as dst:
        for line in src:
            if line.startswith("##FORMAT=<ID=DP,"):
                continue
            if line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            format_fields = fields[8].split(":")
            if "DP" in format_fields:
                dp_index = format_fields.index("DP")
                keep_indices = [
                    index
                    for index in range(len(format_fields))
                    if index != dp_index
                ]
                fields[8] = ":".join(format_fields[index] for index in keep_indices)
                for sample_index in range(9, len(fields)):
                    genotype_fields = fields[sample_index].split(":")
                    fields[sample_index] = ":".join(
                        genotype_fields[index]
                        for index in keep_indices
                        if index < len(genotype_fields)
                    )
            dst.write("\t".join(fields) + "\n")

    return out_vcf


def get_vcf_contig_headers(path):
    """Return contig IDs declared in the VCF header."""
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = []
    pattern = re.compile(r"^##contig=<ID=([^,>]+)")
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                break
            match = pattern.match(line.rstrip())
            if match:
                contigs.append(match.group(1))
    return contigs


def extract_r_function_source(path, function_name):
    """Extract an R function body from an Rmd file by balanced braces."""
    text = Path(path).read_text()
    marker = f"{function_name} <- function"
    start = text.find(marker)
    if start == -1:
        raise AssertionError(f"Function '{function_name}' not found in {path}")

    brace_start = text.find("{", start)
    if brace_start == -1:
        raise AssertionError(f"Could not find opening brace for '{function_name}' in {path}")

    depth = 0
    end = None
    for idx in range(brace_start, len(text)):
        char = text[idx]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                end = idx + 1
                break

    if end is None:
        raise AssertionError(f"Could not find closing brace for '{function_name}' in {path}")

    return text[start:end]


def write_reference_source_config(base_config, out_dir, *, source):
    """Write a config copy with reference.source overridden."""
    text = Path(base_config).read_text()
    pattern = re.compile(
        r'(^reference:\n  name: ".*"\n  source: ").*(".*$)',
        re.MULTILINE,
    )
    if not pattern.search(text):
        raise AssertionError("Expected reference block not found in config")

    out_path = Path(out_dir) / "config_reference_source.yaml"
    out_path.write_text(pattern.sub(rf'\1{source}\2', text, count=1))
    return out_path


def write_gvcf_sample_sheet(out_dir, *, sample_id, gvcf_path):
    """Write a sample sheet for one external gVCF input."""
    out_path = Path(out_dir) / "external_gvcf_samples.csv"
    out_path.write_text(
        "sample_id,input_type,input\n"
        f"{sample_id},gvcf,{gvcf_path}\n"
    )
    return out_path


def write_mixed_bam_gvcf_sample_sheet(out_dir, *, bam_path, gvcf_path):
    """Write a sample sheet with one BAM-backed sample and one external gVCF sample."""
    out_path = Path(out_dir) / "mixed_bam_gvcf_samples.csv"
    out_path.write_text(
        "sample_id,input_type,input\n"
        f"sample_bam,bam,{bam_path}\n"
        f"sample_gvcf,gvcf,{gvcf_path}\n"
    )
    return out_path


def write_fastq_sample_sheet(out_dir, *, sample_id, library_id):
    """Write a one-row FASTQ sample sheet with explicit sample and library IDs."""
    out_path = Path(out_dir) / "numeric_id_fastqs.csv"
    out_path.write_text(
        "sample_id,input_type,input,library_id,mark_duplicates\n"
        f"{sample_id},fastq,tests/data/fastq/sample1_1.fastq.gz;tests/data/fastq/sample1_2.fastq.gz,{library_id},true\n"
    )
    return out_path


@contextmanager
def serve_directory(directory):
    """Serve a directory over HTTP for reference URL tests."""
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        host, port = sock.getsockname()

    server = http.server.ThreadingHTTPServer((host, port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://{host}:{port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


@pytest.mark.dry_run
def test_v1_style_config_exits_with_migration_message(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=CONFIGS_DIR / "legacy_v1_style.yaml",
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "Detected a v1-style snpArcher config" in output
        assert "reference/changelog/#config-key-mapping-v1-to-v2" in output
        assert "sentieon" in output


@pytest.mark.dry_run
def test_v2_config_without_variant_calling_block_is_not_v1(request):
    """`intervals` is a v1 boolean but also a v2 block. A v2 config that sets
    intervals and relies on the variant_calling defaults is not a v1 config."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        text, count = re.subn(
            r"^variant_calling:\n(?:(?:  .*)?\n)*",
            "",
            get_config_file().read_text(),
            flags=re.M,
        )
        assert count == 1 and "\nintervals:\n" in text
        cfg = Path(tmpdir) / "config_no_variant_calling.yaml"
        cfg.write_text(text)

        result = smk.dry_run(
            target="setup",
            configfile=cfg,
            samples=get_samples_file(),
        )

        result.assert_success()
        assert "Detected a v1-style snpArcher config" not in result.stdout + result.stderr


@pytest.mark.dry_run
def test_extra_config_keys_warn_but_do_not_fail(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=CONFIGS_DIR / "local_genome_extra_keys.yaml",
            samples=get_samples_file(),
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert "Ignoring unsupported config key(s):" in output
        assert "reads.sort" in output
        assert "remote_reads" in output
        assert "variant_calling.gatk.min_pruning" in output
        assert "reference.path" in output


@pytest.mark.dry_run
def test_global_mark_duplicates_config_default_is_respected(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=CONFIGS_DIR / "local_genome_global_markdup.yaml",
            samples=SAMPLES_DIR / "local_fastqs_global_markdup.csv",
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "merge_library_level_bams")
        assert not scheduled_rule_present(output, "markdup_library")
        assert not scheduled_rule_present(output, "merge_dedup_libraries")

@pytest.mark.dry_run
def test_setup_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        expected_rules = [
            "prepare_reference",
            "index_reference",
            "picard_intervals",
            "filter_picard_intervals",
            "create_gvcf_intervals",
            "create_db_intervals",
        ]
        output = result.stdout + result.stderr
        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found"


@pytest.mark.dry_run
def test_setup_dry_run_accepts_interval_complexity_config(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_interval_complexity_config(get_config_file(), tmpdir)

        result = smk.dry_run(
            target="setup",
            configfile=cfg,
            samples=get_samples_file(),
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert "filter_picard_intervals" in output


@pytest.mark.dry_run
def test_full_pipeline_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.dry_run
def test_gvcfs_target_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="gvcfs",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.dry_run
def test_gvcfs_target_explains_bcftools_has_no_gvcfs(request):
    """bcftools calls all samples together and writes no per-sample gVCFs.
    The gvcfs target must say so, not fail on missing input files."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "bcftools")

        result = smk.dry_run(
            target="gvcfs",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "writes no per-sample gVCFs, so the gvcfs target is not available" in output


@pytest.mark.dry_run
def test_multirow_same_library_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_multirow_same_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "InputFunctionException" not in output
        assert "merge_library_bams" in output
        assert "markdup_library" in output
        assert "merge_dedup_libraries" in output
        assert "input_unit=u1" in output
        assert "input_unit=u2" in output


@pytest.mark.dry_run
def test_multirow_multi_library_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_multirow_multi_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "InputFunctionException" not in output
        assert "merge_library_bams" in output
        assert "markdup_library" in output
        assert "merge_dedup_libraries" in output
        assert "library=libA" in output
        assert "library=libB" in output


@pytest.mark.dry_run
def test_mixed_srr_and_fastq_same_sample_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=SAMPLES_DIR / "local_fastqs_and_srr_same_library.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "unsupported mixed input_type values" not in output
        assert "download_sra" in output
        assert "while read -r url md5; do" in output
        assert 'curl -fSL "$url"' in output
        # The ENA path has no extractor read count, so its completeness
        # evidence is the md5 ENA publishes per file.
        assert "md5sum -c -" in output
        assert "merge_library_bams" in output
        assert "library=libA" in output
        assert "input_unit=u1" in output
        assert "input_unit=u2" in output


@pytest.mark.dry_run
@pytest.mark.parametrize(
    ("generate_bed_file", "coverage_enabled", "mappability_enabled"),
    [
        (True, True, True),
        (True, False, True),
        (True, True, False),
        (False, True, True),
    ],
)
def test_callable_sites_target_dry_run(
    request,
    generate_bed_file,
    coverage_enabled,
    mappability_enabled,
):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=generate_bed_file,
            coverage_enabled=coverage_enabled,
            mappability_enabled=mappability_enabled,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert (
            scheduled_rule_present(output, "callable_coverage_thresholds")
            == coverage_enabled
        )
        assert scheduled_rule_present(output, "clam_loci") == coverage_enabled
        assert scheduled_rule_present(output, "coverage_bed") == (
            generate_bed_file and coverage_enabled
        )
        assert scheduled_rule_present(output, "mappability_bed") == mappability_enabled
        assert scheduled_rule_present(output, "callable_sites_bed") == (
            generate_bed_file and (coverage_enabled or mappability_enabled)
        )


@pytest.mark.dry_run
def test_callable_sites_disabled_sources_warn_and_skip_final_bed(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=False,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "Skipping results/callable_sites/callable_sites.bed generation" in output
        assert not scheduled_rule_present(output, "callable_sites_bed")
        assert not scheduled_rule_present(output, "clam_loci")
        assert not scheduled_rule_present(output, "mappability_bed")


@pytest.mark.dry_run
def test_callable_sites_mixed_bam_gvcf_warns_and_skips_coverage_bed(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gvcf_path = tmp_path / "external_name.g.vcf.gz"
        samples = write_mixed_bam_gvcf_sample_sheet(
            tmpdir,
            bam_path=FIXTURES_DIR / "results/bams/markdup/sample0.bam",
            gvcf_path=gvcf_path,
        )
        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=True,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=samples,
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "gVCF input sample(s): sample_gvcf" in output
        assert "Coverage statistics will be computed only for BAM-backed samples" in output
        assert "coverage BED generation is disabled" in output
        assert "callable_coverage_thresholds" in output
        assert "clam_loci" in output
        assert "coverage_bed" not in output
        assert "mappability_bed" in output
        assert "callable_sites_bed" in output


@pytest.mark.dry_run
def test_callable_sites_mixed_bam_gvcf_coverage_only_targets_stats(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gvcf_path = tmp_path / "external_name.g.vcf.gz"
        samples = write_mixed_bam_gvcf_sample_sheet(
            tmpdir,
            bam_path=FIXTURES_DIR / "results/bams/markdup/sample0.bam",
            gvcf_path=gvcf_path,
        )
        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=samples,
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "gVCF input sample(s): sample_gvcf" in output
        assert "Coverage statistics will be computed only for BAM-backed samples" in output
        assert "no callable-sites BED sources are enabled for this cohort" in output
        assert "callable_coverage_thresholds" in output
        assert "clam_loci" in output
        assert "coverage_bed" not in output
        assert "mappability_bed" not in output
        assert "callable_sites_bed" not in output


@pytest.mark.dry_run
def test_coverage_bed_rejected_for_mixed_gvcf_inputs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gvcf_path = tmp_path / "external_name.g.vcf.gz"
        samples = write_mixed_bam_gvcf_sample_sheet(
            tmpdir,
            bam_path=FIXTURES_DIR / "results/bams/markdup/sample0.bam",
            gvcf_path=gvcf_path,
        )
        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="results/callable_sites/coverage.bed",
            configfile=cfg,
            samples=samples,
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "Cannot create results/callable_sites/coverage.bed" in output
        assert "contains gVCF inputs" in output


@pytest.mark.dry_run
def test_callable_sites_numeric_thresholds_validate(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
            fraction="0.75",
            min_sample_depth="5",
            max_depth="1.5",
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "coverage_bed" in output


@pytest.mark.dry_run
def test_callable_sites_invalid_fraction_fails_validation(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=False,
            fraction="1.5",
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "callable_sites.coverage.fraction" in output


@pytest.mark.dry_run
def test_clam_loci_uses_population_counts_and_thresholds(request):
    """clam loci runs in count mode on thresholds chosen from the config."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=True,
            mappability_enabled=True,
            fraction="0.5",
            min_sample_depth="3",
            max_depth="2",
        )

        result = smk.dry_run(
            target="results/callable_sites/callable_loci.zarr",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "--per-sample" not in output
        assert "--max-depth 2 " in output
        assert "--min-sample-depth 3 " in output
        assert "--fraction 0.5 " in output
        assert "--mappability-bed results/callable_sites/mappability.bed" in output
        assert '--max-mean-depth "$(threshold clam_max_mean_depth)"' in output
        assert '-d "$(threshold clam_min_proportion)"' in output


@pytest.mark.dry_run
def test_replaced_coverage_keys_block_only_callable_sites(request):
    """min_coverage/max_coverage skip callable sites but not the rest of the run.

    The default target still maps and calls variants, with no callable-sites
    rules. Asking for callable sites explicitly fails with the migration message.
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=True,
            mappability_enabled=True,
        )
        cfg.write_text(
            cfg.read_text().replace(
                "    merge_distance: 100\n  mappability:",
                "    merge_distance: 100\n    max_coverage: 40\n  mappability:",
                1,
            )
        )

        result = smk.dry_run(target="all", configfile=cfg, samples=get_samples_file())
        result.assert_success()
        output = result.stdout + result.stderr
        assert "callable_sites.coverage.max_coverage are no longer supported" in output
        assert "Callable sites will not be generated for this run" in output
        assert scheduled_rule_present(output, "bwa_mem")
        for rule in ("mosdepth", "clam_loci", "genmap_index", "callable_sites_bed"):
            assert not scheduled_rule_present(output, rule), rule

        result = smk.dry_run(
            target="callable_sites", configfile=cfg, samples=get_samples_file()
        )
        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "callable_sites.coverage.max_coverage are no longer supported" in output


@pytest.mark.dry_run
def test_callable_sites_leaves_out_excluded_samples(request):
    """Samples excluded in sample_metadata do not enter the coverage mask."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/callable_sites/depths.zarr",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "exclude_and_outgroup.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "results/callable_sites/depths/sample2.per-base.d4" in output
        assert "results/callable_sites/depths/sample1." not in output


@pytest.mark.dry_run
def test_gatk_without_intervals_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_intervals_config(get_config_file(), tmpdir, enabled=False)

        result = smk.dry_run(
            target="call_variants",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "gatk_haplotypecaller" in output
        assert "joint_genomics_db_import" in output

        source = workflow_source("rules", "variant_calling", "gatk.smk")
        assert_haplotypecaller_uses_explicit_bam_csi_index(source, expected_count=2)


@pytest.mark.full_run
@pytest.mark.parametrize("long_contig", [False, True], ids=["standard", "long_contig"])
def test_gatk_without_intervals_full_run(request, long_contig):
    """GATK joint calling runs end to end with intervals disabled.

    With intervals off, one GenomicsDBImport covers the whole genome. It was
    given the reference .fai as -L, which GATK rejects as an interval file, so
    every such run died after HaplotypeCaller; only dry runs covered this path.
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_intervals_config(get_config_file(), tmpdir, enabled=False)
        if long_contig:
            cfg = write_long_contig_config(cfg, tmpdir, "gatk")
        index = ".csi" if long_contig else ".tbi"

        result = smk.run(target="call_variants", configfile=cfg, samples=get_samples_file())
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/filtered.vcf.gz", f"results/vcfs/filtered.vcf.gz{index}"
        )
        assert list(iter_vcf_records(Path(tmpdir) / "results/vcfs/filtered.vcf.gz"))


@pytest.mark.dry_run
def test_joint_genomicsdb_import_dry_run_uses_profile_heap_and_contig_merge_guard(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_intervals_config(get_config_file(), tmpdir, enabled=False)

        result = smk.dry_run(
            target="call_variants",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "joint_genomics_db_import")
        assert scheduled_rule_present(output, "joint_genotype_gvcfs")
        assert "-Xmx3072m" not in output

        source = workflow_source("rules", "variant_calling", "joint_gvcf.smk")
        assert_gatk_rule_uses_profile_heap(source)
        assert "genomicsdb-merge-contigs-arg" in source
        assert "--threshold {params.merge_contig_threshold}" in source
        assert_profile_resource(
            WORKFLOW_PROFILE_DIR,
            "joint_genomics_db_import",
            "mem_mb_reduced",
            "attempt * 48000",
        )


@pytest.mark.dry_run
def test_interval_genomicsdb_import_dry_run_uses_profile_heap(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
        )

        result = smk.dry_run(target="results/gatk_genomics_db/L0000.tar")
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "gatk_genomics_db_import")
        assert "-Xmx3072m" not in output

        source = workflow_source("rules", "variant_calling", "gatk_intervals.smk")
        assert_gatk_rule_uses_profile_heap(source)
        assert_haplotypecaller_uses_explicit_bam_csi_index(source, expected_count=1)
        assert "genomicsdb-merge-contigs-arg" in source
        assert "--threshold {params.merge_contig_threshold}" in source
        assert_profile_resource(
            WORKFLOW_PROFILE_DIR,
            "gatk_genomics_db_import",
            "mem_mb_reduced",
            "attempt * 24000",
        )


@pytest.mark.dry_run
def test_interval_genomicsdb_import_profile_heap_override_is_authoritative(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        profile_dir = write_profile_with_genomicsdb_heap_override(
            tmpdir, mem_mb_reduced=12345
        )
        smk = SnakemakeRunner(
            Path(tmpdir), use_conda=not no_conda, workflow_profile=profile_dir
        )
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
        )

        result = smk.dry_run(target="results/gatk_genomics_db/L0000.tar")
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "gatk_genomics_db_import")
        assert "-Xmx3072m" not in output
        assert_profile_resource(
            profile_dir,
            "gatk_genomics_db_import",
            "mem_mb_reduced",
            "12345",
        )


@pytest.mark.full_run
@pytest.mark.parametrize("compressed", [False, True])
def test_reference_url_sources(request, compressed):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        source_dir = tmp_path / "reference_server"
        source_dir.mkdir()

        gz_ref = TEST_DATA_DIR / "reference" / "test_genome.fa.gz"
        plain_ref = source_dir / "test_genome.fa"
        plain_ref.write_text(gzip.decompress(gz_ref.read_bytes()).decode())
        shutil.copy2(gz_ref, source_dir / "test_genome.fa.gz")

        filename = "test_genome.fa.gz" if compressed else "test_genome.fa"

        with serve_directory(source_dir) as base_url:
            smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
            cfg = write_reference_source_config(
                get_config_file(),
                tmpdir,
                source=f"{base_url}/{filename}",
            )

            result = smk.run(
                target="setup",
                configfile=cfg,
                samples=get_samples_file(),
            )
            result.assert_success()
            result.assert_output_exists(
                "results/reference/test_genome.fa.gz",
                "results/reference/test_genome.fa.gz.fai",
                "results/reference/test_genome.dict",
            )


@pytest.mark.full_run
@pytest.mark.parametrize("intervals_enabled", [False, True])
@pytest.mark.parametrize("sample_id", ["sample_gvcf", "00123"])
def test_create_db_mapfile_preserves_external_gvcf_sample_id(request, intervals_enabled, sample_id):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gvcf_path = tmp_path / "external_name.g.vcf.gz"
        gvcf_path.write_text("")

        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=False,
            mappability_enabled=True,
        )
        cfg = write_intervals_config(cfg, tmpdir, enabled=intervals_enabled)
        samples = write_gvcf_sample_sheet(
            tmpdir,
            sample_id=sample_id,
            gvcf_path=gvcf_path,
        )

        result = smk.run(
            target="create_db_mapfile",
            configfile=cfg,
            samples=samples,
        )
        result.assert_success()

        mapfile = (tmp_path / "results/genomics_db/mapfile.txt").read_text().strip()
        assert mapfile == f"{sample_id}\t{gvcf_path}"


@pytest.mark.dry_run
def test_fastq_dry_run_accepts_numeric_like_sample_and_library_ids(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        samples = write_fastq_sample_sheet(tmpdir, sample_id="00123", library_id="456")

        result = smk.dry_run(
            target=[
                "results/filtered_fastqs/00123/456/u1_1.fastq.gz",
                "results/filtered_fastqs/00123/456/u1_2.fastq.gz",
            ],
            configfile=get_config_file(),
            samples=samples,
        )

        result.assert_success()


@pytest.mark.full_run
def test_setup(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()
        result.assert_output_exists(
            "results/reference/test_genome.fa.gz",
            "results/reference/test_genome.fa.gz.fai",
            "results/reference/test_genome.dict",
            "results/intervals/gvcf/intervals.txt",
            "results/intervals/db/intervals.txt",
        )


@pytest.mark.full_run
def test_full_pipeline(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()
        result.assert_output_exists(
            # gatk is GATK-lineage and generate_filtered_vcf defaults true, so the
            # hard-filtered VCF is the final call set, with its index.
            "results/vcfs/filtered.vcf.gz",
            "results/vcfs/filtered.vcf.gz.tbi",
            "results/qc_metrics/qc_report.tsv",
            "results/callable_sites/callable_sites.bed",
        )
        # filtered.vcf.gz holds every raw record with FILTER set, so the raw
        # VCF is a temporary intermediate, not a second copy of the call set.
        assert not (Path(tmpdir) / "results/vcfs/raw.vcf.gz").exists()


@pytest.mark.dry_run
def test_bcftools_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "bcftools")

        result = smk.dry_run(
            target="call_variants",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "bcftools_regions" in output
        assert "bcftools_concat_regions" in output
        # bcftools is not a GATK-lineage caller: no GATK hard filtering, so
        # call_variants resolves to the raw VCF and variant_filtration never runs.
        assert "variant_filtration" not in output
        assert "results/vcfs/raw.vcf.gz" in output


@pytest.mark.dry_run
def test_deepvariant_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "deepvariant")

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
            extra_args=["--use-apptainer"],
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "deepvariant_call" in output
        assert "/opt/deepvariant/bin/run_deepvariant" in output
        assert "glnexus_joint" in output

        source = workflow_source("rules", "variant_calling", "glnexus_gvcf.smk")
        assert_glnexus_rule_uses_profile_memory(source)
        assert_profile_resource_scales_by_attempt(
            WORKFLOW_PROFILE_DIR,
            "glnexus_joint",
            "mem_mb",
        )
        assert_profile_resource_scales_by_attempt(
            WORKFLOW_PROFILE_DIR,
            "glnexus_joint",
            "mem_mb_reduced",
        )


@pytest.mark.dry_run
def test_deepvariant_requires_apptainer(request):
    """run_deepvariant exists only in the container image, which Snakemake uses
    only when Apptainer is enabled. Without it the run must stop while the DAG
    is built, not when the first DeepVariant job starts after mapping."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "deepvariant")

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "so Snakemake must run with Apptainer enabled" in output
        assert "--use-apptainer" in output


@pytest.mark.dry_run
def test_deepvariant_removes_intermediate_results(request):
    """DeepVariant's intermediate tfrecords can take tens of GB per sample.
    Each job must write them to its own directory and remove it on exit."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "deepvariant")

        result = smk.dry_run(
            target="results/gvcfs/sample1.g.vcf.gz",
            configfile=cfg,
            samples=get_samples_file(),
            extra_args=["--use-apptainer"],
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert 'WORK=$(mktemp -d "results/deepvariant/sample1.XXXXXX")' in output
        assert "trap 'rm -rf \"$WORK\"' EXIT" in output
        assert '--intermediate_results_dir "$WORK"' in output


@pytest.mark.dry_run
def test_gatk_long_contig_dry_run_uses_work_gvcfs_and_csi_archives(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "gatk",
        )

        result = smk.dry_run(
            target="results/vcfs/raw.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "gatk_genomics_db_import")
        assert scheduled_rule_present(output, "gatk_genotype_gvcfs")
        assert not scheduled_rule_present(output, "glnexus_joint")
        assert "results/gvcfs/work/sample0.g.vcf" in output
        assert "results/gvcfs/sample0.g.vcf.gz.csi" in output

        source = workflow_source("rules", "variant_calling", "gatk_intervals.smk")
        assert "INTERVAL_GVCF_PATTERN" in source
        assert "results/interval_gvcfs/{sample}/{interval}.g.vcf" in source
        assert "gatk IndexFeatureFile -I {output.gvcf}" in source


@pytest.mark.dry_run
def test_gatk_long_contig_without_intervals_uses_work_gvcfs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "gatk",
        )
        cfg = write_intervals_config(cfg, tmpdir, enabled=False)

        result = smk.dry_run(
            target="results/vcfs/raw.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "gatk_haplotypecaller")
        assert scheduled_rule_present(output, "joint_genomics_db_import")
        assert scheduled_rule_present(output, "joint_genotype_gvcfs")
        assert scheduled_rule_present(output, "compress_joint_raw_vcf")
        assert "results/gvcfs/work/sample0.g.vcf" in output
        assert "results/gvcfs/sample0.g.vcf.gz.csi" in output


@pytest.mark.dry_run
def test_bcftools_long_contig_dry_run_uses_csi_indexes(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        smk = SnakemakeRunner(tmpdir, use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/bams/markdup",
        )
        regions = tmpdir / "results/vcfs/regions/regions.tsv"
        regions.parent.mkdir(parents=True, exist_ok=True)
        regions.write_text("L000000\tcontig0\n")
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "bcftools",
        )

        result = smk.dry_run(
            target="results/vcfs/regions/L000000.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "bcftools_call")
        assert "results/vcfs/regions/L000000.vcf.gz.csi" in output
        assert "bcftools index -f -c" in output
        source = workflow_source("rules", "variant_calling", "bcftools.smk")
        assert "tabix -p vcf" not in source


@pytest.mark.dry_run
def test_deepvariant_long_contig_dry_run_archives_compressed_csi_gvcfs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "deepvariant",
        )

        result = smk.dry_run(
            target="results/gvcfs/sample0.g.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
            extra_args=["--use-apptainer"],
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "deepvariant_call")
        assert scheduled_rule_present(output, "archive_deepvariant_gvcf")
        assert "results/deepvariant/sample0.g.vcf" in output
        assert '--intermediate_results_dir "$WORK"' in output
        assert "results/gvcfs/sample0.g.vcf.gz.csi" in output
        assert "bcftools index -f -c" in output


@pytest.mark.dry_run
def test_deepvariant_long_contig_dry_run_glnexus_consumes_csi_gvcfs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "deepvariant",
        )

        result = smk.dry_run(
            target="results/vcfs/raw.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
            extra_args=["--use-apptainer"],
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "archive_deepvariant_gvcf")
        assert scheduled_rule_present(output, "glnexus_joint")
        assert "results/gvcfs/sample0.g.vcf.gz.csi" in output
        assert "results/vcfs/raw.vcf.gz.csi" in output


@pytest.mark.dry_run
def test_long_contig_hard_filters_run_on_uncompressed_work_vcf(request):
    """GATK cannot read a CSI-indexed bgzipped VCF, so hard filtering in
    long-contig mode must consume the uncompressed work VCF and compress after."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "gatk",
        )

        result = smk.dry_run(
            target="results/vcfs/filtered.vcf.gz.csi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "variant_filtration")
        assert scheduled_rule_present(output, "compress_filtered_vcf")
        # GATK reads the plain work VCF + Tribble .idx, never the CSI-indexed one.
        assert "-V results/vcfs/work/raw.vcf " in output
        assert "-V results/vcfs/raw.vcf.gz" not in output
        assert "results/vcfs/work/filtered.vcf.idx" in output
        # GATK writes the .idx for plain-VCF output, so indexing is not disabled.
        assert "--create-output-variant-index false" not in output
        # The published call set is still bgzipped and CSI-indexed.
        assert "results/vcfs/filtered.vcf.gz.csi" in output
        assert "bcftools index -f -c results/vcfs/filtered.vcf.gz" in output


@pytest.mark.dry_run
def test_short_contig_hard_filters_filter_compressed_vcf(request):
    """Without long-contig mode, filtering keeps reading the bgzipped raw VCF."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            "gatk",
            mode="false",
        )

        result = smk.dry_run(
            target="results/vcfs/filtered.vcf.gz.tbi",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert scheduled_rule_present(output, "variant_filtration")
        assert not scheduled_rule_present(output, "compress_filtered_vcf")
        assert "-V results/vcfs/raw.vcf.gz" in output
        assert "--create-output-variant-index false" in output
        assert "bcftools index -f -t results/vcfs/filtered.vcf.gz" in output


@pytest.mark.dry_run
def test_generate_filtered_vcf_false_skips_hard_filtering(request):
    """generate_filtered_vcf: false must skip filtering outright, so downstream
    consumers fall back to the raw VCF instead of depending on a VCF that is
    never produced."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_long_contig_config(get_qc_config(), tmpdir, "gatk")
        text = Path(cfg).read_text()
        text = re.sub(
            r"(^variant_calling:\n)",
            r"\1  generate_filtered_vcf: false\n",
            text,
            count=1,
            flags=re.M,
        )
        assert "generate_filtered_vcf: false" in text
        cfg = Path(tmpdir) / "config_no_filtered_vcf.yaml"
        cfg.write_text(text)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert not scheduled_rule_present(output, "variant_filtration")
        assert not scheduled_rule_present(output, "compress_filtered_vcf")
        assert "results/vcfs/filtered.vcf.gz" not in output
        # QC reads the raw call set instead.
        assert re.search(
            r"rule qc_vcftools_individuals:\n    input: results/vcfs/raw\.vcf\.gz\n",
            output,
        )


@pytest.mark.dry_run
def test_parabricks_dry_run(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "parabricks")

        result = smk.dry_run(
            target=["call_variants", "results/gvcfs/sample1.g.vcf.gz"],
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "parabricks_haplotypecaller" in output
        assert "concat_interval_vcfs" in output


@pytest.mark.dry_run
@pytest.mark.parametrize("tool", ["bcftools", "deepvariant", "parabricks"])
def test_gvcf_input_rejected_for_new_callers(request, tool):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, tool)

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=SAMPLES_DIR / "local_gvcf.csv",
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "does not support samples with input_type='gvcf'" in output


@pytest.mark.dry_run
@pytest.mark.parametrize(
    "tool, ploidy, accepted",
    [
        ("bcftools", 4, False),
        ("parabricks", 3, False),
        ("deepvariant", 1, False),
        ("bcftools", 1, True),
        ("parabricks", 1, True),
        ("deepvariant", 2, True),
        ("gatk", 4, True),
    ],
)
def test_ploidy_validated_against_caller(request, tool, ploidy, accepted):
    """bcftools call and Parabricks HaplotypeCaller support only haploid and
    diploid calling, and DeepVariant always calls diploid genotypes. Any other
    ploidy must stop the run at startup, not fail in the caller or be ignored."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, tool)
        cfg.write_text(cfg.read_text().replace("ploidy: 2", f"ploidy: {ploidy}", 1))

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
            extra_args=["--use-apptainer"] if tool == "deepvariant" else None,
        )

        if accepted:
            result.assert_success()
        else:
            assert not result.succeeded
            output = result.stdout + result.stderr
            assert f"'{tool}' supports variant_calling.ploidy" in output
            assert f"(got {ploidy})" in output


@pytest.mark.dry_run
@pytest.mark.parametrize(
    "tool, key, flag",
    [
        ("deepvariant", "num_shards", "--num_shards"),
        ("parabricks", "num_cpu_threads", "--num-cpu-threads"),
    ],
)
def test_caller_thread_settings_reach_command(request, tool, key, flag):
    """These rules take their threads from the config. A set-threads entry in
    the default workflow profile would override the config key and make it dead."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, tool)
        text, count = re.subn(rf"{key}: \d+", f"{key}: 3", cfg.read_text())
        assert count == 1
        cfg.write_text(text)

        # Snakemake caps a job's threads at --cores; raise it above the setting.
        extra_args = ["--cores", "32"]
        if tool == "deepvariant":
            extra_args.append("--use-apptainer")
        result = smk.dry_run(
            target="results/gvcfs/sample1.g.vcf.gz",
            configfile=cfg,
            samples=get_samples_file(),
            extra_args=extra_args,
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert f"{flag} 3 " in output


@pytest.mark.dry_run
def test_callable_sites_coverage_warns_for_gvcf_only_inputs(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_config_file(),
            tmpdir,
            coverage_enabled=True,
            mappability_enabled=True,
        )

        result = smk.dry_run(
            target="callable_sites",
            configfile=cfg,
            samples=SAMPLES_DIR / "local_gvcf.csv",
        )

        result.assert_success()
        output = result.stdout + result.stderr
        assert "sample sheet contains only gVCF input sample(s): sample_gvcf" in output
        assert "Coverage statistics and coverage BED generation are disabled" in output
        assert not scheduled_rule_present(output, "callable_coverage_thresholds")
        assert not scheduled_rule_present(output, "clam_loci")
        assert not scheduled_rule_present(output, "coverage_bed")
        assert scheduled_rule_present(output, "mappability_bed")
        assert scheduled_rule_present(output, "callable_sites_bed")


@pytest.mark.dry_run
def test_parabricks_requires_container_image(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(
            get_config_file(), tmpdir, "parabricks", parabricks_image=""
        )

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "parabricks.container_image is required" in output


@pytest.mark.dry_run
@pytest.mark.parametrize("tool", ["sentieon", "parabricks"])
def test_long_contig_mode_rejects_unsupported_callers(request, tool):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_long_contig_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            tool,
        )

        result = smk.dry_run(
            target="setup",
            configfile=cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )

        assert not result.succeeded
        output = result.stdout + result.stderr
        assert f"long_contig_mode is not implemented for the {tool}" in output


@pytest.mark.dry_run
def test_long_contig_auto_detects_existing_fai(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        short_cfg = write_auto_long_contig_reference_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            contig_length=1000,
            tool="parabricks",
        )
        short_result = smk.dry_run(
            target="setup",
            configfile=short_cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )
        short_result.assert_success()

        long_cfg = write_auto_long_contig_reference_config(
            FIXTURES_DIR / "config" / "config.yaml",
            tmpdir,
            contig_length=(2**29),
            tool="parabricks",
        )
        long_result = smk.dry_run(
            target="setup",
            configfile=long_cfg,
            samples=FIXTURES_DIR / "config" / "samples.csv",
        )

        assert not long_result.succeeded
        output = long_result.stdout + long_result.stderr
        assert "long_contig_mode is not implemented for the parabricks" in output


@pytest.mark.full_run
def test_full_pipeline_bcftools(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_config_for_tool(get_config_file(), tmpdir, "bcftools")

        result = smk.run(
            target=["all", "call_variants"],
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()
        # bcftools is not GATK-lineage: no hard-filtered VCF is produced, and
        # call_variants resolves to the raw VCF, which must keep its index.
        result.assert_output_exists(
            "results/vcfs/raw.vcf.gz",
            "results/vcfs/raw.vcf.gz.tbi",
            "results/qc_metrics/qc_report.tsv",
        )
        assert not (Path(tmpdir) / "results/vcfs/filtered.vcf.gz").exists()


@pytest.mark.full_run
def test_bcftools_caller_quotes_contig_names(request):
    """SAM allows shell metacharacters in contig names, and NCBI-style
    gi|...| names are common. The per-contig mpileup must get the name quoted."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        contig = "gi|15829254|ref|AB011548.2|"
        fasta = gzip.decompress((TEST_DATA_DIR / "reference" / "test_genome.fa.gz").read_bytes())
        reference = tmp_path / "test_genome.fa"
        reference.write_text(fasta.decode().replace(">AB011548.2 ", f">{contig} ", 1))

        smk = SnakemakeRunner(tmp_path, use_conda=not no_conda)
        cfg = write_reference_source_config(get_config_file(), tmpdir, source=str(reference))
        cfg = write_config_for_tool(cfg, tmpdir, "bcftools")

        result = smk.run(target="call_variants", configfile=cfg, samples=get_samples_file())
        result.assert_success()
        result.assert_output_exists("results/vcfs/raw.vcf.gz", "results/vcfs/raw.vcf.gz.tbi")
        with gzip.open(tmp_path / "results/vcfs/raw.vcf.gz", "rt") as handle:
            chroms = [line.split("\t", 1)[0] for line in handle if not line.startswith("#")]
        assert contig in chroms


@pytest.mark.full_run
def test_incremental(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        # Run setup
        result = smk.run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        # Dry run all - setup rules should not appear
        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "prepare_reference" not in output
        assert "index_reference" not in output

        # Run all
        result = smk.run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.full_run
def test_multistage_interval_concat(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        configfile = get_multistage_config_file()
        samples = get_samples_file()

        result = smk.run(
            target="setup",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()

        result = smk.dry_run(
            target="results/vcfs/raw.vcf.gz",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "concat_interval_gvcfs_stage" in output
        assert "concat_interval_vcfs_stage" in output
        gvcf_stage_match = re.search(r"concat_interval_gvcfs_stage\s+(\d+)", output)
        assert gvcf_stage_match is not None
        assert int(gvcf_stage_match.group(1)) >= 1

        vcf_stage_match = re.search(r"concat_interval_vcfs_stage\s+(\d+)", output)
        assert vcf_stage_match is not None
        assert int(vcf_stage_match.group(1)) > 1

        result = smk.run(
            target="results/vcfs/raw.vcf.gz",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/raw.vcf.gz",
            "results/vcfs/raw.vcf.gz.tbi",
            "logs/concat_interval_gvcfs/staged/sample1/r1/c0.txt",
            "logs/concat_interval_gvcfs/staged/sample2/r1/c0.txt",
            "logs/concat_interval_vcfs/staged/r2/c0.txt",
        )


@pytest.mark.full_run
def test_long_contig_hard_filters_full_run(request):
    """Execute hard filtering with long_contig_mode forced on the small test
    genome. Long-contig mode is a pure code-path switch, so forcing it here
    exercises the uncompressed-work-VCF filtering path without a >512 Mb
    reference. A dry run cannot catch this: GATK only rejects the CSI-indexed
    input at runtime."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        configfile = write_long_contig_config(get_config_file(), tmpdir, "gatk")
        samples = get_samples_file()

        result = smk.run(target="setup", configfile=configfile, samples=samples)
        result.assert_success()

        targets = [
            "results/vcfs/filtered.vcf.gz",
            "results/vcfs/filtered.vcf.gz.csi",
        ]
        result = smk.run(target=targets, configfile=configfile, samples=samples)
        result.assert_success()
        result.assert_output_exists(*targets)

        # No TBI is produced (tabix cannot index long contigs at all).
        assert not (Path(tmpdir) / "results/vcfs/filtered.vcf.gz.tbi").exists()

        log = (Path(tmpdir) / "logs/variant_filtration.txt").read_text()
        assert "An index is required but was not found" not in log
        assert "results/vcfs/work/raw.vcf" in log

        # The published VCF carries the GATK hard-filter FILTER labels.
        with gzip.open(Path(tmpdir) / "results/vcfs/filtered.vcf.gz", "rt") as handle:
            header = [line for line in handle if line.startswith("##FILTER=")]
        declared = {
            match.group(1)
            for line in header
            for match in [re.match(r"##FILTER=<ID=([^,]+)", line)]
            if match
        }
        assert {"RPRS_filter", "FS_SOR_filter", "MQ_filter", "QUAL_filter"} <= declared

        # The raw VCF is not a product when hard filtering runs, so nothing
        # asks for it later. When `all` did, the deleted temp work VCF meant
        # re-running joint genotyping for every interval to rebuild it.
        result = smk.dry_run(target="all", configfile=configfile, samples=samples)
        result.assert_success()
        output = result.stdout + result.stderr
        for rule in ("gatk_genotype_gvcfs", "concat_interval_vcfs", "compress_interval_raw_vcf"):
            assert not scheduled_rule_present(output, rule), rule


@pytest.mark.full_run
def test_long_contig_qc_consumes_filtered_vcf_full_run(request):
    """The QC path was blocked in long-contig mode because its input VCF could
    never be built. Run QC's VCF-processing rules through to the pruned SNP set
    with long_contig_mode forced on. Stops short of plink/admixture/R so the
    test runs on every platform."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        configfile = write_long_contig_config(get_qc_config(), tmpdir, "gatk")
        samples = get_samples_file()

        result = smk.run(target="setup", configfile=configfile, samples=samples)
        result.assert_success()

        targets = [
            "results/qc/pruned.vcf.gz",
            "results/qc/snpqc.txt",
            "results/qc/individuals.idepth",
        ]
        result = smk.run(target=targets, configfile=configfile, samples=samples)
        result.assert_success()
        result.assert_output_exists(*targets)
        result.assert_output_exists("results/vcfs/filtered.vcf.gz.csi")

        samps = Path(tmpdir) / "results" / "qc" / "individuals.samps.txt"
        assert samps.read_text().split(), "QC sample filter kept no samples"


@pytest.mark.full_run
def test_long_contig_multistage_interval_gather(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        configfile = write_long_contig_config(
            get_multistage_config_file(),
            tmpdir,
            "gatk",
        )
        samples = get_samples_file()
        # A small heap makes SortVcf's records-in-RAM differ from Picard's
        # default, so the logs show it is sized from the heap and sample count.
        small_heap = [
            "--set-resources",
            "concat_interval_gvcfs_stage:mem_mb_reduced=100",
            "concat_interval_vcfs_stage:mem_mb_reduced=100",
        ]

        result = smk.run(
            target="setup",
            configfile=configfile,
            samples=samples,
        )
        result.assert_success()

        staged_gvcf_targets = [
            "results/gvcfs/work/staged/sample1/r1/c0.g.vcf",
            "results/gvcfs/work/staged/sample1/r1/c0.g.vcf.idx",
        ]
        result = smk.run(
            target=staged_gvcf_targets,
            configfile=configfile,
            samples=samples,
            extra_args=["--notemp", *small_heap],
        )
        result.assert_success()
        result.assert_output_exists(*staged_gvcf_targets)
        gvcf_run_output = result.stdout + result.stderr

        staged_vcf_targets = [
            "results/vcfs/work/staged/r1/c0.vcf",
            "results/vcfs/work/staged/r1/c0.vcf.idx",
        ]
        result = smk.run(
            target=staged_vcf_targets,
            configfile=configfile,
            samples=samples,
            extra_args=["--notemp", *small_heap],
        )
        result.assert_success()
        result.assert_output_exists(*staged_vcf_targets)

        gvcf_log = result.workdir / "logs/concat_interval_gvcfs/staged/sample1/r1/c0.txt"
        vcf_log = result.workdir / "logs/concat_interval_vcfs/staged/r1/c0.txt"
        gvcf_log_text = gvcf_log.read_text()
        vcf_log_text = vcf_log.read_text()
        assert "IndexFeatureFile" in gvcf_log_text
        assert "IndexFeatureFile" in vcf_log_text

        # Half of a 100 MB heap: one sample per gVCF, two in the cohort VCF.
        assert "-Xmx100m" in gvcf_run_output
        assert "-Xmx100m" in result.stdout + result.stderr
        assert "MAX_RECORDS_IN_RAM=22814" in gvcf_log_text
        assert "MAX_RECORDS_IN_RAM=20576" in vcf_log_text

        source = workflow_source("rules", "variant_calling", "gatk_intervals.smk")
        assert "picard SortVcf" in source
        assert "O={output.gvcf}" in source
        assert "O={output.vcf}" in source
        assert "CREATE_INDEX=false" in source
        assert "gatk GatherVcfs" not in source


# --- Metadata tests ---

@pytest.mark.dry_run
def test_no_metadata_dry_run(request):
    """Pipeline works without sample_metadata configured."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()


@pytest.mark.dry_run
def test_valid_metadata_dry_run(request):
    """Pipeline works with valid sample_metadata."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "valid_metadata.csv")},
        )
        result.assert_success()


@pytest.mark.dry_run
def test_minimal_metadata_dry_run(request):
    """Pipeline works with metadata that only has sample_id and exclude."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "minimal_metadata.csv")},
        )
        result.assert_success()


@pytest.mark.dry_run
def test_partial_metadata_warns(request):
    """Pipeline warns when some samples lack metadata rows."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "partial_metadata.csv")},
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "sample2" in output, "Expected warning about sample2 missing from metadata"


@pytest.mark.dry_run
def test_unknown_sample_in_metadata_fails(request):
    """Pipeline fails when metadata contains sample_ids not in the sample sheet."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="setup",
            configfile=get_config_file(),
            samples=get_samples_file(),
            config_overrides={"sample_metadata": str(METADATA_DIR / "unknown_sample.csv")},
        )
        assert not result.succeeded, "Expected failure for unknown sample in metadata"
        output = result.stdout + result.stderr
        assert "unknown_sample" in output


# --- Postprocess module tests ---

def get_postprocess_config():
    """Config with postprocess module enabled."""
    return CONFIGS_DIR / "local_genome_postprocess.yaml"


@pytest.mark.dry_run
def test_postprocess_dry_run(request):
    """Dry run postprocess target with callable-sites BED generation enabled."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/postprocess/clean_snps.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "callable_sites_bed" in output
        assert "postprocess_filter_individuals" in output, \
            "Expected postprocess_filter_individuals rule in DAG"
        assert "postprocess_basic_filter" in output, \
            "Expected postprocess_basic_filter rule in DAG"
        assert "postprocess_update_bed" in output, \
            "Expected postprocess_update_bed rule in DAG"


@pytest.mark.dry_run
def test_postprocess_disabled_no_rules(request):
    """When postprocess is disabled, no postprocess rules appear."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "postprocess_" not in output, \
            "Postprocess rules should not appear when module is disabled"


@pytest.mark.dry_run
def test_generate_filtered_vcf_auto_disabled_for_non_gatk(request):
    """generate_filtered_vcf: true with a non-GATK caller warns and auto-disables."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        # Explicitly request the filtered VCF with a non-GATK caller.
        text = Path(get_config_file()).read_text()
        text = text.replace(
            'tool: "gatk"',
            'tool: "bcftools"\n  generate_filtered_vcf: true',
            1,
        )
        cfg = Path(tmpdir) / "config_bcftools_generate.yaml"
        cfg.write_text(text)

        # The run still succeeds; the flag is disabled with a warning rather than erroring.
        result = smk.dry_run(target="all", configfile=cfg, samples=get_samples_file())
        result.assert_success()
        output = result.stdout + result.stderr
        assert "generate_filtered_vcf" in output  # warning was emitted
        assert "variant_filtration" not in output  # no hard filtering scheduled


@pytest.mark.dry_run
def test_generate_filtered_vcf_false_disables_filtering_everywhere(request):
    """generate_filtered_vcf: false -> no hard filtering at all, and the raw VCF
    becomes the final call set that call_variants and the modules consume."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        text = Path(get_config_file()).read_text()
        text = text.replace(
            'tool: "gatk"',
            'tool: "gatk"\n  generate_filtered_vcf: false',
            1,
        )
        cfg = Path(tmpdir) / "config_gatk_no_filtered.yaml"
        cfg.write_text(text)

        # Default target: raw only, no hard filtering.
        all_result = smk.dry_run(target="all", configfile=cfg, samples=get_samples_file())
        all_result.assert_success()
        all_output = all_result.stdout + all_result.stderr
        assert "variant_filtration" not in all_output

        # call_variants resolves to the raw VCF rather than building a filtered
        # VCF the user asked not to generate.
        cv_result = smk.dry_run(
            target="call_variants", configfile=cfg, samples=get_samples_file()
        )
        cv_result.assert_success()
        cv_output = cv_result.stdout + cv_result.stderr
        assert "variant_filtration" not in cv_output
        assert "results/vcfs/filtered.vcf.gz" not in cv_output
        assert "results/vcfs/raw.vcf.gz" in cv_output


@pytest.mark.dry_run
def test_gatk_postprocess_consumes_filtered_vcf(request):
    """GATK-lineage: postprocess consumes the shared hard-filtered VCF (built once)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        result = smk.dry_run(
            target="all",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "variant_filtration" in output
        assert "postprocess_basic_filter" in output
        # basic_filter consumes the main-workflow hard-filtered VCF.
        assert "results/vcfs/filtered.vcf.gz" in output


@pytest.mark.dry_run
def test_non_gatk_postprocess_consumes_raw_vcf(request):
    """Non-GATK: no hard filtering; postprocess consumes the raw VCF directly."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        text = Path(get_postprocess_config()).read_text()
        text = text.replace(
            'tool: "gatk"',
            'tool: "bcftools"\n  generate_filtered_vcf: false',
            1,
        )
        cfg = Path(tmpdir) / "postprocess_bcftools.yaml"
        cfg.write_text(text)

        result = smk.dry_run(target="all", configfile=cfg, samples=get_samples_file())
        result.assert_success()
        output = result.stdout + result.stderr
        assert "variant_filtration" not in output
        assert "postprocess_basic_filter" in output
        assert "results/vcfs/raw.vcf.gz" in output


@pytest.mark.dry_run
def test_postprocess_split_by_type_false_skips_subsets(request):
    """split_by_type: false keeps the strict filtered VCF but no SNP/indel split."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        text = Path(get_postprocess_config()).read_text()
        text = text.replace(
            'exclude_scaffolds: "mtDNA,Y"',
            'exclude_scaffolds: "mtDNA,Y"\n      split_by_type: false',
            1,
        )
        cfg = Path(tmpdir) / "postprocess_no_split.yaml"
        cfg.write_text(text)

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "postprocess_strict_filter" in output
        assert "postprocess_subset_snps" not in output
        assert "postprocess_subset_indels" not in output
        assert "postprocess_drop_indel_SNPs" not in output


@pytest.mark.dry_run
def test_postprocess_standalone_config_flag_coercion(request):
    """Standalone module coerces string --config flags (split_by_type=false)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        fixture_dir = TEST_DATA_DIR / "postprocess"
        result = smk.dry_run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                "samples": str(fixture_dir / "samples.csv"),
                "sample_metadata": str(fixture_dir / "sample_metadata.csv"),
                "vcf": str(fixture_dir / "raw.vcf"),
                "ref_fai": str(fixture_dir / "ref.fai"),
                "callable_sites_bed": str(fixture_dir / "callable_sites.bed"),
                # Passed as the string "false"; the module must coerce it to a bool.
                "split_by_type": "false",
            },
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "filtered.vcf.gz" in output
        assert "subset_snps" not in output
        assert "subset_indels" not in output


@pytest.mark.dry_run
def test_postprocess_with_metadata_dry_run(request):
    """Postprocess works with sample metadata (exclude column)."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/postprocess/clean_snps.vcf.gz",
            configfile=get_postprocess_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "exclude_and_outgroup.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "callable_sites_bed" in output
        assert "postprocess_filter_individuals" in output


@pytest.mark.dry_run
def test_postprocess_warns_and_disables_without_callable_bed(request):
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_callable_sites_config(
            get_postprocess_config(),
            tmpdir,
            generate_bed_file=True,
            coverage_enabled=False,
            mappability_enabled=False,
        )

        result = smk.dry_run(
            target="all",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "Disabling postprocess" in output
        assert "postprocess_" not in output


@pytest.mark.full_run
def test_postprocess_standalone_full_run(request):
    """Full execution of postprocess module against dedicated standalone fixtures."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        fixture_dir = TEST_DATA_DIR / "postprocess"
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                "samples": str(fixture_dir / "samples.csv"),
                "sample_metadata": str(fixture_dir / "sample_metadata.csv"),
                "vcf": str(fixture_dir / "raw.vcf"),
                "ref_fai": str(fixture_dir / "ref.fai"),
                "callable_sites_bed": str(fixture_dir / "callable_sites.bed"),
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        result.assert_success()
        result.assert_output_exists(
            "results/postprocess/filtered.vcf.gz",
            "results/postprocess/clean_snps.vcf.gz",
            "results/postprocess/clean_indels.vcf.gz",
        )

        filtered_vcf = Path(tmpdir) / "results" / "postprocess" / "filtered.vcf.gz"
        clean_snps_vcf = Path(tmpdir) / "results" / "postprocess" / "clean_snps.vcf.gz"
        clean_indels_vcf = Path(tmpdir) / "results" / "postprocess" / "clean_indels.vcf.gz"

        assert get_vcf_samples(filtered_vcf) == ["sample1"]
        assert get_vcf_samples(clean_snps_vcf) == ["sample1"]
        assert get_vcf_samples(clean_indels_vcf) == ["sample1"]

        snp_records = list(iter_vcf_records(clean_snps_vcf))
        indel_records = list(iter_vcf_records(clean_indels_vcf))
        assert snp_records, "Expected at least one SNP record in clean_snps.vcf.gz"
        assert indel_records, "Expected at least one indel record in clean_indels.vcf.gz"
        assert all(len(record[3]) == 1 and all(len(alt) == 1 for alt in record[4].split(",")) for record in snp_records)
        assert all(
            len(record[3]) != 1 or any(len(alt) != 1 for alt in record[4].split(","))
            for record in indel_records
        )


def write_postprocess_inputs(out_dir):
    """Small postprocess inputs where each record's fate follows by hand.

    s3 is excluded in the metadata. Callable regions: ctg1:1000-2000 and the
    whole of every other contig. contig_size is 10000 (the default):
    `small` (5000 bp) is dropped, `exact` (10000 bp) is kept.
    """
    out_dir = Path(out_dir)
    contigs = {"ctg1": 100000, "small": 5000, "exact": 10000, "mtDNA": 16000}
    records = [
        ("ctg1", 995, "AGGGGGGGGG", "A", "0/1\t0/0\t0/0"),  # POS not callable; REF reaches into the region
        ("ctg1", 1500, "A", "G", "0/1\t0/0\t0/0"),          # callable SNP
        ("ctg1", 1600, "A", "G,*", "0/2\t0/2\t1/1"),        # only * is left once s3 is dropped
        ("ctg1", 1650, "A", "G,*", "0/1\t0/2\t0/0"),        # G and * both remain: multiallelic
        ("ctg1", 1700, "A", "AT", "0/1\t0/0\t0/0"),         # callable indel
        ("ctg1", 3000, "A", "G", "0/1\t0/0\t0/0"),          # not callable
        ("exact", 100, "A", "G", "0/1\t0/0\t0/0"),          # contig exactly contig_size
        ("mtDNA", 100, "A", "G", "0/1\t0/0\t0/0"),          # excluded scaffold
        ("small", 100, "A", "G", "0/1\t0/0\t0/0"),          # contig shorter than contig_size
    ]
    vcf = out_dir / "pp_raw.vcf"
    with open(vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for name, length in contigs.items():
            handle.write(f"##contig=<ID={name},length={length}>\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n")
        for chrom, pos, ref, alt, gts in records:
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gts}\n")
    fai = out_dir / "pp_ref.fai"
    fai.write_text("".join(f"{name}\t{length}\t0\t60\t61\n" for name, length in contigs.items()))
    bed = out_dir / "pp_callable.bed"
    bed.write_text("ctg1\t1000\t2000\n" + "".join(
        f"{name}\t0\t{length}\n" for name, length in contigs.items() if name != "ctg1"
    ))
    samples = out_dir / "pp_samples.csv"
    samples.write_text("sample_id\ns1\ns2\ns3\n")
    metadata = out_dir / "pp_metadata.csv"
    metadata.write_text("sample_id,exclude\ns1,false\ns2,false\ns3,true\n")
    return {
        "samples": str(samples),
        "sample_metadata": str(metadata),
        "vcf": str(vcf),
        "ref_fai": str(fai),
        "callable_sites_bed": str(bed),
    }


@pytest.mark.full_run
def test_postprocess_filters_by_position_contig_and_scaffold(request):
    """Records are kept by their POS, small contigs and excluded scaffolds
    leave the VCF and the filtered BED alike, and a record whose only ALT
    left is * is dropped, while one that keeps a real ALT beside * is not."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                **write_postprocess_inputs(tmpdir),
                "exclude_scaffolds": "'mtDNA, Y'",
                "keep_basic_filter": "true",
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        result.assert_success()

        pp = Path(tmpdir) / "results" / "postprocess"

        def sites(name):
            return [f"{r[0]}:{r[1]}" for r in iter_vcf_records(pp / name)]

        # basic.vcf.gz has no position, contig or multiallelic filters yet.
        assert sites("basic.vcf.gz") == [
            "ctg1:995", "ctg1:1500", "ctg1:1650", "ctg1:1700", "ctg1:3000",
            "exact:100", "mtDNA:100", "small:100",
        ]
        assert sites("filtered.vcf.gz") == ["ctg1:1500", "ctg1:1700", "exact:100"]
        assert sites("clean_snps.vcf.gz") == ["ctg1:1500", "exact:100"]
        assert sites("clean_indels.vcf.gz") == ["ctg1:1700"]
        assert (pp / "callable_sites_filtered.bed").read_text() == (
            "ctg1\t1000\t2000\nexact\t0\t10000\n"
        )
        logs = Path(tmpdir) / "logs" / "postprocess"
        for rule in ("basic_filter", "update_bed", "strict_filter",
                     "subset_snps", "subset_indels", "drop_indel_snps"):
            assert (logs / f"{rule}.txt").exists(), rule


@pytest.mark.full_run
def test_postprocess_finishes_when_no_snp_survives(request):
    """Only the indel at ctg1:1700 is callable, so clean_snps.vcf.gz is empty.

    drop_indel_SNPs used to pass bcftools a plain-text list of SNP positions,
    which bcftools rejects when empty ("Failed to read the targets").
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        bed = Path(tmpdir) / "indel_only.bed"
        bed.write_text("ctg1\t1690\t1710\n")
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                **write_postprocess_inputs(tmpdir),
                "callable_sites_bed": str(bed),
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        result.assert_success()

        pp = Path(tmpdir) / "results" / "postprocess"
        assert [f"{r[0]}:{r[1]}" for r in iter_vcf_records(pp / "filtered.vcf.gz")] == ["ctg1:1700"]
        assert list(iter_vcf_records(pp / "clean_snps.vcf.gz")) == []


@pytest.mark.full_run
def test_postprocess_stops_when_no_callable_region_remains(request):
    """Every callable region is on an excluded scaffold or a small contig.

    The filtered callable BED is then empty, which bcftools rejects as a
    targets file, so update_bed stops with an explanation instead.
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                **write_postprocess_inputs(tmpdir),
                "exclude_scaffolds": "'ctg1,exact,mtDNA'",
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        assert not result.succeeded
        output = result.stdout + result.stderr
        assert "Error in rule update_bed" in output
        assert "No callable sites remain after removing contigs shorter than" in output


@pytest.mark.full_run
def test_postprocess_include_list_reruns_on_exclusion_change(request):
    """Editing exclude in the sample metadata reruns filter_individuals."""
    with tempfile.TemporaryDirectory() as tmpdir:
        inputs = write_postprocess_inputs(tmpdir)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=False,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        kwargs = dict(
            target="results/postprocess/include_samples.txt",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides=inputs,
        )
        smk.run(**kwargs).assert_success()

        metadata = Path(inputs["sample_metadata"])
        metadata.write_text(metadata.read_text().replace("s2,false", "s2,true"))
        result = smk.dry_run(**kwargs)
        result.assert_success()
        assert scheduled_rule_present(result.stdout + result.stderr, "filter_individuals")


@pytest.mark.dry_run
def test_postprocess_unsupported_filtering_key_cannot_replace_vcf(request):
    """modules.postprocess.filtering.vcf is ignored, as the warning says."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = Path(tmpdir) / "config_pp_vcf.yaml"
        text = get_postprocess_config().read_text()
        assert text.count("    filtering:\n") == 1
        cfg.write_text(text.replace("    filtering:\n", "    filtering:\n      vcf: /elsewhere/other.vcf.gz\n"))
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        result = smk.dry_run(
            target="results/postprocess/filtered.vcf.gz",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()
        output = result.stdout + result.stderr
        assert "Ignoring unsupported config key(s): modules.postprocess.filtering.vcf" in output
        assert "/elsewhere/other.vcf.gz" not in output


@pytest.mark.dry_run
def test_postprocess_maf_above_half_fails_validation(request):
    """maf > 0.5 would make the lower bound exceed the upper one."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        cfg = Path(tmpdir) / "config_pp_maf.yaml"
        text = get_postprocess_config().read_text()
        assert text.count("      maf: 0.01\n") == 1
        cfg.write_text(text.replace("      maf: 0.01\n", "      maf: 0.6\n"))
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        result = smk.dry_run(target="all", configfile=cfg, samples=get_samples_file())
        assert not result.succeeded
        assert "modules.postprocess.filtering.maf" in result.stdout + result.stderr


@pytest.mark.dry_run
def test_postprocess_module_rejects_maf_above_half():
    """A standalone run skips the main workflow's schema, so the module checks
    maf itself; at 0.6 the filter AF<0.6 | AF>0.4 removes every variant."""
    fixture_dir = TEST_DATA_DIR / "postprocess"
    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            [
                "snakemake", "-n",
                "-s", str(WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile"),
                "--directory", tmpdir,
                "--config",
                f"samples={fixture_dir / 'samples.csv'}",
                f"vcf={fixture_dir / 'raw.vcf'}",
                f"ref_fai={fixture_dir / 'ref.fai'}",
                f"callable_sites_bed={fixture_dir / 'callable_sites.bed'}",
                "maf=0.6",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "maf must be between 0 and 0.5, got 0.6" in result.stdout + result.stderr


@pytest.mark.dry_run
def test_postprocess_module_runs_without_a_configfile():
    """The module's documented standalone command works from any directory."""
    fixture_dir = TEST_DATA_DIR / "postprocess"
    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            [
                "snakemake", "-n",
                "-s", str(WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile"),
                "--directory", tmpdir,
                "--config",
                f"samples={fixture_dir / 'samples.csv'}",
                f"vcf={fixture_dir / 'raw.vcf'}",
                f"ref_fai={fixture_dir / 'ref.fai'}",
                f"callable_sites_bed={fixture_dir / 'callable_sites.bed'}",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.full_run
def test_run_postprocess_long_contig_mode_writes_csi(request):
    """run-postprocess.sh --long-contig-mode switches the module to CSI indexes."""
    if request.config.getoption("--no-conda"):
        pytest.skip("run-postprocess.sh always uses conda")
    fixture_dir = TEST_DATA_DIR / "postprocess"
    repo_root = WORKFLOW_DIR.parent
    with tempfile.TemporaryDirectory() as tmpdir:
        env = dict(os.environ, SNAKEMAKE_CONDA_PREFIX=str(repo_root / ".snakemake" / "conda"))
        result = subprocess.run(
            [
                "bash", str(repo_root / "run-postprocess.sh"), str(repo_root),
                str(fixture_dir / "samples.csv"), str(fixture_dir / "raw.vcf"),
                str(fixture_dir / "ref.fai"), str(fixture_dir / "callable_sites.bed"),
                "--long-contig-mode", "--cores", "1",
            ],
            cwd=tmpdir,
            env=env,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        pp = Path(tmpdir) / "results" / "postprocess"
        assert (pp / "filtered.vcf.gz.csi").exists()
        assert not (pp / "filtered.vcf.gz.tbi").exists()


@pytest.mark.full_run
def test_postprocess_standalone_keep_basic_and_no_split(request):
    """keep_basic_filter retains basic.vcf.gz; split_by_type=false skips subsets."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "postprocess" / "Snakefile",
        )
        fixture_dir = TEST_DATA_DIR / "postprocess"
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "postprocess" / "config" / "config.yaml",
            config_overrides={
                "samples": str(fixture_dir / "samples.csv"),
                "sample_metadata": str(fixture_dir / "sample_metadata.csv"),
                "vcf": str(fixture_dir / "raw.vcf"),
                "ref_fai": str(fixture_dir / "ref.fai"),
                "callable_sites_bed": str(fixture_dir / "callable_sites.bed"),
                "keep_basic_filter": "true",
                "split_by_type": "false",
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "bedtools")
        result.assert_success()
        result.assert_output_exists(
            "results/postprocess/filtered.vcf.gz",
            "results/postprocess/basic.vcf.gz",
        )
        # split_by_type=false -> no per-type subsets
        assert not (Path(tmpdir) / "results/postprocess/clean_snps.vcf.gz").exists()
        assert not (Path(tmpdir) / "results/postprocess/clean_indels.vcf.gz").exists()


# --- QC module tests ---

def get_qc_config():
    """Config with QC module enabled."""
    return CONFIGS_DIR / "local_genome_qc.yaml"


def get_qc_underscore_ref_config():
    """Config with QC enabled and a local reference name containing an underscore."""
    return CONFIGS_DIR / "local_genome_ref_name_underscore.yaml"


@pytest.mark.dry_run
def test_qc_dry_run(request):
    """Dry run QC dashboard target — exercises all QC rules."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_contig_map" in output, \
            "Expected qc_contig_map rule in DAG"
        assert "qc_vcftools_individuals" in output, \
            "Expected qc_vcftools_individuals rule in DAG"
        assert "qc_subsample_snps" in output, \
            "Expected qc_subsample_snps rule in DAG"
        assert "qc_prepare_plink_inputs" in output, \
            "Expected qc_prepare_plink_inputs rule in DAG"
        assert "qc_plink" in output, \
            "Expected qc_plink rule in DAG"
        assert "qc_admixture" in output, \
            "Expected qc_admixture rule in DAG"
        assert "qc_qc_dashboard" in output, \
            "Expected qc_qc_dashboard rule in DAG"
        # copy_qc_report should appear since main workflow provides qc_report
        assert "qc_copy_qc_report" in output, \
            "Expected qc_copy_qc_report rule in DAG"
        assert "bcftools query --allow-undef-tags" in output
        # QC indexes are explicitly CSI, not reliant on the bcftools default.
        assert "bcftools index -f -c results/qc/filtered.vcf.gz" in output


@pytest.mark.dry_run
def test_default_target_with_qc_enabled_runs_full_pipeline(request):
    """No explicit target should still resolve to the main workflow's all rule."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target=[],
            configfile=get_qc_underscore_ref_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "rule all:" in output, "Expected the main all rule to be the default target"
        assert "bwa_mem" in output, "Expected core mapping/variant DAG rules in the default target"
        assert "qc_qc_dashboard" in output, "Expected the QC dashboard DAG in the default target"


@pytest.mark.dry_run
def test_qc_disabled_no_rules(request):
    """When QC is disabled, no QC rules appear."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="all",
            configfile=get_config_file(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_contig_map" not in output, \
            "QC rules should not appear when module is disabled"
        assert "qc_plink" not in output, \
            "QC rules should not appear when module is disabled"


def test_qc_indexes_are_explicitly_csi():
    """QC writes CSI indexes explicitly so the module stays long-contig safe
    without depending on the bcftools default index format."""
    source = workflow_source("modules", "qc", "Snakefile")
    assert source.count("bcftools index -f -c") == 2
    assert not re.search(r"bcftools index (?!-f -c)", source)
    # Declared index outputs match what the commands write.
    assert source.count(".vcf.gz.csi") == 3
    assert ".vcf.gz.tbi" not in source


@pytest.mark.dry_run
def test_qc_no_coords_without_metadata(request):
    """QC works without metadata — generate_coords_file should not run."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
        )
        result.assert_success()

        output = result.stdout + result.stderr
        # No metadata → no coords → no generate_coords_file rule
        assert "qc_generate_coords_file" not in output, \
            "generate_coords_file should not run without metadata"


@pytest.mark.dry_run
def test_qc_with_coords_metadata(request):
    """QC generates coords when metadata has lat/long columns."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=get_qc_config(),
            samples=get_samples_file(),
            config_overrides={
                "sample_metadata": str(METADATA_DIR / "valid_metadata.csv"),
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "qc_generate_coords_file" in output, \
            "generate_coords_file should run when metadata has lat/long"


@pytest.mark.dry_run
def test_qc_defaults_exclude_scaffolds_when_omitted(request):
    """QC dry-run should succeed when modules.qc.exclude_scaffolds is omitted."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda)
        cfg = write_qc_config_without_exclude_scaffolds(get_qc_config(), tmpdir)

        result = smk.dry_run(
            target="results/qc/qc_dashboard.html",
            configfile=cfg,
            samples=get_samples_file(),
        )
        result.assert_success()


def test_write_numeric_qc_inputs_rewrites_contigs():
    """The synthetic numeric QC fixture should rewrite both header and record contigs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        numeric_vcf, numeric_fai = write_numeric_qc_inputs(tmpdir)

        headers = get_vcf_contig_headers(numeric_vcf)
        assert "1" in headers
        assert "2" in headers
        assert "JAKDEW010000001.1" not in headers
        assert "JAKDEW010000002.1" not in headers

        records = list(iter_vcf_records(numeric_vcf))
        chroms = {record[0] for record in records}
        assert "1" in chroms
        assert "2" in chroms

        fai_lines = Path(numeric_fai).read_text().strip().splitlines()
        assert fai_lines[:2] == ["1\t1\t3", "2\t1\t3"]


def test_qc_dashboard_helper_preserves_numeric_like_ids():
    """QC dashboard helper should preserve leading zeros for headered and headerless tables."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not available")

    helper_source = extract_r_function_source(
        WORKFLOW_DIR / "modules" / "qc" / "scripts" / "qc_dashboard_interactive.Rmd",
        "read_table_preserve_ids",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        headerless = tmp_path / "dist.id"
        headerless.write_text("0 000123\n0 000456\n")
        headered = tmp_path / "depth.tsv"
        headered.write_text("INDV\tMEAN_DEPTH\n000123\t5.2\n")

        script = tmp_path / "validate_read_table_preserve_ids.R"
        script.write_text(
            "\n".join(
                [
                    "args <- commandArgs(trailingOnly = TRUE)",
                    "headerless <- args[[1]]",
                    "headered <- args[[2]]",
                    helper_source,
                    "headerless_df <- read_table_preserve_ids(headerless, id_cols = c('V1', 'V2'))",
                    "if (!is.character(headerless_df$V1) || !is.character(headerless_df$V2)) stop('Headerless ID columns were not read as character')",
                    "if (!identical(headerless_df$V2[[1]], '000123')) stop('Headerless leading zeros were not preserved')",
                    "headered_df <- read_table_preserve_ids(headered, header = TRUE, sep = '\\t', id_cols = c('INDV'))",
                    "if (!is.character(headered_df$INDV)) stop('Headered ID column was not read as character')",
                    "if (!identical(headered_df$INDV[[1]], '000123')) stop('Headered leading zeros were not preserved')",
                ]
            )
        )

        result = subprocess.run(
            ["Rscript", str(script), str(headerless), str(headered)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, (result.stdout + result.stderr).strip()


def test_qc_dashboard_depth_helpers_handle_nonfinite_values():
    """QC dashboard helpers should sanitize non-finite depth and avoid empty lm fits."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not available")

    rmd_path = WORKFLOW_DIR / "modules" / "qc" / "scripts" / "qc_dashboard_interactive.Rmd"
    helper_source = "\n".join(
        [
            extract_r_function_source(rmd_path, "sanitize_depth_values"),
            extract_r_function_source(rmd_path, "depth_pc_r2"),
        ]
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        script = Path(tmpdir) / "validate_depth_helpers.R"
        script.write_text(
            "\n".join(
                [
                    helper_source,
                    "depth <- sanitize_depth_values(c('5.5', 'NaN', '-nan', 'Inf', '-Inf', 'NA', '.'))",
                    "if (!isTRUE(all.equal(depth[[1]], 5.5))) stop('Finite depth was not preserved')",
                    "if (!all(is.na(depth[2:7]))) stop('Non-finite depth values were not converted to NA')",
                    "empty_depth <- data.frame(MEAN_DEPTH = c(NA_real_, NaN, Inf), value = c(1, 2, 3))",
                    "if (!is.na(depth_pc_r2(empty_depth))) stop('Expected NA R2 for zero finite depth rows')",
                    "mixed_depth <- data.frame(MEAN_DEPTH = c(1, 2, NA, 4), value = c(1, 2, 3, 4))",
                    "if (!is.finite(depth_pc_r2(mixed_depth))) stop('Expected finite R2 for complete finite rows')",
                ]
            )
        )

        result = subprocess.run(
            ["Rscript", str(script)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, (result.stdout + result.stderr).strip()


@pytest.mark.dry_run
def test_qc_plink_filters_sparse_samples_before_pca(request):
    """QC PLINK step should pass a per-sample missingness filter to avoid NaN GRMs."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.dry_run(
            target="results/qc/plink.bed",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
                "vcf": str(TEST_DATA_DIR / "qc" / "raw.vcf.gz"),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
                "max_sample_missingness": "0.49",
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "--mind 0.49" in output
        assert "--vcf-half-call missing" in output


@pytest.mark.full_run
def test_qc_plink_accepts_half_call_genotypes(request):
    """QC PLINK import should treat VCF half-calls as missing genotypes."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        half_call_vcf = write_qc_vcf_with_half_call(tmpdir)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="results/qc/plink.bed",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
                "vcf": str(half_call_vcf),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools", "vcftools", "plink2", "plink")
        result.assert_success()
        result.assert_output_exists(
            "results/qc/plink.bed",
            "results/qc/plink.bim",
            "results/qc/plink.fam",
        )


@pytest.mark.full_run
def test_qc_subsample_snps_allows_missing_gatk_annotations(request):
    """SNP QC export should emit dots when caller-specific annotations are absent."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        missing_annotations_vcf, sample_names = write_qc_vcf_without_gatk_annotations(tmpdir)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        sample_filter = Path(tmpdir) / "results" / "qc" / "individuals.samps.txt"
        sample_filter.parent.mkdir(parents=True, exist_ok=True)
        sample_filter.write_text("\n".join(sample_names) + "\n")
        result = smk.run(
            target="results/qc/snpqc.txt",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "vcf": str(missing_annotations_vcf),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
            },
        )
        skip_if_arm64_packages_unavailable(result, "bcftools")
        result.assert_success()
        result.assert_output_exists("results/qc/snpqc.txt")

        snpqc = Path(tmpdir) / "results" / "qc" / "snpqc.txt"
        rows = [
            line.rstrip("\n").split("\t")
            for line in snpqc.read_text().splitlines()
            if line.strip()
        ]
        assert rows, "Expected SNP QC rows"
        assert all(len(row) == 10 for row in rows)

        with open(TEST_DATA_DIR / "qc" / "ref.fai") as handle:
            contig = handle.readline().split()[0]
        first = rows[0]
        assert first[0] == contig
        assert first[1] == "1"
        # QC recomputes AF over the samples it keeps: 4 hets in 150 samples.
        assert first[3] == "0.0133333"
        assert first[4] == "60"
        assert first[5:8] == [".", ".", "."]
        assert first[8] == "50"
        assert first[9] == "."


@pytest.mark.full_run
def test_qc_ad_only_vcf_computes_finite_depth(request):
    """AD-only VCFs should produce finite SNP depth from summed FORMAT/AD values."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        ad_only_vcf = write_ad_only_qc_vcf(tmpdir)
        config_overrides = {
            "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
            "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
            "vcf": str(ad_only_vcf),
            "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
            "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
        }
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        depth_result = smk.run(
            target="results/qc/individuals.idepth",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides=config_overrides,
        )
        skip_if_arm64_packages_unavailable(depth_result, "bcftools", "vcftools")
        depth_result.assert_success()
        depth_result.assert_output_exists("results/qc/individuals.idepth")

        depth_path = Path(tmpdir) / "results" / "qc" / "individuals.idepth"
        lines = [
            line.rstrip("\n").split("\t")
            for line in depth_path.read_text().splitlines()
            if line.strip()
        ]
        header = lines[0]
        nsites_index = header.index("N_SITES")
        mean_depth_index = header.index("MEAN_DEPTH")
        rows = lines[1:]

        assert rows, "Expected AD-derived depth rows"
        assert all(int(row[nsites_index]) > 0 for row in rows)
        assert all(math.isfinite(float(row[mean_depth_index])) for row in rows)


@pytest.mark.full_run
def test_qc_min_depth_removing_all_samples_fails_clearly(request):
    """QC should stop at the sample filter, naming min_depth, when no sample passes it."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="results/qc/individuals.samps.txt",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "vcf": str(TEST_DATA_DIR / "qc" / "raw.vcf.gz"),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
                "min_depth": "100000",
            },
        )
        skip_if_arm64_packages_unavailable(result, "vcftools")
        assert not result.succeeded, "QC should fail when min_depth removes every sample"

        output = result.stdout + result.stderr
        assert "min_depth filter removed all" in output
        assert "modules.qc.min_depth (100000)" in output


@pytest.mark.full_run
def test_qc_ad_only_vcf_renders_dashboard(request):
    """AD-only VCFs should still render the complete QC dashboard."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        ad_only_vcf = write_ad_only_qc_vcf(tmpdir)
        config_overrides = {
            "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
            "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
            "vcf": str(ad_only_vcf),
            "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
            "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
        }
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        dashboard_result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides=config_overrides,
        )
        skip_if_arm64_packages_unavailable(
            dashboard_result,
            "admixture",
            "bcftools",
            "plink2",
            "plink",
            "r-ggmap",
            "r-ape",
            "bioconductor-ggtree",
        )
        dashboard_result.assert_success()
        dashboard_result.assert_output_exists("results/qc/qc_dashboard.html")


@pytest.mark.full_run
def test_qc_standalone_full_run(request):
    """Full execution of QC module as standalone workflow against test fixtures."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
                "vcf": str(TEST_DATA_DIR / "qc" / "raw.vcf.gz"),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
                "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
            },
        )
        skip_if_arm64_packages_unavailable(
            result,
            "admixture",
            "plink2",
            "plink",
            "r-ggmap",
            "r-ape",
            "bioconductor-ggtree",
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc/qc_dashboard.html",
            "results/qc/contig_map.tsv",
            "results/qc/plink.eigenvec",
            "results/qc/plink.bim_fixed",
            "results/qc/plink.3.Q",
            "results/qc/coords.txt",
            "results/qc/qc_report.tsv",
        )

        contig_map = Path(tmpdir) / "results" / "qc" / "contig_map.tsv"
        contig_lines = contig_map.read_text().strip().splitlines()
        assert contig_lines[0].split("\t") == ["original_contig", "plink_contig", "admixture_id"]
        assert [line.split("\t")[2] for line in contig_lines[1:]] == ["1", "2"]

        # Copy dashboard HTML out before tmpdir cleanup so CI can upload it
        artifacts_dir = Path("test-artifacts")
        artifacts_dir.mkdir(exist_ok=True)
        dashboard = Path(tmpdir) / "results" / "qc" / "qc_dashboard.html"
        if dashboard.exists():
            shutil.copy2(dashboard, artifacts_dir / "qc_dashboard.html")


@pytest.mark.full_run
def test_qc_plink_handles_sex_organelle_contigs_and_small_cohorts(request):
    """PLINK and ADMIXTURE run on chrX/MT contigs with only 5 samples.

    plink2 needs sex information for chrX and rejects MT, so contigs are always
    renamed for PLINK. It also cannot compute 10 PCs from 5 samples, so the
    PC count follows the samples that remain (min(10, n - 1)).
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        vcf, fai, samples, _ = write_qc_cohort(tmpdir, contig_names=["chrX", "MT"], n_samples=5)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target=["results/qc/plink.eigenvec", "results/qc/plink.3.Q"],
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={"samples": str(samples), "vcf": str(vcf), "fai": str(fai)},
        )
        skip_if_arm64_packages_unavailable(result, "admixture", "plink2", "plink")
        result.assert_success()

        qc_dir = Path(tmpdir) / "results" / "qc"
        n_samples = len((qc_dir / "plink.fam").read_text().splitlines())
        header = (qc_dir / "plink.eigenvec").read_text().splitlines()[0].split()
        assert n_samples == 5
        assert [col for col in header if col.startswith("PC")] == [f"PC{i}" for i in range(1, 5)]
        # ADMIXTURE gets integer codes in its own copy; plink.bim keeps the
        # PLINK names that plink.bed and the dashboard were built with.
        bim_contigs = {line.split("\t", 1)[0] for line in (qc_dir / "plink.bim").read_text().splitlines()}
        assert bim_contigs <= {"qcctg1", "qcctg2"}, bim_contigs


@pytest.mark.full_run
def test_qc_dashboard_renders_three_samples(request):
    """Three samples, the fewest QC accepts, with the default 3 clusters.

    R's default k-means needs fewer clusters than points, so with one cluster
    per sample the dashboard failed with "number of cluster centres must lie
    between 1 and nrow(x)".
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(TEST_DATA_DIR / "qc" / "ref.fai") as handle:
            contigs = [line.split()[0] for line in handle if line.strip()]
        vcf, fai, samples, _ = write_qc_cohort(tmpdir, contig_names=contigs, n_samples=3)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(samples),
                "vcf": str(vcf),
                "fai": str(fai),
                "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
            },
        )
        skip_if_arm64_packages_unavailable(
            result, "admixture", "plink2", "plink", "r-ggmap", "r-ape", "bioconductor-ggtree"
        )
        result.assert_success()

        html = (Path(tmpdir) / "results" / "qc" / "qc_dashboard.html").read_text()
        assert "a k-means clustering with 2 clusters is applied" in html


@pytest.mark.full_run
def test_qc_dashboard_counts_pass_records_and_many_clusters(request):
    """The dashboard's numbers, text and map match the data.

    Passing records are labelled PASS, as tools such as bcftools filter and
    GATK's ApplyVQSR write them (snpArcher's own GATK hard filtering leaves
    them "."), and there are 12 clusters (more than the palette's 10 colours). The module runs from a
    copy whose scripts directory is read-only, as in a shared or read-only
    install: the dashboard must render without writing there. One sample is
    in Europe, so the map must not be scoped to North America.
    """
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        module_copy = Path(tmpdir) / "qc_module"
        shutil.copytree(WORKFLOW_DIR / "modules" / "qc", module_copy)
        scripts_dir = module_copy / "scripts"
        with open(TEST_DATA_DIR / "qc" / "ref.fai") as handle:
            contigs = [line.split()[0] for line in handle if line.strip()]
        vcf, fai, samples, n_records = write_qc_cohort(
            tmpdir, contig_names=contigs, n_samples=150, filter_label="PASS"
        )
        metadata = Path(tmpdir) / "metadata_with_europe.csv"
        metadata.write_text(
            (TEST_DATA_DIR / "qc" / "sample_metadata.csv").read_text().replace(
                "test_A01,34.05,-118.24", "test_A01,52.52,13.40"
            )
        )
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=module_copy / "Snakefile",
        )
        scripts_dir.chmod(0o555)
        try:
            result = smk.run(
                target="all",
                configfile=module_copy / "config" / "config.yaml",
                config_overrides={
                    "samples": str(samples),
                    "sample_metadata": str(metadata),
                    "vcf": str(vcf),
                    "fai": str(fai),
                    "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
                    "clusters": "12",
                    "min_depth": "1.5",
                },
            )
        finally:
            scripts_dir.chmod(0o755)
        skip_if_arm64_packages_unavailable(
            result, "admixture", "plink2", "plink", "r-ggmap", "r-ape", "bioconductor-ggtree"
        )
        result.assert_success()

        qc_dir = Path(tmpdir) / "results" / "qc"
        html = (qc_dir / "qc_dashboard.html").read_text()
        n_pruned = len((qc_dir / "plink.bim").read_text().splitlines())
        n_analyzed = len((qc_dir / "plink.fam").read_text().splitlines())
        assert f"In total, {n_records} variants were discovered" in html
        assert "NA variants" not in html
        assert f"(in this report, {n_pruned})" in html
        assert f"Of the 150 samples in the VCF, {n_analyzed} passed" in html
        assert "mean depth of 1.5x or less" in html
        assert '"scope":"world"' in html


def run_qc_site_filter(request, tmpdir, kept, dropped, rows):
    """Run only subsample_snps on a small VCF; return the positions it keeps.

    `rows` maps POS to (kept genotypes, dropped genotypes, INFO/AF). Only the
    `kept` samples are in individuals.samps.txt, as if they alone had passed
    the depth filter.
    """
    tmp_path = Path(tmpdir)
    vcf = tmp_path / "subset_raw.vcf.gz"
    with gzip.open(vcf, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=ctg1,length=10000>\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(kept + dropped) + "\n")
        for pos, (kept_gts, dropped_gts, af) in rows.items():
            handle.write(f"ctg1\t{pos}\t.\tA\tG\t50\t.\tAF={af:.3f}\tGT\t")
            handle.write("\t".join(kept_gts + dropped_gts) + "\n")
    fai = tmp_path / "subset_ref.fai"
    fai.write_text("ctg1\t10000\t6\t60\t61\n")
    samples = tmp_path / "subset_samples.csv"
    samples.write_text("sample_id\n" + "".join(f"{s}\n" for s in kept + dropped))
    samps = tmp_path / "results" / "qc" / "individuals.samps.txt"
    samps.parent.mkdir(parents=True)
    samps.write_text("".join(f"{s}\n" for s in kept))

    smk = SnakemakeRunner(
        tmp_path,
        use_conda=not request.config.getoption("--no-conda"),
        snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
    )
    result = smk.run(
        target="results/qc/snpqc.txt",
        configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
        config_overrides={"samples": str(samples), "vcf": str(vcf), "fai": str(fai)},
        extra_args=["--allowed-rules", "subsample_snps"],
    )
    skip_if_arm64_packages_unavailable(result, "bcftools")
    result.assert_success()
    snpqc = (tmp_path / "results" / "qc" / "snpqc.txt").read_text().splitlines()
    return [int(line.split("\t")[1]) for line in snpqc]


@pytest.mark.full_run
def test_qc_site_filters_use_only_kept_samples(request):
    """AF and missingness filters are evaluated over the samples QC keeps.

    4 samples pass the depth filter; 13 do not. Site 100 is monomorphic in the
    kept samples (its INFO/AF counts the others) and must go. Site 200 is
    polymorphic and fully called in the kept samples, but missing in all 13
    others (13/17 > 0.75), and must stay.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        positions = run_qc_site_filter(
            request,
            tmpdir,
            kept=[f"k{i}" for i in range(4)],
            dropped=[f"d{i}" for i in range(13)],
            rows={
                100: (["0/0"] * 4, ["1/1"] * 13, 26 / 34),
                200: (["0/1", "0/0", "0/0", "0/0"], ["./."] * 13, 1 / 8),
                300: (["0/1", "0/1", "0/0", "0/0"], ["0/0"] * 13, 2 / 34),
            },
        )
    assert positions == [200, 300]


@pytest.mark.full_run
def test_qc_maf_filter_applies_to_either_allele(request):
    """A site whose REF allele is the rare one is filtered like any other.

    With 60 kept samples, site 100 has ALT frequency 119/120 = 0.992, so its
    minor allele (REF) is at 0.008 < 0.01 and the site must go. Site 200 has
    ALT frequency 118/120 = 0.983 (minor 0.017) and stays.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        positions = run_qc_site_filter(
            request,
            tmpdir,
            kept=[f"k{i}" for i in range(60)],
            dropped=[],
            rows={
                100: (["0/1"] + ["1/1"] * 59, [], 119 / 120),
                200: (["0/1"] * 2 + ["1/1"] * 58, [], 118 / 120),
            },
        )
    assert positions == [200]


@pytest.mark.full_run
def test_qc_coords_rebuild_when_metadata_changes(request):
    """Editing sample coordinates reruns generate_coords_file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        metadata = Path(tmpdir) / "metadata.csv"
        metadata.write_text((TEST_DATA_DIR / "qc" / "sample_metadata.csv").read_text())
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=False,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        kwargs = dict(
            target="results/qc/coords.txt",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(metadata),
            },
        )
        smk.run(**kwargs).assert_success()

        time.sleep(1.1)
        metadata.write_text(metadata.read_text().replace("34.05,-118.24", "35.00,-118.00"))
        result = smk.dry_run(**kwargs)
        result.assert_success()
        assert scheduled_rule_present(result.stdout + result.stderr, "generate_coords_file")


@pytest.mark.dry_run
def test_qc_module_runs_without_a_configfile():
    """The module's documented standalone command works from any directory.

    With a `configfile: "config/config.yaml"` directive it failed unless the
    working directory happened to contain that file (and then silently read
    it, e.g. a project's main config).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            [
                "snakemake", "-n",
                "-s", str(WORKFLOW_DIR / "modules" / "qc" / "Snakefile"),
                "--directory", tmpdir,
                "--config",
                f"samples={TEST_DATA_DIR / 'qc' / 'samples.csv'}",
                f"vcf={TEST_DATA_DIR / 'qc' / 'raw.vcf.gz'}",
                f"fai={TEST_DATA_DIR / 'qc' / 'ref.fai'}",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.dry_run
def test_qc_subsample_snps_command(request):
    """AF is recomputed after subsetting, pruning is seeded, and scaffolds are quoted."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.dry_run(
            target="results/qc/pruned.vcf.gz",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "vcf": str(TEST_DATA_DIR / "qc" / "raw.vcf.gz"),
                "fai": str(TEST_DATA_DIR / "qc" / "ref.fai"),
                "exclude_scaffolds": "'mtDNA, Y'",
            },
        )
        result.assert_success()

        output = result.stdout + result.stderr
        assert "bcftools +fill-tags -Ou -- -t AF" in output
        assert "--random-seed 42" in output
        assert "-t '^mtDNA,Y'" in output


@pytest.mark.full_run
def test_qc_numeric_contigs_full_run(request):
    """QC should normalize numeric-only contig labels and still run end-to-end."""
    no_conda = request.config.getoption("--no-conda")
    with tempfile.TemporaryDirectory() as tmpdir:
        numeric_vcf, numeric_fai = write_numeric_qc_inputs(tmpdir)
        smk = SnakemakeRunner(
            Path(tmpdir),
            use_conda=not no_conda,
            snakefile=WORKFLOW_DIR / "modules" / "qc" / "Snakefile",
        )
        result = smk.run(
            target="all",
            configfile=WORKFLOW_DIR / "modules" / "qc" / "config" / "config.yaml",
            config_overrides={
                "samples": str(TEST_DATA_DIR / "qc" / "samples.csv"),
                "sample_metadata": str(TEST_DATA_DIR / "qc" / "sample_metadata.csv"),
                "vcf": str(numeric_vcf),
                "fai": str(numeric_fai),
                "qc_report": str(TEST_DATA_DIR / "qc" / "qc_report.tsv"),
            },
        )
        skip_if_arm64_packages_unavailable(
            result,
            "admixture",
            "plink2",
            "plink",
            "r-ggmap",
            "r-ape",
            "bioconductor-ggtree",
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc/contig_map.tsv",
            "results/qc/plink_input.vcf.gz",
            "results/qc/plink.bim_fixed",
            "results/qc/qc_dashboard.html",
        )

        contig_map = Path(tmpdir) / "results" / "qc" / "contig_map.tsv"
        rows = [line.split("\t") for line in contig_map.read_text().strip().splitlines()[1:]]
        assert [row[1] for row in rows] == ["qcctg1", "qcctg2"]
        assert [row[2] for row in rows] == ["1", "2"]

        bim_fixed = Path(tmpdir) / "results" / "qc" / "plink.bim_fixed"
        chroms = []
        for line in bim_fixed.read_text().strip().splitlines():
            chrom = line.split("\t", 1)[0]
            if chrom not in chroms:
                chroms.append(chrom)
        assert chroms[:2] == ["1", "2"]

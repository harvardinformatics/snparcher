import gzip
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
import yaml

from conftest import SnakemakeRunner, FIXTURES_DIR

SAMPLES = FIXTURES_DIR / "config" / "samples.csv"


@pytest.mark.unit
def test_setup(request):
    """Setup from scratch: reference prep + indexing + intervals."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data")

        result = smk.run(target="setup", samples=SAMPLES)
        result.assert_success()
        result.assert_output_exists(
            "results/reference/my_organism.fa.gz",
            "results/reference/my_organism.fa.gz.fai",
            "results/reference/my_organism.dict",
            "results/intervals/gvcf/intervals.txt",
            "results/intervals/db/intervals.txt",
        )


@pytest.mark.unit
def test_fastp(request):
    """Just fastp for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data", "results/reference", "results/intervals")

        result = smk.run(
            target=[
                "results/filtered_fastqs/sample0/sample0/u1_1.fastq.gz",
                "results/filtered_fastqs/sample0/sample0/u1_2.fastq.gz",
            ],
            samples=SAMPLES,
        )
        result.print_log()
        result.assert_success()
        result.assert_output_exists(
            "results/filtered_fastqs/sample0/sample0/u1_1.fastq.gz",
        )


@pytest.mark.unit
def test_bwa_mem(request):
    """Just bwa mem for one sample — fastp output provided."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/filtered_fastqs/sample0",
        )

        result = smk.run(
            target=[
                "results/bams/raw/sample0/sample0/u1.bam",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/raw/sample0/sample0/u1.bam",
        )


@pytest.mark.unit
def test_markdup(request):
    """Just merge + markdup for one sample — raw BAMs provided."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
        )
        # The raw BAMs are temp() inputs that Snakemake deletes after use.
        smk.copy_fixtures("results/bams/raw/sample0")

        result = smk.run(
            target=[
                "results/bams/markdup/sample0.bam",
                "results/bams/markdup/sample0.bam.csi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/markdup/sample0.bam",
            "results/bams/markdup/sample0.bam.csi",
        )
        assert (FIXTURES_DIR / "results/bams/raw/sample0/sample0/u1.bam").exists()


@pytest.mark.unit
def test_external_bam_staging_and_indexing(request):
    """External BAM inputs are staged and indexed with workflow-managed CSI."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        samples = tmpdir / "external_bam_samples.csv"
        external_bam = FIXTURES_DIR / "results/bams/markdup/sample0.bam"
        samples.write_text(
            "sample_id,input_type,input\n"
            f"sample_ext,bam,{external_bam}\n"
        )

        smk = SnakemakeRunner(tmpdir, use_conda=not no_conda, conda_prefix=conda_prefix)

        result = smk.run(
            target="results/bams/input/sample_ext.bam.csi",
            samples=samples,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/bams/input/sample_ext.bam",
            "results/bams/input/sample_ext.bam.csi",
        )


@pytest.mark.unit
def test_gvcfs(request):
    """Per-sample gVCF generation for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target=[
                "results/gvcfs/sample0.g.vcf.gz",
                "results/gvcfs/sample0.g.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/gvcfs/sample0.g.vcf.gz",
            "results/gvcfs/sample0.g.vcf.gz.tbi",
        )


@pytest.mark.unit
def test_qc_metrics(request):
    """QC metrics for one sample + combine."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
            "results/fastp",
        )

        result = smk.run(
            target="results/qc_metrics/bam/sample0.json",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc_metrics/bam/sample0.json",
        )


@pytest.mark.unit
def test_callable_sites_mappability(request):
    """Mappability only — no per-sample work needed."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
        )

        result = smk.run(
            target="results/callable_sites/mappability.bed",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/mappability.bed",
        )


@pytest.mark.unit
def test_callable_sites_coverage(request):
    """Coverage for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target="results/callable_sites/depths/sample0.per-base.d4",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/depths/sample0.per-base.d4",
        )


@pytest.mark.unit
def test_haplotypecaller(request):
    """HaplotypeCaller for one sample, one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target=[
                "results/interval_gvcfs/sample0/0000.g.vcf.gz",
                "results/interval_gvcfs/sample0/0000.g.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/interval_gvcfs/sample0/0000.g.vcf.gz",
        )


@pytest.mark.unit
def test_genomics_db_import(request):
    """GenomicsDB import for one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
        )

        result = smk.run(
            target="results/gatk_genomics_db/L0000.tar",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/gatk_genomics_db/L0000.tar",
        )


@pytest.mark.unit
def test_genotype_gvcfs(request):
    """GenotypeGVCFs for one interval."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/gvcfs",
            "results/genomics_db",
            "results/gatk_genomics_db",
        )

        result = smk.run(
            target=[
                "results/vcfs/intervals/L0000.vcf.gz",
                "results/vcfs/intervals/L0000.vcf.gz.tbi",
            ],
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/vcfs/intervals/L0000.vcf.gz",
        )


@pytest.mark.unit
def test_collect_fastp_stats(request):
    """Collect fastp stats for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/fastp",
        )

        result = smk.run(
            target="results/qc_metrics/fastp/sample0.json",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/qc_metrics/fastp/sample0.json",
        )


@pytest.mark.unit
def test_clam_collect_keeps_sample_ids_with_dots(request):
    """depths.zarr names its columns by sample_id, dots included.

    clam names samples by file name up to the first dot unless it is given
    the names, so ind.1 and ind.2 both became "ind" and clam_collect failed.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )
        result = smk.run(
            target=[
                f"results/callable_sites/depths/{s}.per-base.d4" for s in ("sample0", "sample1")
            ],
            samples=SAMPLES,
        )
        result.assert_success()

        # Rename the samples, and their depth files, to IDs with dots.
        depths_dir = Path(tmpdir) / "results" / "callable_sites" / "depths"
        dotted = {"sample0": "ind.1", "sample1": "ind.2"}
        for old, new in dotted.items():
            (depths_dir / f"{old}.per-base.d4").rename(depths_dir / f"{new}.per-base.d4")
        sheet_lines = SAMPLES.read_text().splitlines()
        samples = Path(tmpdir) / "dotted_samples.csv"
        samples.write_text("\n".join(
            [sheet_lines[0]]
            + [
                line.replace(old, new, 1)
                for line in sheet_lines[1:]
                for old, new in dotted.items()
                if line.startswith(f"{old},")
            ]
        ) + "\n")

        result = smk.run(
            target="results/callable_sites/depths.zarr",
            samples=samples,
            extra_args=["--allowed-rules", "clam_collect"],
        )
        result.print_log()
        result.assert_success()

        metadata = json.loads(
            (Path(tmpdir) / "results/callable_sites/depths.zarr/zarr.json").read_text()
        )["attributes"]["clam_metadata"]
        assert metadata["column_names"] == ["ind.1", "ind.2"]


@pytest.mark.unit
def test_mosdepth(request):
    """Mosdepth for one sample."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/bams/markdup",
        )

        result = smk.run(
            target="results/callable_sites/depths/sample0.mosdepth.summary.txt",
            samples=SAMPLES,
        )
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/depths/sample0.mosdepth.summary.txt",
        )
@pytest.mark.unit
def test_coverage_bed(request):
    """Coverage BED generation from callable loci zarr."""
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures(
            "config",
            "data",
            "results/reference",
            "results/intervals",
            "results/callable_sites/callable_loci.zarr",
        )

        result = smk.run(
            target="results/callable_sites/coverage.bed",
            samples=SAMPLES,
        )
        result.print_log()
        result.assert_success()
        result.assert_output_exists(
            "results/callable_sites/coverage.bed",
        )

        # Verify BED file is non-empty and well-formed
        bed_path = Path(tmpdir) / "results" / "callable_sites" / "coverage.bed"
        lines = bed_path.read_text().strip().splitlines()
        assert len(lines) > 0, "coverage.bed should not be empty"
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 3, f"BED line should have 3 fields: {line}"
            assert fields[0] == "chr2l", f"Expected contig chr2l: {line}"
            assert int(fields[1]) < int(fields[2]), f"Start should be < end: {line}"


@pytest.mark.unit
def test_callable_sites_bed_is_intersection(request):
    """The final callable BED keeps only sites passing coverage AND mappability.

    The sources are written directly, and --allowed-rules stops Snakemake from
    rebuilding them, so only callable_sites_bed runs. A union of the two
    sources would yield chr2l:0-300 rather than the two overlapping pieces.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data", "results/reference")

        callable_dir = Path(tmpdir) / "results" / "callable_sites"
        callable_dir.mkdir(parents=True)
        # Coverage BED is written in reference order, not sorted.
        (callable_dir / "coverage.bed").write_text("chr2l\t200\t300\nchr2l\t0\t100\n")
        (callable_dir / "mappability.bed").write_text("chr2l\t50\t250\n")

        result = smk.run(
            target="results/callable_sites/callable_sites.bed",
            samples=SAMPLES,
            extra_args=["--allowed-rules", "callable_sites_bed"],
        )
        result.print_log()
        result.assert_success()

        bed = (callable_dir / "callable_sites.bed").read_text()
        assert bed == "chr2l\t50\t100\nchr2l\t200\t250\n", bed


def write_clam_depth_zarr(zarr_dir, contig, depths, chunk_rows=1000):
    """Write a depth store laid out the way `clam collect` writes one.

    `depths` is a (positions, samples) array. clam compresses its chunks; this
    writes them uncompressed (bytes codec only), which clam reads the same way.
    """
    length, n_samples = depths.shape
    samples = [f"s{i}" for i in range(n_samples)]
    zarr_dir.mkdir(parents=True)
    (zarr_dir / "zarr.json").write_text(json.dumps({
        "zarr_format": 3,
        "node_type": "group",
        "attributes": {
            "clam_metadata": {
                "contigs": [{"name": contig, "length": length}],
                "column_names": samples,
                "chunk_size": chunk_rows,
                "populations": [{"name": "default", "samples": samples}],
            }
        },
    }))

    array_dir = zarr_dir / contig
    array_dir.mkdir()
    (array_dir / "zarr.json").write_text(json.dumps({
        "zarr_format": 3,
        "node_type": "array",
        "shape": [length, n_samples],
        "data_type": "uint32",
        "chunk_grid": {
            "name": "regular",
            "configuration": {"chunk_shape": [chunk_rows, n_samples]},
        },
        "chunk_key_encoding": {"name": "default", "configuration": {"separator": "/"}},
        "fill_value": 0,
        "codecs": [{"name": "bytes", "configuration": {"endian": "little"}}],
    }))

    for index, start in enumerate(range(0, length, chunk_rows)):
        chunk = np.zeros((chunk_rows, n_samples), dtype="<u4")
        block = depths[start : start + chunk_rows]
        chunk[: len(block)] = block
        chunk_path = array_dir / "c" / str(index) / "0"
        chunk_path.parent.mkdir(parents=True)
        chunk_path.write_bytes(chunk.tobytes())


def run_coverage_mask(request, tmpdir, depths, mappability_bed, **coverage):
    """Run thresholds -> clam loci -> coverage.bed on a hand-written depth store.

    Returns the BED lines and the thresholds TSV as a dict.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
    smk.link_fixtures("config", "data")

    callable_dir = Path(tmpdir) / "results" / "callable_sites"
    write_clam_depth_zarr(callable_dir / "depths.zarr", "chr2l", depths)
    if mappability_bed is not None:
        (callable_dir / "mappability.bed").write_text(mappability_bed)

    config = yaml.safe_load((FIXTURES_DIR / "config" / "config.yaml").read_text())
    config["callable_sites"]["coverage"].update(merge_distance=0, **coverage)
    config["callable_sites"]["mappability"]["enabled"] = mappability_bed is not None
    configfile = Path(tmpdir) / "coverage_mask_config.yaml"
    configfile.write_text(yaml.safe_dump(config))

    # --allowed-rules keeps Snakemake from rebuilding the inputs with mosdepth,
    # clam collect and genmap.
    result = smk.run(
        target="results/callable_sites/coverage.bed",
        configfile=configfile,
        samples=SAMPLES,
        extra_args=[
            "--allowed-rules", "callable_coverage_thresholds", "clam_loci", "coverage_bed",
        ],
    )
    result.print_log()
    result.assert_success()

    thresholds = dict(
        line.split("\t")
        for line in (callable_dir / "coverage_thresholds.tsv").read_text().splitlines()
    )
    return (callable_dir / "coverage.bed").read_text().splitlines(), thresholds


def coverage_mask_depths():
    """5 samples on one contig, built so each region's fate follows by hand.

    Single-copy sites sum to exactly 50 (per-sample depths 8-12). Over the
    mappable sites with data the median is 50 and the robust spread is zero,
    so the spread floor (0.5 on the square-root scale) sets the auto cutoff:
    (sqrt(50) + 4 * 0.5) ** 2 = 82.3, i.e. summed depth <= 82.
    """
    length, n_samples = 5000, 5
    positions = np.arange(length)[:, None]
    depths = (10 + (positions + np.arange(n_samples)) % 5 - 2).astype(np.uint32)
    depths[500:600] = 20                    # collapsed duplication: summed 100
    depths[700:800] = 16                    # summed 80: within auto, above 1.5x
    depths[900:1000] = 0
    depths[900:1000, 0] = 10                # data in 1 of 5 samples: below 0.25
    depths[1000:1100] = 0
    depths[1000:1100, :2] = 25              # data in 2 of 5 samples: passes
    depths[1100:1200] = 0                   # no data anywhere
    # Unmappable and very deep, and 60% of all sites: if it entered the depth
    # distribution, the median would be 1000 and nothing would be excluded.
    depths[2000:] = 200
    return depths


@pytest.mark.unit
@pytest.mark.parametrize(
    "max_depth, expected",
    [
        ("auto", ["chr2l\t0\t500", "chr2l\t600\t900", "chr2l\t1000\t1100", "chr2l\t1200\t2000"]),
        # 1.5 x 50 = 75: summed 80 is now excluded too.
        (1.5, ["chr2l\t0\t500", "chr2l\t600\t700", "chr2l\t800\t900",
               "chr2l\t1000\t1100", "chr2l\t1200\t2000"]),
        # 2 x 50 = 100: the boundary is inclusive, so summed 100 is kept.
        (2, ["chr2l\t0\t900", "chr2l\t1000\t1100", "chr2l\t1200\t2000"]),
    ],
)
def test_coverage_mask(request, max_depth, expected):
    """Summed-depth upper bound, presence fraction, and mappable-only anchor.

    Runs the real threshold script and clam loci. The unmappable region
    (2000-5000) lies outside mappability.bed and fails the depth bound
    anyway; it is in the test to show that it does not move the median.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        lines, thresholds = run_coverage_mask(
            request,
            tmpdir,
            coverage_mask_depths(),
            mappability_bed="chr2l\t0\t2000\n",
            max_depth=max_depth,
        )

    assert lines == expected
    assert float(thresholds["median_summed_depth"]) == 50
    assert int(thresholds["min_samples_with_data"]) == 2
    if max_depth == "auto":
        assert int(thresholds["summed_depth_cutoff"]) == 82


@pytest.mark.unit
def test_coverage_mask_fraction_boundary_is_exact(request):
    """fraction 0.55 of 200 samples keeps a site where exactly 110 have data.

    0.55 * 200 is 110.00000000000001 in floating point, and clam's -d drops
    a site at exactly that proportion. The thresholds step computes the count
    exactly and hands clam a proportion half a sample lower.
    """
    depths = np.ones((400, 200), dtype=np.uint32)
    for start, with_data in ((300, 109), (310, 110), (320, 111)):
        depths[start : start + 10, with_data:] = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        lines, thresholds = run_coverage_mask(
            request, tmpdir, depths, mappability_bed=None, fraction=0.55
        )

    assert int(thresholds["min_samples_with_data"]) == 110
    assert lines == ["chr2l\t0\t300", "chr2l\t310\t400"]


@pytest.mark.unit
def test_coverage_mask_max_depth_boundary_is_exact(request):
    """max_depth 2.3 of a median of 50 keeps a site whose summed depth is 115.

    2.3 * 50 is 114.99999999999999 in floating point, which floors to 114.
    """
    depths = np.full((400, 5), 10, dtype=np.uint32)  # summed 50
    depths[300:310] = 23  # summed 115: on the bound
    depths[310:320] = (23, 23, 23, 23, 24)  # summed 116: above it

    with tempfile.TemporaryDirectory() as tmpdir:
        lines, thresholds = run_coverage_mask(
            request, tmpdir, depths, mappability_bed=None, max_depth=2.3
        )

    assert int(thresholds["summed_depth_cutoff"]) == 115
    assert lines == ["chr2l\t0\t310", "chr2l\t320\t400"]


@pytest.mark.unit
def test_genmap_index_leaves_reference_fasta_alone(request):
    """genmap_index must not touch results/reference/{name}.fa.

    It used to decompress the reference to exactly that path and delete it on
    exit, destroying any FASTA already there.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data", "results/reference/my_organism.fa.gz")

        fasta = Path(tmpdir) / "results" / "reference" / "my_organism.fa"
        fasta.write_text(">user_owned\nACGT\n")

        result = smk.run(
            target="results/callable_sites/genmap_index",
            samples=SAMPLES,
            extra_args=["--allowed-rules", "genmap_index"],
        )
        result.print_log()
        result.assert_success()

        assert fasta.read_text() == ">user_owned\nACGT\n"
        assert any((Path(tmpdir) / "results" / "callable_sites" / "genmap_index").iterdir())


@pytest.mark.unit
def test_download_sra_completeness(request):
    """download_sra emits every read the archive holds.

    Exercises the completeness assertion in the rule: the extractor's output is
    counted and compared against the archive's own spot count, rather than the
    rule merely checking that the files exist. A truncated extractor, or a pigz
    that compresses nothing and exits 0, fails here.

    The output directory starts with a stale single-read file, as a retry finds
    it: the raw FASTQs are temp(), but that file was never declared. The rule
    must neither trip over it nor count its reads.

    Requires network access: prefetch pulls SRR000001 (~312 MB) from NCBI.
    """
    no_conda = request.config.getoption("--no-conda")
    conda_prefix = request.config.getoption("--conda-prefix")
    samples = Path(__file__).parent / "sample_sheets" / "local_fastqs_and_srr_same_library.csv"

    with tempfile.TemporaryDirectory() as tmpdir:
        smk = SnakemakeRunner(Path(tmpdir), use_conda=not no_conda, conda_prefix=conda_prefix)
        smk.link_fixtures("config", "data")

        r1 = "results/fastqs/sample_mixed_reads/libA/u2/SRR000001_1.fastq.gz"
        r2 = "results/fastqs/sample_mixed_reads/libA/u2/SRR000001_2.fastq.gz"
        leftover = Path(tmpdir) / "results/fastqs/sample_mixed_reads/libA/u2/SRR000001.fastq.gz"
        leftover.parent.mkdir(parents=True)
        with gzip.open(leftover, "wt") as handle:
            handle.write("@stale\nACGT\n+\nIIII\n" * 3)

        result = smk.run(target=[r1, r2], samples=samples)
        result.print_log()
        result.assert_success()
        result.assert_output_exists(r1, r2)

        # The rule logs what it counted; assert it actually ran the check
        # rather than passing because the block was skipped.
        log = (
            Path(tmpdir)
            / "logs"
            / "download_sra"
            / "sample_mixed_reads"
            / "libA"
            / "u2"
            / "SRR000001.txt"
        ).read_text()
        assert "completeness:" in log, f"completeness check did not run:\n{log}"

        emitted, expected = re.search(
            r"completeness: emitted=(\d+) expected=(\d+)", log
        ).groups()
        assert emitted == expected, f"emitted {emitted} != expected {expected}"


@pytest.mark.unit
def test_csra_expected_reads():
    """The cSRA census predicts what sam-dump | collate | fastq emits.

    Each spot is one vdb-dump tab line: READ_LEN, PRIMARY_ALIGNMENT_ID and
    ALIGNMENT_COUNT, one comma-separated value per read. The layouts are the
    ones found in ficAlb's ERR2167036 and ERR2167049, where the prediction
    matched the extracted counts exactly.
    """
    script = Path(__file__).parents[1] / "workflow" / "scripts" / "csra_expected_reads.awk"
    spots = [
        "151, 151\t12, 9\t1, 3",  # both primary: a pair
        "151, 151\t5, 0\t1, 0",  # half-aligned: the unaligned mate is emitted
        "151, 151\t0, 0\t0, 0",  # fully unaligned: a pair
        "151, 151\t0, 0\t0, 2",  # secondary-only beside unaligned: still visited
        "151, 151\t7, 0\t1, 2",  # secondary-only mate of a primary: a singleton
        "151, 151\t0, 0\t2, 1",  # both secondary-only: nothing emitted
        "151, 0\t8, 0\t1, 0",  # mate-less BAM record: unpaired
        "0, 151\t0, 0\t0, 0",  # mate-less and unaligned: unpaired
        "151\t3\t1",  # one-read spot: not modelled
    ]
    result = subprocess.run(
        ["awk", "-f", str(script)],
        input="\n".join(spots) + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    census = dict(line.split("=") for line in result.stdout.split())
    assert census == {
        "spots": "9",
        "spots_unmodelled": "1",
        "lines_malformed": "0",
        "reads_primary": "5",
        "reads_zero_length": "2",
        "reads_lost_secondary_only": "3",
        "pairs": "4",
        "singletons": "1",
        "unpaired": "2",
    }

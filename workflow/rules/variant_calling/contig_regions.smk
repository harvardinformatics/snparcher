# Per-contig calling regions, shared by the bcftools and repadapt callers.
# Region IDs are L{idx:06d} in .fai order, so sorted IDs give .fai order.

from pathlib import Path


def bcftools_call_input(wc):
    bams = [get_final_bam(s) for s in SAMPLES_WITH_BAM]
    return {
        "bams": bams,
        "bam_indexes": [get_bam_index(bam) for bam in bams],
        **REF_FILES,
    }


checkpoint bcftools_regions:
    input:
        ref_fai=REF_FILES["ref_fai"],
    output:
        tsv="results/vcfs/regions/regions.tsv",
    run:
        Path("results/vcfs/regions").mkdir(parents=True, exist_ok=True)
        with open(input.ref_fai) as fin, open(output.tsv, "w") as fout:
            for idx, line in enumerate(fin):
                contig = line.split("\t", 1)[0].strip()
                if contig:
                    fout.write(f"L{idx:06d}\t{contig}\n")


def _read_bcftools_regions(regions_tsv):
    regions = {}
    with open(regions_tsv) as f:
        for line in f:
            if not line.strip():
                continue
            region_id, contig = line.rstrip("\n").split("\t", 1)
            regions[region_id] = contig
    return regions


def _get_bcftools_regions_file(wc):
    regions_file = "results/vcfs/regions/regions.tsv"
    if exists(regions_file):
        return regions_file
    return checkpoints.bcftools_regions.get(**wc).output.tsv


def get_bcftools_region_ids(wc):
    regions = _read_bcftools_regions(_get_bcftools_regions_file(wc))
    return sorted(regions.keys())


def get_bcftools_region_name(region_id):
    regions_file = "results/vcfs/regions/regions.tsv"
    if not exists(regions_file):
        raise ValueError(
            "Region map not available yet. "
            "Run through checkpoint bcftools_regions first."
        )
    regions = _read_bcftools_regions(regions_file)
    if region_id not in regions:
        raise ValueError(f"Unknown region id: {region_id}")
    return regions[region_id]

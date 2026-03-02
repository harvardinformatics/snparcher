from pathlib import Path


def bcftools_call_input(wc):
    bams = [get_final_bam(s) for s in SAMPLES_WITH_BAM]
    return {
        "bams": bams,
        "bais": [f"{bam}.bai" for bam in bams],
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
            "Region map not available yet for bcftools caller. "
            "Run through checkpoint bcftools_regions first."
        )
    regions = _read_bcftools_regions(regions_file)
    if region_id not in regions:
        raise ValueError(f"Unknown bcftools region id: {region_id}")
    return regions[region_id]


def get_bcftools_region_vcfs(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand("results/vcfs/regions/{region_id}.vcf.gz", region_id=region_ids)


def get_bcftools_region_vcf_tbis(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand("results/vcfs/regions/{region_id}.vcf.gz.tbi", region_id=region_ids)


rule bcftools_call:
    input:
        unpack(bcftools_call_input),
        regions_tsv="results/vcfs/regions/regions.tsv",
    output:
        vcf=temp("results/vcfs/regions/{region_id}.vcf.gz"),
        tbi=temp("results/vcfs/regions/{region_id}.vcf.gz.tbi"),
    params:
        min_mapq=config["variant_calling"]["bcftools"]["min_mapq"],
        min_baseq=config["variant_calling"]["bcftools"]["min_baseq"],
        max_depth=config["variant_calling"]["bcftools"]["max_depth"],
        ploidy=config["variant_calling"]["ploidy"],
        contig=lambda wc: get_bcftools_region_name(wc.region_id),
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_call/{region_id}.txt"
    log:
        "logs/bcftools_call/{region_id}.txt"
    shell:
        """
        bcftools mpileup \
            -f {input.ref} \
            -q {params.min_mapq} \
            -Q {params.min_baseq} \
            -d {params.max_depth} \
            -r {params.contig} \
            -Ou {input.bams} 2> {log} \
        | bcftools call \
            -m \
            --ploidy {params.ploidy} \
            -v \
            -Oz \
            -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule bcftools_concat_regions:
    input:
        vcfs=get_bcftools_region_vcfs,
        tbis=get_bcftools_region_vcf_tbis,
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_concat_regions.txt"
    log:
        "logs/bcftools_concat_regions.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir} -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

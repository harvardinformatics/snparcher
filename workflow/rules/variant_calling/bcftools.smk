BCFTOOLS_REGIONS_TSV = "results/vcfs/regions/regions.tsv"


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
        tsv=BCFTOOLS_REGIONS_TSV,
    run:
        with open(input.ref_fai) as fin, open(output.tsv, "w") as fout:
            for idx, line in enumerate(fin):
                contig = line.split("\t", 1)[0].strip()
                if contig:
                    fout.write(f"L{idx:06d}\t{contig}\n")


@cache
def _read_bcftools_regions(regions_tsv, mtime_ns):
    """Region ID -> contig. Each bcftools_call job asks for its own contig, and
    the file has a line per contig, so it is parsed once per version."""
    regions = {}
    with open(regions_tsv) as f:
        for line in f:
            if not line.strip():
                continue
            region_id, contig = line.rstrip("\n").split("\t", 1)
            regions[region_id] = contig
    return regions


def _bcftools_regions(regions_tsv):
    return _read_bcftools_regions(regions_tsv, os.stat(regions_tsv).st_mtime_ns)


def get_bcftools_region_ids(wc):
    if exists(BCFTOOLS_REGIONS_TSV):
        regions_tsv = BCFTOOLS_REGIONS_TSV
    else:
        regions_tsv = checkpoints.bcftools_regions.get(**wc).output.tsv
    return sorted(_bcftools_regions(regions_tsv))


def get_bcftools_region_name(region_id):
    if not exists(BCFTOOLS_REGIONS_TSV):
        raise ValueError(
            "Region map not available yet for bcftools caller. "
            "Run through checkpoint bcftools_regions first."
        )
    regions = _bcftools_regions(BCFTOOLS_REGIONS_TSV)
    if region_id not in regions:
        raise ValueError(f"Unknown bcftools region id: {region_id}")
    return regions[region_id]


def get_bcftools_region_vcfs(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand("results/vcfs/regions/{region_id}.vcf.gz", region_id=region_ids)


def get_bcftools_region_vcf_tbis(wc):
    region_ids = get_bcftools_region_ids(wc)
    return expand(
        get_compressed_vcf_index("results/vcfs/regions/{region_id}.vcf.gz"),
        region_id=region_ids,
    )


rule bcftools_call:
    input:
        unpack(bcftools_call_input),
        regions_tsv=BCFTOOLS_REGIONS_TSV,
    output:
        vcf=temp("results/vcfs/regions/{region_id}.vcf.gz"),
        idx=temp(get_compressed_vcf_index("results/vcfs/regions/{region_id}.vcf.gz")),
    params:
        min_mapq=config["variant_calling"]["bcftools"]["min_mapq"],
        min_baseq=config["variant_calling"]["bcftools"]["min_baseq"],
        max_depth=config["variant_calling"]["bcftools"]["max_depth"],
        ploidy=config["variant_calling"]["ploidy"],
        contig=lambda wc: get_bcftools_region_name(wc.region_id),
        # "<bam> <sample_id>" pairs for mpileup -G.
        bam_samples=" ".join(f"{get_final_bam(s)} {s}" for s in SAMPLES_WITH_BAM),
        index_args=BCFTOOLS_INDEX_ARGS,
    threads: 1
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_call/{region_id}.txt"
    log:
        "logs/bcftools_call/{region_id}.txt"
    shell:
        """
        # mpileup names samples by the BAMs' read-group SM, and merges BAMs
        # that share one. -G assigns every read in each BAM to its sample_id.
        READ_GROUPS=$(mktemp "{resources.tmpdir}/bcftools_read_groups_XXXXXX")
        trap 'rm -f "$READ_GROUPS"' EXIT
        printf '*\t%s\t%s\n' {params.bam_samples} > "$READ_GROUPS"
        bcftools mpileup \
            -f {input.ref} \
            -G "$READ_GROUPS" \
            -q {params.min_mapq} \
            -Q {params.min_baseq} \
            -d {params.max_depth} \
            -r {params.contig:q} \
            --threads {threads} \
            -Ou {input.bams} 2> {log} \
        | bcftools call \
            -m \
            --ploidy {params.ploidy} \
            --threads {threads} \
            -v \
            -Oz \
            -o {output.vcf} - 2>> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """


rule bcftools_concat_regions:
    input:
        vcfs=get_bcftools_region_vcfs,
        tbis=get_bcftools_region_vcf_tbis,
    output:
        vcf=raw_vcf_output(RAW_VCF),
        idx=raw_vcf_output(RAW_VCF_INDEX),
    params:
        index_args=BCFTOOLS_INDEX_ARGS,
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/bcftools_concat_regions.txt"
    log:
        "logs/bcftools_concat_regions.txt"
    shell:
        """
        # Each region VCF holds one contig, so nothing overlaps and concat can
        # read the files one at a time; -a would open them all at once, which
        # exceeds the open-file limit on fragmented assemblies.
        LIST=$(mktemp "{resources.tmpdir}/bcftools_concat_regions_XXXXXX")
        trap 'rm -f "$LIST"' EXIT
        printf '%s\\n' {input.vcfs} > "$LIST"
        bcftools concat -f "$LIST" -Ou 2> {log} \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} - 2>> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """

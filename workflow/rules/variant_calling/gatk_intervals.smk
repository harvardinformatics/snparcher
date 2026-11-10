wildcard_constraints:
    sample="[^/]+",


INTERVAL_GVCF_PATTERN = (
    "results/interval_gvcfs/{sample}/{interval}.g.vcf"
    if LONG_CONTIG_MODE
    else "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz"
)


STAGED_GVCF_PATTERN = (
    "results/gvcfs/work/staged/{sample}/r{round}/c{chunk}.g.vcf"
    if LONG_CONTIG_MODE
    else "results/gvcfs/staged/{sample}/r{round}/c{chunk}.g.vcf.gz"
)


def haplotype_caller_input(wildcards):
    if sample_has_input_type(wildcards.sample, "gvcf"):
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller")
    
    return {
        **get_indexed_final_bam_input(wildcards.sample),
        "interval": f"results/intervals/gvcf/{wildcards.interval}-scattered.interval_list",
        **REF_FILES,
    }


def get_interval_gvcfs(wc):
    """Get interval gVCF files for a sample."""
    intervals = interval_list_ids("gvcf", wc)
    return expand(INTERVAL_GVCF_PATTERN, sample=wc.sample, interval=intervals)


def staged_gvcf_path(wc, round_idx, chunk_idx):
    return STAGED_GVCF_PATTERN.format(
        sample=wc.sample,
        round=round_idx,
        chunk=chunk_idx,
    )


def get_interval_gvcf_stage_inputs(wc):
    return get_stage_inputs(
        base_files=get_interval_gvcfs(wc),
        round_idx=int(wc.round),
        chunk_idx=int(wc.chunk),
        wc=wc,
        path_builder=staged_gvcf_path,
    )


def get_interval_gvcf_stage_tbis(wc):
    stage_inputs = get_interval_gvcf_stage_inputs(wc)
    return [get_vcf_index(gvcf) for gvcf in stage_inputs]


def get_interval_gvcf_stage_picard_inputs(wc):
    return format_picard_vcf_inputs(get_interval_gvcf_stage_inputs(wc))


def get_final_interval_gvcf_stage_file(wc):
    return get_final_stage_file(
        base_files=get_interval_gvcfs(wc),
        wc=wc,
        path_builder=staged_gvcf_path,
    )


def get_final_interval_gvcf_stage_tbi(wc):
    return get_vcf_index(get_final_interval_gvcf_stage_file(wc))


rule gatk_haplotypecaller_interval:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf=temp(INTERVAL_GVCF_PATTERN),
        idx=temp(get_vcf_index(INTERVAL_GVCF_PATTERN)),
    params:
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=HC_MIN_PRUNING,
        min_dangling=HC_MIN_DANGLING,
    threads: 1
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_haplotypecaller/{sample}/{interval}.benchmark.txt"
    log:
        "logs/gatk_haplotypecaller/{sample}/{interval}.log.txt",
    shell:
        """
        gatk HaplotypeCaller \
        --java-options '-Xmx{resources.mem_mb_reduced}m' \
        -R {input.ref} \
        -I {input.bam} \
        --read-index {input.bam_index} \
        -O {output.gvcf} \
        -L {input.interval} \
        -ploidy {params.ploidy} \
        --native-pair-hmm-threads {threads} \
        --emit-ref-confidence GVCF \
        --min-pruning {params.min_pruning} \
        --min-dangling-branch-length {params.min_dangling} &> {log}
        """


rule concat_interval_gvcfs_stage:
    input:
        gvcfs=get_interval_gvcf_stage_inputs,
        tbis=get_interval_gvcf_stage_tbis,
    output:
        gvcf=temp(STAGED_GVCF_PATTERN),
        idx=temp(get_vcf_index(STAGED_GVCF_PATTERN)),
    params:
        index_args=BCFTOOLS_INDEX_ARGS,
        picard_inputs=get_interval_gvcf_stage_picard_inputs,
        sortvcf_record_bytes=sortvcf_record_bytes(1),
        long_mode=LONG_CONTIG_MODE,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/concat_interval_gvcfs/staged/{sample}/r{round}/c{chunk}.txt"
    log:
        "logs/concat_interval_gvcfs/staged/{sample}/r{round}/c{chunk}.txt"
    shell:
        """
        if [ "{params.long_mode}" = "True" ]; then
            # Keep half the heap's worth of records in RAM, at most Picard's
            # default of 500,000 (see SORTVCF_RECORD_BYTES).
            max_records=$(( {resources.mem_mb_reduced} * 1048576 / 2 / {params.sortvcf_record_bytes} ))
            max_records=$(( max_records > 500000 ? 500000 : (max_records < 1 ? 1 : max_records) ))
            picard SortVcf \
                -Xmx{resources.mem_mb_reduced}m \
                {params.picard_inputs} \
                O={output.gvcf} \
                TMP_DIR={resources.tmpdir} \
                MAX_RECORDS_IN_RAM=$max_records \
                CREATE_INDEX=false \
                > {log} 2>&1
            gatk IndexFeatureFile -I {output.gvcf} >> {log} 2>&1
        else
            bcftools concat -D -a -Ou {input.gvcfs} 2> {log} \
                | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.gvcf} - 2>> {log}
            bcftools index {params.index_args} {output.gvcf} 2>> {log}
        fi
        """


rule concat_interval_gvcfs:
    input:
        gvcf=get_final_interval_gvcf_stage_file,
        tbi=get_final_interval_gvcf_stage_tbi,
    output:
        **GATK_SAMPLE_GVCF,
    benchmark:
        "benchmarks/concat_interval_gvcfs/{sample}.txt"
    log:
        "logs/concat_interval_gvcfs/{sample}.txt"
    shell:
        """
        mv {input.gvcf} {output.gvcf} 2> {log}
        mv {input.tbi} {output.idx} 2>> {log}
        touch {output.gvcf} {output.idx}
        """

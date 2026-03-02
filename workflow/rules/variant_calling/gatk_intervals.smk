localrules: create_db_mapfile

wildcard_constraints:
    sample="[^/]+",
    round=r"\d+",
    chunk=r"\d+"

def haplotype_caller_input(wildcards):
    input_type = get_sample_input_type(wildcards.sample)
    
    if input_type == "gvcf":
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotype_caller")
    
    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        "interval": f"results/intervals/gvcf/{wildcards.interval}-scattered.interval_list",
        **REF_FILES,
    }


def get_interval_gvcfs(wc):
    """Get interval gVCF files for a sample."""
    intervals_file = "results/intervals/gvcf/intervals.txt"
    
    if exists(intervals_file):
        # Read directly from file if it exists
        # We do this because checkpoint.get() doesn't seem to pick up the checkpoint output if it was created in previous run
        # I.e if user does setup target rule, then rule all, the workflow fails 
        # Related GH issue: https://github.com/snakemake/snakemake/issues/3879
        with open(intervals_file) as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        # Fall back to checkpoint mechanism to trigger creation
        checkpoint_output = checkpoints.create_gvcf_intervals.get().output[0]
        with open(checkpoint_output) as f:
            lines = [l.strip() for l in f.readlines()]
    
    list_files = [os.path.basename(x) for x in lines]
    intervals = [f.replace("-scattered.interval_list", "") for f in list_files]
    return expand(
        "results/interval_gvcfs/{sample}/{interval}.g.vcf.gz",
        sample=wc.sample,
        interval=intervals,
    )


def get_db_intervals(wc):
    """Get DB interval IDs from checkpoint output."""
    intervals_file = "results/intervals/db/intervals.txt"
    
    if exists(intervals_file):
        with open(intervals_file) as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        checkpoint_output = checkpoints.create_db_intervals.get().output[0]
        with open(checkpoint_output) as f:
            lines = [l.strip() for l in f.readlines()]
    
    list_files = [os.path.basename(x) for x in lines]
    return [f.replace("-scattered.interval_list", "") for f in list_files]


def get_interval_vcfs(wc):
    """Get unfiltered interval VCF files."""
    intervals = get_db_intervals(wc)
    return expand(
        "results/vcfs/intervals/L{interval}.vcf.gz",
        interval=intervals,
    )


def get_gvcfs_for_db(wc):
    gvcfs = [get_final_gvcf(s) for s in SAMPLES_ALL]
    return {
        "gvcfs": gvcfs,
        "tbis": [g + ".tbi" for g in gvcfs],
        "interval": f"results/intervals/db/{wc.interval}-scattered.interval_list",
        "db_mapfile": "results/genomics_db/mapfile.txt",
    }


def get_concat_batch_size():
    """Get batch size for staged bcftools concat operations."""
    size = int(config["variant_calling"]["gatk"]["concat_batch_size"])
    if size < 2:
        raise ValueError("variant_calling.gatk.concat_batch_size must be >= 2")
    return size


def get_concat_max_rounds():
    """Get max allowed rounds for staged concat operations."""
    rounds = int(config["variant_calling"]["gatk"]["concat_max_rounds"])
    if rounds < 1:
        raise ValueError("variant_calling.gatk.concat_max_rounds must be >= 1")
    return rounds


def _ceil_div(a, b):
    return (a + b - 1) // b


def get_stage_chunk_counts(num_files):
    """Compute chunk counts for each staged concat round."""
    if num_files < 1:
        raise ValueError("Staged concat requires at least one input file")

    batch_size = get_concat_batch_size()
    max_rounds = get_concat_max_rounds()

    counts = []
    current = num_files
    while current > 1:
        if len(counts) >= max_rounds:
            raise ValueError(
                f"Staged concat exceeded variant_calling.gatk.concat_max_rounds={max_rounds}. "
                "Increase concat_batch_size or concat_max_rounds."
            )
        n_chunks = _ceil_div(current, batch_size)
        counts.append(n_chunks)
        current = n_chunks
    return counts


def staged_vcf_path(wc, round_idx, chunk_idx):
    return f"results/vcfs/staged/r{round_idx}/c{chunk_idx}.vcf.gz"


def staged_gvcf_path(wc, round_idx, chunk_idx):
    return f"results/gvcfs/staged/{wc.sample}/r{round_idx}/c{chunk_idx}.g.vcf.gz"


def get_stage_inputs(base_files, round_idx, chunk_idx, wc, path_builder):
    """Return input files for one staged concat chunk."""
    stage_counts = get_stage_chunk_counts(len(base_files))
    if not stage_counts:
        raise ValueError("Requested staged concat inputs with only one base file")

    if round_idx < 1 or round_idx > len(stage_counts):
        raise ValueError(
            f"Invalid staged concat round {round_idx}. Valid rounds are 1..{len(stage_counts)}"
        )

    n_chunks = stage_counts[round_idx - 1]
    if chunk_idx < 0 or chunk_idx >= n_chunks:
        raise ValueError(
            f"Invalid staged concat chunk {chunk_idx} for round {round_idx}. "
            f"Valid chunks are 0..{n_chunks - 1}"
        )

    if round_idx == 1:
        round_inputs = list(base_files)
    else:
        prev_chunks = stage_counts[round_idx - 2]
        round_inputs = [path_builder(wc, round_idx - 1, i) for i in range(prev_chunks)]

    batch_size = get_concat_batch_size()
    start = chunk_idx * batch_size
    end = min(start + batch_size, len(round_inputs))
    selected = round_inputs[start:end]

    if not selected:
        raise ValueError(
            f"No files selected for staged concat round={round_idx}, chunk={chunk_idx}"
        )
    return selected


def get_final_stage_file(base_files, wc, path_builder):
    """Return final staged output path (or original file if no staging is needed)."""
    if not base_files:
        raise ValueError("No files available for concat")

    stage_counts = get_stage_chunk_counts(len(base_files))
    if not stage_counts:
        return base_files[0]

    final_round = len(stage_counts)
    return path_builder(wc, final_round, 0)


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
    return [f"{gvcf}.tbi" for gvcf in stage_inputs]


def get_final_interval_gvcf_stage_file(wc):
    return get_final_stage_file(
        base_files=get_interval_gvcfs(wc),
        wc=wc,
        path_builder=staged_gvcf_path,
    )


def get_final_interval_gvcf_stage_tbi(wc):
    return get_final_interval_gvcf_stage_file(wc) + ".tbi"


def get_interval_vcf_stage_inputs(wc):
    return get_stage_inputs(
        base_files=get_interval_vcfs(wc),
        round_idx=int(wc.round),
        chunk_idx=int(wc.chunk),
        wc=wc,
        path_builder=staged_vcf_path,
    )


def get_interval_vcf_stage_tbis(wc):
    stage_inputs = get_interval_vcf_stage_inputs(wc)
    return [f"{vcf}.tbi" for vcf in stage_inputs]


def get_final_interval_vcf_stage_file(wc):
    return get_final_stage_file(
        base_files=get_interval_vcfs(wc),
        wc=wc,
        path_builder=staged_vcf_path,
    )


def get_final_interval_vcf_stage_tbi(wc):
    return get_final_interval_vcf_stage_file(wc) + ".tbi"

rule gatk_haplotypecaller_interval:
    input:
        unpack(haplotype_caller_input),
    output:
        gvcf=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz"),
        tbi=temp("results/interval_gvcfs/{sample}/{interval}.g.vcf.gz.tbi"),
    params:
        java_mem=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        ploidy=config["variant_calling"]["ploidy"],
        min_pruning=1 if config["variant_calling"]["expected_coverage"] == "low" else 2,
        min_dangling=1 if config["variant_calling"]["expected_coverage"] == "low" else 4,
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_haplotypecaller/{sample}/{interval}.benchmark.txt"
    log:
        "logs/gatk_haplotypecaller/{sample}/{interval}.log.txt",
    shell:
        """
        gatk HaplotypeCaller \
        --java-options {params.java_mem} \
        -R {input.ref} \
        -I {input.bam} \
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
        gvcf=temp("results/gvcfs/staged/{sample}/r{round}/c{chunk}.g.vcf.gz"),
        tbi=temp("results/gvcfs/staged/{sample}/r{round}/c{chunk}.g.vcf.gz.tbi"),
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/concat_interval_gvcfs/staged/{sample}/r{round}/c{chunk}.txt"
    log:
        "logs/concat_interval_gvcfs/staged/{sample}/r{round}/c{chunk}.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.gvcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.gvcf} - 2>> {log}
        tabix -p vcf {output.gvcf} 2>> {log}
        """


rule concat_interval_gvcfs:
    input:
        gvcf=get_final_interval_gvcf_stage_file,
        tbi=get_final_interval_gvcf_stage_tbi,
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    benchmark:
        "benchmarks/concat_interval_gvcfs/{sample}.txt"
    log:
        "logs/concat_interval_gvcfs/{sample}.txt"
    shell:
        """
        mv {input.gvcf} {output.gvcf} 2> {log}
        mv {input.tbi} {output.tbi} 2>> {log}
        """

rule create_db_mapfile:
    input:
        gvcfs=[get_final_gvcf(s) for s in SAMPLES_ALL],
    output:
        mapfile="results/genomics_db/mapfile.txt",
    run:
        with open(output.mapfile, "w") as f:
            for gvcf in input.gvcfs:
                sample = os.path.basename(gvcf).replace(".g.vcf.gz", "")
                print(sample, gvcf, sep="\t", file=f)


rule gatk_genomics_db_import:
    input:
        unpack(get_gvcfs_for_db),
    output:
        db=temp(directory("results/gatk_genomics_db/L{interval}")),
        tar="results/gatk_genomics_db/L{interval}.tar",
    params:
        java_mem=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genomics_db_import/{interval}.txt"
    log:
        "logs/gatk_genomics_db_import/{interval}.txt"
    shell:
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '{params.java_mem}' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals \
            --reader-threads {threads} \
            -L {input.interval} \
            --tmp-dir {resources.tmpdir} \
            --sample-name-map {input.db_mapfile} \
            &> {log}
        tar -cf {output.tar} {output.db} &>> {log}
        """

rule gatk_genotype_gvcfs:
    input:
        db="results/gatk_genomics_db/L{interval}.tar",
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/intervals/L{interval}.vcf.gz"),
        tbi=temp("results/vcfs/intervals/L{interval}.vcf.gz.tbi"),
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{int(resources.mem_mb * 0.9)}m",
        het_prior=config["variant_calling"]["gatk"]["het_prior"],
        db=subpath(input.db, strip_suffix=".tar"),
    resources:
        mem_mb=4096,
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/gatk_genotype_gvcfs/{interval}.txt"
    log:
        "logs/gatk_genotype_gvcfs/{interval}.txt"
    shell:
        """
        tar -xf {input.db}
        gatk GenotypeGVCFs \
            --java-options '{params.java_opts}' \
            -R {input.ref} \
            --heterozygosity {params.het_prior} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} \
            &> {log}
        """

rule concat_interval_vcfs_stage:
    input:
        vcfs=get_interval_vcf_stage_inputs,
        tbis=get_interval_vcf_stage_tbis,
    output:
        vcf=temp("results/vcfs/staged/r{round}/c{chunk}.vcf.gz"),
        tbi=temp("results/vcfs/staged/r{round}/c{chunk}.vcf.gz.tbi"),
    conda:
        "../../envs/bcftools.yaml"
    benchmark:
        "benchmarks/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    log:
        "logs/concat_interval_vcfs/staged/r{round}/c{chunk}.txt"
    shell:
        """
        bcftools concat -D -a -Ou {input.vcfs} 2> {log} \
            | bcftools sort -T {resources.tmpdir}/ -Oz -o {output.vcf} - 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_interval_vcfs:
    input:
        vcf=get_final_interval_vcf_stage_file,
        tbi=get_final_interval_vcf_stage_tbi,
    output:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
    benchmark:
        "benchmarks/concat_interval_vcfs/benchmark.txt"
    log:
        "logs/concat_interval_vcfs/log.txt"
    shell:
        """
        mv {input.vcf} {output.vcf} 2> {log}
        mv {input.tbi} {output.tbi} 2>> {log}
        """

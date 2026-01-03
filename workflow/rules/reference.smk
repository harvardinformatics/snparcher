"""
Reference genome processing rules for snpArcher v2.

Rules:
- reference_download: Download reference from NCBI if accession specified
- reference_link: Symlink local reference to results directory
- reference_index: Index reference with BWA, samtools, and create dict

Reference files must be bgzip-compressed (.fna.gz). Regular gzip is not supported.
"""

localrules:
    reference_link,


if REF_ACCESSION:

    rule reference_download:
        """Download reference genome from NCBI and bgzip compress."""
        output:
            ref=f"results/reference/{REF_FILE}",
        params:
            accession=REF_ACCESSION,
            dataset=f"results/reference/{REF_NAME}_dataset.zip",
            outdir=f"results/reference/{REF_NAME}",
            uncompressed=f"results/reference/{REF_NAME}.fna",
        conda:
            "../envs/fastq2bam.yml"
        log:
            "logs/reference_download.txt"
        benchmark:
            "benchmarks/reference_download.txt"
        shell:
            """
            mkdir -p {params.outdir}
            datasets download genome accession {params.accession} \
                --include genome \
                --filename {params.dataset} \
            && (7z x {params.dataset} -aoa -o{params.outdir} || unzip -o {params.dataset} -d {params.outdir}) \
            && cat {params.outdir}/ncbi_dataset/data/{params.accession}/*.fna > {params.uncompressed} \
            && bgzip -f {params.uncompressed} \
            2> {log}
            """


if REF_PATH:

    rule reference_link:
        """Symlink user-provided reference to results directory.

        Reference must be bgzip-compressed (.fna.gz). Regular gzip is not supported
        and will cause downstream tools to fail.
        """
        input:
            ref=REF_PATH,
        output:
            ref=f"results/reference/{REF_FILE}",
        log:
            "logs/reference_link.txt"
        shell:
            "ln -sf $(realpath {input.ref}) {output.ref} 2> {log}"


rule reference_index:
    """Index reference genome with BWA and samtools."""
    input:
        ref=f"results/reference/{REF_FILE}",
    output:
        indexes=expand(
            f"results/reference/{REF_FILE}.{{ext}}",
            ext=["sa", "pac", "bwt", "ann", "amb"],
        ),
        fai=f"results/reference/{REF_FILE}.fai",
        dictf=f"results/reference/{REF_NAME}.dict",
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/reference_index.txt"
    benchmark:
        "benchmarks/reference_index.txt"
    threads: 4
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} --output {output.fai} >> {log} 2>&1
        samtools dict {input.ref} -o {output.dictf} >> {log} 2>&1
        """

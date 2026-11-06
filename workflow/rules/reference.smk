def get_reference_type():
    """Determine reference source type."""
    source = config["reference"]["source"]
    if source.startswith(("http://", "https://", "ftp://")):
        return "url"
    if Path(source).exists():
        return "local"
    if re.fullmatch(r"GC[AF]_\d+(\.\d+)?", source):
        return "accession"
    if Path(REF_FILES["ref"]).exists():
        # Already prepared, so the source is never read again.
        return "local"
    return "invalid"


REF_TYPE = get_reference_type()


if REF_TYPE == "invalid":

    def invalid_reference_source(wildcards):
        raise WorkflowError(
            f"reference.source '{config['reference']['source']}' is not an existing file, "
            "a URL, or an NCBI assembly accession (GCA_ or GCF_, e.g. GCF_000002315.6). "
            "Relative paths are resolved from the working directory."
        )

    # Fails while the job graph is built, and only for runs that need the
    # reference; targets such as staging an external BAM still run.
    rule prepare_reference:
        input:
            invalid_reference_source,
        output:
            ref="results/reference/{name}.fa.gz",


elif REF_TYPE == "local":

    rule prepare_reference:
        output:
            ref="results/reference/{name}.fa.gz",
        params:
            source=config["reference"]["source"],
        conda:
            "../envs/reference.yaml"
        log:
            "logs/prepare_reference/{name}.txt"
        shell:
            """
            if [[ "{params.source}" == *.gz ]]; then
                gunzip -c "{params.source}" 2> {log} | bgzip -c > {output.ref} 2>> {log}
            else
                bgzip -c "{params.source}" > {output.ref} 2> {log}
            fi
            """


elif REF_TYPE == "url":

    rule prepare_reference:
        output:
            ref="results/reference/{name}.fa.gz",
        params:
            source=config["reference"]["source"],
        conda:
            "../envs/reference.yaml"
        log:
            "logs/prepare_reference/{name}.txt"
        shell:
            """
            if [[ "{params.source}" == *.gz ]]; then
                curl -fSL "{params.source}" 2> {log} | gunzip -c 2>> {log} | bgzip -c > {output.ref} 2>> {log}
            else
                curl -fSL "{params.source}" 2> {log} | bgzip -c > {output.ref} 2>> {log}
            fi
            """


else:  # accession

    rule prepare_reference:
        output:
            ref="results/reference/{name}.fa.gz",
        params:
            source=config["reference"]["source"],
        conda:
            "../envs/reference.yaml"
        log:
            "logs/prepare_reference/{name}.txt"
        shell:
            """
            TMPDIR=$(mktemp -d)
            datasets download genome accession "{params.source}" \
                --include genome \
                --filename "$TMPDIR/dataset.zip" &> {log}
            unzip -o "$TMPDIR/dataset.zip" -d "$TMPDIR" >> {log} 2>&1
            cat "$TMPDIR"/ncbi_dataset/data/{params.source}/*.fna | bgzip -c > {output.ref} 2>> {log}
            rm -rf "$TMPDIR"
            """


rule index_reference:
    input:
        ref="results/reference/{name}.fa.gz",
    output:
        gzi="results/reference/{name}.fa.gz.gzi",
        fai="results/reference/{name}.fa.gz.fai",
        dict="results/reference/{name}.dict",
        bwa_idx=multiext("results/reference/{name}.fa.gz", ".sa", ".pac", ".bwt", ".ann", ".amb"),
    conda:
        "../envs/reference.yaml"
    log:
        "logs/index_reference/{name}.txt"
    shell:
        """
        samtools faidx {input.ref} 2> {log}
        samtools dict {input.ref} -o {output.dict} >> {log} 2>&1
        bwa index {input.ref} >> {log} 2>&1
        """

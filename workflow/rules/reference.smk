def get_reference_type():
    """Determine reference source type."""
    source = config["reference"]["source"]
    if source.startswith(("http://", "https://", "ftp://")):
        return "url"
    elif Path(source).exists() or source.startswith("/"):
        return "local"
    else:
        return "accession"


REF_TYPE = get_reference_type()


if REF_TYPE == "local":

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

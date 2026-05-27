def glnexus_joint_input(wc):
    return {
        "gvcfs": get_joint_gvcf_paths(),
        "indexes": get_glnexus_joint_gvcf_indexes(),
    }


rule glnexus_joint:
    input:
        unpack(glnexus_joint_input),
    output:
        vcf=temp(RAW_VCF),
        idx=temp(RAW_VCF_INDEX),
    params:
        db="results/vcfs/GLnexus.DB",
        config="DeepVariant",
        index_args=BCFTOOLS_INDEX_ARGS,
    threads: 1
    conda:
        "../../envs/glnexus.yaml"
    benchmark:
        "benchmarks/glnexus_joint.txt"
    log:
        "logs/glnexus_joint.txt"
    shell:
        """
        glnexus_db="{params.db}"
        rm -rf "$glnexus_db"
        trap 'rm -rf "$glnexus_db"' EXIT

        glnexus_mem_gbytes=$(( {resources.mem_mb_reduced} / 1024 ))
        if [ "$glnexus_mem_gbytes" -lt 1 ]; then
            glnexus_mem_gbytes=1
        fi

        glnexus_cli --dir "$glnexus_db" --mem-gbytes "$glnexus_mem_gbytes" \
            --config {params.config} --threads {threads} {input.gvcfs} 2> {log} \
            | bcftools view -Oz -o {output.vcf} - 2>> {log}
        bcftools index {params.index_args} {output.vcf} 2>> {log}
        """

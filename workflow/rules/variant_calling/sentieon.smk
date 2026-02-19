def sentieon_haplotyper_input(wildcards):
    input_type = get_sample_input_type(wildcards.sample)
    
    if input_type == "gvcf":
        raise ValueError(f"Sample {wildcards.sample} has input_type 'gvcf', should not call haplotyper")
    
    bam = get_final_bam(wildcards.sample)
    return {
        "bam": bam,
        "bai": bam + ".bai",
        **REF_FILES,
    }


def sentieon_combine_gvcf_input(wc):
    return {
        "gvcfs": [get_final_gvcf(s) for s in SAMPLES_ALL],
        "tbis": [get_final_gvcf(s) + ".tbi" for s in SAMPLES_ALL],
    }


rule sentieon_haplotyper:
    input:
        unpack(sentieon_haplotyper_input),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        lic=config["variant_calling"]["sentieon"]["license"],
        ploidy=config["variant_calling"]["ploidy"],
    threads: 8
    conda:
        "../../envs/sentieon.yaml"
    benchmark:
        "benchmarks/sentieon_haplotyper/{sample}.txt"
    log:
        "logs/sentieon_haplotyper/{sample}.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver \
            -r {input.ref} \
            -t {threads} \
            -i {input.bam} \
            --algo Haplotyper \
            --genotype_model multinomial \
            --emit_mode gvcf \
            --emit_conf 30 \
            --call_conf 30 \
            --ploidy {params.ploidy} \
            {output.gvcf} \
            2> {log}
        """

rule sentieon_combine_gvcf:
    input:
        unpack(sentieon_combine_gvcf_input),
        **REF_FILES,
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        lic=config["variant_calling"]["sentieon"]["license"],
        gvcf_args=lambda wc, input: " ".join([f"-v {g}" for g in input.gvcfs]),
    threads: 8
    conda:
        "../../envs/sentieon.yaml"
    benchmark:
        "benchmarks/sentieon_combine_gvcf.txt"
    log:
        "logs/sentieon_combine_gvcf.txt"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver \
            -r {input.ref} \
            -t {threads} \
            --algo GVCFtyper \
            --emit_mode VARIANT \
            {output.vcf} \
            {params.gvcf_args} \
            2> {log}
        """


rule variant_filtration:
    input:
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
        **REF_FILES,
    output:
        vcf="results/vcfs/filtered.vcf.gz",
        tbi="results/vcfs/filtered.vcf.gz.tbi",
    conda:
        "../../envs/gatk.yaml"
    benchmark:
        "benchmarks/variant_filtration.txt"
    log:
        "logs/variant_filtration.txt"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            --output {output.vcf} \
            --filter-name "RPRS_filter" \
            --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" \
            --filter-name "FS_SOR_filter" \
            --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \
            --filter-name "MQ_filter" \
            --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
            --filter-name "QUAL_filter" \
            --filter-expression "QUAL < 30.0" \
            --create-output-variant-index \
            --invalidate-previous-filters true \
            &> {log}
        """

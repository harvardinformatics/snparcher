"""
Sentieon variant calling rules for snpArcher v2.

Requires Sentieon license configured in config['sentieon_lic'].
"""


def _sentieon_gvcf_inputs(wildcards):
    gvcfs = expand("results/gvcfs/{sample}.g.vcf.gz", sample=SAMPLE_IDS)
    tbis = expand("results/gvcfs/{sample}.g.vcf.gz.tbi", sample=SAMPLE_IDS)
    return {"gvcfs": gvcfs, "tbis": tbis}


def _sentieon_gvcf_args(wildcards):
    gvcfs = expand("results/gvcfs/{sample}.g.vcf.gz", sample=SAMPLE_IDS)
    return " ".join(["-v " + gvcf for gvcf in gvcfs])


rule call_sentieon_haplotyper:
    """Generate GVCF for a single sample using Sentieon Haplotyper."""
    wildcard_constraints:
        sample="|".join(BAM_REQUIRED_SAMPLES) if BAM_REQUIRED_SAMPLES else "NOMATCH",
    input:
        unpack(get_final_bam),
        unpack(get_ref_bundle),
    output:
        gvcf="results/gvcfs/{sample}.g.vcf.gz",
        tbi="results/gvcfs/{sample}.g.vcf.gz.tbi",
    params:
        lic=config["variant_calling"]["sentieon"]["license"],
        ploidy=config["variant_calling"]["gatk"]["ploidy"],
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/call_sentieon_haplotyper/{sample}.txt"
    benchmark:
        "benchmarks/call_sentieon_haplotyper/{sample}.txt"
    threads: 8
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} -t {threads} -i {input.bam} \
            --algo Haplotyper \
            --genotype_model multinomial \
            --emit_mode gvcf \
            --emit_conf 30 \
            --call_conf 30 \
            --ploidy {params.ploidy} \
            {output.gvcf} 2> {log}
        """


rule call_sentieon_gvcftyper:
    input:
        unpack(_sentieon_gvcf_inputs),
        unpack(get_ref_bundle),
    output:
        vcf=temp("results/vcfs/raw.vcf.gz"),
        tbi=temp("results/vcfs/raw.vcf.gz.tbi"),
    params:
        glist=_sentieon_gvcf_args,
        lic=config["variant_calling"]["sentieon"]["license"],
    conda:
        "../../envs/sentieon.yml"
    log:
        "logs/call_sentieon_gvcftyper.txt"
    benchmark:
        "benchmarks/call_sentieon_gvcftyper.txt"
    threads: 8
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} -t {threads} \
            --algo GVCFtyper --emit_mode VARIANT \
            {output.vcf} {params.glist} 2> {log}
        """


rule call_sentieon_filter:
    """Apply GATK recommended hard filters to VCF using parallel filtering."""
    input:
        unpack(get_ref_bundle),
        vcf="results/vcfs/raw.vcf.gz",
        tbi="results/vcfs/raw.vcf.gz.tbi",
    output:
        vcf=f"results/{REF_NAME}_raw.vcf.gz",
        tbi=f"results/{REF_NAME}_raw.vcf.gz.tbi",
    conda:
        "../../envs/bam2vcf.yml"
    log:
        "logs/call_sentieon_filter.txt"
    benchmark:
        "benchmarks/call_sentieon_filter.txt"
    threads: 8
    shadow:
        "minimal"
    shell:
        """
        # Get the contig names from the .fai index
        contigs=$(cut -f1 {input.fai})
        
        # Create a function that will be passed to gnu parallel
        filter_contig() {{
            contig=$1
            echo $contig

            gatk --java-options "-Xmx4g" VariantFiltration \
                -R {input.ref} \
                -L ${{contig}} \
                -V {input.vcf} \
                --output filter_${{contig}}.vcf.gz \
                --filter-name "RPRS_filter" \
                --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" \
                --filter-name "FS_SOR_filter" \
                --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \
                --filter-name "MQ_filter" \
                --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
                --filter-name "QUAL_filter" \
                --filter-expression "QUAL < 30.0" \
                --invalidate-previous-filters true
        }}
        
        export -f filter_contig
        
        # Pass each contig to gnu parallel
        parallel -j {threads} filter_contig ::: ${{contigs}} 2> {log}
        
        bcftools concat filter_*.vcf.gz --threads {threads} -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

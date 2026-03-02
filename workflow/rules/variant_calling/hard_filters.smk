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

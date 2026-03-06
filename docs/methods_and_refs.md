--------------------------------------------------------------------------------
VERSION 1: STANDARD GATK PIPELINE (sentieon=false)
--------------------------------------------------------------------------------
Raw sequencing reads were quality-filtered using fastp v0.20.1 (Chen et al. 2018) with automatic adapter detection for paired-end reads. Reads were aligned to the reference genome using BWA v0.7.17 MEM (Li 2013) with the -M flag to mark shorter split hits as secondary alignments. BAM files were sorted and merged with samtools v1.14 (Li et al. 2009), and duplicate reads were marked using sambamba v0.8.0 (Tarasov et al. 2015). Per-sample variant calling was performed with GATK v4.1.8.0 HaplotypeCaller (McKenna et al. 2010) in GVCF mode (--emit-ref-confidence GVCF) with ploidy set to 2. GVCFs were consolidated using GenomicsDBImport with POSIX filesystem optimizations enabled and a batch size of 25. Joint genotyping was performed with GenotypeGVCFs using a heterozygosity prior of 0.005. Variants were filtered using GATK VariantFiltration following best practices (Van der Auwera et al. 2013) with the following criteria: ReadPosRankSum < -8.0 for SNPs and < -20.0 for indels, QD < 2.0, FS > 60.0 or SOR > 3.0 for SNPs, FS > 200.0 or SOR > 10.0 for indels, MQ < 40.0, MQRankSum < -12.5, and QUAL < 30.0. Filtered variants were annotated with bcftools v1.10 (Danecek et al. 2021). Callable sites were determined by computing per-base depth using mosdepth v0.3.10 (Pedersen & Quinlan 2018) in D4 format (Hou et al. 2021) and k-mer mappability (k=150, exact matches only) using genmap v1.3.0 (Pockrandt et al. 2020), then combined using bedtools v2.30.0 (Quinlan & Hall 2010) and clam v0.2.2 (Mirchandani et al. 2025). As part of the quality control and postprocessing of VCFs, snpArcher uses Plink2 v2.00a2.3 (Chang et al. 2015) and ADMIXTURE (Alexander et al. 2009). The workflow was orchestrated using Snakemake (Mölder et al. 2021).

If using low coverage mode, please add “HaplotypeCaller was run using parameters optimized for low-coverage data (--min-pruning 1, --min-dangling-branch-length 1).

--------------------------------------------------------------------------------
VERSION 2: SENTIEON PIPELINE (sentieon=true)
--------------------------------------------------------------------------------

Raw sequencing reads were quality-filtered using fastp v0.20.1 (Chen et al. 2018) with automatic adapter detection for paired-end reads. Reads were aligned to the reference genome using Sentieon's optimized BWA implementation (Kendig et al. 2019) with the -M flag and batch processing size of 10,000,000 (-K 10000000), followed by sorting with maximum BAM compression. BAM files were merged with samtools (Li et al. 2009) and duplicate reads were identified using Sentieon's two-pass algorithm (LocusCollector for scoring followed by Dedup for marking). Per-sample variant calling was performed with Sentieon Haplotyper in GVCF mode using the multinomial genotype model (--genotype_model multinomial), emit and call confidence thresholds of 30, and ploidy set to 2. Joint genotyping was executed using Sentieon GVCFtyper with variant-only emission mode (--emit_mode VARIANT). Variants were filtered using GATK v4.1.8.0 VariantFiltration (McKenna et al. 2010) following best practices (Van der Auwera et al. 2013) with the following criteria: ReadPosRankSum < -8.0 for SNPs and < -20.0 for indels, QD < 2.0, FS > 60.0 or SOR > 3.0 for SNPs, FS > 200.0 or SOR > 10.0 for indels, MQ < 40.0, MQRankSum < -12.5, and QUAL < 30.0. Callable sites were determined by computing per-base depth using mosdepth v0.3.10 (Pedersen & Quinlan 2018) in D4 format (Hou et al. 2021) and k-mer mappability (k=150, exact matches only) using genmap v1.3.0 (Pockrandt et al. 2020), then combined using bedtools v2.30.0 (Quinlan & Hall 2010) and clam v0.2.2 (Mirchandani et al. 2025). As part of the quality control and postprocessing of VCFs, snpArcher uses Plink2 v2.00a2.3 (Chang et al. 2015) and ADMIXTURE (Alexander et al. 2009). The workflow was orchestrated using Snakemake (Mölder et al. 2021).

================================================================================
REFERENCES
================================================================================

Alexander DH, Novembre J, Lange K. Fast model-based estimation of ancestry in unrelated individuals. Genome Research. 2009;19(9):1655–1664.

Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience. 2015;4(1):7.

Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884–i890.

Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, et al. The variant call format and VCFtools. Bioinformatics. 2011;27(15):2156–2158.

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008.

Hou H, Pedersen B, Quinlan A. Balancing efficient analysis and storage of quantitative genomics data with the D4 format and d4tools. Nature Computational Science. 2021;1(6):441–447.

Kendig KI, Baheti S, Bockol MA, Drucker TM, Hart SN, Heldenbrand JR, Hernaez M, Hudson ME, Kalmbach MT, Klee EW, et al. Sentieon DNASeq variant calling workflow demonstrates strong computational performance and accuracy. Frontiers in Genetics. 2019;10:736.

Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv. 2013. [preprint] arXiv:1303.3997.

Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009;25(14):1754–1760.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. The sequence alignment/map format and SAMtools. Bioinformatics. 2009;25(16):2078–2079.

McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research. 2010;20(9):1297–1303.

Mirchandani C, Enbody E, Sackton TB, Corbett-Detig R. Efficient Estimation of Nucleotide Diversity and Divergence Using Callable Loci (and More). Molecular Biology and Evolution. 2025;42(12):msaf282.

Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, et al. Sustainable data analysis with Snakemake. F1000Research. 2021;10:33.

Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867–868.

Pockrandt C, Alzamel M, Iliopoulos CS, Reinert K. GenMap: ultra-fast computation of genome mappability. Bioinformatics. 2020;36(12):3687–3692.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics. 2007;81(3):559–575.

Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841–842.

Tarasov A, Vilella AJ, Cuppen E, Nijman IJ, Prins P. Sambamba: fast processing of NGS alignment formats. Bioinformatics. 2015;31(12):2032–2034.

Van der Auwera GA, Carneiro MO, Hartl C, Poplin R, Del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, et al. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics. 2013;43(1):11.10.1–11.10.33.


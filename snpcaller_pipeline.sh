#!/bin/bash
### Pipeline for bam file processing post-bwa mem in order to call SNPs

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load bwa/0.7.17.r1188
module load picardtools/2.7.1
module load GATK/3.6
module load snpeff/4.0

picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa
snpeff=/share/apps/snpEff-4.0/snpEff.jar
sample="$1"

#Call Variants
java -Xmx60g -jar $GATK \
     -T HaplotypeCaller \
     -R $genome \
     -I "$sample"_removedup_realigned.bam \
     -o "$sample"_raw_variants.vcf


#Should I go the GVCF?
#java -Xmx60g -jar $GATK \
#     -T HaplotypeCaller \
#     -R $genome \
#     -I "$sample"_removedup_realigned.bam \
#     -o "$sample"_raw_variants.g.vcf \
#     -ERC GVCF \
#     --variant_index_type LINEAR \
#     --variant_index_parameter 128000


#java -Xmx60g -jar $GATK \
#     -T GenotypeGVCFs \
#     -R $genome \
#     -V "$sample"_raw_variants.g.vcf \
#     -o "$sample"_raw_varinats.vcf



#Extract SNPs & Indels
java -Xmx60g -jar $GATK \
     -T SelectVariants \
     -R $genome \
     -V "$sample"_raw_variants.vcf \
     -selectType SNP \
     -o "$sample"_raw_snps.vcf

java -Xmx60g -jar $GATK \
     -T SelectVariants \
     -R $genome \
     -V "$sample"_raw_variants.vcf \
     -selectType INDEL \
     -o "$sample"_raw_indels.vcf

#Filter SNPs
java -Xmx60g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V "$sample"_raw_snps.vcf \
     --filterExpression 'QD < 2.0 || FS > 60.0 || MG < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
     --filterName "basic_snp_filter" \
     -o "$sample"_filtered_snps.vcf

#Filter Indels
java -Xmx60g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V "$sample"_raw_indels.vcf \
     --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
     --filterName "basic_indel_filter" \
     -o "$sample"_filtered_indels.vcf

#Base Quality Score Recalibration (BQSR #1)
java -Xmx60g -jar $GATK \
     -T BaseRecalibrator \
     -R $genome \
     -I "$sample"_removedup_realigned.bam \
     -knownSites "$sample"_filtered_snps.vcf \
     -knownSites "$sample"_filtered_indels.vcf \
     -o "$sample"_recal_data.table

# BQSR #2
java -Xmx60g -jar $GATK \
     -T BaseRecalibrator \
     -R $genome \
     -I "$sample"_removedup_realigned.bam \
     -knownSites "$sample"_filtered_snps.vcf \
     -knownSites "$sample"_filtered_indels.vcf \
     -BQSR "$sample"_recal_data.table \
     -o "$sample"_post_recal_data.table

#Analyze Covariates
java -Xmx60g -jar $GATK \
     -T AnalyzeCovariates \
     -R $genome \
     -before "$sample"_recal_data.table \
     -after "$sample"_post_recal_data.table \
     -plots "$sample"_recalibration_plots.pdf

#Apply BQSR
java -Xmx60g -jar $GATK \
     -T PrintReads \
     -R $genome \
     -I "$sample"_removedup_realigned.bam \
     -BQSR "$sample"_recal_data.table \
     -o "$sample"_recal_reads.bam

#Call Variants
java -Xmx60g -jar $GATK \
     -T HaplotypeCaller \
     -R $genome \
     -I "$sample"_recal_reads.bam \
     -o "$sample"_raw_variants_recal.vcf

#Extract SNPS & Indels
java -Xmx60g -jar $GATK \
     -T SelectVariants \
     -R $genome \
     -V "$sample"_raw_variants_recal.vcf \
     -selectType SNP \
     -o "$sample"_raw_snps_recal.vcf

java -Xmx60g -jar $GATK\
     -T SelectVariants \
     -R $genome\
     -V "$sample"_raw_variants_recal.vcf \
     -selectType INDEL \
     -o "$sample"_raw_indels_recal.vcf

#Filter SNPs & Indels
java -Xmx60g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V "$sample"_raw_snps_recal.vcf \
     --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
     --filterName "basic_snp_filter" \
     -o "$sample"_filtered_snps_final.vcf

java -Xmx60g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V "$sample"_raw_indels_recal.vcf \
     --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' \
     --filterName "basic_indel_filter" \
     -o "$sample"_filtered_indels_final.vcf

#Annotate SNPs 
#java -Xmx40g -jar $snpeff -v snpeff_db "$sample"_filtered_snps_final.vcf > "$sample"_filtered_snps_final.ann.vcf
     

#Compute coverage statistics
bedtools genomecov -bga -ibam "$sample"_recal_reads.bam > "$sample"_genomecov.bedgraph

#Compile statistics
#./parse_metrics.sh 







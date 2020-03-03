
1;95;0c#!/bin/bash
### Pipeline for bam file processing post-bwa mem in order to call SNPs

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load picardtools/2.7.1
module load GATK/3.6

#mkdir /scratch/sodell
#cp /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa /scratch/sodell/Zea_mays.AGPv3.dna.toplevel.fa
#cp /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.* /scratch/sodell/
#cp /group/jrigrp6/RILAB_data/N2/alignment/Tvid_AGPv4/Tvid_removedup_realigned.bam /scratch/sodell/Tvid_removedup_realigned.bam
#cp /group/jrigrp6/RILAB_data/N2/alignment/Tvic_AGPv4/Tvic_removedup_realigned.bam /scratch/sodell/Tvic_removedup_realigned.bam
#cp /group/jrigrp6/RILAB_data/N2/alignment/Lo8_AGPv4/Lo8-AGPv4_removedup_realigned.bam /scratch/sodell/Lo8_removedup_realigned.bam
#cp /group/jrigrp6/RILAB_data/N2/alignment/Tvid_AGPv4/Tvid_removedup_realigned.bai /scratch/sodell/Tvid_removedup_realigned.bai
#cp /group/jrigrp6/RILAB_data/N2/alignment/Tvic_AGPv4/Tvic_removedup_realigned.bai /scratch/sodell/Tvic_removedup_realigned.bai
#cp /group/jrigrp6/RILAB_data/N2/alignment/Lo8_AGPv4/Lo8-AGPv4_removedup_realigned.bai /scratch/sodell/Lo8_removedup_realigned.bai


picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/scratch/sodell/Zea_mays.AGPv3.dna.toplevel.fa
home=/scratch/sodell/

# add a Read group
#java -Xmx20g -jar $picard AddOrReplaceReadGroups \
#     I="$sample".bam \
#     O="$sample"_with_RG.bam \
#     SORT_ORDER=coordinate \
#     RGID=$sample \
#     RGLB=$sample \
#     RGPL=illumina \
#     RGPU=$sample \
#     RGSM=$sample
# remove duplicate reads
#java -Xmx20g -jar $picard MarkDuplicates \
#     VALIDATION_STRINGENCY=LENIENT \
#     INPUT="$sample"_with_RG.bam \
#     OUTPUT="$sample"_sorted_markdup.bam \
#     METRICS_FILE="$sample"_metrics.txt
# Realignment Around Indels
#java -Xmx20g -jar $picard BuildBamIndex \
#     INPUT="$sample"_sorted_markdup.bam
#java -Xmx20g -jar $GATK \
#     -T RealignerTargetCreator \
#     -R $genome \
#     -I "$sample"_sorted_markdup.bam \
#     -o "$sample"_sorted_markdup.intervals
#java -Xmx20g -jar $GATK \
#     -T IndelRealigner \
#     -R $genome \
#     -I "$sample"_sorted_markdup.bam \
#     -targetIntervals "$sample"_sorted_markdup.intervals \
#     -o "$sample"_removedup_realigned.bam
# HaplotypeCaller
#for sample in Tvic Tvid; do
    #    java -Xmx128g -jar $GATK \
#	 -T HaplotypeCaller \
#	 -R $genome \
#	 -I $home/"$sample"_removedup_realigned.bam \
#	 --genotyping_mode DISCOVERY \
#	 -o $home/"$sample"_removedup_realigned.g.vcf \
#	 -ERC GVCF
#done
# once you have all your gvcf (from the last command line). You can call your snps together with the GenotypeGVCFs below
java -Xmx128g -jar $GATK \
     -T GenotypeGVCFs \
     -R $genome \
     --variant $home/Lo8_removedup_realigned.g.vcf\
     --variant $home/Tvic_removedup_realigned_all.g.vcf\
     --variant $home/Tvid_removedup_realigned_all.g.vcf\
     -o $home/AGPv3_raw_variants.vcf
     #this is the part were you extract and filter your snps

java -Xmx128g -jar $GATK \
     -T SelectVariants \
     -R $genome \
     -V $home/AGPv3_raw_variants.vcf \
     -selectType SNP \
     -o $home/AGPv3_raw_snps.vcf
#This creates a VCF file called raw_snps.vcf, containing just the SNPs from the original file of raw variants and the VariantFiltration apply the filters to the snps file
java -Xmx128g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V $home/AGPv3_raw_snps.vcf \
     --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
     --filterName "my_snp_filter" \
     -o $home/AGPv3_filtered_snps.vcf

#you can do the same with the indels: Extract the Indels from the call set
java -Xmx128g -jar $GATK \
     -T SelectVariants \
     -R $genome \
     -V $home/AGPv3_raw_variants.vcf \
     -selectType INDEL \
     -o $home/AGPv3_raw_indels.vcf
#Apply the filter to the Indel call set
java -Xmx128g -jar $GATK \
     -T VariantFiltration \
     -R $genome \
     -V $home/AGPv3_raw_indels.vcf \
     --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
     --filterName "my_indel_filter" \
     -o $home/AGPv3_filtered_indels.vcf

#mv $home/AGPv3_raw_variants.vcf /group/jripgr6/RILAB_data/N2/alignment/AGPv3_raw_variants.vcf
#mv $home/AGPv3_raw_snps.vcf /group/jripgr6/RILAB_data/N2/alignment/AGPv3_raw_snps.vcf
#mv $home/AGPv3_raw_indels.vcf /group/jripgr6/RILAB_data/N2/alignment/AGPv3_raw_indels.vcf
#mv $home/AGPv3_filtered_snps.vcf /group/jrigrp6/RILAB_data/N2/alignment/AGPv3_filtered_snps.vcf
#mv $home/AGPv3_filtered_indels.vcf /group/jrigrp6/RILAB_data/N2/alignment/AGPv3_filtered_indels.vcf


















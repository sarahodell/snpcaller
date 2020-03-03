#!/bin/bash
### Pipeline for creating alternate reference for sample using FastaAlternateReferenceMaker

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load bwa/0.7.17.r1188
module load picardtools/2.7.1
module load GATK/3.6
module load snpeff/4.0

chr=$SLURM_ARRAY_TASK_ID
picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/scratch/sodell/Zea_mays.AGPv3.dna.toplevel.fa
home=/scratch/sodell/
sample="Tvid"

#Call Variants
java -Xmx60g -jar $GATK \
     -T FastaAlternateReferenceMaker \
     -R $genome \
     -o "$sample"_AGPv3_pseudoref.fa \
     -V "$sample"_filtered_snps_final.vcf \
     -V "$sample"_filtered_indels_final.vcf 

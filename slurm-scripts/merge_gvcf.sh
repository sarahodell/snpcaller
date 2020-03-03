#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J callvariants
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 600:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load picardtools/2.7.1
module load GATK/3.6

#chr=$SLURM_ARRAY_TASK_ID
picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/scratch/sodell/Zea_mays.AGPv3.dna.toplevel.fa
home=/scratch/sodell/
sample="Lo8"

java -Xms64g -jar $GATK \
     -T CombineGVCFs \
     -R $genome \
     --variant $home/"$sample"_removedup_realigned_chr1.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr2.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr3.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr4.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr5.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr6.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr7.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr8.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr9.g.vcf \
     --variant $home/"$sample"_removedup_realigned_chr10.g.vcf \
     -o $home/"$sample"_removedup_realigned_all.g.vcf 



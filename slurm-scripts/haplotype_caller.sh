#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J callvariants
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 96:00:00
#SBATCH --ntasks=8
#SBATCH --mem 64G

module load jdk/1.8
module load R
module load java/1.8
module load maven
module load picardtools/2.7.1
module load GATK/3.6

chr=2
#chr=$SLURM_ARRAY_TASK_ID
picard=/share/apps/picard-tools-2.7.1/picard.jar
GATK=/share/apps/GATK-3.6/GenomeAnalysisTK.jar
genome=/scratch/sodell/Zea_mays.AGPv3.dna.toplevel.fa
home=/scratch/sodell/
sample="Lo8"

java -Xms64g -jar $GATK \
     -T HaplotypeCaller \
     -R $genome \
     -I $home/"$sample"_removedup_realigned.bam \
     -L $chr \
     --genotyping_mode DISCOVERY \
     -o $home/"$sample"_removedup_realigned_chr$chr.g.vcf \
     -ERC GVCF


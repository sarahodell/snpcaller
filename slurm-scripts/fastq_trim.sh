#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J fqtrim
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load python/2.7.6
module load cutadapt

cutadapt -m 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTAT Lo8_R1_filtered.fastq.gz -o Lo8_R1_filtered_cut.fastq.gz

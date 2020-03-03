#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J fastqc
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load perlbrew
module load fastx
module load python
module load cutadapt

fastq_quality_filter -v -q 20 -p 75 -i Lo8_R2.fastq -o Lo8_R2_filtered.fastq

cutadapt -m 20 -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG Lo8_R2_filtered.fastq -o Lo8_R2_filtered_cut.fastq

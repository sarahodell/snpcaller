#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J fastqc
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load python
module load cutadapt


cutadapt -m 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG Tvic_R2_filtered.fastq -o Tvic_R2_filtered_cut.fastq

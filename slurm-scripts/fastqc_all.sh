#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J fastqc
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load perlbrew
module load fastqc

fastqc Tvic_R1.fastq -o Tvic_qc

fastqc Tvic_R2.fastq -o Tvic_qc

fastqc Tvid_R1.fastq -o Tvid_qc

fastqc Tvid_R2.fastq -o Tvid_qc

#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J fastqc
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 24:00:00

module load perlbrew
module load fastqc

gunzip -ck /group/jrigrp6/RILAB_data/N2/raw_data/Lo8_R1.fastq.gz > Lo8_R1.fastq
fastqc Lo8_R1.fastq -o Lo8_qc

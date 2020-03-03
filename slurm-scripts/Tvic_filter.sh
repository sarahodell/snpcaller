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

fastq_quality_filter -v -q 20 -p 75 -i Tvic_R1.fastq -o Tvic_R1_filtered.fastq
fastq_quality_filter -v -q 20 -p 75 -i Tvic_R2.fastq -o Tvic_R2_filtered.fastq

cutadapt -m 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGC Tvic_R1_filtered.fastq -o Tvic_R1_filtered_cut.fastq

cutadapt -m 20 -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG Tvic_R2_filtered.fastq -o Tvic_R2_filtered_cut.fastq

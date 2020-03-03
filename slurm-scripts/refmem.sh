#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J refmem
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 336:00:00
#SBATCH --ntasks=16

module load bwa/0.7.17.r1188

bwa mem -M -t 16 /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa Lo8_qc/Lo8_R1.fastq Lo8_qc/Lo8_R2.fastq > Lo8_qc/Lo8_bwa_mem.sam







#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J alignR1
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 200:00:00

module load bwa

bwa aln /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa Lo8_qc/Lo8_R1_filtered_cut.fastq -f Lo8_qc/Lo8_R1_alignment.sai







#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J TvicR2
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 200:00:00

module load bwa

bwa aln /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa Tvic_qc/Tvic_R2_filtered_cut.fastq -f Tvic_qc/Tvic_R2_alignment.sai





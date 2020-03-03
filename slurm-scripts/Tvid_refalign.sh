#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J TvidR1
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 200:00:00

module load bwa

bwa aln /home/sodell/projects/resources/Zea_maysAGPv3/Zea_mays.AGPv3.dna.toplevel.fa Tvid_qc/Tvid_R1_filtered_cut.fastq -f Tvid_qc/Tvid_R1_alignment.sai






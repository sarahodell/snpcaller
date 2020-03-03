#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J Tvidbam
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 336:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

#./sam_to_bam.sh Tvid_qc/Tvid

#./bam_pipeline.sh Tvid_qc/Tvid 

./bam_finish.sh Tvid_qc/Tvid








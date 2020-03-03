#!/bin/bash -l
#SBATCH -D /home/sodell/projects/N2
#SBATCH -J callvariants
#SBATCH -o /home/sodell/projects/N2/slurm-logs/out-%j.txt
#SBATCH -e /home/sodell/projects/N2/slurm-logs/error-%j.txt
#SBATCH -t 600:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

#./sam_to_bam.sh Lo8_qc/Lo8

#./bam_pipeline.sh Lo8_qc/Lo8

#./bam_finish.sh Lo8_qc/Lo8

#./snpcaller_pipeline.sh Lo8_qc/Lo8

./anne_workflow.sh








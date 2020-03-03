#!/bin/bash

date=$(data "+Y_%m_%d")
file="slurm-logs/${data}_test.log"

snakemake --jobs 1 \
	  --rerun-incomplete \
	  --latency-wait 60 \
	  --cluster-config cluster.json \
	  --cluster "sbatch -p {cluster.p} -o {cluster.o} -e {cluster.e} \
	  --cpus-per-task {cluter.cpus-per-task}" \
	  -p &>> ${file}


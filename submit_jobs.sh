#!/bin/bash
#SBATCH --job-name=shotgun_preprocess
#SBATCH --time=5-00:00:00
#SBATCH --partition=medium
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

snakemake \
    --use-conda \
    --slurm-logdir slurm_logs \
    --conda-prefix=~/snakemake_conda_prefix \
    --executor slurm \
    --jobs 3000 \
    --rerun-incomplete \
    --configfile config/config.yaml \
    --latency-wait 60 \
    --keep-going \
    --printshellcmds \
    2>&1 | tee snakemake_$(date +%Y%m%d_%H%M%S).log
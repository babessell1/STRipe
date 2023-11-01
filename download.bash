#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=download_bams
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=11GB
#SBATCH --ntasks=11
#SBATCH --time=6:00:00
#SBATCH --output=logs/stripe_download.out
#SBATCH --error=logs/stripe_download.err

# 16 cores
#salloc -A remills99 -p standard -N 1 -n 11 -t 6:00:00 --mem=11G

#eval "$(conda shell.bash hook)"
#conda activate tricolor
snakemake -s download.smk --unlock
snakemake -s download.smk \
    --cores 11 \
    --rerun-incomplete \
    --use-conda
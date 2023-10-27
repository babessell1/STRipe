#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=stripe_download
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=1GB
#SBATCH --time=6:00:00
#SBATCH --output=logs/stripe_download.out
#SBATCH --error=logs/stripe_download.err

#eval "$(conda shell.bash hook)"
#conda activate tricolor
snakemake -s download.smk \
    -n -r \
    --dag \
    --verbose \
    --cores 1 \
    --rerun-incomplete \
    --use-conda

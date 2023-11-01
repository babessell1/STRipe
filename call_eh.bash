#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=call_eh
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=11
#SBATCH --time=6:00:00
#SBATCH --output=logs/stripe_download.out
#SBATCH --error=logs/stripe_download.err

#eval "$(conda shell.bash hook)"
#conda activate tricolor

#salloc -A remills99 -p standard -N 1 -n 11 -t 6:00:00 --mem=22G

snakemake -s call_eh.smk --unlock
snakemake -s call_eh.smk \
    --cores 32 \
    --rerun-incomplete \
    --use-conda

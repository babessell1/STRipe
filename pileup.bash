#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=get_pileups
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=44GB
#SBATCH --ntasks=11
#SBATCH --time=6:00:00
#SBATCH --output=logs/stripe_download.out
#SBATCH --error=logs/stripe_download.err

#eval "$(conda shell.bash hook)"
#conda activate tricolor

#salloc -A bioinf593f23_class -p standard -N 1 -n 11 -t 6:00:00 --mem=44G

snakemake -s get_pileup.smk --unlock
snakemake -s get_pileup.smk \
    --cores 11 \
    --rerun-incomplete \
    --use-conda

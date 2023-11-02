#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=call_eh
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=11
#SBATCH --time=7:59:00
#SBATCH --output=logs/eh.out
#SBATCH --error=logs/eh.err

#eval "$(conda shell.bash hook)"
#conda activate tricolor

#salloc -A bioinf593f23_class -p standard -N 1 -n 11 -t 6:00:00 --mem=32G

snakemake -s call_eh.smk --unlock
snakemake -s call_eh.smk \
    --cores 32 \
    --rerun-incomplete \
    --use-conda

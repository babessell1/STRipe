#!/bin/bash

#SBATCH --account=
#SBATCH --job-name=stromboli
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --time=1:00:00
#SBATCH --output=logs/stripe.out
#SBATCH --error=logs/stripe.err

#eval "$(conda shell.bash hook)"
conda activate tricolor
snakemake

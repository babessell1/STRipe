#!/bin/bash

#SBATCH --job-name=stromboli
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --time=1:00:00
#SBATCH --output=logs/stromboli.out
#SBATCH --error=logs/stromboli.err

#eval "$(conda shell.bash hook)"
conda activate tricolor
snakemake

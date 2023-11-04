#!/bin/bash

#SBATCH --account=bioinf593f23_class 
#SBATCH --job-name=stripe_download
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=11
#SBATCH --time=6:00:00
#SBATCH --output=logs/hifi.out
#SBATCH --error=logs/hifi.err




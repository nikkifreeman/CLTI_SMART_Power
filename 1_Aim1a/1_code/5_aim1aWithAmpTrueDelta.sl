#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 03-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=nlbf@live.unc.edu

Rscript 5_aim1aWithAmpTrueDelta.R

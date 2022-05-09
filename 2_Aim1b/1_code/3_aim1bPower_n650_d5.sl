#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 07-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --mail-type=BEGIN, END
#SBATCH --mail-user=nlbf@live.unc.edu

Rscript 3_aim1bPower_n650_d5.R

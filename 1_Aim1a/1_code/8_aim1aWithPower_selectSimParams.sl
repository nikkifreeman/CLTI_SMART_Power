#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1

studySize=$1
p=$2

Rscript 8_aim1aWithAmpPower_selectSimParams.R $studySize $p

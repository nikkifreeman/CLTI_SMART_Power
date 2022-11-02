#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 10-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1

studySize=$1
amp0=$2
amp1=$3
gridRow=$4

Rscript 6_aim1aWithAmpPower.R $studySize $amp0 $amp1 $gridRow

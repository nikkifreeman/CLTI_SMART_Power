#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 07-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=nlbf@live.unc.edu

studySize=$1
dominantRegime=$2
amp0=$3
amp1=$4
gridRow=$5

Rscript 3_aim1bPower.R $studySize $dominantRegime $amp0 $amp1 $gridRow

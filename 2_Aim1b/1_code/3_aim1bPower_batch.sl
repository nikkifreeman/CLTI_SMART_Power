#!/bin/bash

for studySize in {650,750}
do
	for dominantRegime in {1..6}
	do
		for amp0 in {0.2,0.25}
		do
			for amp1 in {0.2,0.25}
			do
			  for gridRow in {1..60}
			  do
				  sbatch 3_aim1bPower.sl $studySize $dominantRegime $amp0 $amp1 $gridRow
				done
			done
		done
	done
done

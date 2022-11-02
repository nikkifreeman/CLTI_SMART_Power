#!/bin/bash

for studySize in {600,650,700}
do
	for amp0 in {0.2,0.25}
	do
		for amp1 in {0.2,0.25}
		do
		  for gridRow in {1..60}
		  do

		sbatch "6_aim1aWithPower.sl" $studySize $amp0 $amp1 $gridRow
		  done
		done
	done
done

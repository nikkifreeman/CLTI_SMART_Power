#!/bin/bash

for studySize in {650, 700}
do
	for amp0 in {0.2, 0.25}
	do
		for amp1 in {0.2, 0.25}
		do

		sbatch "6_aim1aWithAmpPower.sl" $studySize $amp0 $amp1
		done
	done
done

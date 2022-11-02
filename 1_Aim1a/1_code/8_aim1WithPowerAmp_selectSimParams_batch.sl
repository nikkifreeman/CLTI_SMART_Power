#!/bin/bash

for studySize in {600,650,700}
do
	for p in {1..327}
	do
		sbatch "8_aim1aWithPower_selectSimParams.sl" $studySize $p
	done
done
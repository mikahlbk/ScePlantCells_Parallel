#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe smp 1-64
#$ -q  long				# Specify queue
#$ -N  run_Animate9	 	# Specify job name

mkdir Animate9
./program Animate9

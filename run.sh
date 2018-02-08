#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe smp 1-64
#$ -q  long				# Specify queue
#$ -N  run_Animate35	 	# Specify job name

mkdir Animate35
./program Animate35


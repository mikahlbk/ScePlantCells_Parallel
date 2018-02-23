#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe smp 1-64
#$ -q  long			# Specify queue
#$ -N  run_Animate1 	# Specify job name

mkdir Animate1
./program Animate1


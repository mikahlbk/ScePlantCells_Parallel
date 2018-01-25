#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -q  long 				# Specify queue
#$ -N  run_DEC11_Animate5   		 	# Specify job name

mkdir Animate5
./program Animate5

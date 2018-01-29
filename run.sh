#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe mpi-* 24 
#$ -q  debug 				# Specify queue
#$ -N  run_Jan28_Animate15   	 	# Specify job name

mkdir Animate15
./program Animate15

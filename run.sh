#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe mpi-* 24 
#$ -q  debug 				# Specify queue
#$ -N  run_Animate4   	 	# Specify job name

mkdir Animate4
./program Animate4

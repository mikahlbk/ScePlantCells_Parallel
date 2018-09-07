#!/bin/csh
#$ -M  mbkuhn@math.ucr.edu
#$ -m  ab
#$ -pe smp 12
#$ -q  long# Specify queue
#$ -N  run_Animate_CK

setenv OMP_NUM_THREADS 16
mkdir Animate_CK
mkdir Nematic_CK
./program Animate_CK Nematic_CK
		
		

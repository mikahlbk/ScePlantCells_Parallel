#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe smp 12
#$ -q  long			# Specify queue
#$ -N  run_Animate46

setenv OMP_NUM_THREADS 16
mkdir Animate46
./program Animate46



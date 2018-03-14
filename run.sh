#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m ab
#$ -pe smp 12
#$ -q  long		# Specify queue
#$ -N  run_Animate5


setenv OMP_NUM_THREADS 16
mkdir Animate5
./program Animate5

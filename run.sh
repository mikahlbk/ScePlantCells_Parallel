#!/bin/csh
#$ -M mbkuhn@math.ucr.edu
#$ -m ab
#$ -pe smp 12
#$ -q  debug # Specify queue
#$ -N  run_Animate3


setenv OMP_NUM_THREADS 16
mkdir Animate3
mkdir Nematic3
./program Animate3 Nematic3

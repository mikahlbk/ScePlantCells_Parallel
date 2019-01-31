#!/bin/csh
#$ -M  mbkuhn@math.ucr.edu
#$ -m  ab
#$ -pe smp 12
#$ -q  long# Specify queue
#$ -N  test_div_new_files_1

setenv OMP_NUM_THREADS 16
mkdir Animate = Animate_test_1
mkdir Nem = Nematic_test_1
mkdir Loc = Locations_test_1
./program Animate Nem Loc
		
		

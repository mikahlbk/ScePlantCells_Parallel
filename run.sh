#!/bin/csh
#$ -M  mbkuhn@math.ucr.edu
#$ -m  ab
#$ -pe smp 12
#$ -q  debug# Specify queue
#$ -N  test_div_new_files_2

setenv OMP_NUM_THREADS 16
mkdir Animate_test_2
mkdir Nematic_test_2
mkdir Locations_test_2
./program Animate_test_2 Nematic_test_2 Locations_test_2
		
		

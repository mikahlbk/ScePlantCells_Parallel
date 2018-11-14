#!/bin/csh
#$ -M  mbkuhn@math.ucr.edu
#$ -m  ab
#$ -pe smp 12
#$ -q  debug# Specify queue
#$ -N  test_div_new_files_1

setenv OMP_NUM_THREADS 16
mkdir Animate_test_1
mkdir Nematic_test_1
mkdir Locations_test_1
./program Animate_test_1 Nematic_test_1 Locations_test_1
		
		

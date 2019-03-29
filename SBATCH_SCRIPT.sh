#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --output=my22.stdout
#SBATCH --job-name="test_22"
#SBATCH -p batch

set OMP_NUM_THREADS=48
mkdir Animate_test_22
mkdir Nematic_test_22     
mkdir Locations_test_20        
./program Animate_test_22 Locations_test_5 Nematic_test_22
                      


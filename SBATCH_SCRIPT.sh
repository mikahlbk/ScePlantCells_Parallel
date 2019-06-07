#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=16:15:00
#SBATCH --output=my35lowtest.stdout
#SBATCH --job-name="test35lowtest"
#SBATCH -p batch

set OMP_NUM_THREADS=12
mkdir Animate_test35lowtest
mkdir Nematic_test_1    
mkdir Locations_test_1       
./program Animate_test35lowtest Locations_test_1 Nematic_test_1
                    

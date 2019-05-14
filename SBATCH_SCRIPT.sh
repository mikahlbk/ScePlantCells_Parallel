#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=my13.stdout
#SBATCH --job-name="test_13"
#SBATCH -p short

set OMP_NUM_THREADS=48
mkdir Animate_test_13
mkdir Nematic_test_1    
mkdir Locations_test_1       
./program Animate_test_13 Locations_test_1 Nematic_test_1
                    

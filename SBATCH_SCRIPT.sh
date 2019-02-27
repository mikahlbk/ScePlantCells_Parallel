#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=my21.stdout
#SBATCH --job-name="test_21"
#SBATCH -p batch

export OMP_NUM_THREADS 48
mkdir Animate_test_21
mkdir Nematic_test_1        
mkdir Locations_test_1         
./program Animate_test_21 Locations_test_1 Nematic_test_1                                  
                      


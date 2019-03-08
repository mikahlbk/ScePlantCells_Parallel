#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=my6.stdout
#SBATCH --job-name="test_6"
#SBATCH -p batch

export OMP_NUM_THREADS 48
mkdir Animate_test_6
mkdir Nematic_test_1       
mkdir Locations_test_1         
./program Animate_test_6 Locations_test_1 Nematic_test_1                                  
                      


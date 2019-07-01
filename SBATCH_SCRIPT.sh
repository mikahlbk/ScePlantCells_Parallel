#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:02:00
#SBATCH --output=my1.stdout
#SBATCH --job-name="test1"
#SBATCH -p short

set OMP_NUM_THREADS=48
mkdir Animate_test_1
mkdir Nematic_test_1    
mkdir Locations_test_1       
mkdir Animate_No_Cyt_1
./program Animate_test_1 Locations_test_1 Nematic_test_1 Animate_No_Cyt_1
                    

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:10:00
#SBATCH --output=test_locations1.stdout
#SBATCH --job-name="test_locations1"
#SBATCH -p batch

set OMP_NUM_THREADS=12
mkdir Animate_test_locations_1
mkdir Locations_no_cyt_test_locations_1
mkdir Locations_cyt_test_locations_1
./program Animate_test_locations_1 Locations_no_cyt_test_locations_1 Locations_cyt_test_locations_1
                    

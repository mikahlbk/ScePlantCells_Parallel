#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=test_mech_div2.stdout
#SBATCH --job-name="test_mech_div2"
#SBATCH -p batch

set OMP_NUM_THREADS=24
mkdir Animate_cyt_mech_div_2
mkdir Locations_test_mech_div
mkdir Locations_test_mech_div_2
mkdir Animate_No_cyt_mech_div_2
./program Animate_cyt_mech_div_2 Locations_test_mech_div Locations_test_mech_div_2 Animate_No_cyt_mech_div_2
                    

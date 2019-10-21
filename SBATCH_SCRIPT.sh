#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --output=div_test_4.stdout
#SBATCH --job-name="div_test_4"
#SBATCH -p short

set OMP_NUM_THREADS=12
mkdir Animate_div_4
mkdir Locations_div_4
mkdir Locations_NC_div_4
mkdir Animate_NC_div_4
./program Animate_div_4 Animate_NC_div_4 Locations_div_3 Locations_NC_div_3
                   

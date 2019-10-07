#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
<<<<<<< HEAD
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
=======
#SBATCH --time=32:00:00
#SBATCH --output=test_chemical_tissue_params1.stdout
#SBATCH --job-name="test_chemical_tissue_params1"
#SBATCH -p batch

set OMP_NUM_THREADS=12
mkdir Animate_chemical_tissue_params1
mkdir Locations_chemical_tissue_params1
mkdir Locations_chemical_NC_tissue_params1
mkdir Animate_chemical_NC_tissue_params1
./program Animate_chemical_tissue_params1 Animate_chemical_NC_tissue_params1 Locations_chemical_tissue_params1 Locations_chemical_NC_tissue_params1
>>>>>>> origin/chemical_only
                    

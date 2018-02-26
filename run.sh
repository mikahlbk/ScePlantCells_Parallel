#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -pe smp 1-64
#$ -q  long			# Specify queue
#$ -N  run_Animate67


mkdir Animate67
./program Animate67


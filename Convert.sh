#!/bin/bash
qsub <<- BODY
#PBS -l select=1:ncpus=24:mem=192g,walltime=71:0:0 -q dl560g10q@vm-pbs2
#PBS -N $1
python Convert.py
BODY
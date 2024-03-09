#!/bin/sh

#SBATCH --job-name=making_cont    
#SBATCH --error=random_contacts-%j.err        
#SBATCH --output=random_contacts-%j.log        
#SBATCH --time=71:00:00 
#SBATCH --partition=operation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G



source /home/popov/miniconda3/bin/activate
conda activate cooltools

python control_for_loop_prediction.py

#python tads_control.py


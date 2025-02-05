#!/bin/bash
#SBATCH --job-name=cleanup_job
#SBATCH --output=slurm_log.log
#SBATCH --error=slurm_log.log
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

N=4000
B=1756

python3 generate_matrix.py $N $B
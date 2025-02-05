#!/bin/bash
#SBATCH --job-name=cleanup_job
#SBATCH --output=/dev/null
#SBATCH --error=err_gen.log
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

N=4000
B=1756

python3 generate_matrix.py $N $B
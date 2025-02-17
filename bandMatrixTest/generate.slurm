#!/bin/bash
#SBATCH --job-name=cleanup_job
#SBATCH --output=/dev/null
#SBATCH --error=err_gen.log
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

N=2000
B=1200

module load python/3.8
python3 generator.py $N $B
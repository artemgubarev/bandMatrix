#!/bin/bash
#SBATCH --job-name=cleanup_job
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=00:01:00
#SBATCH --ntasks=1

rm -f err* out* bandMatrix comparator solution.txt msolution.txt

echo -e "\e[32mrm SLURM finished$(date)\e[0m"
#!/bin/bash
#SBATCH --job-name=matlab_only
#SBATCH --output=matlab_out.log
#SBATCH --error=matlab_err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00

export INPUT_MATRIX_FILE="testData/matrix2000x1600.txt"
matlab -nodisplay -nosplash -r "run('solve_slae.m'); exit;" > matlab_run.log 2>&1
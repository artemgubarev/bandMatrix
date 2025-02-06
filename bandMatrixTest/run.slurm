#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

export INPUT_MATRIX_FILE="testData/matrix4x4x1.txt"

# module load compiler/gcc/9

g++ -o bandMatrix bandMatrixTest.cpp

srun matlab -nodisplay -nosplash -nodesktop -r "run('solve_slae.m'); exit"
srun ./bandMatrix

echo -e "\e[32mrun SLURM finished$(date)\e[0m"
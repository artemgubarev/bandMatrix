#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
module load mpi/impi/2019.4.243
export INPUT_MATRIX_FILE="testData/matrix2000.txt"
mpicxx -o bandMatrix bandMatrixTest.cpp
srun --mpi=pmix ./bandMatrix
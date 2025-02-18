#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=2              
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

export INPUT_MATRIX_FILE="testData/matrix2000x1200.txt"
export MODE=2
export OMP_NUM_THREADS=4

if [ "$MODE" -eq 0 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp
elif [ "$MODE" -eq 1 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp -fopenmp
elif [ "$MODE" -eq 2 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp
elif [ "$MODE" -eq 3 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp -fopenmp
else
    echo "MODE value incorrect"
    exit 1
fi

srun --mpi=pmi2 -n 4 ./bandMatrix
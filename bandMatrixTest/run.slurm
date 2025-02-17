#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1              
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=2     
#SBATCH --time=01:00:00

module load mpi/impi/2021.4.0
module load mpi/openmpi/1.10.5/gcc/6.3.0

export INPUT_MATRIX_FILE="testData/matrix2000x1200.txt"
export MODE=0
export OMP_NUM_THREADS=2

if [ "$MODE" -eq 0 ]; then
    g++ -o bandMatrix bandMatrixTest.cpp
elif [ "$MODE" -eq 1 ]; then
    g++ -o bandMatrix bandMatrixTest.cpp -fopenmp
elif [ "$MODE" -eq 2 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp
elif [ "$MODE" -eq 3 ]; then
    mpicxx -o bandMatrix bandMatrixTest.cpp -fopenmp
else
    echo "MODE value incorrect"
    exit 1
fi

srun ./bandMatrix
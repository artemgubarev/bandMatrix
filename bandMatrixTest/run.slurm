#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1              
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=2     
#SBATCH --time=01:00:00

export INPUT_MATRIX_FILE="testData/matrix2000x1600.txt"
export MODE=2

export OMP_NUM_THREADS=2

module load mpi

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

if [ "$MODE" -eq 0 ]; then
    srun --cpu-bind=cores ./bandMatrix
elif [ "$MODE" -eq 1 ]; then
    srun --cpu-bind=cores ./bandMatrix
elif [ "$MODE" -eq 2 ]; then
    srun --ntasks=4 --cpu-bind=cores mpirun -np 4 ./bandMatrix
elif [ "$MODE" -eq 3 ]; then
    export OMP_NUM_THREADS=2
    srun --ntasks=2 --cpus-per-task=2 --cpu-bind=cores mpirun -np 2 ./bandMatrix
fi
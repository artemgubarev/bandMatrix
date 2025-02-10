#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=1              
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=2     
#SBATCH --time=01:00:00

export INPUT_MATRIX_FILE="testData/matrix2000x1600.txt"
# 0 - Serial version
# 1 - OpenMP
export MODE=1

export OMP_NUM_THREADS=2   
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic,4"
export OMP_STACKSIZE=512M
export OMP_WAIT_POLICY=passive

matlab -nodisplay -nosplash -r "run('solve_slae.m'); exit;" > /dev/null 2>&1

g++ -o bandMatrix bandMatrixTest.cpp -fopenmp

srun --cpu-bind=cores ./bandMatrix
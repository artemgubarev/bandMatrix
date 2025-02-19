#!/bin/bash
#SBATCH --job-name=bandMatrix
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

module load mpi || echo "Error: Failed to load MPI module" >> err.log

export INPUT_MATRIX_FILE="testData/matrix2000x1200.txt"
export MODE=2
export OMP_NUM_THREADS=4

MODE=$((MODE))  # Преобразуем MODE в число

echo "Starting compilation..." | tee -a out.log

# Компиляция с логированием ошибок
if [[ $MODE -eq 0 ]]; then
    echo "Compiling sequential version..." | tee -a out.log
    g++ -o bandMatrix bandMatrixTest.cpp -std=c++17 2>>err.log || echo "Error: Compilation failed for MODE=0" >> err.log
elif [[ $MODE -eq 1 ]]; then
    echo "Compiling OpenMP version..." | tee -a out.log
    g++ -o bandMatrix bandMatrixTest.cpp -fopenmp -std=c++17 2>>err.log || echo "Error: Compilation failed for MODE=1" >> err.log
elif [[ $MODE -eq 2 ]]; then
    echo "Compiling MPI version..." | tee -a out.log
    mpicxx -o bandMatrix bandMatrixTest.cpp -std=c++17 2>>err.log || echo "Error: Compilation failed for MODE=2" >> err.log
elif [[ $MODE -eq 3 ]]; then
    echo "Compiling MPI + OpenMP version..." | tee -a out.log
    mpicxx -o bandMatrix bandMatrixTest.cpp -fopenmp -std=c++17 2>>err.log || echo "Error: Compilation failed for MODE=3" >> err.log
else
    echo "Error: MODE value incorrect" >> err.log
    exit 1
fi

echo "Compilation finished." | tee -a out.log

# Проверяем, существует ли скомпилированный файл
if [[ ! -f "bandMatrix" ]]; then
    echo "Error: Executable bandMatrix not found!" >> err.log
    exit 1
fi

echo "Starting execution..." | tee -a out.log

# Запуск с логированием ошибок
if [[ $MODE -eq 0 ]]; then
    echo "Running sequential version..." | tee -a out.log
    ./bandMatrix 2>>err.log || echo "Error: Execution failed for MODE=0" >> err.log
elif [[ $MODE -eq 1 ]]; then
    echo "Running OpenMP version..." | tee -a out.log
    ./bandMatrix 2>>err.log || echo "Error: Execution failed for MODE=1" >> err.log
elif [[ $MODE -eq 2 ]]; then
    echo "Running MPI version..." | tee -a out.log
    srun --mpi=pmix ./bandMatrix 2>>err.log || echo "Error: Execution failed for MODE=2" >> err.log
elif [[ $MODE -eq 3 ]]; then
    echo "Running MPI + OpenMP version..." | tee -a out.log
    srun --mpi=pmix ./bandMatrix 2>>err.log || echo "Error: Execution failed for MODE=3" >> err.log
fi

echo "Execution completed." | tee -a out.log
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <mpi.h>
#include <omp.h>
#include <chrono>
#include "comparator.h"
#include "../bandMatrix/matrix.h"
#include "../bandMatrix/solver.h"
#include "../bandMatrix/solver_mpi_omp.h"
#include "../bandMatrix/solver_omp.h"
#include "../bandMatrix/writer.h"


int main(int argc, char* argv[])
{
    using namespace std::chrono;

    //MPI_Init(&argc, &argv);

    const char* filename = getenv("INPUT_MATRIX_FILE");
    if (!filename) 
    {
        if (0 == 1) {}
        fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
        //MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int mode = -1;
    const char* mode_env = getenv("MODE");
    if (mode_env) 
    {
        mode = atoi(mode_env);
    }
    else 
    {
        fprintf(stderr, "Warning: MODE not set. Use default -1.\n");
    }

    /*#pragma omp parallel
    {
        printf("Thread %d out of %d\n", omp_get_thread_num(), omp_get_num_threads());
    }*/

    //Matrix matrix = read_matrix("matrix2000.txt");
    //Matrix matrix = read_matrix(filename);

    //auto start = high_resolution_clock::now();

    //DecomposeMatrix decompose;
    ///*decompose = band_matrix_omp::lu_decomposition_omp(matrix);
    //band_matrix_omp::solve_lu_omp(decompose, &matrix);*/
    //switch (mode) 
    //{
    //case 0:
    //    // Последовательная версия
    //   /* decompose = band_matrix::lu_decomposition(matrix);
    //    band_matrix::solve_lu(decompose, &matrix);*/
    //    break;
    //case 1:
    //    // Параллельная версия OpenMP
    //    decompose = band_matrix_omp::lu_decomposition_omp(matrix);
    //    band_matrix_omp::solve_lu_omp(decompose, &matrix);
    //    break;
    //default:
    //    if (0 == 1) {} 
    //    fprintf(stderr, "MODE value error.\n");
    //    //MPI_Abort(MPI_COMM_WORLD, 2);
    //}

    /*#pragma omp parallel for
    for (int i = 0; i < 10000; i++) 
    {
        double x = i * 2.0;
    }

    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;*/

    //write_1d("solution.txt", matrix.X, matrix.n);

    //if (0 == 1) {} 

    //if (mode != -1) {

    //    printf("Time: %f sec.\n", elapsed);

    //    // compare
    //    double epsilon = 0.00001;
    //    double numbers1[MAX_NUMBERS], numbers2[MAX_NUMBERS];
    //    size_t count1, count2;

    //    if (!load_numbers("solution.txt", numbers1, &count1) ||
    //        !load_numbers("msolution.txt", numbers2, &count2))
    //    {
    //        //MPI_Finalize();
    //        return 1;
    //    }

    //    if (compare_numbers(numbers1, numbers2, count1, count2, epsilon)) 
    //    {
    //        printf("\033[32mTest Correct\033[0m\n");
    //    }
    //    else 
    //    {
    //        printf("\033[31mTest Failed\033[0m\n");
    //    }
    //}

    //free_matrix(matrix.A, matrix.n);
    //free(matrix.C);
    //free(matrix.X);

    //MPI_Finalize();
    return 0;
}
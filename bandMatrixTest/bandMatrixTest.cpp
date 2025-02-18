//#define HAVE_STRUCT_TIMESPEC
//#define _NO_DEBUG_HEAP 1
//
//#include <cstdlib>
//#include <stdio.h>
//#include <stdlib.h>
//#include <cstring>
//#include <time.h>
//#include <mpi.h>
//#include <omp.h>
//#include "comparator.h"
//#include "../bandMatrix/matrix.h"
//#include "../bandMatrix/solver_serial.h"
//#include "../bandMatrix/solver_mpi.h"
//#include "../bandMatrix/solver_omp.h"
//#include "../bandMatrix/solver_mpi_omp.h"
//#include "../bandMatrix/solver_pthreads.h"
//#include "../bandMatrix/writer.h"
//
//#ifdef _WIN32
//#include <windows.h>
//#else
//#include <sys/time.h>
//#endif
//
//double get_time() {
//#ifdef _WIN32
//    LARGE_INTEGER frequency, start;
//    QueryPerformanceFrequency(&frequency);
//    QueryPerformanceCounter(&start);
//    return (double)start.QuadPart / frequency.QuadPart;
//#else
//    struct timespec t;
//    clock_gettime(CLOCK_MONOTONIC, &t);
//    return t.tv_sec + t.tv_nsec * 1e-9;
//#endif
//}
//
//void get_output_filename(const char* input_file, char* output_filename, size_t size) 
//{
//    const char* filename = strrchr(input_file, '/');
//    filename = (filename) ? filename + 1 : input_file;
//
//    const char* name_part = strstr(filename, "matrix");
//    if (name_part) {
//        name_part += 6;
//    }
//    else {
//        name_part = filename;
//    }
//
//    char name_only[256];
//    strncpy(name_only, name_part, sizeof(name_only) - 1);
//    name_only[sizeof(name_only) - 1] = '\0';
//
//    char* dot = strrchr(name_only, '.');
//    if (dot) {
//        *dot = '\0';
//    }
//
//    snprintf(output_filename, size, "matlabSolutions/msolution%s.txt", name_only);
//}
//
//int main(int argc, char* argv[])
//{
//    const char* filename = getenv("INPUT_MATRIX_FILE");
//    if (!filename) 
//    {
//        if (0 == 1) {}
//        fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
//    }
//
//    int mode = -1;
//    const char* mode_env = getenv("MODE");
//    if (mode_env) 
//    {
//        mode = atoi(mode_env);
//    }
//    else 
//    {
//        fprintf(stderr, "Warning: MODE not set. Use default -1.\n");
//    }
//
//  /*  const char* filename = "matrix2000.txt";
//    int mode = 2;*/
//
//    if (mode == 2 || mode == 3)
//    {
//        MPI_Init(&argc, &argv);
//    }
//
//    int rank = 0, size = 1;
//    if (mode == 2 || mode == 3)
//    {
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//    }
//
//    // ✅ Теперь можно безопасно замерять MPI-время
//    double start_time = get_time();
//    double start_time_mpi = (mode == 2 || mode == 3) ? MPI_Wtime() : 0.0;
//
//    DecomposeMatrix decompose;
//
//    if (mode == 0 || mode == 1)
//    {
//        Matrix matrix = read_matrix(filename);
//        if (mode == 0) // Serial
//        {
//            decompose = band_matrix_serial::lu_decomposition(matrix);
//            band_matrix_serial::solve_lu(decompose, &matrix);
//        }
//        if (mode == 1)// OpenMP
//        {
//            decompose = band_matrix_omp::lu_decomposition_omp(matrix);
//            band_matrix_omp::solve_lu_omp(decompose, &matrix);
//        }
//        write_1d("solution.txt", matrix.X, matrix.n);
//
//        free_matrix(matrix.A, matrix.n);
//        free(matrix.C);
//        free(matrix.X);
//    }
//    else if (mode == 2 || mode == 3)
//    {
//        Matrix matrix;
//        if (rank == 0)
//        {
//            matrix = read_matrix(filename);
//        }
//        else
//        {
//            matrix.n = 0;
//            matrix.b = 0;
//            matrix.A = NULL;
//            matrix.C = NULL;
//            matrix.X = NULL;
//        }
//
//        if (mode == 2) // MPI
//        {
//            decompose = band_matrix_mpi::lu_decomposition_mpi(matrix);
//            band_matrix_mpi::solve_lu_mpi(decompose, &matrix);
//        }
//        else  // MPI + OpenMP
//        {
//            decompose = band_matrix_mpi_omp::lu_decomposition_mpi_omp(matrix);
//            band_matrix_mpi_omp::solve_lu_mpi_omp(decompose, &matrix);
//        }
//
//        if (rank == 0)
//        {
//            write_1d("solution.txt", matrix.X, matrix.n);
//        }
//
//        free_matrix(matrix.A, matrix.n);
//        free(matrix.C);
//        free(matrix.X);
//
//       
//    }
//
//    double end_time_mpi = (mode == 2 || mode == 3) ? MPI_Wtime() : 0.0;
//    double end_time = get_time();
//
//    if (mode != -1) {
//
//        if (rank == 0) {  // ✅ Выводим только на Rank 0
//            if (mode == 2 || mode == 3)
//            {
//                printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
//            }
//            else
//            {
//                printf("Time: %f sec.\n", end_time - start_time);
//            }
//
//            // compare
//            double epsilon = 0.00001;
//            double* numbers1 = new double[MAX_NUMBERS];
//            double* numbers2 = new double[MAX_NUMBERS];
//            size_t count1, count2;
//
//            char output_filename[512];
//            get_output_filename(filename, output_filename, sizeof(output_filename));
//
//            if (!load_numbers("solution.txt", numbers1, &count1) ||
//                !load_numbers(output_filename, numbers2, &count2))
//            {
//                return 1;
//            }
//
//            if (compare_numbers(numbers1, numbers2, count1, count2, epsilon))
//            {
//                printf("\033[32mTest Correct\033[0m\n");
//            }
//            else
//            {
//                printf("\033[31mTest Failed\033[0m\n");
//            }
//
//            delete[] numbers1;
//            delete[] numbers2;
//        }
//    }
//
//    if (mode == 2 || mode == 3)
//    {
//        MPI_Finalize();
//    }
//
//    return 0;
//}

#define HAVE_STRUCT_TIMESPEC
#define _NO_DEBUG_HEAP 1

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include "comparator.h"
#include "../bandMatrix/matrix.h"
#include "../bandMatrix/solver_mpi.h"  // новый заголовочный файл с параллельными функциями
#include "../bandMatrix/writer.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

// Функция замера времени (используем MPI_Wtime для согласованности)
double get_time() {
#ifdef _WIN32
    LARGE_INTEGER frequency, start;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);
    return (double)start.QuadPart / frequency.QuadPart;
#else
    return MPI_Wtime();
#endif
}

void get_output_filename(const char* input_file, char* output_filename, size_t size)
{
    const char* filename = strrchr(input_file, '/');
    filename = (filename) ? filename + 1 : input_file;

    const char* name_part = strstr(filename, "matrix");
    if (name_part) {
        name_part += 6;
    }
    else {
        name_part = filename;
    }

    char name_only[256];
    strncpy(name_only, name_part, sizeof(name_only) - 1);
    name_only[sizeof(name_only) - 1] = '\0';

    char* dot = strrchr(name_only, '.');
    if (dot) {
        *dot = '\0';
    }

    snprintf(output_filename, size, "matlabSolutions/msolution%s.txt", name_only);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const char* filename = getenv("INPUT_MATRIX_FILE");
    if (!filename)
    {
        if (rank == 0)
            fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
        MPI_Finalize();
        return 1;
    }

    double start_time = get_time();

    // Rank 0 читает всю матрицу
    Matrix full_matrix;
    if (rank == 0) {
        full_matrix = read_matrix(filename);
    }

    // Передаём размеры матрицы всем процессам
    int n, b;
    if (rank == 0) {
        n = full_matrix.n;
        b = full_matrix.b;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Определяем, сколько строк получит каждый процесс
    int base = n / nprocs;
    int rem = n % nprocs;
    int local_n = (rank < rem) ? base + 1 : base;
    int global_start = (rank < rem) ? rank * (base + 1) : rem * (base + 1) + (rank - rem) * base;

    // Распределяем строки матрицы A по процессам.
    // На Rank 0 матрица хранится как массив указателей на строки, скопируем её в один непрерывный буфер.
    double* local_A = (double*)malloc(local_n * n * sizeof(double));
    double* sendbuf = NULL;
    int* sendcounts = NULL;
    int* displs = NULL;
    if (rank == 0) {
        sendbuf = (double*)malloc(n * n * sizeof(double));
        for (int i = 0; i < n; i++) {
            memcpy(&sendbuf[i * n], full_matrix.A[i], n * sizeof(double));
        }
        sendcounts = new int[nprocs];
        displs = new int[nprocs];
        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int rows = (p < rem) ? base + 1 : base;
            sendcounts[p] = rows * n;
            displs[p] = offset;
            offset += rows * n;
        }
    }
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        free(sendbuf);
        delete[] sendcounts;
        delete[] displs;
    }

    // Распространяем вектор C аналогичным образом
    double* local_C = (double*)malloc(local_n * sizeof(double));
    double* sendC = NULL;
    int* sendcounts_C = NULL;
    int* displs_C = NULL;
    if (rank == 0) {
        sendC = full_matrix.C; // C уже хранится непрерывно
        sendcounts_C = new int[nprocs];
        displs_C = new int[nprocs];
        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int rows = (p < rem) ? base + 1 : base;
            sendcounts_C[p] = rows;
            displs_C[p] = offset;
            offset += rows;
        }
    }
    MPI_Scatterv(sendC, sendcounts_C, displs_C, MPI_DOUBLE, local_C, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        delete[] sendcounts_C;
        delete[] displs_C;
    }

    // Формируем локальную структуру матрицы для параллельных вычислений
    LocalMatrix localMat;
    localMat.local_n = local_n;
    localMat.global_start = global_start;
    localMat.n = n;
    localMat.b = b;
    localMat.A = local_A;
    localMat.C = local_C;
    localMat.X = (double*)malloc(local_n * sizeof(double)); // сюда запишется решение для локальных строк

    // Выполняем параллельное LU‑разложение и решаем систему
    band_matrix_mpi::lu_decomposition_mpi(localMat);
    band_matrix_mpi::solve_lu_mpi(localMat);

    // Собираем решение X с разных процессов в глобальный вектор на Rank 0
    double* global_X = NULL;
    int* recvcounts = NULL;
    int* recvdispls = NULL;
    if (rank == 0) {
        global_X = (double*)malloc(n * sizeof(double));
        recvcounts = new int[nprocs];
        recvdispls = new int[nprocs];
        int offset = 0;
        for (int p = 0; p < nprocs; p++) {
            int rows = (p < rem) ? base + 1 : base;
            recvcounts[p] = rows;
            recvdispls[p] = offset;
            offset += rows;
        }
    }
    MPI_Gatherv(localMat.X, local_n, MPI_DOUBLE, global_X, recvcounts, recvdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        write_1d("solution.txt", global_X, n);
        free(global_X);
        delete[] recvcounts;
        delete[] recvdispls;
    }

    free(local_A);
    free(local_C);
    free(localMat.X);
    if (rank == 0) {
        free_matrix(full_matrix.A, full_matrix.n);
        free(full_matrix.C);
        free(full_matrix.X);
    }

    double end_time = get_time();
    if (rank == 0) {
        printf("Time: %f sec.\n", end_time - start_time);

        // Сравнение решения (как в оригинале)
        double epsilon = 0.00001;
        double* numbers1 = new double[MAX_NUMBERS];
        double* numbers2 = new double[MAX_NUMBERS];
        size_t count1, count2;
        char output_filename[512];
        get_output_filename(filename, output_filename, sizeof(output_filename));

        if (!load_numbers("solution.txt", numbers1, &count1) ||
            !load_numbers(output_filename, numbers2, &count2))
        {
            return 1;
        }
        if (compare_numbers(numbers1, numbers2, count1, count2, epsilon))
        {
            printf("\033[32mTest Correct\033[0m\n");
        }
        else
        {
            printf("\033[31mTest Failed\033[0m\n");
        }
        delete[] numbers1;
        delete[] numbers2;
    }

    MPI_Finalize();
    return 0;
}

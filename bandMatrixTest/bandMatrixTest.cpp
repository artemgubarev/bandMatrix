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
#include "../bandMatrix/solver_serial.h"
#include "../bandMatrix/solver_mpi.h"
#include "../bandMatrix/solver_omp.h"
#include "../bandMatrix/solver_mpi_omp.h"
#include "../bandMatrix/solver_pthreads.h"
#include "../bandMatrix/writer.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

double get_time() {
#ifdef _WIN32
	LARGE_INTEGER frequency, start;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	return (double)start.QuadPart / frequency.QuadPart;
#else
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + t.tv_nsec * 1e-9;
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

double** allocate_matrix_mpi(size_t n) {
	double** matrix = (double**)malloc(n * sizeof(double*));
	if (!matrix) {
		perror("Error allocating matrix pointers");
		exit(EXIT_FAILURE);
	}
	matrix[0] = (double*)malloc(n * n * sizeof(double));
	if (!matrix[0]) {
		perror("Error allocating matrix memory");
		exit(EXIT_FAILURE);
	}
	for (size_t i = 1; i < n; i++) {
		matrix[i] = matrix[0] + i * n;
	}
	return matrix;
}

Matrix read_matrix_mpi(const char* filename) {
	Matrix mat;
	FILE* file = fopen(filename, "r");
	if (!file) {
		perror("Error opening file");
		exit(EXIT_FAILURE);
	}

	if (fscanf(file, "%zu", &mat.n) != 1 || fscanf(file, "%zu", &mat.b) != 1) {
		perror("Error reading matrix dimensions");
		exit(EXIT_FAILURE);
	}

	mat.A = allocate_matrix(mat.n);
	mat.C = (double*)malloc(mat.n * sizeof(double));
	if (!mat.C) {
		perror("Error allocating vector C");
		exit(EXIT_FAILURE);
	}

	for (size_t i = 0; i < mat.n; i++) {
		for (size_t j = 0; j < mat.n; j++) {
			if (fscanf(file, "%lf", &mat.A[i][j]) != 1) {
				perror("Error reading matrix data");
				exit(EXIT_FAILURE);
			}
		}
	}

	for (size_t i = 0; i < mat.n; i++) {
		if (fscanf(file, "%lf", &mat.C[i]) != 1) {
			perror("Error reading vector C");
			exit(EXIT_FAILURE);
		}
	}

	fclose(file);
	return mat;
}

int main(int argc, char* argv[])
{
	const char* filename = getenv("INPUT_MATRIX_FILE");
	if (!filename) 
	{
		fprintf(stderr, "Error: environment variable INPUT_MATRIX_FILE not set.\n");
		return 1;
	}

	//int mode = -1;
	//const char* mode_env = getenv("MODE");
	//if (mode_env) 
	//{
	//	mode = atoi(mode_env);
	//}

	//int rank, size;
	//if (mode == 2 || mode == 3) 
	//{
	//	MPI_Init(&argc, &argv);
	//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//	printf("Process %d of %d started\n", rank, size);
	//	fflush(stdout);  // Важно для немедленного вывода
	//}

	//double start_time = get_time();
	//double start_time_mpi = (mode == 2 || mode == 3) ? MPI_Wtime() : 0.0;

	//DecomposeMatrix decompose;

	//if (mode == 0 || mode == 1) {
	//	// Последовательная или OpenMP версия - без изменений
	//	Matrix matrix = read_matrix(filename);
	//	if (mode == 0) {
	//		decompose = band_matrix_serial::lu_decomposition(matrix);
	//		band_matrix_serial::solve_lu(decompose, &matrix);
	//	}
	//	if (mode == 1) {
	//		decompose = band_matrix_omp::lu_decomposition_omp(matrix);
	//		band_matrix_omp::solve_lu_omp(decompose, &matrix);
	//	}
	//	write_1d("solution.txt", matrix.X, matrix.n);

	//	free_matrix(matrix.A, matrix.n);
	//	free(matrix.C);
	//	free(matrix.X);
	//}
	//else if (mode == 2 || mode == 3) {
	//	Matrix matrix;
	//	size_t n, b;

	//	// Чтение матрицы только в корневом процессе
	//	if (rank == 0) {
	//		matrix = read_matrix(filename);
	//		n = matrix.n;
	//		b = matrix.b;
	//	}

	//	// Broadcast размеров матрицы всем процессам
	//	MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	//	MPI_Bcast(&b, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

	//	// Выделение памяти в некорневых процессах
	//	if (rank != 0) {
	//		matrix.n = n;
	//		matrix.b = b;
	//		matrix.A = (double**)malloc(n * sizeof(double*));
	//		for (size_t i = 0; i < n; i++) {
	//			matrix.A[i] = (double*)malloc(n * sizeof(double));
	//		}
	//		matrix.C = (double*)malloc(n * sizeof(double));
	//		matrix.X = NULL;
	//	}

	//	// Распределение данных между процессами
	//	for (size_t i = 0; i < n; i++) {
	//		MPI_Bcast(matrix.A[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//	}
	//	MPI_Bcast(matrix.C, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//	if (mode == 2) {
	//		decompose = band_matrix_mpi::lu_decomposition_mpi(matrix, rank, size);
	//		band_matrix_mpi::solve_lu_mpi(decompose, &matrix, rank, size);
	//	}
	//	else {
	//		decompose = band_matrix_mpi_omp::lu_decomposition_mpi_omp(matrix);
	//		band_matrix_mpi_omp::solve_lu_mpi_omp(decompose, &matrix);
	//	}

	//	// Сбор результатов в корневом процессе
	//	if (rank == 0) {
	//		write_1d("solution.txt", matrix.X, matrix.n);
	//	}

	//	// Освобождение памяти
	//	for (size_t i = 0; i < n; i++) {
	//		free(matrix.A[i]);
	//	}
	//	free(matrix.A);
	//	free(matrix.C);
	//	if (matrix.X != NULL) {
	//		free(matrix.X);
	//	}
	//}

	//double end_time_mpi = (mode == 2 || mode == 3) ? MPI_Wtime() : 0.0;
	//double end_time = get_time();

	//// Вывод результатов и сравнение только в корневом процессе
	//if (mode != -1 && rank == 0) {
	//	if (mode == 2 || mode == 3) {
	//		printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
	//	}
	//	else {
	//		printf("Time: %f sec.\n", end_time - start_time);
	//	}

	//	// Сравнение результатов
	//	double epsilon = 0.00001;
	//	double* numbers1 = new double[MAX_NUMBERS];
	//	double* numbers2 = new double[MAX_NUMBERS];
	//	size_t count1, count2;

	//	char output_filename[512];
	//	get_output_filename(filename, output_filename, sizeof(output_filename));

	//	if (!load_numbers("solution.txt", numbers1, &count1) ||
	//		!load_numbers(output_filename, numbers2, &count2)) {
	//		return 1;
	//	}

	//	if (compare_numbers(numbers1, numbers2, count1, count2, epsilon)) {
	//		printf("\033[32mTest Correct\033[0m\n");
	//	}
	//	else {
	//		printf("\033[31mTest Failed\033[0m\n");
	//	}

	//	delete[] numbers1;
	//	delete[] numbers2;
	//}

	//if (mode == 2 || mode == 3) {
	//	MPI_Finalize();
	//}


	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	Matrix matrix;
	size_t n = 0, b = 0;

	if (rank == 0) {
		matrix = read_matrix_mpi(filename);
		//matrix = read_matrix_mpi("C:\\Users\\artem\\source\\repos\\BandMatrixPython\\matrix2000.txt");
		n = matrix.n;
		b = matrix.b;
	}

	// Рассылаем параметры
	MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, comm);
	MPI_Bcast(&b, 1, MPI_UNSIGNED_LONG, 0, comm);

	if (rank != 0) {
		matrix.n = n;
		matrix.b = b;
		matrix.A = allocate_matrix_mpi(n);
		matrix.C = (double*)malloc(n * sizeof(double));
		if (!matrix.C) {
			perror("Error allocating vector C on non-root process");
			exit(EXIT_FAILURE);
		}
	}

	// Рассылаем данные только если память выделена
	MPI_Bcast(matrix.A[0], n * n, MPI_DOUBLE, 0, comm);
	MPI_Bcast(matrix.C, n, MPI_DOUBLE, 0, comm);

	double start_time_mpi = MPI_Wtime();

	DecomposeMatrix decomp = band_matrix_mpi::lu_decomposition(matrix, comm);
	band_matrix_mpi::solve_lu(decomp, &matrix, comm);

	double end_time_mpi = MPI_Wtime();
	if (rank == 0) 
	{
		printf("Time: %f sec.\n", end_time_mpi - start_time_mpi);
	}

	free(matrix.A[0]);
	free(matrix.A);
	free(matrix.C);
	free(matrix.X);

	for (size_t i = 0; i < n; i++) {
		free(decomp.l[i]);
		free(decomp.u[i]);
	}
	free(decomp.l);
	free(decomp.u);

	MPI_Finalize();

	return 0;
}


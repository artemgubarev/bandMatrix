#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_mpi_omp
{
	void reverse_array(double* array, int n) {
		for (int i = 0; i < n / 2; i++) {
			double tmp = array[i];
			array[i] = array[n - 1 - i];
			array[n - 1 - i] = tmp;
		}
	}

    DecomposeMatrix lu_decomposition(Matrix matrix) {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        size_t n = matrix.n;
        size_t b = matrix.b;

        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (size_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

#pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)n; i++) {
            for (int j = 0; j < (int)n; j++) {
                result.u[i][j] = matrix.A[i][j];
            }
        }

#pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)n; i++) {
            result.l[i][i] = 1.0;
        }

        size_t start_row = (n * rank) / size;
        size_t end_row = (n * (rank + 1)) / size;

        for (size_t k = 0; k < n - 1; k++) {
            size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

            int owner = (k * size) / n;

            int segment_length = (int)(upper_bound - k);
            double* u_row_segment = (double*)malloc(segment_length * sizeof(double));
            if (rank == owner) {
                for (size_t j = k; j < upper_bound; j++) {
                    u_row_segment[j - k] = result.u[k][j];
                }
            }
            if (u_row_segment[0] == 0.0) {
                if (rank == 0) fprintf(stderr, "Zero pivot encountered at row %zu\n", k);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            MPI_Request req;
            MPI_Ibcast(u_row_segment, segment_length, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);

            size_t local_start = (k + 1 > start_row) ? k + 1 : start_row;
            int loop_end = (end_row < upper_bound) ? (int)end_row : (int)upper_bound;
            int local_start_int = (int)local_start;
#pragma omp parallel for schedule(static)
            for (int i = local_start_int; i < loop_end; i++) {
                result.l[i][k] = result.u[i][k] / u_row_segment[0];
                for (int j = (int)k; j < (int)upper_bound; j++) {
                    result.u[i][j] -= result.l[i][k] * u_row_segment[j - (int)k];
                }
            }
            free(u_row_segment);
        }
        return result;
    }

    void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix) {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        size_t n = matrix->n;

        size_t start_row = (n * rank) / size;
        size_t end_row = (n * (rank + 1)) / size;

        double* y = (double*)malloc(n * sizeof(double));

        for (size_t i = 0; i < n; i++) {
            if (i >= start_row && i < end_row) {
                double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
                for (int j = 0; j < (int)i; j++) {
                    s += decompose_matrix.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }
            int owner = (i * size) / n;
            MPI_Request req;
            MPI_Ibcast(&y[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }

        if (matrix->X == NULL) {
            matrix->X = (double*)malloc(n * sizeof(double));
        }
        for (int i = n - 1; i >= 0; i--) {
            if (i >= (int)start_row && i < (int)end_row) {
                double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
                for (int j = i + 1; j < (int)n; j++) {
                    s += decompose_matrix.u[i][j] * matrix->X[j];
                }
                if (decompose_matrix.u[i][i] == 0.0) {
                    if (rank == 0) fprintf(stderr, "Zero pivot in back substitution at row %d\n", i);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
            }
            int owner = (i * size) / n;
            MPI_Request req;
            MPI_Ibcast(&matrix->X[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, MPI_STATUS_IGNORE);
        }

        free(y);
    }
}
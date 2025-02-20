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

	/* ������������ LU-���������� ��������� ������� � ��������� MPI+OpenMP ���������������.
	   ������ �������������� ������ ����� ���������� MPI, � ���������� ����� ���������� ������������������ OpenMP. */
	DecomposeMatrix lu_decomposition(Matrix matrix) {
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		int n = matrix.n;
		int b = matrix.b;

		DecomposeMatrix result;
		result.l = (double**)malloc(n * sizeof(double*));
		result.u = (double**)malloc(n * sizeof(double*));
		for (int i = 0; i < n; i++) {
			result.l[i] = (double*)calloc(n, sizeof(double));
			result.u[i] = (double*)malloc(n * sizeof(double));
		}

		/* �������� �������� ������� A � U.
		   ���� ���� ����� �������������� �� ������� � OpenMP. */
#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				result.u[i][j] = matrix.A[i][j];
			}
		}

		/* ������������ �������� L ����� 1 */
#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++) {
			result.l[i][i] = 1.0;
		}

		// ���������� �������� �����, �������������� ������ ��������� MPI
		int start_row = (n * rank) / size;
		int end_row = (n * (rank + 1)) / size;

		// �������� ���� ������������
		for (int k = 0; k < n - 1; k++) {
			// ���������� ����������� ������ � �������� �����: �� k �� upper_bound
			int upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

			// ���������� �������-��������� ������ k
			int owner = (k * size) / n;

			// ����� ��� �������� ����� ������ k (�� ������� k �� upper_bound)
			int segment_length = upper_bound - k;
			double* u_row_segment = (double*)malloc(segment_length * sizeof(double));
			if (rank == owner) {
				for (int j = k; j < upper_bound; j++) {
					u_row_segment[j - k] = result.u[k][j];
				}
			}
			// ���������� ������������� Bcast ��� �������� ��������
			MPI_Request req;
			MPI_Ibcast(u_row_segment, segment_length, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);

			// ������ ������� ��������� ���� ������, ���������� � �������� [max(k+1, start_row), end_row)
			int local_start = (k + 1 > start_row) ? k + 1 : start_row;
			int loop_end = (end_row < upper_bound) ? (int)end_row : (int)upper_bound;
			int local_start_int = (int)local_start;
#pragma omp parallel for schedule(static)
			for (int i = local_start_int; i < loop_end; i++) {
				// ...
				result.l[i][k] = result.u[i][k] / u_row_segment[0];
				for (int j = k; j < upper_bound; j++) {
					result.u[i][j] -= result.l[i][k] * u_row_segment[j - k];
				}
			}
			free(u_row_segment);
			// ����� ������ �� ����� � ������������� ���������� ��� MPI_Ibcast.
		}
		return result;
	}

	/* ������������ ������� ������� ������� LU-���������� � ��������� ���������������.
	   ������ � �������� ���� ������������������ �� ������� MPI, � ���������� ������������ ���������� OpenMP.
	*/
	void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix) {
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		int n = matrix->n;

		// ���������� ������������� ����� �� ��������� MPI
		int start_row = (n * rank) / size;
		int end_row = (n * (rank + 1)) / size;

		double* y = (double*)malloc(n * sizeof(double));

		/* ������ ���: ������ L*y = C.
		   ��� ������ ������ i:
			 - ���� ������ ����������� �������� ��������, ��������� y[i] = C[i] - sum(L[i][j]*y[j])
			   � ������������ �������������.
			 - ����� � ������� �������������� Bcast ������� y[i] ���� ���������.
		*/
		for (int i = 0; i < n; i++) {
			if (i >= start_row && i < end_row) {
				double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
				for (int j = 0; j < i; j++) {
					s += decompose_matrix.l[i][j] * y[j];
				}
				y[i] = matrix->C[i] - s;
			}
			int owner = (i * size) / n;  // �������-�������� ������ i
			MPI_Request req;
			MPI_Ibcast(&y[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}

		/* �������� ���: ������ U*x = y.
		   �������� ��� �� ��������� ������ � ������.
		   ��� ������ ������ i:
			 - ���� ������ ����������� ��������, ��������� x[i] = (y[i] - sum(U[i][j]*x[j])) / U[i][i]
			   � ������������ �������������.
			 - ����� ��������� x[i] ���� ���������.
		*/
		if (matrix->X == NULL) {
			matrix->X = (double*)malloc(n * sizeof(double));
		}
		for (int i = n - 1; i >= 0; i--) {
			if (i >= start_row && i < end_row) {
				double s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
				for (int j = i + 1; j < n; j++) {
					s += decompose_matrix.u[i][j] * matrix->X[j];
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
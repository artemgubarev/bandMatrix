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
    //=============================================================================
    // ��������� ������� LU-���������� ��������� (���������) ������� � �������������� MPI + OpenMP.
    // ������:
    //  1) Rank 0 ������ � ������ ������ ������� A.
    //  2) Rank 0 ������ ������������ LU, �� ������ ������������ ���������������� ����� � ������� OpenMP.
    //  3) �� ��������� � ��� ������� L � U ����������� (MPI_Bcast) ��������� ���������.
    //=============================================================================
    DecomposeMatrix lu_decomposition_mpi_omp(const Matrix& matrix)
    {
        // ���������� rank � size
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix.n;
        int b = matrix.b;

        // �� ������ �������� �������� ������ ��� L � U (n x n).
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // ====== Rank 0 ������ ��� ������ �� ������������, ������������ OpenMP ======
        if (rank == 0)
        {
            // 1) ������������� U ������ A, � L � ��������� ���������
            // ���������� ������������� ������ � ������� OpenMP
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result.u[i][j] = matrix.A[i][j];
                }
                result.l[i][i] = 1.0;
            }

            // 2) ���������� LU-������������ � �������������� OpenMP
            //    (��� ������ �������� ��������).
            for (int k = 0; k < n - 1; k++) {
                double pivotVal = result.u[k][k];

                // ��������� ������ i = k+1..k+b
                // ������������� ���� �� i
                #pragma omp parallel for schedule(static)
                for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
                    result.l[i][k] = result.u[i][k] / pivotVal;
                    for (int j = k; j < MIN(k + b + 1, n); j++) {
                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
                    }
                }
            }
        }

        // ====== ��������� ���������� L � U ���� ��������� ======
        // ����� ��� �������� (������� Rank 0) ��������� � Bcast,
        // ����� � ���� ��� ���������� L,U.
        for (int i = 0; i < n; i++) {
            MPI_Bcast(result.l[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(result.u[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        return result;
    }

    //=============================================================================
    // ������� ������� A*x = C (��� A = L*U) � �������������� MPI + OpenMP.
    // ������:
    //  1) Rank 0 ������ ������/�������� ���.
    //     - ������ ��� (L*y = C) ����� �������� �������������� ���������� ���� ������������.
    //     - �������� ��� (U*x = y) ���� �������� ��������������.
    //  2) ��������� (MPI_Bcast) ��������� ������ X ���� ��������� ���������.
    //=============================================================================
    void solve_lu_mpi_omp(const DecomposeMatrix& decomp, Matrix* matrix)
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix->n;

        // �������� ������ ��� ��������������� ������ y
        double* y = (double*)malloc(n * sizeof(double));

        // Rank 0 � ��������� ������ � �������� ���
        if (rank == 0)
        {
            // ---- ������ ��� L*y = C ----
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                // ���������� ���� ������������ ����� ��������������
                #pragma omp parallel for reduction(+:s)
                for (int j = 0; j < i; j++) {
                    s += decomp.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }

            // ---- �������� ��� U*x = y ----
            matrix->X = (double*)malloc(n * sizeof(double));

            for (int i = n - 1; i >= 0; i--) {
                double s = 0.0;
                // ����� ������������� ���������� ���� ������������
                #pragma omp parallel for reduction(+:s)
                for (int j = i + 1; j < n; j++) {
                    s += decomp.u[i][j] * matrix->X[j];
                }
                matrix->X[i] = (y[i] - s) / decomp.u[i][i];
            }
        }
        else {
            // �� ������ ������ ������ ������� X, ����� ����� ���������� ���������
            matrix->X = (double*)calloc(n, sizeof(double));
        }

        // ��������� (bcast) ������� ������ X �� Rank 0 �� ����
        MPI_Bcast(matrix->X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // ����������� ��������� �����
        free(y);
    }

}

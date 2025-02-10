//#include <stdio.h>
//#include <stdlib.h>
//#include <mpi.h>
//#include <omp.h>
//#include "matrix.h"
//
//#ifndef MIN
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#endif
//
//namespace band_matrix_mpi_omp
//{
//    static void reverse_array(double* array, size_t n)
//    {
//        if (!array || n == 0) return;
//        for (size_t i = 0; i < n / 2; i++) {
//            double tmp = array[i];
//            array[i] = array[n - 1 - i];
//            array[n - 1 - i] = tmp;
//        }
//    }
//
//    /**
//     * ������ ������������� LU-���������� ��� ��������� �������
//     * � ������� MPI + OpenMP.
//     *
//     * ���������: ������ ������ ������� A �� ������ ��������.
//     */
//    DecomposeMatrix lu_decomposition_mpi_omp(const Matrix matrix)
//    {
//        int rank, size;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        size_t n = matrix.n;
//        size_t b = matrix.b;
//
//        // ������� L � U
//        DecomposeMatrix result;
//        result.l = (double**)malloc(n * sizeof(double*));
//        result.u = (double**)malloc(n * sizeof(double*));
//        for (size_t i = 0; i < n; i++) {
//            result.l[i] = (double*)calloc(n, sizeof(double));
//            result.u[i] = (double*)malloc(n * sizeof(double));
//        }
//
//        // �������������
//        // U = A, L = ��������� ���������
//        for (size_t i = 0; i < n; i++) {
//            for (size_t j = 0; j < n; j++) {
//                result.u[i][j] = matrix.A[i][j];
//                // L[i][j] ��� ��������� ������
//            }
//            result.l[i][i] = 1.0;
//        }
//
//        // ������� ���� LU-����������
//        for (size_t k = 0; k < n - 1; k++) {
//            // ����������, "����" ����������� ������ k.
//            // ���������� ������� ������:
//            size_t blockSize = n / size;
//            size_t start = blockSize * rank;
//            size_t end = (rank == size - 1) ? (n - 1) : (blockSize * (rank + 1) - 1);
//
//            // ��������, ����� �� k � ��������� ����� ������� ��������
//            int owner = (int)(k / blockSize);
//            if (owner >= size) owner = size - 1; // �� ������, ���� n �� ������� �����
//
//            // ���� � � �������� pivot-������, � � �������� � ����� ��������
//            if (rank == owner) {
//                double pivotVal = result.u[k][k];
//                // ��������� ������� ������ �� ������������� �������� (������ ������, 
//                // ��� ������������� LU ������ �� ����� ��� ������, � ��������� L[i][k].
//                // �� � ��������� ����������� ��������� ��������� ��� ������������� ������� ������)
//                // ����� ����� ���������� ����������, �.�. ������������ ����� ������ pivot � U[k][k].
//                // ��� ��������������� � ����� ���� ����� ������������ ������ ����� L[i][k] = U[i][k]/U[k][k].
//                // �� ���� ����� ����� ������������� �������:
//                // for (size_t j = k; j < n; j++) {
//                //     result.u[k][j] /= pivotVal;
//                // }
//            }
//
//            // ������ �������-�������� ��������� pivot-������ k ����
//            MPI_Bcast(result.u[k], (int)n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
//
//            // ��������� ����� ��������������� ����� L? ���, L[k][k] = 1. ��������� �� �����.
//            // ��� �� ������� �� k+1 �� k+b, ��������� ��
//            // ������ ������� ������ ��� ��� "�����" �����
//#pragma omp parallel for
//            for (size_t i = start; i <= end; i++) {
//                if (i > k && i < MIN(k + b + 1, n)) {
//                    // ��������� L[i][k]
//                    result.l[i][k] = result.u[i][k] / result.u[k][k];
//                    // ��������� U[i][j]
//                    for (size_t j = k; j < MIN(k + b + 1, n); j++) {
//                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
//                    }
//                }
//            }
//
//            // ����� ���������������� U (���� ������ ������� ������� ������ �����), 
//            // ����� ����� ����� ������� U ��� ��������� ���������� ������. 
//            // � ���������� ������� ����� �������� ��� ������ i= k+1..k+b 
//            // (�� ����� ���� ����� ���� �� �� ������� ������ U �� ������ ��������, 
//            // �� ��� ����������� � ������� Allgather).
//            for (size_t i = k + 1; i < MIN(k + b + 1, n); i++) {
//                MPI_Bcast(result.u[i], (int)n, MPI_DOUBLE, (int)(i / blockSize), MPI_COMM_WORLD);
//                MPI_Bcast(result.l[i], (int)n, MPI_DOUBLE, (int)(i / blockSize), MPI_COMM_WORLD);
//            }
//        }
//
//        return result;
//    }
//
//    void solve_lu_mpi_omp(const DecomposeMatrix decompose_matrix, Matrix* matrix)
//    {
//        // ���������� ������, � ������� �������� �� ������� �� ������ ��������
//        // � ������ ������/�������� ��� ��������. � ���������� ����� ���������
//        // ������ ������������� ������/�������� ���.
//        int rank, size;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        size_t n = matrix->n;
//
//        double* y = (double*)malloc(n * sizeof(double));
//        // ������ ���: L * y = C
//        for (size_t i = 0; i < n; i++) {
//            double s = 0.0;
//            for (size_t j = 0; j < i; j++) {
//                s += decompose_matrix.l[i][j] * y[j];
//            }
//            y[i] = matrix->C[i] - s;
//        }
//
//        // �������� ���: U * x = y
//        matrix->X = (double*)malloc(n * sizeof(double));
//        for (int i = (int)n - 1, k = 0; i >= 0; --i, k++) {
//            double s = 0.0;
//            for (int j = (int)n - 1; j > i; --j) {
//                s += decompose_matrix.u[i][j] * matrix->X[n - 1 - j];
//            }
//            matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
//        }
//        reverse_array(matrix->X, n);
//
//        free(y);
//    }
//}
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_mpi
{
    struct DecomposeMatrix lu_decomposition_mpi(struct Matrix matrix, int rank, int size)
    {
        struct DecomposeMatrix result;
        size_t n = matrix.n;
        size_t b = matrix.b;

        // Вычисление локального диапазона строк
        size_t rows_per_proc = n / size;
        size_t start_row = rank * rows_per_proc;
        size_t end_row = (rank == size - 1) ? n : start_row + rows_per_proc;

        // Выделение памяти
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));

        for (size_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));

            // Копирование только своей части данных
            if (i >= start_row && i < end_row) {
                for (size_t j = 0; j < n; j++) {
                    result.u[i][j] = matrix.A[i][j];
                }
            }
        }

        // Инициализация диагональных элементов L
        for (size_t i = start_row; i < end_row; i++) {
            result.l[i][i] = 1.0;
        }

        // Основной цикл LU-разложения
        for (size_t k = 0; k < n - 1; k++) {
            // Broadcast k-й строки
            MPI_Bcast(result.u[k], n, MPI_DOUBLE, k / rows_per_proc, MPI_COMM_WORLD);

            size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

            // Обработка только своих строк
            for (size_t i = k + 1; i < upper_bound; i++) {
                if (i >= start_row && i < end_row) {
                    result.l[i][k] = result.u[i][k] / result.u[k][k];
                    for (size_t j = k; j < upper_bound; j++) {
                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
                    }
                }
            }

            // Синхронизация результатов
            for (size_t i = k + 1; i < upper_bound; i++) {
                int owner = i / rows_per_proc;
                MPI_Bcast(result.u[i], n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
                MPI_Bcast(result.l[i], n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
            }
        }

        return result;
    }

    void reverse_array_mpi(double* array, size_t n)
    {
        size_t i;
        double temp;
        for (i = 0; i < n / 2; ++i)
        {
            temp = array[i];
            array[i] = array[n - 1 - i];
            array[n - 1 - i] = temp;
        }
    }

    void solve_lu_mpi(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix, int rank, int size)
    {
        size_t n = matrix->n;
        size_t rows_per_proc = n / size;
        size_t start_row = rank * rows_per_proc;
        size_t end_row = (rank == size - 1) ? n : start_row + rows_per_proc;

        double* y = (double*)malloc(n * sizeof(double));

        // Прямой ход (решение Ly = b)
        for (size_t i = 0; i < n; i++)
        {
            if (i >= start_row && i < end_row)
            {
                double s = 0.0;
                for (size_t j = 0; j < i; j++)
                {
                    s += decompose_matrix.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }
            // Broadcast промежуточных результатов
            MPI_Bcast(&y[i], 1, MPI_DOUBLE, i / rows_per_proc, MPI_COMM_WORLD);
        }

        // Обратный ход (решение Ux = y)
        if (rank == 0)
        {
            matrix->X = (double*)malloc(n * sizeof(double));
        }

        for (int i = n - 1, k = 0; i >= 0; --i, k++)
        {
            if (i >= start_row && i < end_row)
            {
                double s = 0.0;
                for (int j = n - 1; j > i; --j)
                {
                    s += decompose_matrix.u[i][j] * matrix->X[n - j - 1];
                }
                matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
            }
            // Broadcast результатов
            if (matrix->X != NULL)
            {
                MPI_Bcast(&matrix->X[k], 1, MPI_DOUBLE, i / rows_per_proc, MPI_COMM_WORLD);
            }
        }

        if (rank == 0)
        {
            reverse_array_mpi(matrix->X, n);
        }

        free(y);
    }
}
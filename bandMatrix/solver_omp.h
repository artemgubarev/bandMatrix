#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_omp
{
    static void reverse_array(double* array, size_t n)
    {
        if (!array || n == 0) return;
        for (size_t i = 0; i < n / 2; i++) {
            double tmp = array[i];
            array[i] = array[n - 1 - i];
            array[n - 1 - i] = tmp;
        }
    }

    DecomposeMatrix lu_decomposition_omp(const Matrix matrix)
    {
        size_t n = matrix.n;
        size_t b = matrix.b;

        // Выделение памяти для L и U
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (size_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // Инициализация L и U
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                result.u[i][j] = matrix.A[i][j];
            }
            result.l[i][i] = 1.0;
        }

        // Главный цикл LU-разложения
        for (size_t k = 0; k < n - 1; k++) {
            double pivotVal = result.u[k][k];

            // Обновление строк от k+1 до k+b параллельно
#pragma omp parallel for
            for (size_t i = k + 1; i < MIN(k + b + 1, n); i++) {
                result.l[i][k] = result.u[i][k] / pivotVal;
                for (size_t j = k; j < MIN(k + b + 1, n); j++) {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }

        return result;
    }

    void solve_lu_omp(const DecomposeMatrix decompose_matrix, Matrix* matrix)
    {
        size_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        // Прямой ход: L * y = C
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
#pragma omp parallel for reduction(+:s)
            for (size_t j = 0; j < i; j++) {
                s += decompose_matrix.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        // Обратный ход: U * x = y
        matrix->X = (double*)malloc(n * sizeof(double));
        for (int i = (int)n - 1, k = 0; i >= 0; --i, k++) {
            double s = 0.0;
#pragma omp parallel for reduction(+:s)
            for (int j = (int)n - 1; j > i; --j) {
                s += decompose_matrix.u[i][j] * matrix->X[n - 1 - j];
            }
            matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
        }
        reverse_array(matrix->X, n);
        free(y);
    }
}
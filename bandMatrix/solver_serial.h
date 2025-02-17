#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include "matrix.h"

namespace band_matrix_serial
{
    struct DecomposeMatrix lu_decomposition(struct Matrix matrix)
    {
        struct DecomposeMatrix result;

        size_t n = matrix.n;
        size_t b = matrix.b;

        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));

        for (size_t i = 0; i < n; i++)
        {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                result.u[i][j] = matrix.A[i][j];
            }
        }

        for (size_t i = 0; i < n; i++)
        {
            result.l[i][i] = 1.0;
        }

        for (size_t k = 0; k < n - 1; ++k)
        {
            size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
            for (size_t i = k + 1; i < upper_bound; ++i)
            {
                result.l[i][k] = result.u[i][k] / result.u[k][k];
                for (size_t j = k; j < upper_bound; ++j)
                {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }
        return result;
    }

    void reverse_array(double* array, size_t n)
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

    void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix)
    {
        size_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        for (size_t i = 0; i < n; i++)
        {
            double s = 0.0;
            for (size_t j = 0; j < i; j++)
            {
                s += decompose_matrix.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        matrix->X = (double*)malloc(n * sizeof(double));

        for (int i = n - 1, k = 0; i >= 0; --i, k++)
        {
            double s = 0.0;
            for (int j = n - 1; j > i; --j)
            {
                s += decompose_matrix.u[i][j] * matrix->X[n - j - 1];
            }
            matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
        }

        reverse_array(matrix->X, n);
        free(y);
    }
}
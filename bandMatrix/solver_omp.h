#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"
#include <utility>

namespace band_matrix_omp
{
    void reverse_array(double* array, size_t n)
    {
        double* reversed = new double[n];

        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            reversed[i] = array[n - 1 - i];
        }

        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            array[i] = reversed[i];
        }

        delete[] reversed;
    }

    void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix)
    {
        size_t n = matrix->n;

        double* y = (double*)malloc(n * sizeof(double));

        for (size_t i = 0; i < n; i++) 
        {
            if (i == 0) 
            {
                y[i] = matrix->C[i];
            }
            else {
                double s = 0;
                #pragma omp parallel for reduction(+:s)
                for (int j = 0; j < static_cast<int>(i); j++) 
                {
                    s += decompose_matrix.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }
        }

        matrix->X = (double*)malloc(n * sizeof(double));

        for (int i = n - 1, k = 0; i >= 0; --i, k++) 
        {
            if (i == n - 1) {
                matrix->X[k] = y[i] / decompose_matrix.u[i][i];
            }
            else {
                double s = 0;
                #pragma omp parallel for reduction(+:s)
                for (int j = n - 1; j > i; --j) {
                    s += decompose_matrix.u[i][j] * matrix->X[n - j - 1];
                }
                matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
            }
        }

        reverse_array(matrix->X, n);
        free(y);
    }

    struct DecomposeMatrix lu_decomposition(Matrix matrix)
    {
        struct DecomposeMatrix result;
        size_t n = matrix.n;
        size_t b = matrix.b;

        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));

        for (size_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                result.u[i][j] = matrix.A[i][j];
            }
            result.l[i][i] = 1.0;
        }

        for (int k = 0; k < static_cast<int>(n) - 1; ++k) 
        {
            #pragma omp parallel for
            for (int i = k + 1; i < std::min(k + static_cast<int>(b) + 1, static_cast<int>(n)); ++i) 
            {
                result.l[i][k] = result.u[i][k] / result.u[k][k];

                #pragma omp parallel for
                for (int j = k; j < std::min(k + static_cast<int>(b) + 1, static_cast<int>(n)); ++j) 
                {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }

        return result;
    }
}
//#pragma once
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <utility>
//#include "matrix.h"
//
//namespace band_matrix
//{
//    struct DecomposeMatrix lu_decomposition(Matrix matrix)
//    {
//        struct DecomposeMatrix result;
//
//        size_t n = matrix.n;
//        size_t b = matrix.b;
//
//        result.l = (double**)malloc(n * sizeof(double*));
//        result.u = (double**)malloc(n * sizeof(double*));
//
//        for (size_t i = 0; i < n; i++)
//        {
//            result.l[i] = (double*)calloc(n, sizeof(double));
//            result.u[i] = (double*)malloc(n * sizeof(double));
//        }
//
//        for (size_t i = 0; i < n; i++)
//        {
//            for (size_t j = 0; j < n; j++)
//            {
//                result.u[i][j] = matrix.A[i][j];
//            }
//        }
//
//        for (size_t i = 0; i < n; i++)
//        {
//            result.l[i][i] = 1.0;
//        }
//
//        for (size_t k = 0; k < n - 1; ++k)
//        {
//            for (size_t i = k + 1; i < std::min(k + b + 1, n); ++i)
//            {
//                result.l[i][k] = result.u[i][k] / result.u[k][k];
//                for (size_t j = k; j < std::min(k + b + 1, n); ++j)
//                {
//                    result.u[i][j] = result.u[i][j] - result.l[i][k] * result.u[k][j];
//                }
//            }
//        }
//        return result;
//    }
//
//    void reverse_array(double* array, size_t n)
//    {
//        double* reversed = new double[n];
//
//        for (int i = 0; i < n; ++i) {
//            reversed[i] = array[n - 1 - i];
//        }
//
//        for (int i = 0; i < n; ++i) {
//            array[i] = reversed[i];
//        }
//
//        delete[] reversed;
//    }
//
//    void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix)
//    {
//        size_t n = matrix->n;
//
//        double* y = (double*)malloc(n * sizeof(double));
//
//        for (size_t i = 0; i < n; i++)
//        {
//            if (i == 0)
//            {
//                y[i] = matrix->C[i];
//            }
//            else
//            {
//                double s = 0;
//                for (size_t j = 0; j < i; j++)
//                {
//                    s += decompose_matrix.l[i][j] * y[j];
//                }
//                y[i] = matrix->C[i] - s;
//            }
//        }
//
//        matrix->X = (double*)malloc(n * sizeof(double));
//
//        for (int i = n - 1, k = 0; i >= 0; --i, k++)
//        {
//            if (i == n - 1)
//            {
//                matrix->X[k] = y[i] / decompose_matrix.u[i][i];
//            }
//            else
//            {
//                double s = 0;
//                for (int j = n - 1; j > i; --j)
//                {
//                    s += decompose_matrix.u[i][j] * matrix->X[n - j - 1];
//                }
//                matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
//            }
//        }
//
//        reverse_array(matrix->X, n);
//        free(y);
//    }
//}
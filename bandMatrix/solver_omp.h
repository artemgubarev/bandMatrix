//#include <stdio.h>
//#include <stdlib.h>
//#include <omp.h>
//#include "matrix.h"
//
//#ifndef MIN
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#endif
//
//namespace band_matrix_omp
//{
//    static void reverse_array(double* array, int n)
//    {
//        if (!array || n == 0) return;
//        for (int i = 0; i < n / 2; i++) {
//            double tmp = array[i];
//            array[i] = array[n - 1 - i];
//            array[n - 1 - i] = tmp;
//        }
//    }
//
//    DecomposeMatrix lu_decomposition_omp(const Matrix matrix)
//    {
//        int n = matrix.n;
//        int b = matrix.b;
//
//        // Выделение памяти для L и U
//        DecomposeMatrix result;
//        result.l = (double**)malloc(n * sizeof(double*));
//        result.u = (double**)malloc(n * sizeof(double*));
//        for (int i = 0; i < n; i++) {
//            result.l[i] = (double*)calloc(n, sizeof(double));
//            result.u[i] = (double*)malloc(n * sizeof(double));
//        }
//
//        // Инициализация L и U
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                result.u[i][j] = matrix.A[i][j];
//            }
//            result.l[i][i] = 1.0;
//        }
//
//        // Главный цикл LU-разложения
//        for (int k = 0; k < n - 1; k++) {
//            double pivotVal = result.u[k][k];
//
//            // Обновление строк от k+1 до k+b параллельно
//            #pragma omp parallel for
//            for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
//                result.l[i][k] = result.u[i][k] / pivotVal;
//                for (int j = k; j < MIN(k + b + 1, n); j++) {
//                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
//                }
//            }
//        }
//
//        return result;
//    }
//
//    void solve_lu_omp(const DecomposeMatrix decompose_matrix, Matrix* matrix)
//    {
//        int n = matrix->n;
//        double* y = (double*)malloc(n * sizeof(double));
//
//        // Прямой ход: L * y = C
//        for (int i = 0; i < n; i++) {
//            double s = 0.0;
//            #pragma omp parallel for reduction(+:s)
//            for (int j = 0; j < i; j++) {
//                s += decompose_matrix.l[i][j] * y[j];
//            }
//            y[i] = matrix->C[i] - s;
//        }
//
//        // Обратный ход: U * x = y
//        matrix->X = (double*)malloc(n * sizeof(double));
//        for (int i = (int)n - 1, k = 0; i >= 0; --i, k++) {
//            double s = 0.0;
//            #pragma omp parallel for reduction(+:s)
//            for (int j = (int)n - 1; j > i; --j) {
//                s += decompose_matrix.u[i][j] * matrix->X[n - 1 - j];
//            }
//            matrix->X[k] = (y[i] - s) / decompose_matrix.u[i][i];
//        }
//        reverse_array(matrix->X, n);
//        free(y);
//    }
//}


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

namespace band_matrix_omp
{
    // LU-разложение для полосатой матрицы с использованием OpenMP.
    DecomposeMatrix lu_decomposition_omp(const Matrix matrix)
    {
        int n = matrix.n;
        int b = matrix.b;

        // Выделение памяти для L и U.
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // Параллельная инициализация матриц L и U.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.u[i][j] = matrix.A[i][j];
            }
            result.l[i][i] = 1.0;
        }

        // Главный цикл LU-разложения.
        for (int k = 0; k < n - 1; k++) {
            double pivotVal = result.u[k][k];

            // Обновление строк от k+1 до k+b (учитывая границу n).
            #pragma omp parallel for schedule(static)
            for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
                result.l[i][k] = result.u[i][k] / pivotVal;
                for (int j = k; j < MIN(k + b + 1, n); j++) {
                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }
        }

        return result;
    }

    // Функция решения системы методом LU-разложения.
    void solve_lu_omp(const DecomposeMatrix decomp, Matrix* matrix)
    {
        int n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        // Прямой ход: решение системы L*y = C.
        // Данный цикл вынужден быть последовательным, т.к. y[i] зависит от y[0..i-1].
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            // Внутренний цикл суммирования можно распараллелить, если число элементов достаточно велико.
            #pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = 0; j < i; j++) {
                s += decomp.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        // Обратный ход: решение системы U*x = y.
        matrix->X = (double*)malloc(n * sizeof(double));
        for (int i = (int)n - 1; i >= 0; i--) {
            double s = 0.0;
            #pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = i + 1; j < (int)n; j++) {
                s += decomp.u[i][j] * matrix->X[j];
            }
            matrix->X[i] = (y[i] - s) / decomp.u[i][i];
        }

        free(y);
    }
}

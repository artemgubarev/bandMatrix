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
    // Гибридная функция LU-разложения полосатой (ленточной) матрицы с использованием MPI + OpenMP.
    // Логика:
    //  1) Rank 0 читает и хранит полную матрицу A.
    //  2) Rank 0 делает факторизацию LU, но внтури факторизации распараллеливает циклы с помощью OpenMP.
    //  3) По окончании — все матрицы L и U рассылаются (MPI_Bcast) остальным процессам.
    //=============================================================================
    DecomposeMatrix lu_decomposition_mpi_omp(const Matrix& matrix)
    {
        // Определяем rank и size
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix.n;
        int b = matrix.b;

        // На каждом процессе выделяем память для L и U (n x n).
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // ====== Rank 0 делает всю работу по факторизации, параллелизуя OpenMP ======
        if (rank == 0)
        {
            // 1) Инициализация U копией A, и L — единичная диагональ
            // Параллелим инициализацию циклов с помощью OpenMP
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result.u[i][j] = matrix.A[i][j];
                }
                result.l[i][i] = 1.0;
            }

            // 2) Выполнение LU-факторизации с использованием OpenMP
            //    (без выбора главного элемента).
            for (int k = 0; k < n - 1; k++) {
                double pivotVal = result.u[k][k];

                // Обновляем строки i = k+1..k+b
                // Распараллелим цикл по i
                #pragma omp parallel for schedule(static)
                for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
                    result.l[i][k] = result.u[i][k] / pivotVal;
                    for (int j = k; j < MIN(k + b + 1, n); j++) {
                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
                    }
                }
            }
        }

        // ====== Рассылаем результаты L и U всем процессам ======
        // Здесь все процессы (включая Rank 0) участвуют в Bcast,
        // чтобы у всех был одинаковый L,U.
        for (int i = 0; i < n; i++) {
            MPI_Bcast(result.l[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(result.u[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        return result;
    }

    //=============================================================================
    // Решение системы A*x = C (где A = L*U) с использованием MPI + OpenMP.
    // Логика:
    //  1) Rank 0 делает прямой/обратный ход.
    //     - Прямой ход (L*y = C) можно частично распараллелить внутренний цикл суммирования.
    //     - Обратный ход (U*x = y) тоже частично распараллелить.
    //  2) Рассылает (MPI_Bcast) найденный вектор X всем остальным процессам.
    //=============================================================================
    void solve_lu_mpi_omp(const DecomposeMatrix& decomp, Matrix* matrix)
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix->n;

        // Выделяем память под вспомогательный вектор y
        double* y = (double*)malloc(n * sizeof(double));

        // Rank 0 — выполняет прямой и обратный ход
        if (rank == 0)
        {
            // ---- Прямой ход L*y = C ----
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                // Внутренний цикл суммирования можно распараллелить
                #pragma omp parallel for reduction(+:s)
                for (int j = 0; j < i; j++) {
                    s += decomp.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }

            // ---- Обратный ход U*x = y ----
            matrix->X = (double*)malloc(n * sizeof(double));

            for (int i = n - 1; i >= 0; i--) {
                double s = 0.0;
                // Снова распараллелим внутренний цикл суммирования
                #pragma omp parallel for reduction(+:s)
                for (int j = i + 1; j < n; j++) {
                    s += decomp.u[i][j] * matrix->X[j];
                }
                matrix->X[i] = (y[i] - s) / decomp.u[i][i];
            }
        }
        else {
            // На других рангах просто выделим X, чтобы иметь корректный указатель
            matrix->X = (double*)calloc(n, sizeof(double));
        }

        // Рассылаем (bcast) готовый вектор X от Rank 0 ко всем
        MPI_Bcast(matrix->X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Освобождаем временный буфер
        free(y);
    }

}

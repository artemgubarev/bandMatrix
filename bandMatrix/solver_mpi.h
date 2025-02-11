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

    //=============================================================================
    // Функция LU-разложения полосатой (ленточной) матрицы с использованием *только* MPI.
    // В данном упрощённом варианте Rank 0 делает всё вычисление локально (как было в OpenMP-версии),
    // а затем рассылает (broadcast) полученные матрицы L и U всем остальным.
    //=============================================================================
    DecomposeMatrix lu_decomposition_mpi(const Matrix& matrix)
    {
        // Инициализируем MPI
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix.n;
        int b = matrix.b;

        // На *каждом* процессе выделяем память для L и U
        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // ====== Rank 0 делает всю работу по факторизации ======
        if (rank == 0)
        {
            // 1) Инициализация U копией A, и L — единичная диагональ
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result.u[i][j] = matrix.A[i][j];
                }
                result.l[i][i] = 1.0;
            }

            // 2) Выполнение LU-факторизации (без выбора главного элемента)
            for (int k = 0; k < n - 1; k++) {
                double pivotVal = result.u[k][k];

                // Обновление строк i = k+1..k+b (но не выходя за границу n)
                for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
                    result.l[i][k] = result.u[i][k] / pivotVal;
                    for (int j = k; j < MIN(k + b + 1, n); j++) {
                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
                    }
                }
            }
        }

        // ====== Рассылаем результаты L и U всем процессам ======
        // Так как L и U имеют размер n x n, рассылаем построчно.
        // Можно (и лучше) делать MPI_Bcast покомпонентно для каждого row.
        // Здесь — простейший пример.
        for (int i = 0; i < n; i++) {
            MPI_Bcast(result.l[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(result.u[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        return result;
    }

    //=============================================================================
    // Решение системы A*x = C (где A = L*U) с использованием *только* MPI.
    // Снова, в упрощённом виде: Rank 0 делает прямой и обратный ход целиком,
    // потом рассылаем (broadcast) вектор X всем.
    //=============================================================================
    void solve_lu_mpi(const DecomposeMatrix& decomp, Matrix* matrix)
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int n = matrix->n;

        // Выделяем память под вспомогательный вектор y, если нужно
        double* y = (double*)malloc(n * sizeof(double));

        // Rank 0 — выполняет прямой и обратный ход
        if (rank == 0)
        {
            // ---- Прямой ход L*y = C ----
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < i; j++) {
                    s += decomp.l[i][j] * y[j];
                }
                y[i] = matrix->C[i] - s;
            }

            // ---- Обратный ход U*x = y ----
            matrix->X = (double*)malloc(n * sizeof(double));
            for (int i = n - 1; i >= 0; i--) {
                double s = 0.0;
                for (int j = i + 1; j < n; j++) {
                    s += decomp.u[i][j] * matrix->X[j];
                }
                matrix->X[i] = (y[i] - s) / decomp.u[i][i];
            }
        }
        else {
            // На других рангах пока просто выделим X,
            // чтобы был корректный указатель при Bcast.
            matrix->X = (double*)calloc(n, sizeof(double));
        }

        // Рассылаем (bcast) готовый вектор X от Rank 0 ко всем
        MPI_Bcast(matrix->X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Освобождаем временный буфер
        free(y);
    }

} // namespace band_matrix_mpi

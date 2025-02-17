//#include <stdio.h>
//#include <stdlib.h>
//#include <mpi.h>
//#include <math.h>
//#include "matrix.h"
//
//#ifndef MIN
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#endif
//
//namespace band_matrix_mpi
//{
//    //=============================================================================
//    // Функция LU-разложения полосатой (ленточной) матрицы с использованием *только* MPI.
//    // В данном упрощённом варианте Rank 0 делает всё вычисление локально (как было в OpenMP-версии),
//    // а затем рассылает (broadcast) полученные матрицы L и U всем остальным.
//    //=============================================================================
//    DecomposeMatrix lu_decomposition_mpi(const Matrix& matrix)
//    {
//        // Инициализируем MPI
//        int rank, size;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        int n = matrix.n;
//        int b = matrix.b;
//
//        // На *каждом* процессе выделяем память для L и U
//        DecomposeMatrix result;
//        result.l = (double**)malloc(n * sizeof(double*));
//        result.u = (double**)malloc(n * sizeof(double*));
//        for (int i = 0; i < n; i++) {
//            result.l[i] = (double*)calloc(n, sizeof(double));
//            result.u[i] = (double*)malloc(n * sizeof(double));
//        }
//
//        // ====== Rank 0 делает всю работу по факторизации ======
//        if (rank == 0)
//        {
//            // 1) Инициализация U копией A, и L — единичная диагональ
//            for (int i = 0; i < n; i++) {
//                for (int j = 0; j < n; j++) {
//                    result.u[i][j] = matrix.A[i][j];
//                }
//                result.l[i][i] = 1.0;
//            }
//
//            // 2) Выполнение LU-факторизации (без выбора главного элемента)
//            for (int k = 0; k < n - 1; k++) {
//                double pivotVal = result.u[k][k];
//
//                // Обновление строк i = k+1..k+b (но не выходя за границу n)
//                for (int i = k + 1; i < MIN(k + b + 1, n); i++) {
//                    result.l[i][k] = result.u[i][k] / pivotVal;
//                    for (int j = k; j < MIN(k + b + 1, n); j++) {
//                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
//                    }
//                }
//            }
//        }
//
//        // ====== Рассылаем результаты L и U всем процессам ======
//        // Так как L и U имеют размер n x n, рассылаем построчно.
//        // Можно (и лучше) делать MPI_Bcast покомпонентно для каждого row.
//        // Здесь — простейший пример.
//        for (int i = 0; i < n; i++) {
//            MPI_Bcast(result.l[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Bcast(result.u[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        }
//
//        return result;
//    }
//
//    //=============================================================================
//    // Решение системы A*x = C (где A = L*U) с использованием *только* MPI.
//    // Снова, в упрощённом виде: Rank 0 делает прямой и обратный ход целиком,
//    // потом рассылаем (broadcast) вектор X всем.
//    //=============================================================================
//    void solve_lu_mpi(const DecomposeMatrix& decomp, Matrix* matrix)
//    {
//        int rank, size;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        int n = matrix->n;
//        double* y = (double*)malloc(n * sizeof(double));
//
//        if (rank == 0)
//        {
//            for (int i = 0; i < n; i++) {
//                double s = 0.0;
//                for (int j = 0; j < i; j++) {
//                    s += decomp.l[i][j] * y[j];
//                }
//                y[i] = matrix->C[i] - s;
//            }
//
//            matrix->X = (double*)malloc(n * sizeof(double));
//            for (int i = n - 1; i >= 0; i--) {
//                double s = 0.0;
//                for (int j = i + 1; j < n; j++) {
//                    s += decomp.u[i][j] * matrix->X[j];
//                }
//                matrix->X[i] = (y[i] - s) / decomp.u[i][i];
//            }
//        }
//        else {
//            matrix->X = (double*)calloc(n, sizeof(double));
//        }
//
//        MPI_Bcast(matrix->X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//        free(y);
//    }
//}


#include <mpi.h>
#include <stdlib.h>
#include <math.h>

// Структура для локальной части матрицы
typedef struct {
    int local_n;      // число строк, принадлежащих данному процессу
    int global_start; // глобальный номер первой строки в локальном блоке
    int n;            // общий размер матрицы
    int b;            // ширина ленты
    double* A;        // локальная часть матрицы A (хранится в виде одномерного массива, строки подряд)
    double* C;        // локальная часть правой части
    double* X;        // локальная часть решения
} LocalMatrix;

namespace band_matrix_mpi {

    // Параллельное LU-разложение для полосатой матрицы с распределением по строкам.
    // Факторизация выполняется "на месте": в A сохраняется U, а ниже диагонали – множители L (с диагональю, равной 1).
    void lu_decomposition_mpi_parallel(LocalMatrix& localMat) {
        int rank, nprocs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        int n = localMat.n;
        int b = localMat.b;
        // Для каждого шага k (пивотной строки)
        for (int k = 0; k < n - 1; k++) {
            // Определяем, какой процесс владеет строкой с номером k.
            // Для простоты вычисляем распределение заново:
            int base = n / nprocs;
            int rem = n % nprocs;
            int owner;
            int cum = 0;
            for (int p = 0; p < nprocs; p++) {
                int rows = (p < rem) ? base + 1 : base;
                if (k < cum + rows) {
                    owner = p;
                    break;
                }
                cum += rows;
            }
            // Подготавливаем буфер для пивотной строки (все процессы выделяют память под массив размера n)
            double* pivot_row = new double[n];

            if (rank == owner) {
                int local_index = k - localMat.global_start;
                // Копируем всю строку (можно ограничиться столбцами от k до min(n, k+b+1), но для простоты копируем полностью)
                for (int j = 0; j < n; j++) {
                    pivot_row[j] = localMat.A[local_index * n + j];
                }
            }
            // Рассылаем пивотную строку всем процессам
            MPI_Bcast(pivot_row, n, MPI_DOUBLE, owner, MPI_COMM_WORLD);

            // Каждый процесс обновляет свои строки, для которых глобальный номер i удовлетворяет: k < i < min(n, k+b+1)
            for (int i = 0; i < localMat.local_n; i++) {
                int global_row = localMat.global_start + i;
                if (global_row > k && global_row < ((k + b + 1) < n ? (k + b + 1) : n)) {
                    double factor = localMat.A[i * n + k] / pivot_row[k];
                    localMat.A[i * n + k] = factor; // сохраняем множитель L
                    int j_end = (k + b + 1 < n) ? (k + b + 1) : n;
                    for (int j = k; j < j_end; j++) {
                        localMat.A[i * n + j] -= factor * pivot_row[j];
                    }
                }
            }
            // Синхронизируем процессы перед следующим шагом
            MPI_Barrier(MPI_COMM_WORLD);
            delete[] pivot_row;
        }
    }

    // Параллельное решение системы A*x = C.
    // Сначала решается L*y = C (прямой ход), затем U*x = y (обратный ход).
    // Для каждого шага владеющий процесс вычисляет элемент, затем транслирует его всем.
    void solve_lu_mpi_parallel(LocalMatrix& localMat) {
        int rank, nprocs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        int n = localMat.n;

        // Решение L*y = C.
        // Создаём глобальный вектор y (на всех процессах он будет храниться полностью).
        double* y = (double*)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            int base = n / nprocs;
            int rem = n % nprocs;
            int owner;
            int cum = 0;
            for (int p = 0; p < nprocs; p++) {
                int rows = (p < rem) ? base + 1 : base;
                if (i < cum + rows) {
                    owner = p;
                    break;
                }
                cum += rows;
            }
            double yi;
            if (rank == owner) {
                int local_index = i - localMat.global_start;
                yi = localMat.C[local_index];
                // Вычитаем сумму L(i,j)*y[j] для j от 0 до i-1.
                for (int j = 0; j < i; j++) {
                    yi -= localMat.A[local_index * n + j] * y[j];
                }
            }
            // Транслируем вычисленное значение y[i] от процесса-владельца
            MPI_Bcast(&yi, 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);
            y[i] = yi;
        }

        // Решение U*x = y.
        double* x = (double*)malloc(n * sizeof(double));
        for (int i = n - 1; i >= 0; i--) {
            int base = n / nprocs;
            int rem = n % nprocs;
            int owner;
            int cum = 0;
            for (int p = 0; p < nprocs; p++) {
                int rows = (p < rem) ? base + 1 : base;
                if (i < cum + rows) {
                    owner = p;
                    break;
                }
                cum += rows;
            }
            double xi;
            if (rank == owner) {
                int local_index = i - localMat.global_start;
                xi = y[i];
                // Вычитаем сумму U(i,j)*x[j] для j от i+1 до n-1.
                for (int j = i + 1; j < n; j++) {
                    xi -= localMat.A[local_index * n + j] * x[j];
                }
                xi = xi / localMat.A[local_index * n + i];
            }
            MPI_Bcast(&xi, 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);
            x[i] = xi;
        }
        // Копируем глобальное решение x в локальный вектор решения (соответствующие строки)
        for (int i = 0; i < localMat.local_n; i++) {
            int global_row = localMat.global_start + i;
            localMat.X[i] = x[global_row];
        }
        free(y);
        free(x);
    }

} // end namespace band_matrix_mpi

#endif

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
//     * Пример параллельного LU-разложения для ленточной матрицы
//     * с помощью MPI + OpenMP.
//     *
//     * Упрощение: храним полную матрицу A на каждом процессе.
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
//        // Выделим L и U
//        DecomposeMatrix result;
//        result.l = (double**)malloc(n * sizeof(double*));
//        result.u = (double**)malloc(n * sizeof(double*));
//        for (size_t i = 0; i < n; i++) {
//            result.l[i] = (double*)calloc(n, sizeof(double));
//            result.u[i] = (double*)malloc(n * sizeof(double));
//        }
//
//        // Инициализация
//        // U = A, L = единичная диагональ
//        for (size_t i = 0; i < n; i++) {
//            for (size_t j = 0; j < n; j++) {
//                result.u[i][j] = matrix.A[i][j];
//                // L[i][j] уже заполнена нулями
//            }
//            result.l[i][i] = 1.0;
//        }
//
//        // Главный цикл LU-разложения
//        for (size_t k = 0; k < n - 1; k++) {
//            // Определяем, "кому" принадлежит строка k.
//            // Простейший блочный подход:
//            size_t blockSize = n / size;
//            size_t start = blockSize * rank;
//            size_t end = (rank == size - 1) ? (n - 1) : (blockSize * (rank + 1) - 1);
//
//            // Выясняем, лежит ли k в диапазоне строк данного процесса
//            int owner = (int)(k / blockSize);
//            if (owner >= size) owner = size - 1; // на случай, если n не делится точно
//
//            // Если я — владелец pivot-строки, я её нормирую и потом рассылаю
//            if (rank == owner) {
//                double pivotVal = result.u[k][k];
//                // Нормируем ведущую строку по диагональному элементу (строго говоря, 
//                // для классического LU обычно не делим всю строку, а вычисляем L[i][k].
//                // Но в некоторых реализациях допустимо сохранять уже нормированную ведущую строку)
//                // Здесь можно пропустить нормировку, т.к. классическая схема хранит pivot в U[k][k].
//                // Для согласованности с кодом выше будем использовать именно форму L[i][k] = U[i][k]/U[k][k].
//                // Но если хотим более классического подхода:
//                // for (size_t j = k; j < n; j++) {
//                //     result.u[k][j] /= pivotVal;
//                // }
//            }
//
//            // Теперь процесс-владелец рассылает pivot-строку k всем
//            MPI_Bcast(result.u[k], (int)n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
//
//            // Рассылаем также соответствующую часть L? Нет, L[k][k] = 1. Остальное не нужно.
//            // Идём по строкам от k+1 до k+b, обновляем их
//            // Каждый процесс делает это для "своих" строк
//#pragma omp parallel for
//            for (size_t i = start; i <= end; i++) {
//                if (i > k && i < MIN(k + b + 1, n)) {
//                    // Вычисляем L[i][k]
//                    result.l[i][k] = result.u[i][k] / result.u[k][k];
//                    // Обновляем U[i][j]
//                    for (size_t j = k; j < MIN(k + b + 1, n); j++) {
//                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
//                    }
//                }
//            }
//
//            // Чтобы синхронизировать U (ведь каждый процесс обновил только часть), 
//            // проще всего снова собрать U или переслать обновлённые строки. 
//            // В упрощённом примере будем собирать все строки i= k+1..k+b 
//            // (на самом деле лучше было бы не держать полный U на каждом процессе, 
//            // но для наглядности — сделаем Allgather).
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
//        // Упрощённая версия, в которой собираем всё решение на каждом процессе
//        // и делаем прямой/обратный ход локально. В реальности нужно аккуратно
//        // делать распределённый прямой/обратный ход.
//        int rank, size;
//        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        size_t n = matrix->n;
//
//        double* y = (double*)malloc(n * sizeof(double));
//        // Прямой ход: L * y = C
//        for (size_t i = 0; i < n; i++) {
//            double s = 0.0;
//            for (size_t j = 0; j < i; j++) {
//                s += decompose_matrix.l[i][j] * y[j];
//            }
//            y[i] = matrix->C[i] - s;
//        }
//
//        // Обратный ход: U * x = y
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
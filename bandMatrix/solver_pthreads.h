#define HAVE_STRUCT_TIMESPEC

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "solver_serial.h"
#include "matrix.h"

namespace band_matrix_pthreads
{
    typedef struct {
        int tid;
        int num_threads;
        size_t n;
        size_t b;
        struct DecomposeMatrix* result;
        pthread_barrier_t* barrier;
    } LUContext;

    void* lu_decomposition_worker(void* arg) {
        LUContext* ctx = (LUContext*)arg;
        size_t n = ctx->n;
        size_t b = ctx->b;
        int num_threads = ctx->num_threads;
        int tid = ctx->tid;

        for (size_t k = 0; k < n - 1; ++k) {
            size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;
            size_t total_rows = (upper_bound > k + 1) ? (upper_bound - (k + 1)) : 0;

            if (total_rows > 0) {
                size_t base = k + 1;
                size_t chunk = total_rows / num_threads;
                size_t remainder = total_rows % num_threads;
                size_t start = base + tid * chunk + (tid < remainder ? tid : remainder);
                size_t count = chunk + (tid < remainder ? 1 : 0);
                size_t end = start + count;

                for (size_t i = start; i < end; ++i) {
                    ctx->result->l[i][k] = ctx->result->u[i][k] / ctx->result->u[k][k];
                    for (size_t j = k; j < upper_bound; ++j) {
                        ctx->result->u[i][j] -= ctx->result->l[i][k] * ctx->result->u[k][j];
                    }
                }
            }
            pthread_barrier_wait(ctx->barrier);
        }
        return NULL;
    }

    struct DecomposeMatrix lu_decomposition(struct Matrix matrix, int num_threads) {
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

        if (n < 1000 || num_threads <= 1) {
            return band_matrix_serial::lu_decomposition(matrix);
        }

        pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
        LUContext* contexts = (LUContext*)malloc(num_threads * sizeof(LUContext));
        pthread_barrier_t barrier;
        pthread_barrier_init(&barrier, NULL, num_threads);

        for (int t = 0; t < num_threads; t++) {
            contexts[t].tid = t;
            contexts[t].num_threads = num_threads;
            contexts[t].n = n;
            contexts[t].b = b;
            contexts[t].result = &result;
            contexts[t].barrier = &barrier;
            pthread_create(&threads[t], NULL, lu_decomposition_worker, &contexts[t]);
        }

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }
        pthread_barrier_destroy(&barrier);
        free(threads);
        free(contexts);

        return result;
    }

    void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix) {
        size_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t j = 0; j < i; j++) {
                s += decompose_matrix.l[i][j] * y[j];
            }
            y[i] = matrix->C[i] - s;
        }

        matrix->X = (double*)malloc(n * sizeof(double));

        for (int i = n - 1; i >= 0; i--) {
            double s = 0.0;
            for (size_t j = i + 1; j < n; j++) {
                s += decompose_matrix.u[i][j] * matrix->X[j];
            }
            matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
        }

        free(y);
    }

    void solve_lu(struct DecomposeMatrix decompose_matrix, struct Matrix* matrix, int num_threads) {
        solve_lu(decompose_matrix, matrix);
    }

    void free_decompose_matrix(struct DecomposeMatrix* dm, size_t n) {
        for (size_t i = 0; i < n; i++) {
            free(dm->l[i]);
            free(dm->u[i]);
        }
        free(dm->l);
        free(dm->u);
    }
}


//// Структура контекста для потоков LU-разложения
//typedef struct {
//    int tid;                // номер потока (от 0 до num_threads-1)
//    int num_threads;        // общее число потоков
//    size_t n;               // размер матрицы
//    size_t b;               // ширина ленточной матрицы
//    DecomposeMatrix* result; // указатель на структуру с результатом разложения
//    pthread_barrier_t* barrier; // указатель на барьер синхронизации
//} LUContext;
//
//// Функция-воркер для LU-разложения (работают потоки из пула)
//void* lu_decomposition_worker(void* arg) {
//    LUContext* ctx = (LUContext*)arg;
//    size_t n = ctx->n;
//    size_t b = ctx->b;
//    int num_threads = ctx->num_threads;
//    int tid = ctx->tid;
//
//    // Внешний цикл по индексам k
//    for (size_t k = 0; k < n - 1; ++k) {
//        // Вычисляем, до какого столбца нужно обновлять (учитывая ленточную структуру)
//        size_t r = MIN(k + b + 1, n);
//        // Строки для обновления: от k+1 до r-1 (если таковых нет – total_rows==0)
//        size_t total_rows = (r > k + 1) ? (r - (k + 1)) : 0;
//
//        if (total_rows > 0) {
//            // Разбиваем диапазон [k+1, r) на num_threads частей.
//            size_t base = k + 1;
//            size_t chunk = total_rows / num_threads;
//            size_t remainder = total_rows % num_threads;
//            size_t start = base + tid * chunk + (tid < remainder ? tid : remainder);
//            size_t count = chunk + (tid < remainder ? 1 : 0);
//            size_t end = start + count;
//
//            for (size_t i = start; i < end; ++i) {
//                // Вычисляем множитель для LU-разложения
//                ctx->result->l[i][k] = ctx->result->u[i][k] / ctx->result->u[k][k];
//                // Обновляем элементы строки i в матрице U
//                for (size_t j = k; j < r; ++j) {
//                    ctx->result->u[i][j] -= ctx->result->l[i][k] * ctx->result->u[k][j];
//                }
//            }
//        }
//        // Синхронизация: ждём, пока все потоки завершат итерацию k
//        pthread_barrier_wait(ctx->barrier);
//    }
//    return NULL;
//}
//
//// Параллельное LU-разложение с использованием пула потоков
//DecomposeMatrix lu_decomposition_parallel(Matrix matrix, int num_threads) {
//    DecomposeMatrix result;
//    size_t n = matrix.n;
//    size_t b = matrix.b;
//
//    // Выделяем память для матриц L и U
//    result.l = (double**)malloc(n * sizeof(double*));
//    result.u = (double**)malloc(n * sizeof(double*));
//    for (size_t i = 0; i < n; i++) {
//        result.l[i] = (double*)calloc(n, sizeof(double));
//        result.u[i] = (double*)malloc(n * sizeof(double));
//    }
//    // Инициализируем: копируем A в U и задаём единичную диагональ в L
//    for (size_t i = 0; i < n; i++) {
//        for (size_t j = 0; j < n; j++) {
//            result.u[i][j] = matrix.A[i][j];
//        }
//        result.l[i][i] = 1.0;
//    }
//
//    // Если матрица слишком мала или число потоков <= 1 – можно использовать последовательное разложение.
//    if (n < 1000 || num_threads <= 1) {
//        // Здесь можно вызвать вашу последовательную функцию lu_decomposition(), если она уже есть.
//        band_matrix::lu_decomposition(matrix);
//        return result;
//    }
//
//    // Создаём пул потоков и барьер для синхронизации.
//    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
//    LUContext* contexts = (LUContext*)malloc(num_threads * sizeof(LUContext));
//    pthread_barrier_t barrier;
//    pthread_barrier_init(&barrier, NULL, num_threads);
//
//    // Инициализируем контексты и запускаем потоки.
//    for (int t = 0; t < num_threads; t++) {
//        contexts[t].tid = t;
//        contexts[t].num_threads = num_threads;
//        contexts[t].n = n;
//        contexts[t].b = b;
//        contexts[t].result = &result;
//        contexts[t].barrier = &barrier;
//        pthread_create(&threads[t], NULL, lu_decomposition_worker, &contexts[t]);
//    }
//
//    // Ожидаем завершения работы всех потоков.
//    for (int t = 0; t < num_threads; t++) {
//        pthread_join(threads[t], NULL);
//    }
//    pthread_barrier_destroy(&barrier);
//    free(threads);
//    free(contexts);
//
//    return result;
//}
//
//// Последовательное решение системы методом подстановок
//void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix) {
//    size_t n = matrix->n;
//    double* y = (double*)malloc(n * sizeof(double));
//
//    // Прямой ход: решаем L*y = C
//    for (size_t i = 0; i < n; i++) {
//        double s = 0.0;
//        for (size_t j = 0; j < i; j++) {
//            s += decompose_matrix.l[i][j] * y[j];
//        }
//        y[i] = matrix->C[i] - s;
//    }
//
//    // Обратный ход: решаем U*X = y
//    matrix->X = (double*)malloc(n * sizeof(double));
//    for (int i = n - 1; i >= 0; i--) {
//        double s = 0.0;
//        for (size_t j = i + 1; j < n; j++) {
//            s += decompose_matrix.u[i][j] * matrix->X[j];
//        }
//        matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
//    }
//
//    free(y);
//}
//
//// Для простоты оставляем решение системы последовательным.
//// Если же требуется параллельное решение, то его можно реализовать с использованием более сложных схем синхронизации.
//void solve_lu_parallel(DecomposeMatrix decompose_matrix, Matrix* matrix, int num_threads) {
//    solve_lu(decompose_matrix, matrix);
//}
//
//// Функция освобождения памяти, выделенной под разложение
//void free_decompose_matrix(DecomposeMatrix* dm, size_t n) {
//    for (size_t i = 0; i < n; i++) {
//        free(dm->l[i]);
//        free(dm->u[i]);
//    }
//    free(dm->l);
//    free(dm->u);
//}
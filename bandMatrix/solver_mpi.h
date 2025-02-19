#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

namespace band_matrix_mpi
{

    // Функция разворота массива (для обратной подстановки)
    void reverse_array(double* array, size_t n) {
        for (size_t i = 0; i < n / 2; ++i) {
            double temp = array[i];
            array[i] = array[n - 1 - i];
            array[n - 1 - i] = temp;
        }
    }

    // Модифицированная функция LU-разложения
    DecomposeMatrix lu_decomposition(Matrix matrix, MPI_Comm comm) {
        int rank, size;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        size_t n = matrix.n;
        size_t b = matrix.b;

        DecomposeMatrix result;
        result.l = (double**)malloc(n * sizeof(double*));
        result.u = (double**)malloc(n * sizeof(double*));
        for (size_t i = 0; i < n; i++) {
            result.l[i] = (double*)calloc(n, sizeof(double));
            result.u[i] = (double*)malloc(n * sizeof(double));
        }

        // Копируем исходную матрицу в U
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                result.u[i][j] = matrix.A[i][j];

        // Заполняем диагональ L единицами
        for (size_t i = 0; i < n; i++)
            result.l[i][i] = 1.0;

        // Основной цикл LU-разложения по столбцам
        for (size_t k = 0; k < n - 1; k++) {
            size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

            // Неблокирующий бродкаст строки k
            int pivot_owner = k % size;
            MPI_Request req_bcast;
            MPI_Ibcast(&(result.u[k][k]), upper_bound - k, MPI_DOUBLE, pivot_owner, comm, &req_bcast);
            MPI_Wait(&req_bcast, MPI_STATUS_IGNORE);

            // Каждый процесс обновляет свою часть строк
            for (size_t i = k + 1; i < upper_bound; i++) {
                if (((i - (k + 1)) % size) == rank) {
                    result.l[i][k] = result.u[i][k] / result.u[k][k];
                    for (size_t j = k; j < upper_bound; j++)
                        result.u[i][j] -= result.l[i][k] * result.u[k][j];
                }
            }

            // Подготовка буфера для обмена обновлениями строк
            size_t num_rows = upper_bound - (k + 1);
            size_t u_width = upper_bound - k; // число элементов в строке для U
            double* send_buffer = (double*)malloc(num_rows * (u_width + 1) * sizeof(double));
            for (size_t i = 0; i < num_rows; i++) {
                size_t global_row = k + 1 + i;
                if (((global_row - (k + 1)) % size) == rank) {
                    memcpy(&send_buffer[i * (u_width + 1)], &(result.u[global_row][k]), u_width * sizeof(double));
                    send_buffer[i * (u_width + 1) + u_width] = result.l[global_row][k];
                }
                else {
                    memset(&send_buffer[i * (u_width + 1)], 0, (u_width + 1) * sizeof(double));
                }
            }

            // Используем неблокирующий Allgather для обмена обновлениями
            double* recv_buffer = (double*)malloc(num_rows * (u_width + 1) * size * sizeof(double));
            MPI_Request req_allgather;
            MPI_Iallgather(send_buffer, num_rows * (u_width + 1), MPI_DOUBLE,
                recv_buffer, num_rows * (u_width + 1), MPI_DOUBLE, comm, &req_allgather);
            MPI_Wait(&req_allgather, MPI_STATUS_IGNORE);

            // Обновляем локальные копии строк
            for (size_t i = 0; i < num_rows; i++) {
                int owner = (i % size);
                size_t index = owner * num_rows * (u_width + 1) + i * (u_width + 1);
                size_t global_row = k + 1 + i;
                memcpy(&(result.u[global_row][k]), &recv_buffer[index], u_width * sizeof(double));
                result.l[global_row][k] = recv_buffer[index + u_width];
            }
            free(send_buffer);
            free(recv_buffer);
        }
        return result;
    }

    // Модифицированная функция решения системы методом LU-разложения
    void solve_lu(DecomposeMatrix decomp, Matrix* matrix, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        size_t n = matrix->n;
        double* y = (double*)malloc(n * sizeof(double));

        if (rank == 0) {
            // Прямая подстановка: решаем L*y = C
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = 0; j < i; j++)
                    s += decomp.l[i][j] * y[j];
                y[i] = matrix->C[i] - s;
            }
            matrix->X = (double*)malloc(n * sizeof(double));
            // Обратная подстановка: решаем U*X = y (используется обратный порядок с последующим разворотом)
            for (int i = n - 1, k = 0; i >= 0; --i, k++) {
                double s = 0.0;
                for (int j = n - 1; j > i; --j)
                    s += decomp.u[i][j] * matrix->X[n - j - 1];
                matrix->X[k] = (y[i] - s) / decomp.u[i][i];
            }
            reverse_array(matrix->X, n);
        }
        if (rank != 0)
            matrix->X = (double*)malloc(n * sizeof(double));

        // Неблокирующий бродкаст решения
        MPI_Request req;
        MPI_Ibcast(matrix->X, n, MPI_DOUBLE, 0, comm, &req);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        free(y);
    }
}
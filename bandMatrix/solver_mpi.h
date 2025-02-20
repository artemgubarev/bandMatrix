#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

namespace band_matrix_mpi
{
	//// Функция разворота массива (для обратной подстановки)
	//void reverse_array(double* array, size_t n) {
	//    for (size_t i = 0; i < n / 2; ++i) {
	//        double temp = array[i];
	//        array[i] = array[n - 1 - i];
	//        array[n - 1 - i] = temp;
	//    }
	//}

	//// Модифицированная функция LU-разложения
	//DecomposeMatrix lu_decomposition(Matrix matrix, MPI_Comm comm) {
	//    int rank, size;
	//    MPI_Comm_rank(comm, &rank);
	//    MPI_Comm_size(comm, &size);

	//    size_t n = matrix.n;
	//    size_t b = matrix.b;

	//    DecomposeMatrix result;
	//    result.l = (double**)malloc(n * sizeof(double*));
	//    result.u = (double**)malloc(n * sizeof(double*));
	//    for (size_t i = 0; i < n; i++) {
	//        result.l[i] = (double*)calloc(n, sizeof(double));
	//        result.u[i] = (double*)malloc(n * sizeof(double));
	//    }

	//    // Копируем исходную матрицу в U
	//    for (size_t i = 0; i < n; i++)
	//        for (size_t j = 0; j < n; j++)
	//            result.u[i][j] = matrix.A[i][j];

	//    // Заполняем диагональ L единицами
	//    for (size_t i = 0; i < n; i++)
	//        result.l[i][i] = 1.0;

	//    // Основной цикл LU-разложения по столбцам
	//    for (size_t k = 0; k < n - 1; k++) {
	//        size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

	//        // Неблокирующий бродкаст строки k
	//        int pivot_owner = k % size;
	//        MPI_Request req_bcast;
	//        MPI_Ibcast(&(result.u[k][k]), upper_bound - k, MPI_DOUBLE, pivot_owner, comm, &req_bcast);
	//        MPI_Wait(&req_bcast, MPI_STATUS_IGNORE);

	//        // Каждый процесс обновляет свою часть строк
	//        for (size_t i = k + 1; i < upper_bound; i++) {
	//            if (((i - (k + 1)) % size) == rank) {
	//                result.l[i][k] = result.u[i][k] / result.u[k][k];
	//                for (size_t j = k; j < upper_bound; j++)
	//                    result.u[i][j] -= result.l[i][k] * result.u[k][j];
	//            }
	//        }

	//        // Подготовка буфера для обмена обновлениями строк
	//        size_t num_rows = upper_bound - (k + 1);
	//        size_t u_width = upper_bound - k; // число элементов в строке для U
	//        double* send_buffer = (double*)malloc(num_rows * (u_width + 1) * sizeof(double));
	//        for (size_t i = 0; i < num_rows; i++) {
	//            size_t global_row = k + 1 + i;
	//            if (((global_row - (k + 1)) % size) == rank) {
	//                memcpy(&send_buffer[i * (u_width + 1)], &(result.u[global_row][k]), u_width * sizeof(double));
	//                send_buffer[i * (u_width + 1) + u_width] = result.l[global_row][k];
	//            }
	//            else {
	//                memset(&send_buffer[i * (u_width + 1)], 0, (u_width + 1) * sizeof(double));
	//            }
	//        }

	//        // Используем неблокирующий Allgather для обмена обновлениями
	//        double* recv_buffer = (double*)malloc(num_rows * (u_width + 1) * size * sizeof(double));
	//        MPI_Request req_allgather;
	//        MPI_Iallgather(send_buffer, num_rows * (u_width + 1), MPI_DOUBLE,
	//            recv_buffer, num_rows * (u_width + 1), MPI_DOUBLE, comm, &req_allgather);
	//        MPI_Wait(&req_allgather, MPI_STATUS_IGNORE);

	//        // Обновляем локальные копии строк
	//        for (size_t i = 0; i < num_rows; i++) {
	//            int owner = (i % size);
	//            size_t index = owner * num_rows * (u_width + 1) + i * (u_width + 1);
	//            size_t global_row = k + 1 + i;
	//            memcpy(&(result.u[global_row][k]), &recv_buffer[index], u_width * sizeof(double));
	//            result.l[global_row][k] = recv_buffer[index + u_width];
	//        }
	//        free(send_buffer);
	//        free(recv_buffer);
	//    }
	//    return result;
	//}

	//// Модифицированная функция решения системы методом LU-разложения
	//void solve_lu(DecomposeMatrix decomp, Matrix* matrix, MPI_Comm comm) {
	//    int rank;
	//    MPI_Comm_rank(comm, &rank);
	//    size_t n = matrix->n;
	//    double* y = (double*)malloc(n * sizeof(double));

	//    if (rank == 0) {
	//        // Прямая подстановка: решаем L*y = C
	//        for (size_t i = 0; i < n; i++) {
	//            double s = 0.0;
	//            for (size_t j = 0; j < i; j++)
	//                s += decomp.l[i][j] * y[j];
	//            y[i] = matrix->C[i] - s;
	//        }
	//        matrix->X = (double*)malloc(n * sizeof(double));
	//        // Обратная подстановка: решаем U*X = y (используется обратный порядок с последующим разворотом)
	//        for (int i = n - 1, k = 0; i >= 0; --i, k++) {
	//            double s = 0.0;
	//            for (int j = n - 1; j > i; --j)
	//                s += decomp.u[i][j] * matrix->X[n - j - 1];
	//            matrix->X[k] = (y[i] - s) / decomp.u[i][i];
	//        }
	//        reverse_array(matrix->X, n);
	//    }
	//    if (rank != 0)
	//        matrix->X = (double*)malloc(n * sizeof(double));

	//    // Неблокирующий бродкаст решения
	//    MPI_Request req;
	//    MPI_Ibcast(matrix->X, n, MPI_DOUBLE, 0, comm, &req);
	//    MPI_Wait(&req, MPI_STATUS_IGNORE);
	//    free(y);
	//}

	void reverse_array(double* array, size_t n) {
		for (size_t i = 0; i < n / 2; i++) {
			double tmp = array[i];
			array[i] = array[n - 1 - i];
			array[n - 1 - i] = tmp;
		}
	}

	/* Параллельное LU-разложение ленточной матрицы.
	   Строки распределяются блочно между процессами. */
	DecomposeMatrix lu_decomposition(Matrix matrix) {
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		size_t n = matrix.n;
		size_t b = matrix.b;

		DecomposeMatrix result;
		result.l = (double**)malloc(n * sizeof(double*));
		result.u = (double**)malloc(n * sizeof(double*));
		for (size_t i = 0; i < n; i++) {
			result.l[i] = (double*)calloc(n, sizeof(double));
			result.u[i] = (double*)malloc(n * sizeof(double));
		}

		// Копируем исходную матрицу A в U
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				result.u[i][j] = matrix.A[i][j];
			}
		}

		// Диагональные элементы L равны 1
		for (size_t i = 0; i < n; i++) {
			result.l[i][i] = 1.0;
		}

		// Определяем диапазон строк, обрабатываемых данным процессом
		size_t start_row = (n * rank) / size;
		size_t end_row = (n * (rank + 1)) / size;

		// Основной цикл факторизации
		for (size_t k = 0; k < n - 1; k++) {
			// Обновление выполняется только в пределах ленты: от k до upper_bound
			size_t upper_bound = (k + b + 1 < n) ? (k + b + 1) : n;

			// Определяем процесс-владельца строки k
			int owner = (k * size) / n;

			// Буфер для передачи части строки k (от столбца k до upper_bound)
			int segment_length = upper_bound - k;
			double* u_row_segment = (double*)malloc(segment_length * sizeof(double));
			if (rank == owner) {
				for (size_t j = k; j < upper_bound; j++) {
					u_row_segment[j - k] = result.u[k][j];
				}
			}
			// Используем неблокирующий Bcast для передачи сегмента
			MPI_Request req;
			MPI_Ibcast(u_row_segment, segment_length, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);

			// Каждый процесс обновляет свои строки, попадающие в [max(k+1, start_row), end_row)
			size_t local_start = (k + 1 > start_row) ? k + 1 : start_row;
			for (size_t i = local_start; i < end_row && i < upper_bound; i++) {
				result.l[i][k] = result.u[i][k] / u_row_segment[0];
				for (size_t j = k; j < upper_bound; j++) {
					result.u[i][j] -= result.l[i][k] * u_row_segment[j - k];
				}
			}
			free(u_row_segment);
			/* Убираем явный MPI_Barrier – синхронизация происходит при MPI_Ibcast */
		}
		return result;
	}

	/* Параллельное решение системы методом LU-разложения.
	   Здесь прямой и обратный ходы распараллеливаются по строкам:
		 - Каждая строка принадлежит определённому процессу,
		   который вычисляет соответствующую компоненту,
		   после чего значение передаётся всем с помощью неблокирующего Bcast.
	*/
	void solve_lu(DecomposeMatrix decompose_matrix, Matrix* matrix) {
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		size_t n = matrix->n;

		// Определяем распределение строк
		size_t start_row = (n * rank) / size;
		size_t end_row = (n * (rank + 1)) / size;

		double* y = (double*)malloc(n * sizeof(double));

		/* Прямой ход: решаем L*y = C.
		   Для каждой строки i:
			 - Если строка принадлежит текущему процессу, вычисляем y[i] = C[i] - sum(L[i][j]*y[j])
			 - Затем с помощью неблокирующего Bcast передаём y[i] всем.
			 Таким образом зависимость сохраняется, а вычисления распределяются.
		*/
		for (size_t i = 0; i < n; i++) {
			if (i >= start_row && i < end_row) {
				double s = 0.0;
				for (size_t j = 0; j < i; j++) {
					s += decompose_matrix.l[i][j] * y[j];
				}
				y[i] = matrix->C[i] - s;
			}
			int owner = (i * size) / n;  // процесс, которому принадлежит строка i
			MPI_Request req;
			MPI_Ibcast(&y[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}

		/* Обратный ход: решаем U*x = y.
		   Итерация идёт от последней строки к первой.
		   Для каждой строки i:
			 - Если строка принадлежит процессу, вычисляем x[i] = (y[i] - sum(U[i][j]*x[j])) / U[i][i]
			 - Рассылаем x[i] всем процессам.
		*/
		if (matrix->X == NULL) {
			matrix->X = (double*)malloc(n * sizeof(double));
		}
		for (int i = n - 1; i >= 0; i--) {
			if (i >= start_row && i < end_row) {
				double s = 0.0;
				for (int j = i + 1; j < n; j++) {
					s += decompose_matrix.u[i][j] * matrix->X[j];
				}
				matrix->X[i] = (y[i] - s) / decompose_matrix.u[i][i];
			}
			int owner = (i * size) / n;
			MPI_Request req;
			MPI_Ibcast(&matrix->X[i], 1, MPI_DOUBLE, owner, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, MPI_STATUS_IGNORE);
		}

		free(y);
	}
}
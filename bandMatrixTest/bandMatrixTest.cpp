#include <iostream>
#include <cstdlib>
#include <time.h>
#include "comparator.h"
#include "../bandMatrix/solver.h"
#include "../bandMatrix/matrix.h"
#include "../bandMatrix/writer.h"

int main()
{
    const char* filename = getenv("INPUT_MATRIX_FILE");
    if (!filename)
    {
        std::cerr << "Error: environment variable INPUT_MATRIX_FILE not set.\n";
        return 1;
    }

    clock_t start, end;
    double cpu_time_used;

    Matrix matrix = read_matrix(filename);

    start = clock();

    DecomposeMatrix decompose = band_matrix::lu_decomposition(matrix);
    band_matrix::solve_lu(decompose, &matrix);

    end = clock();

    write_1d("solution.txt", matrix.X, matrix.n);

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time: %f sec.\n", cpu_time_used);

    double epsilon = 0.00001;
    double numbers1[MAX_NUMBERS], numbers2[MAX_NUMBERS];
    size_t count1, count2;

    if (!load_numbers("solution.txt", numbers1, &count1) || !load_numbers("msolution.txt", numbers2, &count2))
    {
        return 1;
    }

    if (compare_numbers(numbers1, numbers2, count1, count2, epsilon))
    {
        printf("\033[32mTest Correct\033[0m\n");
    }
    else
    {
        printf("\033[31mTest Failed (Red)\033[0m\n");
    }

    return 0;      
}
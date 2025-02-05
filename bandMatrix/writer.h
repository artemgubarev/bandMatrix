#include <stdio.h>

void write_1d(const char* filename, double* array, size_t n)
{
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    for (size_t i = 0; i < n; i++)
    {
        fprintf(file, "%lf ", array[i]);
    }

    fclose(file);
}

void write_2d(const char* filename, double** array, size_t n)
{
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(file, "%lf ", array[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}
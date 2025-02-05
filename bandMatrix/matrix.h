#pragma once

#include <stdlib.h>
#include <stdio.h>

struct DecomposeMatrix
{
    double** l;
    double** u;
};

struct Matrix
{
    size_t n;
    size_t b;
    double* C;
    double** A;
    double* X;
};

double** allocate_matrix(size_t n) 
{
    double** matrix = (double**)malloc(n * sizeof(double*));
    for (size_t i = 0; i < n; i++) 
    {
        matrix[i] = (double*)malloc(n * sizeof(double));
    }
    return matrix;
}

void free_matrix(double** matrix, size_t n) 
{
    for (size_t i = 0; i < n; i++) 
    {
        free(matrix[i]);
    }
    free(matrix);
}

struct Matrix read_matrix(const char* filename) 
{
    struct Matrix mat;
    FILE* file = fopen(filename, "r");
    if (!file) 
    {
        perror("Error openning file.\n");
        perror(filename);


        exit(EXIT_FAILURE);
    }

    fscanf(file, "%zu", &mat.n);
    fscanf(file, "%zu", &mat.b);

    mat.A = allocate_matrix(mat.n);
    mat.C = (double*)malloc(mat.n * sizeof(double));

    for (size_t i = 0; i < mat.n; i++) 
    {
        for (size_t j = 0; j < mat.n; j++) 
        {
            fscanf(file, "%lf", &mat.A[i][j]);
        }
    }

    for (size_t i = 0; i < mat.n; i++) 
    {
        fscanf(file, "%lf", &mat.C[i]);
    }
    fclose(file);
    return mat;
}

void print_matrix(Matrix mat) 
{
    printf("n: %zu\n", mat.n);
    printf("b: %zu\n", mat.b);
    printf("A:\n");
    for (size_t i = 0; i < mat.n; i++) 
    {
        for (size_t j = 0; j < mat.n; j++) 
        {
            printf("%lf ", mat.A[i][j]);
        }
        printf("\n");
    }
    printf("C:\n");
    for (size_t i = 0; i < mat.n; i++) 
    {
        printf("%lf ", mat.C[i]);
    }
    printf("\n");
}
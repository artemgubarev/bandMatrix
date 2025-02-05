#pragma once
#include <stdio.h>

void print_1d(double* array, size_t n)
{
    for (size_t i = 0; i < n; i++)
    {
        printf("%lf ", array[i]);
    }
}

void print_2d(double** array, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            printf("%lf ", array[i][j]);
        }
        printf("\n");
    }
}
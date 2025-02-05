#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NUMBERS 1000000

int load_numbers(const char* filename, double* numbers, size_t* count)
{
    FILE* file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return 0;
    }

    *count = 0;
    while (fscanf(file, "%lf", &numbers[*count]) == 1)
    {
        (*count)++;
        if (*count >= MAX_NUMBERS)
        {
            break;
        }
    }

    fclose(file);
    return 1;
}

int compare_numbers(const double* vec1, const double* vec2, size_t count1, size_t count2, double epsilon)
{
    if (count1 != count2)
    {
        fprintf(stderr, "Error: Files contain different number of values!\n");
        return 0;
    }

    for (size_t i = 0; i < count1; i++)
    {
        if (fabs(vec1[i] - vec2[i]) > epsilon) {
            printf("Mismatch at index %zu: %lf vs %lf (diff = %lf)\n", i, vec1[i], vec2[i], fabs(vec1[i] - vec2[i]));
            return 0;
        }
    }
    return 1;
}
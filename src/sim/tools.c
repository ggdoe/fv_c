#include "sim.h"

void linspace(real* dest, real min, real max, u32 n)
{
    real dx = (max - min) / (n-1);

    for (int i = 0; i < n; i++)
       dest[i] = min + dx * i;
}

void print_mat(real* mat, u32 n, u32 m)
{
    for (u32 j = 0; j < n; j++)
    {
        for (u32 i = 0; i < m; i++)
            printf("%6.3f ", mat[j * m + i]);
        printf("\n");
    }
}


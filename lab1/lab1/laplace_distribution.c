#include "laplace_distribution.h"

double LaplaceIntervalP(int a, int b)
{
    return LaplaceFn(b) - LaplaceFn(a);
}

double LaplaceFn(int x)
{
    int mu = 5;
    int b = 4;
    if (x < mu) return exp( (x - mu) / b ) / 2;
    return 1 - exp( (-1) * ((double)x - (double)mu) / (double)b ) / 2;
}

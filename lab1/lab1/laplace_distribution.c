#include "laplace_distribution.h"

int LaplaceIntervalP(int a, int b)
{
    return LaplaceFn(b) - LaplaceFn(a);
}

int LaplaceFn(int x)
{
    int mu = 5;
    int b = 4;
    if (x < mu) return PositiveLaplace(x, mu, b);
    return NegativeLaplace(x, mu, b);
}

int PositiveLaplace(int x, int mu, int b)
{
    return 1/2 * exp( (x - mu) / b );
}

int NegativeLaplace(int x, int mu, int b)
{
    return 1/2 * exp( (-1) * (x - mu) / b );
}

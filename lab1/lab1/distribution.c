// файл distribution.c

#include "distribution.h"

// параметри для обчислень
#define MU_PARAM 6
#define SIGMA_PARAM 7
#define D1_PARAM 3
#define D2_PARAM 2

#define MAXIT 200
#define EPS 3.0e-7
#define FPMIN 1.0e-30

// P(a, b) = F(b) - F(a)
double ZFisherDistP(int a, int b)
{
    return ZFisherDistFunc(b) - ZFisherDistFunc(a);
}

// F = Ix*
double ZFisherDistFunc(int x)
{
    // обчислення x*
    double z = ((double)x - (double)MU_PARAM) / (double)SIGMA_PARAM;
    double xi = CountXI(z, D1_PARAM, D2_PARAM);
    // обчислення I
    double incompleteBeta = IncompleteBetaFunction((double)(D1_PARAM / 2), (double)(D2_PARAM / 2), xi);
    double beta = BetaFunction((double)(D1_PARAM / 2), (double)(D2_PARAM / 2));
    return incompleteBeta / beta;
}

// допоміжні функції для обчислення функції розподілу

double BetaFunction(double a, double b)
{
    return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

double IncompleteBetaFunction(double a, double b, double x)
{
    if (x < 0 || x > 1) return 1.0 / 0;
    double beta;
    if (x == 0 || x == 1) beta = 0;
    beta = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return beta * IncompleteBetaContinuousFraction(a, b, x) / a;
    else return 1 - beta * IncompleteBetaContinuousFraction(b, a, 1 - x) / b;
}

double IncompleteBetaContinuousFraction(double a, double b, double x)
{
    int m, m2;
    float aa, c, d, del, h, qab, qam, qap;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;

    c = 1.0;
    d = 1.0 - qab * x / qap;

    if (fabs(d) < FPMIN) d=FPMIN;
    d = 1.0 / d;
    h = d;
    for (m = 1; m <= MAXIT; m++)
    {
        m2 = 2 * m;
        aa = m * (b - m)*x / ((qam + m2)*(a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m)*(qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c = 1.0 + aa/c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;

        if (fabs(del - 1.0) < EPS) break;
    }
    if (m > MAXIT) return 1.0 / 0;
    return h;
}

// x*
double CountXI(double z, double d1, double d2)
{
    double a = d2 * exp(2 * z);
    double b = d1 + d2 * exp(2 * z);
    return a / b;
}

// файл - intervals.c

#include "intervals.h"
#include "distribution.h"

#define MIN_INTERVAL_SIZE 2

// розділення на можливі інтервали
int DivideIntoIntervalIndexes(int n, int m, int possibleIntervals[n * n * 2][m])
{
    int intervalSizes[m];
    for (int i = 0; i < m; i++) intervalSizes[i] = i == 0 ? n - (m - 1) * MIN_INTERVAL_SIZE : MIN_INTERVAL_SIZE;
    return GetPossibleIntervalIndexes(n, m, possibleIntervals, intervalSizes, 0, 0);
}


int GetPossibleIntervalIndexes(int n, int m, int possibleIntervals[n * n * 2][m], int intervalSizes[], int sizeIndex, int counter)
{
    int localIntervalSizes[m];
    for (int i = 0; i < m; i++) localIntervalSizes[i] = intervalSizes[i];
    while(localIntervalSizes[sizeIndex] >= MIN_INTERVAL_SIZE)
    {
        if (sizeIndex + 1 < m - 1)
        {
            counter = GetPossibleIntervalIndexes(n, m, possibleIntervals, localIntervalSizes, sizeIndex + 1, counter);
        }
        if (sizeIndex + 1 == m-1)
        {
            int sumSize = 0;
            for (int i = 0; i < m; i++)
            {
                possibleIntervals[counter][i] = sumSize;
                sumSize += localIntervalSizes[i];
            }
            counter += 1;
        }
        localIntervalSizes[sizeIndex] -= 1;
        localIntervalSizes[sizeIndex + 1] += 1;
    }
    return counter;
}

// вибір підходящої комбінації інтервалів
void GetMostProbableInterval(int m, int intervals[m][2], int n, int possibleIntervalIndexes[n * n * 2][m], int numbers[], int count)
{
    double currentProbability = 0.0;
    double highestProbability = 0.0;
    for (int i = 0; i < count; i++)
    {
        currentProbability = 0.0;
        for (int j = 0; j < m; j++)
        {
            int a = numbers[possibleIntervalIndexes[i][j]];
            int b = (j == m - 1 ? numbers[n - 1] : numbers[possibleIntervalIndexes[i][j + 1] - 1]);
            currentProbability += ZFisherDistP(a, b);
        }
        if (currentProbability > highestProbability && currentProbability <= 1) {
            highestProbability = currentProbability;
            for (int j = 0; j < m; j++)
            {
                intervals[j][0] = numbers[possibleIntervalIndexes[i][j]];
                intervals[j][1] = (j == m - 1 ? numbers[n - 1] : numbers[possibleIntervalIndexes[i][j + 1] - 1]);
            }
            if (currentProbability == 1.0) return;
        }
    }
}


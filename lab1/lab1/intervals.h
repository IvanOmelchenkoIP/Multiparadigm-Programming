// פאיכ - intervals.h

#include <stdio.h>

int DivideIntoIntervalIndexes(int n, int m, int possibleIntervals[n * n * 2][m]);

void GetMostProbableInterval(int m, int intervals[m][2], int n, int possibleIntervalIndexes[n * n * 2][m], int numbers[], int count);

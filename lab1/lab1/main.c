// ���� - main.c

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "merge_sort.h"
#include "utils.h"
#include "alphabet.h"
#include "intervals.h"

#include "distribution.h"

#define ALPHABET_POWER 4

int main()
{
    // �������� ��� � ������
    int numbers[] = { 13, 1, 2, 15, -10, 20, 3, 12, 4, 5, 6, 26, 7, 8 };
    char * alphabet[ALPHABET_POWER] = { "A", "B", "C", "D" };

    // ������� ���������� ���� ��������� ��������
    clock_t start = clock();

    // ���������� ���������� ������
    int numbersCount = sizeof(numbers) / sizeof(numbers[0]);
    int sortedNumbers[numbersCount];

    int possibleIntervalIndexes[numbersCount * numbersCount * 2][ALPHABET_POWER];
    int intervals[ALPHABET_POWER][2];

    // ��������� � ���������� ��ﳿ ��������� ����
    CopyIntArray(numbers, sortedNumbers, numbersCount);
    MergeSort(sortedNumbers, 0, numbersCount - 1);

    // �������� �� ������ ���������
    int count = DivideIntoIntervalIndexes(numbersCount, ALPHABET_POWER, possibleIntervalIndexes);
    // ���� ���������� ���������'
    GetMostProbableInterval(ALPHABET_POWER, intervals, numbersCount, possibleIntervalIndexes, sortedNumbers, count);

    // ������������ �� �������� ���, ��������� �� �����
    char * lingChain[numbersCount];
    Translate(numbers, alphabet, numbersCount, ALPHABET_POWER, intervals, lingChain);
    PrintLingChain(lingChain, numbersCount);

    // ����� ��������� ��������
    clock_t end = clock() - start;
    double seconds = ((double) end) / CLOCKS_PER_SEC;
    printf("Program finished in %f seconds", seconds);

    char key[1];
    printf("\n\nPress any key to exit... ");
    scanf("%c", &key);
}

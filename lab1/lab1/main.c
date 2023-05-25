// файл - main.c

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
    // цисловий ряд і алфавіт
    int numbers[] = { 13, 1, 2, 15, -10, 20, 3, 12, 4, 5, 6, 26, 7, 8 };
    char * alphabet[ALPHABET_POWER] = { "A", "B", "C", "D" };

    // початок вимірювання часу виконання програми
    clock_t start = clock();

    // декларація нербхідних масивів
    int numbersCount = sizeof(numbers) / sizeof(numbers[0]);
    int sortedNumbers[numbersCount];

    int possibleIntervalIndexes[numbersCount * numbersCount * 2][ALPHABET_POWER];
    int intervals[ALPHABET_POWER][2];

    // сопіювання і сортування копії числового ряду
    CopyIntArray(numbers, sortedNumbers, numbersCount);
    MergeSort(sortedNumbers, 0, numbersCount - 1);

    // розбиття на можливі інтервали
    int count = DivideIntoIntervalIndexes(numbersCount, ALPHABET_POWER, possibleIntervalIndexes);
    // вибір підходящого інтервалу'
    GetMostProbableInterval(ALPHABET_POWER, intervals, numbersCount, possibleIntervalIndexes, sortedNumbers, count);

    // перетворення на числовий ряд, виведення на екран
    char * lingChain[numbersCount];
    Translate(numbers, alphabet, numbersCount, ALPHABET_POWER, intervals, lingChain);
    PrintLingChain(lingChain, numbersCount);

    // кінець виконання програми
    clock_t end = clock() - start;
    double seconds = ((double) end) / CLOCKS_PER_SEC;
    printf("Program finished in %f seconds", seconds);

    char key[1];
    printf("\n\nPress any key to exit... ");
    scanf("%c", &key);
}

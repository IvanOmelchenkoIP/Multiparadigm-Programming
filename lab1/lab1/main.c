#include <stdio.h>
#include <stdlib.h>

#include "merge_sort.h"
#include "laplace_distribution.h"

#define ALPHABET_POWER 6

int main()
{
    char * alphabet[ALPHABET_POWER] = { "A", "B", "C", "D", "E", "F" };
    int numbers[] = {13, 1, 2, 15, 3, 12, 4, 5, 6, 7, 8, 14, 9, 10};
    int interval[ALPHABET_POWER][2];
    int numbersToIntervals[sizeof(numbers) / sizeof(numbers[0])];
    char * lingChain[sizeof(numbers) / sizeof(numbers[0])];
    MergeSort(numbers, 0, sizeof(numbers) / sizeof(numbers[0]) - 1);
    printf("Hello world!\n");
    for (int i = 0; i < sizeof(numbers) / sizeof(numbers[0]); i++) {
        printf("%d %f\n", numbers[i], LaplaceFn(numbers[i]));
    }


    return 0;
}
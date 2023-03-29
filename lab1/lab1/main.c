#include <stdio.h>
#include <stdlib.h>

#include "quick_sort.h"

int main()
{
    int numbers[] = {13, 1, 2, 15, 3, 12, 4, 5, 6, 7, 8, 14, 9, 10};

    MergeSort(numbers, 0, sizeof(numbers) / sizeof(numbers[0]) - 1);

    printf("Hello world!\n");
    for (int i = 0; i < sizeof(numbers) / sizeof(numbers[0]); i++) {
        printf("%d ", numbers[i]);
    }
    return 0;
}

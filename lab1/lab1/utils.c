// файл - utils.c

#include "utils.h"

// копіювання масиву цілочисельних значень
void CopyIntArray(int source[], int destination[], int elements)
{
    for (int i = 0; i < elements; i++) destination[i] = source[i];
}


#include "utils.h"

int intArrLength(int arr[])
{
    return sizeof(arr) / sizeof(arr[0]);
}

int charArrLength(char * arr[])
{
    return sizeof(arr) / sizeof(arr[0])
}

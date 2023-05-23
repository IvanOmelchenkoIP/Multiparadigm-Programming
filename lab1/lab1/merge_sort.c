// файл - merge_sort.c

#include "merge_sort.h"

// сортування масиву за алгоритмом merge sort
void MergeSort(int arr[], int leftInd, int rightInd)
{
    if (leftInd < rightInd)
    {
        int middleInd = leftInd + (rightInd - leftInd) / 2;
        MergeSort(arr, leftInd, middleInd);
        MergeSort(arr, middleInd + 1, rightInd);
        Merge(arr, leftInd, middleInd, rightInd);
    }
}

void Merge(int arr[], int leftInd, int middleInd, int rightInd)
{
    int leftN = middleInd - leftInd + 1;
    int rightN = rightInd - middleInd;
    int leftArr[leftN], rightArr[rightN];
    for (int i = 0; i < leftN; i++) leftArr[i] = arr[leftInd + i];
    for (int j = 0; j < rightN; j++) rightArr[j] = arr[middleInd + 1 + j];

    int leftArrInd = 0;
    int rightArrInd = 0;
    int arrInd = leftInd;
    while (leftArrInd < leftN && rightArrInd < rightN)
    {
        if (leftArr[leftArrInd] <= rightArr[rightArrInd]) arr[arrInd++] = leftArr[leftArrInd++];
        else arr[arrInd++] = rightArr[rightArrInd++];
    }
    while (leftArrInd < leftN) arr[arrInd++] = leftArr[leftArrInd++];
    while (rightArrInd < rightN) arr[arrInd++] = rightArr[rightArrInd++];
}


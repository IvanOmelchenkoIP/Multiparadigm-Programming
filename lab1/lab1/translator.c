#include "translator.h"

void translate(int numIntervals[], char * alphabet[], char * lingChain[])
{
    for (int i = 0; i < sizeof(numIntervals) / sizeof(numIntervals[0]); i++)
    {
        lingChain[i] = alphabet[numIntervals[i]];
    }
}

void printLingChain(char * lingChain[])
{
    printf("Linguistic Chain: \n");
    for (i = 0; i < sizeof(lingChain) / sizeof(lingChain[0]); i++)
    {
        printf("%s ", lingchain[i]);
    }
    printf("\n");
}


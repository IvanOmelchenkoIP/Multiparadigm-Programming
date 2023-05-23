// файл - alphabet.c

#include "alphabet.h"

// переклад чисельного ряду
void Translate(int numbers[], char * alphabet[], int n, int m, int intervals[m][2], char * lingChain[])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (intervals[j][0] <= numbers[i] && numbers[i] <= intervals[j][1])
            {
                lingChain[i] = alphabet[j];
                break;
            }
        }
    }
}

// виведення лінгвістичного ряду
void PrintLingChain(char * lingChain[], int n)
{
    printf("Linguistic Chain: \n");
    for (int i = 0; i < n; i++) printf("%s ", lingChain[i]);
    printf("\n");
}



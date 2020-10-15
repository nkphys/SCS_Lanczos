#include "Base_to_Decimal.h"
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;



int val(char c)
{
    if (c >= '0' && c <= '9')
        return (int)c - '0';
    else
        return (int)c - 'A' + 10;
}

bool shouldSwap(char str[], int start, int curr)
{
    for (int i = start; i < curr; i++){
        if (str[i] == str[curr])
            return 0;
        }
    return 1;
}

char reVal(int num)
{
    if (num >= 0 && num <= 9)
        return (char)(num + '0');
    else
        return (char)(num - 10 + 'A');
}

// Utility function to reverse a string
void strev(char *str)
{
    int len = strlen(str);
    int i;
    for (i = 0; i < len/2; i++)
    {
        char temp = str[i];
        str[i] = str[len-i-1];
        str[len-i-1] = temp;
    }
}

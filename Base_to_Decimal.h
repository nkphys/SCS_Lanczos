#include <assert.h>
#include <stdio.h> 
#include <string.h>

int val(char c); 
char reVal(int num); 
// Utility function to reverse a string 
void strev(char *str); 
// Function to convert a number from given base 'b' 
// to decimal 
int toDeci(char *str, int base);
// Function to convert a given decimal number 
// to a base 'base' and 
char* fromDeci(char res[], int base, int inputNum);
int value_at_pos(int n, int pos, int base);
int Updated_decimal_with_value_at_pos(int n, int pos, int base, int new_value);
//ADD functions analogous to following:
/*
int bit_value(int n, int pos);
*/

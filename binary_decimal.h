/*
*/
#include <complex>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "tensor_type.h"
#include <stdlib.h>


int Act_Translation_assuming_PBC(int dec_state, int l1_, int l2_ );
int Act_Inv_Translation_assuming_PBC(int dec_state, int l1_, int l2_ );
void print_binary_of_decimal(int n);
int bit_value(int n, int pos);
Mat_1_int decimal_to_binary(int n);
int countCommonBits(int a,int b) ;
int one_bits_in_bw(int n, int m, int Decimal);

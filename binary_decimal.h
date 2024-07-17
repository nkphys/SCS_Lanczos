/*
*/
#include <complex>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "tensor_type.h"
#include <stdlib.h>

#ifndef binary_decimal
#define binary_decimal

int Act_Translation_2D_alongX_assuming_PBC(int dec_state, int Lx, int Ly, int iy);
int Act_Translation_2D_alongY_assuming_PBC(int dec_state, int Lx, int Ly, int ix);
int Act_Translation_assuming_PBC(int dec_state, int l1_, int l2_ );
int Act_Inv_Translation_assuming_PBC(int dec_state, int l1_, int l2_ );
void print_binary_of_decimal(int n);

Mat_1_int decimal_to_binary(int n);
int countCommonBits(int a,int b) ;
int one_bits_in_bw(int n, int m, int Decimal);

template <typename T>
int bit_value(T n, int pos){
    T i;
    i= (T)1<<pos;
    if((n & i) !=0){
        return 1;
    }
    else{
        return 0;
    }
}

#endif

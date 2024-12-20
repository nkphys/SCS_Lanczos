#include <assert.h>
#include <stdio.h> 
#include <string.h>
#include <bits/stdc++.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;
#include "tensor_type.h"

//bool shouldSwap(char str[], int start, int curr);

//template <typename T>
//void findPermutations(char str[], int index, int n, vector<T> & dec_vec, int base);

//int val(char c);

//char reVal(int num);

//void strev(char *str);

//template <typename T>
//T toDeci(char *str, int base);

//template <typename T>
//char* fromDeci(char res[], int base, T inputNum);

//template <typename T>
//int value_at_pos(T n, int pos, int base);

//template <typename T>
//int Sum_of_Values(T inputNum, int base);

//template <typename T>
//T Updated_decimal_with_value_at_pos(T n, int pos, int base, int new_value);




// Returns true if str[curr] does not matches with any of the
// characters after str[start]

#ifndef Base_to_Decimal
#define Base_to_Decimal


int val(char c);

bool shouldSwap(char str[], int start, int curr);

void from_deci_type2_to_n_array(int dec, int base, Mat_1_int &n_array);

int bit_val_at_site_n_array(int site, int Length, Mat_1_int n_array);

void Convert_n_array_to_bit_array(int Length, Mat_1_int n_array, Mat_1_int &bit_array);

int FromBitArray_toDeci_type2(Mat_1_int bit_array, int base);

int Count_bits_in_bw(int site_i,int  site_j, Mat_1_int bit_array);

int Count_bits_in_bw(int site_i, string type_i,int site_j, string type_j, Mat_1_int bit_array);
// Function to convert a number from given base 'b'
// to decimal
template <typename T> T toDeci(char *str, int base, T)
{
    int len = strlen(str);
    T power = 1; // Initialize power of base
    T num = 0;  // Initialize result
    int i;

    // Decimal equivalent is str[len-1]*1 +
    // str[len-1]*base + str[len-1]*(base^2) + ...
    for (i = len - 1; i >= 0; i--)
    {
        // A digit in input number must be
        // less than number's base
        if (val(str[i]) >= base)
        {
            printf("Invalid Number");
            return -1;
        }

        num += val(str[i]) * power;
        power = power * base;
    }

    return num;
}



template <typename T> T toDeci_type2(char *str, int base, T)
{
    int len = strlen(str);
    int power = 1; // Initialize power of base
    T num = 0;  // Initialize result
    int i;

    // Decimal equivalent is str[len-1]*1 +
    // str[len-1]*base + str[len-1]*(base^2) + ...

   Mat_1_int n_vals;
   n_vals.clear();
   int n_=0;

    for(int i=0;i<len;i++){
        if(val(str[i])==1){
           n_vals.push_back(n_+1);
           n_=0;
        }
        else{
        n_ +=1;
        }
    }

    for (i = 0; i < n_vals.size(); i++)
    {
        // A digit in input number must be
        // less than number's base
        if ( n_vals[i] >= base)
        {
            printf("Invalid Number");
            assert(false);
            return -1;
        }

        num += n_vals[i] * power;
        power = power * base;
    }

    return num;
}


pair_int toDeci_type3(char *str, int base1_, int base2_);

// Prints all distinct permutations in str[0..n-1]
template <typename T>
void findPermutations(char str[], int index, int n, vector<T> & dec_vec, int base)
{
    T dec_temp;

    if (index >= n) {
        dec_temp=toDeci(str, base, T());
        dec_vec.push_back(dec_temp);
        //cout << str << endl;
        return;
    }

    for (int i = index; i < n; i++) {

        // Proceed further for str[i] only if it
        // doesn't match with any of the characters
        // after str[index]
        bool check = shouldSwap(str, index, i);
        if (check) {
            swap(str[index], str[i]);
            findPermutations(str, index + 1, n, dec_vec, base);
            swap(str[index], str[i]);
        }
    }
}

template <typename T>
void findPermutations_type2(char str[], int index, int n, vector<T> & dec_vec, int base)
{
    T dec_temp;

    if (index >= n) {
        dec_temp=toDeci_type2(str, base, T());
        dec_vec.push_back(dec_temp);
        //cout << str << endl;
        return;
    }

    for (int i = index; i < n; i++) {

        // Proceed further for str[i] only if it
        // doesn't match with any of the characters
        // after str[index]
        bool check = shouldSwap(str, index, i);
        if (check) {
            swap(str[index], str[i]);
            findPermutations_type2(str, index+1, n, dec_vec, base);
            swap(str[index], str[i]);
        }
    }
}

void findPermutations_type3(char str[], int index, int n, Mat_1_intpair &dec_vec, int base1_, int base2_);

char reVal(int num);

// Utility function to reverse a string
void strev(char *str);





// Function to convert a given decimal number
// to a base 'base' and
template <typename T>
char* fromDeci(char res[], int base, T inputNum, int Length)
{
    int index = 0;  // Initialize index of result

    if(inputNum==0){
        res[0]='0';
        for(int i=1;i<Length;i++){
        res[i] = '0';
        }
        res[Length] = '\0';
        strev(res);
    }
    else{
        // Convert input number is given base by repeatedly
        // dividing it by base and taking remainder
        while (inputNum > 0)
        {
            res[index++] = reVal(inputNum % base);
            inputNum /= base;
        }
        for(int i=index;i<Length;i++){
            res[i]='0';
        }
        res[Length] = '\0';

        // Reverse the result
        strev(res);
    }

    return res;

}



// Function to convert a given decimal number
// to a base 'base' and
template <typename T>
char* fromDeci(char res[], int base, T inputNum)
{
    int index = 0;  // Initialize index of result

    if(inputNum==0){
        res[0]='0';
        res[1] = '\0';
        strev(res);
    }
    else{
        // Convert input number is given base by repeatedly
        // dividing it by base and taking remainder
        while (inputNum > 0)
        {
            res[index++] = reVal(inputNum % base);
            inputNum /= base;
        }
        res[index] = '\0';

        // Reverse the result
        strev(res);
    }

    return res;

}



template <typename T>
void fromVecint_to_Deci(Mat_1_int &Vec_out, int base, T &outputNum, int sites)
{
    ulli base_ulli = (ulli) base;
    assert(Vec_out.size() == sites);

    ulli power=1;
    outputNum=0;

    for(int i=0;i<sites;i++){
    outputNum += Vec_out[i] * power;
    power = power * base_ulli;
    }

}

template <typename T>
void fromDeci_to_Vecint(Mat_1_int &Vec_out, int base, T inputNum, int sites)
{

    ulli base_ulli= (ulli) base;
    assert( inputNum <= (pow(base,sites) - 1));
    Vec_out.clear();
    Vec_out.resize(sites);
    for(int i=0;i<sites;i++){
        Vec_out[i]=0;
    }

    int index = 0;  // Initialize index of result


    // Convert input number is given base by repeatedly
    // dividing it by base and taking remainder
    while (inputNum > 0)
    {
        Vec_out[index] = inputNum % base_ulli;
        inputNum /= base_ulli;
        index +=1;
    }

    assert(index<=sites);

}

template <typename T>
int Sum_of_Values(T inputNum, int base){

    int sum;

    sum=0;
    if(inputNum==0){
        sum=0;
    }
    else{
        while (inputNum > 0)
        {
            sum += inputNum % base;
            inputNum /= base;
        }
    }

    return sum;
}

template <typename T>
int value_at_pos(T n, int pos, int base){

    assert(base<=10);
    int value;
    char temp;
    char res[100];
    fromDeci<T>(res, base, n);
    int temp_int;

    if(n==0){
        temp_int=0;
        //cout<<"--------------------"<<endl;
        //cout<<"pos ="<<pos<<", strlen(res)-1 = "<<strlen(res)-1<<endl;
        //cout<<"--------------------"<<endl;
    }

    /*
cout<<"-----------------------"<<endl;
cout<<"length of res = "<<strlen(res)<<endl;
for(int j=0;j<strlen(res);j++){
temp = res[j];
val_temp = temp - '0';
cout<<val_temp;
}
cout<<endl;
cout<<"-----------------------"<<endl;
*/


    if(pos>strlen(res)-1){
        value=0;
        //cout<<"pos ="<<pos<<", strlen(res)-1 = "<<strlen(res)-1<<endl;
    }
    else{
        temp = res[strlen(res) - pos -1];
        int ia = temp - '0';
        value=ia;
    }

    return value;
}


template <typename T>
T Updated_decimal_with_value_at_pos(T n, int pos, int base, int new_value){

    assert(base<=10);
    assert(new_value<base);
    int value;
    T new_n;
    char temp;
    char res[100];
    char res_return[100];
    fromDeci<T>(res, base, n);
    strev(res);


    if(pos>strlen(res)-1){
        char res_new[100];
        for(int index=0;index<strlen(res);index++){
            res_new[index]=res[index];
        }

        for(int index=strlen(res);index<pos;index++){
            res_new[index]='0';
        }

        res_new[pos]=(char)(new_value + '0');
        res_new[pos+1]='\0';

        //cout<<"pos ="<<pos<<", strlen(res)-1 = "<<strlen(res)-1<<endl;

        for(int j=0;j<(pos+1);j++){
            res_return[j]=res_new[j];
        }
        res_return[pos+1]='\0';

    }
    else{

        for(int j=0;j<strlen(res);j++){
            res_return[j]=res[j];
        }
        res_return[pos]=(char)(new_value + '0');
        res_return[strlen(res)]='\0';
    }

    strev(res_return);
    new_n=toDeci(res_return, base, T());

    return new_n;
}

#endif

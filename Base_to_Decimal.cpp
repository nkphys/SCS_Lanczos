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


int Count_bits_in_bw(int site_i,int  site_j, Mat_1_int bit_array){

    int bits_in_bw=0;

    if(site_i==site_j){
      bits_in_bw=0;
    }
    else if(site_j>site_i){
        for(int l=site_i+1;l<=site_j-1;l++){
          bits_in_bw += bit_array[l];
        }
    }
    else{
        assert(site_i>site_j);
        for(int l=site_j+1;l<=site_i-1;l++){
          bits_in_bw += bit_array[l];
        }
    }

    return bits_in_bw;
}

void from_deci_type2_to_n_array(int dec, int base, Mat_1_int &n_array){


n_array.clear();
int dec_temp=dec;
int n_;
while(dec_temp>0){
    n_=dec_temp%base;
    dec_temp = dec_temp/base;
    n_array.push_back(n_);
}


}


int FromBitArray_toDeci_type2(Mat_1_int bit_array, int base)
{
    int len = bit_array.size();
    int power = 1; // Initialize power of base
    int num = 0;  // Initialize result
    int i;

    // Decimal equivalent is str[len-1]*1 +
    // str[len-1]*base + str[len-1]*(base^2) + ...

   Mat_1_int n_vals;
   n_vals.clear();
   int n_=0;

    for(int i=0;i<len;i++){
        if((bit_array[i])==1){
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


void Convert_n_array_to_bit_array(int Length, Mat_1_int n_array, Mat_1_int &bit_array){

bit_array.clear();
//bit_array.resize(Length);

//int n_offset1=0;
//int n_offset2=0;
for(int n=0;n<n_array.size();n++){
    //n_offset2=n_array[n]+n_offset1;

    for(int i=0;i<n_array[n]-1;i++){
    bit_array.push_back(0);
    }
    bit_array.push_back(1);

   // n_offset1=n_offset1+n_offset2+1;
}

for(int i=bit_array.size();i<Length;i++){
bit_array.push_back(0);
}

assert(bit_array.size()==Length);


}

int bit_val_at_site_n_array(int site, int Length, Mat_1_int n_array){

assert(n_array.size()<=Length);

int bit_val=0;
int n_offset=-1;

if(n_array.size()>0){
    assert(n_array[0]>0);
}

for(int i=0;i<n_array.size();i++){
n_offset = n_offset + n_array[0];
if(site==n_offset){
 bit_val=1;
 break;
}
}

return bit_val;

}

pair_int toDeci_type3(char *str, int base1_, int base2_)
{
    pair_int dec_pair;
    int len = strlen(str);
    int power; // Initialize power of base
    int num;  // Initialize result
    int i;

    // Decimal equivalent is str[len-1]*1 +
    // str[len-1]*base + str[len-1]*(base^2) + ...

   Mat_1_int n1_vals;
   n1_vals.clear();
   int n_;

    n_=0;
    for(int i=0;i<len;i++){
        if(val(str[i])==1){
           n1_vals.push_back(n_+1);
           n_=0;
        }
        else{
        n_ +=1;
        }
    }

    Mat_1_int n2_vals;
    n2_vals.clear();

     n_=0;
     for(int i=0;i<len;i++){
         if(val(str[i])==2){
            n2_vals.push_back(n_+1);
            n_=0;
         }
         else{
         n_ +=1;
         }
     }



    power = 1; // Initialize power of base
    num = 0;  // Initialize result
    for (i = 0; i < n1_vals.size(); i++)
    {
        // A digit in input number must be
        // less than number's base
        if ( n1_vals[i] >= base1_)
        {
            printf("Invalid Number");
            assert(false);
            //return -1;
        }

        num += n1_vals[i] * power;
        power = power * base1_;
    }
    dec_pair.first = num;


    power = 1; // Initialize power of base
    num = 0;  // Initialize result
    for (i = 0; i < n2_vals.size(); i++)
    {
        // A digit in input number must be
        // less than number's base
        if ( n2_vals[i] >= base2_)
        {
            printf("Invalid Number");
            assert(false);
            //return -1;
        }

        num += n2_vals[i] * power;
        power = power * base2_;
    }
    dec_pair.second = num;

    return dec_pair;
}


void findPermutations_type3(char str[], int index, int n, Mat_1_intpair &dec_vec, int base1_, int base2_){


    pair_int dec_temp;

    if (index >= n) {
        dec_temp=toDeci_type3(str, base1_, base2_);
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
            findPermutations_type3(str, index+1, n, dec_vec, base1_, base2_);
            swap(str[index], str[i]);
        }
    }

}

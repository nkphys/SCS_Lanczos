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


// Function to convert a number from given base 'b' 
// to decimal 
int toDeci(char *str, int base) 
{ 
    int len = strlen(str); 
    int power = 1; // Initialize power of base 
    int num = 0;  // Initialize result 
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

// Function to convert a given decimal number 
// to a base 'base' and 
char* fromDeci(char res[], int base, int inputNum) 
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

int value_at_pos(int n, int pos, int base){

assert(base<=10);
int value;
int val_temp;
char temp;
char res[100];
fromDeci(res, base, n);
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



int Updated_decimal_with_value_at_pos(int n, int pos, int base, int new_value){

assert(base<=10);
assert(new_value<base);
int value;
int new_n;
char temp;
char res[100];
char res_return[100];
fromDeci(res, base, n);
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
new_n=toDeci(res_return, base);

return new_n;
} 

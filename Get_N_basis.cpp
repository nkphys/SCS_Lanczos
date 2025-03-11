#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include "tensor_type.h"

ulli factorial(int I){

ulli val;

if(I==0){
val=1;
}
else if (I>0){
val=I*(factorial(I-1));
}
else {
cout<<"negative integer not allowed"<<endl;
assert(false);
}

return val;
}

int main(){

int Sz=0;
int N=10;


ulli val=0;
for(int i=0;i<=int((N/2)+0.5);i++){
val+=factorial(N)/(factorial(N-(2*i))*factorial(i)*factorial(i));
//cout<<"i : "<<factorial(N)/(factorial(N-2*i)*factorial(i)*factorial(i))<<endl;
cout<<"i : "<<val<<endl;
}

cout<<val<<endl;


}

#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "tensor_type.h"

using namespace std;

void DirectProduct_IndicesSet_by_recursion(Mat_1_int &TempState, Mat_1_int &VecSizes, int & Total_no_elements ,Mat_2_int &IndicesSet){

Mat_1_int Diff;
Total_no_elements +=1;
int size = TempState.size();

/*
for(int i=0;i<TempState.size();i++){
cout<<TempState[i]<<" ";
}
cout<<endl;
*/

IndicesSet.push_back(TempState);

//TempState[size-1] +=1;

Diff.resize(TempState.size());
for(int i=0;i<TempState.size();i++){
Diff[i] = (VecSizes[i]-1)-TempState[i];
}

int Sum_diff=0;
int n_to_grow;
for(int n=(size-1);n>=0;n--){
Sum_diff +=Diff[n];
if(Sum_diff!=0){
n_to_grow=n;
break;
}
}

bool proceed=(Sum_diff!=0);
for(int n=0;n<size;n++){
proceed = proceed && (TempState[n]<VecSizes[n]); 
}

if(proceed){
TempState[n_to_grow] +=1;

for(int n=n_to_grow+1; n<size;n++){
TempState[n] = 0;
}

DirectProduct_IndicesSet_by_recursion(TempState, VecSizes, Total_no_elements, IndicesSet);
}

}


int main (){

cout<<"Direct product by recursion"<<endl;

Mat_1_int A, B;
A.resize(20);B.resize(12);

for(int i=0;i<A.size();i++){
A[i]=i;
}
for(int i=0;i<B.size();i++){
B[i]=i;
}

Mat_1_int C;
//C.resize(A.size()*B.size());

for(int i=0;i<A.size();i++){
for(int j=0;j<B.size();j++){
C.push_back(A[i]*B[j]);
}
}

//cout<<"no. of elements in C = A (X) B : "<<C.size()<<endl;


int No_of_vecs=4;
Mat_2_int IndicesSet;
Mat_1_int tempSet;
Mat_1_int VecSizes;
tempSet.resize(No_of_vecs);
VecSizes.resize(No_of_vecs);

for(int i=0;i<No_of_vecs;i++){
tempSet[i]=0;
}

VecSizes[0]=10;
VecSizes[1]=21;
VecSizes[2]=3;
VecSizes[3]=35;

int Total_no_elements=0;
DirectProduct_IndicesSet_by_recursion(tempSet, VecSizes, Total_no_elements, IndicesSet);

for(int i=0;i<(No_of_vecs-1);i++){
cout<<"("<<VecSizes[i]<<") (X) ";
}
cout<<"("<<VecSizes[No_of_vecs-1]<<") = "<<Total_no_elements<<endl;

cout<<"Set.size() = "<<IndicesSet.size()<<endl;

return 0;
}

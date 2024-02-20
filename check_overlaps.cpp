#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


#include "tensor_type.h"

using namespace std;
int main(){

string file01_str = "overlap_01.txt";
string file10_str = "overlap_10.txt";

ifstream file01_in(file01_str.c_str());
ifstream file10_in(file10_str.c_str());

int size=210;

Mat_2_Complex_doub overlap_01, overlap_10;
overlap_01.resize(size); overlap_10.resize(size);

for(int i=0;i<size;i++){
overlap_01[i].resize(size);
overlap_10[i].resize(size);
}

string temp_str;
int temp_int,i_ind, j_ind;
complex<double> val;


for(int line=0;line<(size*size);line++){
file01_in >>temp_str>>temp_str>>temp_int>>temp_int;
file01_in>>i_ind>>j_ind>>val;
overlap_01[i_ind][j_ind]=val;
//cout<<val<<endl;
}

for(int line=0;line<(size*size);line++){
file10_in >>temp_str>>temp_str>>temp_int>>temp_int;
file10_in>>i_ind>>j_ind>>val;
overlap_10[i_ind][j_ind]=val;
//cout<<val<<endl;
}


string file_01_out_str="overlap_01_ordered.txt";
string file_10_out_str="overlap_10_ordered.txt";

ofstream file_01_out(file_01_out_str.c_str());
ofstream file_10_out(file_10_out_str.c_str());

complex<double> val_01_total, val_10_total;
val_01_total=0.0;
file_01_out<<scientific<<setprecision(20);
for(int i=0;i<size;i++){
for(int j=0;j<size;j++){
val_01_total +=overlap_01[i][j];
file_01_out<<i<<"  "<<j<<"  "<<overlap_01[i][j]<<"  "<<val_01_total.real()<<"  "<<val_01_total.imag()<<endl;
}
}

file_10_out<<scientific<<setprecision(20);
val_10_total=0.0;
for(int j=0;j<size;j++){
for(int i=0;i<size;i++){
val_10_total +=overlap_10[i][j];
file_10_out<<i<<"  "<<j<<"  "<<overlap_10[i][j]<<"  "<<val_10_total.real()<<"  "<<val_10_total.imag()<<endl;
}
}


}

#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){

string file_in_name = "Overlap_with_SingleSiteStates.txt";
ifstream file_in(file_in_name.c_str());

string file_out_name = "Ordered_Overlap_with_SingleSiteStates.txt";
ofstream file_out(file_out_name.c_str());


Mat_2_Complex_doub Overlaps;
Mat_2_doub Overlaps_abs;
int row, col;
complex<double> val;
int size=150;
Overlaps.resize(size);
Overlaps_abs.resize(size);
for (int i=0;i<size;i++){
Overlaps[i].resize(size);Overlaps_abs[i].resize(size);
//Overlaps[i].resize(1);Overlaps_abs[i].resize(1);
}

string line;
if(file_in.is_open())
{
while(!file_in.eof())
{
getline(file_in,line);
stringstream line_stream;
line_stream<<line;
line_stream >> row >> col >> val;
//line_stream >> row >> val;col=0;
Overlaps[row][col] = val;
Overlaps_abs[row][col] = abs(val); 
}
}


double val_max=0.0;
int row_max, col_max;
double percentage_so_far=0.0;
for(int index=0;index<Overlaps.size()*Overlaps[0].size();index++){
val_max=-100;
for(int i=0;i<Overlaps.size();i++){
for(int j=0;j<Overlaps[i].size();j++){

if((Overlaps_abs[i][j])>=val_max){
val_max = abs(Overlaps[i][j]);
row_max = i;
col_max = j;
}

}
}

percentage_so_far +=(Overlaps[row_max][col_max].real()*Overlaps[row_max][col_max].real() +  Overlaps[row_max][col_max].imag()*Overlaps[row_max][col_max].imag())*100;
file_out<<row_max<<"  "<<col_max<<"  "<<Overlaps[row_max][col_max].real()<<"  "<<Overlaps[row_max][col_max].imag()<<"  "<<sqrt(Overlaps[row_max][col_max].real()*Overlaps[row_max][col_max].real() +  Overlaps[row_max][col_max].imag()*Overlaps[row_max][col_max].imag())<<"  "<<percentage_so_far<<endl;
Overlaps_abs[row_max][col_max]=-1000;

}


return 0;
}

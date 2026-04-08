#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include "tensor_type.h"
#include <algorithm>
using namespace std;
#include <vector>
#include <complex>

typedef vector< double >  Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;

typedef vector< int >  Mat_1_int;

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;

int main(){


int N_atoms=2;
int N_orbs=1;

int Length=12;


complex<double> one = complex<double>(1.0,0.0);
complex<double> iota = complex<double>(0.0,1.0);

Mat_2_Complex_doub Sigmax, Sigmay, Iden;

Sigmax.resize(2);Sigmay.resize(2);Iden.resize(2);

for(int i=0;i<2;i++){
Sigmax[i].resize(2);
Sigmay[i].resize(2);
Iden[i].resize(2);
for(int j=0;j<2;j++){
Sigmax[i][j]=0.0;
Sigmay[i][j]=0.0;
Iden[i][j]=0.0;
}
}

Iden[0][0]=1.0;Iden[1][1]=1.0;
Sigmax[0][1]=1.0;Sigmax[1][0]=1.0;
Sigmay[0][1]=complex<double>(0.0,-1.0);Sigmay[1][0]=complex<double>(0.0,1.0);


//e1=(sqrt(3)/2,1/2)
//e2=(-sqrt(3)/2,1/2)
//e3=(0,-1)



//atom + ix*(2) + iy*(2*Lx) + spin*(2*Lx*Ly)

int neigh_x, neigh_y;
int neigh_x_bare, neigh_y_bare;


string Out_file_str = "4pt_Cd_C_Cd_C.txt" ;
ofstream Out_file(Out_file_str.c_str());


Out_file<<"#No.of Oprs c1_comp/real CDag_p1,s1 C_p2,s2 CDag_p3,s3 C_p4,s4  c2_comp/real CDag_p1,s1 C_p2,s2 CDag_p3,s3 C_p4,s4 ..."<<endl;


Mat_1_intpair PairReference;
Mat_1_doub PairReferenceSign;

Mat_1_intpair Pair;
Mat_1_doub PairSign;

pair_int tmp_pair_int;

tmp_pair_int.first=6;  //Cdag
tmp_pair_int.second=3;  //C
PairReference.push_back(tmp_pair_int);
PairReferenceSign.push_back(1.0);


Mat_1_int FirstSite=  {3, 7, 0, 10, 5, 9, 2, 6, 1, 11, 4, 8, 1, 9, 4, 6, 3, 11, 0, 8, 5 , 7,  2 , 10}; //Cdag
Mat_1_int SecondSite= {0, 0, 1, 1 , 2, 2, 3, 3, 4, 4 , 5, 5, 6, 6, 7, 7, 8, 8 , 9, 9, 10, 10, 11, 11}; //C

Mat_1_doub Signs = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};



//Out_file<<"#No.of Oprs c1_comp/real CDag_p1,s1 C_p2,s2 CDag_p3,s3 C_p4,s4  c2_comp/real CDag_p1,s1 C_p2,s2 CDag_p3,s3 C_p4,s4 ..."<<endl;

complex<double> Coeff;

for(int bond=0;bond<FirstSite.size();bond++){
Out_file<<16;

for(int spin=0;spin<2;spin++){
for(int spinp=0;spinp<2;spinp++){
Coeff = one_comp*pow(-1.0,spin+spinp)*Signs[bond];
Out_file<<" "<<-1.0*Coeff<<" "<<FirstSite[bond]<<" "<<spin<<" "<<SecondSite[bond]<<" "<<spin<<" "<<PairReference[0].first<<" "<<spinp<<" "<<PairReference[0].second<<" "<<spinp;

Out_file<<" "<<1.0*Coeff<<" "<<FirstSite[bond]<<" "<<spin<<" "<<SecondSite[bond]<<" "<<spin<<" "<<PairReference[0].second<<" "<<spinp<<" "<<PairReference[0].first<<" "<<spinp;

Out_file<<" "<<1.0*Coeff<<" "<<SecondSite[bond]<<" "<<spin<<" "<<FirstSite[bond]<<" "<<spin<<" "<<PairReference[0].first<<" "<<spinp<<" "<<PairReference[0].second<<" "<<spinp;

Out_file<<" "<<-1.0*Coeff<<" "<<SecondSite[bond]<<" "<<spin<<" "<<FirstSite[bond]<<" "<<spin<<" "<<PairReference[0].second<<" "<<spinp<<" "<<PairReference[0].first<<" "<<spinp;

}
}
Out_file<<endl;
}


return 0;
}

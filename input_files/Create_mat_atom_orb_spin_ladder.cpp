#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
//#include "tensor.h"
#include <algorithm>
using namespace std;
#include <vector>
#include <complex>

typedef vector< double >  Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;

typedef vector< int >  Mat_1_int;

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;

#define PI acos(-1.0)


int main(){


bool PBCX, PBCY;
PBCX=true;
PBCY=true;

bool t1_PBCX, t1_PBCY;
t1_PBCX=true;
t1_PBCY=true;



int Lx=6;
int Ly=1;

/*
1,----1
'     '
'     '
'     '
0-----0
*/


complex<double> const one(1,0);
complex<double> const iota(0,1);

int N_atoms=2;
int N_orbs=1;

double t1_parameter=-0.5;
double t1p_parameter=-0.2;
double t2_parameter=-0.8;

double phi=0.0*PI;


//a + ix*(2) + spin*(2*Lx)
int neigh_x, neigh_y;
int neigh_x_bare, neigh_y_bare;
string Hopping_file_str = "Hopping.txt" ;
ofstream Hopping_file(Hopping_file_str.c_str());

string DenDen_file_str = "DenDen.txt" ;
ofstream DenDen_file(DenDen_file_str.c_str());

int dof1, dof2;

Mat_2_Complex_doub HopMat;
HopMat.resize(4*Lx*Ly);
for(int i=0;i<4*Lx*Ly;i++){
HopMat[i].resize(4*Lx*Ly);
for(int j=0;j<4*Lx*Ly;j++){
HopMat[i][j]=0.0;
}
}



Mat_2_doub DenDenMat;
DenDenMat.resize(2*Lx*Ly);
for(int i=0;i<2*Lx*Ly;i++){
DenDenMat[i].resize(2*Lx*Ly);
for(int j=0;j<2*Lx*Ly;j++){
DenDenMat[i][j]=0.0;
}
}


complex<double> val=0.0;

//t0[2][1] c_{2}^{dag}c_{1}
//bond-1
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==1) && (atom2==0) && (spin1==spin2) ){
val=t1_parameter;

for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){

neigh_x=ix;
neigh_y=iy;

dof1 = atom1 + ix*(2) + spin1*(2*Lx);
dof2 = atom2 + neigh_x*(2) + spin2*(2*Lx);

HopMat[dof2][dof1] = val;
HopMat[dof1][dof2] =conj(HopMat[dof2][dof1]);


}
}

}

}}}

}}}



//t1_plus_a1[2][1] c_{2}^{dag}c_{1}
//site---->neigh
//neigh
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

//site
for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==1) && (spin1==spin2) ){
val=t2_parameter;

for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){
neigh_x_bare=(ix+1);
neigh_y_bare=iy;

if( ((neigh_x_bare<Lx && neigh_x_bare>=0) || PBCX)
     ){
neigh_x = (neigh_x_bare +  Lx)%Lx;
neigh_y = (neigh_y_bare +  Ly)%Ly;

dof1 = atom1 + ix*(2) + spin1*(2*Lx);
dof2 = atom2 + neigh_x*(2) + spin2*(2*Lx);

HopMat[dof2][dof1] = val;
HopMat[dof1][dof2] =conj(HopMat[dof2][dof1]);
}
}
}
}

if((atom1==1) && (atom2==0) && (spin1==spin2) ){
val=-t2_parameter;

for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){
neigh_x_bare=(ix+1);
neigh_y_bare=iy;

if( ((neigh_x_bare<Lx && neigh_x_bare>=0) || PBCX)
     ){
neigh_x = (neigh_x_bare +  Lx)%Lx;
neigh_y = (neigh_y_bare +  Ly)%Ly;

dof1 = atom1 + ix*(2) + spin1*(2*Lx);
dof2 = atom2 + neigh_x*(2) + spin2*(2*Lx);

HopMat[dof2][dof1] = val;
HopMat[dof1][dof2] =conj(HopMat[dof2][dof1]);
}
}
}
}

if((atom1==1) && (atom2==1) && (spin1==spin2) ){
val=t1p_parameter;

for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){
neigh_x_bare=(ix+1);
neigh_y_bare=iy;

if( ((neigh_x_bare<Lx && neigh_x_bare>=0) || PBCX)
     ){
neigh_x = (neigh_x_bare +  Lx)%Lx;
neigh_y = (neigh_y_bare +  Ly)%Ly;

dof1 = atom1 + ix*(2) + spin1*(2*Lx*Ly);
dof2 = atom2 + neigh_x*(2) + spin2*(2*Lx*Ly);

HopMat[dof2][dof1] = val;
HopMat[dof1][dof2] =conj(HopMat[dof2][dof1]);
}

}
}

}

if((atom1==0) && (atom2==0) && (spin1==spin2) ){

val=t1p_parameter;
for(int ix=0;ix<Lx;ix++){
for(int iy=0;iy<Ly;iy++){
neigh_x_bare=(ix+1);
neigh_y_bare=iy;

if( ((neigh_x_bare<Lx && neigh_x_bare>=0) || PBCX)
     ){
neigh_x = (neigh_x_bare +  Lx)%Lx;
neigh_y = (neigh_y_bare +  Ly)%Ly;

dof1 = atom1 + ix*(2) + spin1*(2*Lx*Ly);
dof2 = atom2 + neigh_x*(2) + spin2*(2*Lx*Ly);
val = 
HopMat[dof2][dof1] = val;
HopMat[dof1][dof2] =conj(HopMat[dof2][dof1]);

}
}
}
}



}}}

}}}





for(int i=0;i<4*Lx*Ly;i++){
for(int j=0;j<4*Lx*Ly;j++){
Hopping_file<<HopMat[i][j]<<" ";
}
Hopping_file<<endl;
}



double temp_doub;
for(int i=0;i<2*Lx*Ly;i++){
for(int j=0;j<2*Lx*Ly;j++){
if(j>=i){
temp_doub = DenDenMat[i][j] + DenDenMat[j][i];
}
else{
temp_doub =0.0;
}
DenDen_file<<temp_doub<<" ";
}
DenDen_file<<endl;
}


return 0;
}

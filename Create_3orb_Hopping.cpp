#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){



double thop=1.0;

int OX, OY, OZ;
OX=0;OY=1;OZ=2;

int L1=2;
int L2=2;
Mat_2_doub HopMat;
HopMat.resize(L1*L2*3);
for(int i=0;i<L1*L2*3;i++){
HopMat[i].resize(L1*L2*3);
}

string Hopp_file_str = "Hopp.dat";
ofstream Hopp_file_stream(Hopp_file_str.c_str());

int ix_neigh, iy_neigh, site_neigh;
int site;
for(int ix=0;ix<L1;ix++){
for(int iy=0;iy<L2;iy++){
site=ix + iy*L1;

//+X
ix_neigh=ix+1;
iy_neigh=iy;
site_neigh= ix_neigh + iy_neigh*L1;
if(ix_neigh<L1){
HopMat[site_neigh + OY*L1*L2][site + OY*L1*L2]=thop;
HopMat[site_neigh + OZ*L1*L2][site + OZ*L1*L2]=thop;
HopMat[site + OY*L1*L2][site_neigh + OY*L1*L2]=thop;
HopMat[site + OZ*L1*L2][site_neigh + OZ*L1*L2]=thop;
}


//+Y
ix_neigh=ix;
iy_neigh=iy+1;
site_neigh= ix_neigh + iy_neigh*L1;
if(iy_neigh<L2){
HopMat[site_neigh + OX*L1*L2][site + OX*L1*L2]=thop;
HopMat[site_neigh + OZ*L1*L2][site + OZ*L1*L2]=thop;
HopMat[site + OX*L1*L2][site_neigh + OX*L1*L2]=thop;
HopMat[site + OZ*L1*L2][site_neigh + OZ*L1*L2]=thop;
}

}
}



double eps_=0.000001;

for(int site=0;site<3*L1*L2;site++){
for(int sitep=0;sitep<3*L1*L2;sitep++){
if(sitep>site){
Hopp_file_stream<<HopMat[site][sitep]<<" ";
}
else{
Hopp_file_stream<<0<<" ";
}
}
Hopp_file_stream<<endl;
}


return 0;
}

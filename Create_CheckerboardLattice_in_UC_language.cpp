#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){

double J1_exc=1.0;
double J2_exc=2.0;

int L1=3;
int L2=4;

int B_=0;
int C_=1;

string Hopp_file_str = "J1file.dat";
ofstream Hopp_file_stream(Hopp_file_str.c_str());


string Translation_x_file_str = "Translation_x.txt";
ofstream Translation_x_file_stream(Translation_x_file_str.c_str());


string Translation_y_file_str = "Translation_y.txt";
ofstream Translation_y_file_stream(Translation_y_file_str.c_str());

//site= site1 + site_2*L1;
int site1_new, site2_new, s_new;
int site, site_new;

for(int site2=0;site2<L2;site2++){
for(int site1=0;site1<L1;site1++){
for(int s=0;s<2;s++){ //s=0=B, s=1=C
site = s + (2*site1) + site2*2*L1;

site1_new = (site1 + 1 + L1)%L1;
site2_new = site2;
s_new = s;
site_new = s + (2*site1_new) + site2_new*2*L1;
Translation_x_file_stream<<site<<"  "<<site_new<<endl;

site2_new = (site2 + 1 + L2)%L2; 
site1_new = site1;
s_new = s;
site_new = s + (2*site1_new) + site2_new*2*L1;
Translation_y_file_stream<<site<<"  "<<site_new<<endl;
}
}
}














Mat_2_doub J1_;
J1_.resize(2*L1*L2);
for(int site=0;site<2*L1*L2;site++){
J1_[site].resize(2*L1*L2);
for(int sitep=0;sitep<2*L1*L2;sitep++){
J1_[site][sitep]=0.0;
}
}


//site= site1 + site_2*L1;
//int site1_new, site2_new, s_new;
//int site, site_new;

for(int site2=0;site2<L2;site2++){
for(int site1=0;site1<L1;site1++){
for(int s=0;s<2;s++){ //s=0=B, s=1=C
site = s + (2*site1) + site2*2*L1;


//Intra-unit cell------------------------------------
if(s==B_){
site1_new=site1;
site2_new=site2;
s_new=C_;
site_new = s_new + (2*site1_new) + site2_new*2*L1; 
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;
}
//----------------------------------------------------


//+x---------------------------------------------------
site1_new = (site1 + 1 + L1)%L1;
site2_new = site2;

//B<--->B,C
if(s==B_){
s_new=B_;
site_new = s_new + (2*site1_new) + site2_new*2*L1;
J1_[site_new][site]=J1_exc;
J1_[site][site_new]=J1_exc;

s_new=C_;
site_new = s_new + (2*site1_new) + site2_new*2*L1;
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;
}

//-------------------------------------------------------


//+y---------------------------------------------------
site2_new = (site2 + 1 + L2)%L2;
site1_new = site1;

//C<--->B,C
if(s==C_){
s_new=B_;
site_new = s_new + (2*site1_new) + site2_new*2*L1;
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;

s_new=C_;
site_new = s_new + (2*site1_new) + site2_new*2*L1;
J1_[site_new][site]=J1_exc;
J1_[site][site_new]=J1_exc;
}

//-------------------------------------------------------



// (-x+y)  ---------------------------------------------------
site2_new = (site2 + 1 + L2)%L2;
site1_new = (site1 - 1 + L1)%L1 ;

//C<--->B
if(s==C_){
s_new=B_;
site_new = s_new + (2*site1_new) + site2_new*2*L1;
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;
}

//-------------------------------------------------------




}

}}

double eps_=0.000001;
for(int site=0;site<2*L1*L2;site++){
for(int sitep=0;sitep<2*L1*L2;sitep++){
if(sitep>=site && abs(J1_[site][sitep])>eps_){

Hopp_file_stream<<site <<"   "<<sitep<<"   "<<J1_[site][sitep]<<endl;

}
}
}


return 0;
}

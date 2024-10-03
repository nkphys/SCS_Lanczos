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
double J2_exc=1.5;

int L1=4;
int L2=4;


string Hopp_file_str = "J1file.dat";
ofstream Hopp_file_stream(Hopp_file_str.c_str());


Mat_2_doub J1_;
J1_.resize(L1*L2);
for(int site=0;site<L1*L2;site++){
J1_[site].resize(L1*L2);
for(int sitep=0;sitep<L1*L2;sitep++){
J1_[site][sitep]=0.0;
}
}


//site= site1 + site_2*L1;
int site1_new, site2_new;
int site, site_new;
for(int site2=0;site2<L2;site2++){
for(int site1=0;site1<L1;site1++){
site=site1 + site2*L1;

//+x
site1_new = (site1 + 1 + L1)%L1;
site2_new = site2;
site_new = site1_new + site2_new*L1;
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;



//+y
site1_new = site1;
site2_new = (site2 + 1 + L2)%L2;
site_new = site1_new + site2_new*L1;
J1_[site_new][site]=J2_exc;
J1_[site][site_new]=J2_exc;



if( ((site1 + site2)%2) == 0){
//+x+y
site1_new = (site1 + 1 +L1)%L1;
site2_new = (site2 + 1 + L2)%L2;
site_new = site1_new + site2_new*L1;
J1_[site_new][site]=J1_exc;
J1_[site][site_new]=J1_exc;

}

if( ((site1 + site2)%2) == 1){
//-x+y
site1_new = (site1 - 1 +L1)%L1;
site2_new = (site2 + 1 + L2)%L2;
site_new = site1_new + site2_new*L1;
J1_[site_new][site]=J1_exc;
J1_[site][site_new]=J1_exc;

}

}}

double eps_=0.000001;
for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
if(sitep>=site && abs(J1_[site][sitep])>eps_){

Hopp_file_stream<<site <<"   "<<sitep<<"   "<<J1_[site][sitep]<<endl;

}
}
}


return 0;
}

#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){

double t1_hop=1.0;	
double J1_FM=-0.0;
double U1=200.0;

int L1=4;
int L2=4;


string Hopp_file_str = "Hopping.txt";
ofstream Hopp_file_stream(Hopp_file_str.c_str());

string JFM_file_str = "DirectExchange.txt";
ofstream JFM_file_stream(JFM_file_str.c_str());

string Interactions_file_str = "Interactions.txt";
ofstream Interactions_file_stream(Interactions_file_str.c_str());

string PairHopping_file_str = "PairHopping.txt";
ofstream PairHopping_file_stream(PairHopping_file_str.c_str());

string InteractionAssistedHopping_file_str = "InteractionAssistedHopping.txt";
ofstream InteractionAssistedHopping_file_stream(InteractionAssistedHopping_file_str.c_str());


Mat_2_doub Hopp_, J1_, Int_;
Hopp_.resize(L1*L2); J1_.resize(L1*L2);
Int_.resize(L1*L2);

for(int site=0;site<L1*L2;site++){
Hopp_[site].resize(L1*L2);
J1_[site].resize(L1*L2);
Int_[site].resize(L1*L2);
for(int sitep=0;sitep<L1*L2;sitep++){
Hopp_[site][sitep]=0.0;
J1_[site][sitep]=0.0;
Int_[site][sitep]=0.0;
}
}


//site= site1 + site_2*L1;
int site1_new, site2_new;
int site, site_new;
for(int site2=0;site2<L2;site2++){
for(int site1=0;site1<L1;site1++){
site=site1 + site2*L1;

//neigh +e1
site1_new = (site1 + 1 + L1)%L1;
site2_new = site2;
site_new = site1_new + site2_new*L1;

Hopp_[site_new][site]=t1_hop;
Hopp_[site][site_new]=t1_hop;
J1_[site_new][site]=J1_FM;
J1_[site][site_new]=J1_FM;
Int_[site][site_new]=U1;
Int_[site_new][site]=U1;


//neigh +e2
site1_new = site1;
site2_new = (site2+1+L2)%L2;
site_new = site1_new + site2_new*L1;

Hopp_[site_new][site]=t1_hop;
Hopp_[site][site_new]=t1_hop;
J1_[site_new][site]=J1_FM;
J1_[site][site_new]=J1_FM;
Int_[site][site_new]=U1;
Int_[site_new][site]=U1;

//neigh -e1+e2
site1_new = (site1 -1 +L1)%L1;
site2_new = (site2+1+L2)%L2;
site_new = site1_new + site2_new*L1;

Hopp_[site_new][site]=t1_hop;
Hopp_[site][site_new]=t1_hop;
J1_[site_new][site]=J1_FM;
J1_[site][site_new]=J1_FM;
Int_[site][site_new]=U1;
Int_[site_new][site]=U1;


}}


for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
if(sitep>=site){
Hopp_file_stream<<Hopp_[site][sitep]<<"  ";
}
else{
Hopp_file_stream<<0<<"  ";
}
}
Hopp_file_stream<<endl;
}




double sum_=0;
for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
if(sitep>=site){
JFM_file_stream<<J1_[site][sitep]<<"  ";
}
else{
JFM_file_stream<<0<<"  ";
}
}
JFM_file_stream<<endl;
}


for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
if(sitep>=site){
Interactions_file_stream<<Int_[site][sitep]<<"  ";
}
else{
Interactions_file_stream<<0<<"  ";
}
}
Interactions_file_stream<<endl;
}


for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
PairHopping_file_stream<<0.0<<"  ";
InteractionAssistedHopping_file_stream<<0.0<<"  ";
}
PairHopping_file_stream<<endl;
InteractionAssistedHopping_file_stream<<endl;
}





return 0;
}

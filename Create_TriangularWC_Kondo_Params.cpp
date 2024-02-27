#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){


double t_hop=1.0;	
double LS_Jexc_val=-1.000;
double KondoExchange_val=0.5;
double KondoHopping_val=0.5;

int L1=9;
int L2=9;


string Hopp_file_str = "Hopping.txt";
ofstream Hopp_file_stream(Hopp_file_str.c_str());


string LS_Jexc_file_str = "LocalSpinsExchange.txt";
ofstream LS_Jexc_file_stream(LS_Jexc_file_str.c_str());


string KondoExchange_file_str = "KondoExchange.txt";
ofstream KondoExchange_file_stream(KondoExchange_file_str.c_str());


string KondoHopping_file_str = "KondoHopping.txt";
ofstream KondoHopping_file_stream(KondoHopping_file_str.c_str());



Mat_2_doub Hopp_, LS_Jexc_, JK_exc_;
Hopp_.resize(L1*L2); LS_Jexc_.resize(L1*L2);
JK_exc_.resize(L1*L2);

for(int site=0;site<L1*L2;site++){
Hopp_[site].resize(L1*L2);
LS_Jexc_[site].resize(L1*L2);
JK_exc_[site].resize(L1*L2);
for(int sitep=0;sitep<L1*L2;sitep++){
Hopp_[site][sitep]=0.0;
LS_Jexc_[site][sitep]=0.0;
JK_exc_[site][sitep]=0.0;
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
Hopp_[site_new][site]=t_hop;
Hopp_[site][site_new]=t_hop;
LS_Jexc_[site_new][site]=LS_Jexc_val;
LS_Jexc_[site][site_new]=LS_Jexc_val;


//neigh +e2
site1_new = site1;
site2_new = (site2+1+L2)%L2;
site_new = site1_new + site2_new*L1;
Hopp_[site_new][site]=t_hop;
Hopp_[site][site_new]=t_hop;
LS_Jexc_[site_new][site]=LS_Jexc_val;
LS_Jexc_[site][site_new]=LS_Jexc_val;

//neigh -e1+e2
site1_new = (site1 -1 +L1)%L1;
site2_new = (site2+1+L2)%L2;
site_new = site1_new + site2_new*L1;
Hopp_[site_new][site]=t_hop;
Hopp_[site][site_new]=t_hop;
LS_Jexc_[site_new][site]=LS_Jexc_val;
LS_Jexc_[site][site_new]=LS_Jexc_val;



//Spin coupling Localized spins and fermions 
//site--> Localized spin
//site_new--> fermion

//same sublattice
site_new = site;
JK_exc_[site_new][site]=KondoExchange_val;


//sublattice +e1
site1_new = (site1 + 1 +L1)%L1;
site2_new = site2;
site_new = site1_new + site2_new*L1;
JK_exc_[site_new][site]=KondoExchange_val;


//sublattice +e2
site1_new = site1;
site2_new = (site2 + 1 + L2)%L2;
site_new = site1_new + site2_new*L1;
JK_exc_[site_new][site]=KondoExchange_val;


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
//LS_Jexc_file_stream
for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
if(sitep>=site){
LS_Jexc_file_stream<<LS_Jexc_[site][sitep]<<"  ";
//sum_ +=LS_Jexc_[site][sitep];
}
else{
LS_Jexc_file_stream<<0<<"  ";
}
}
LS_Jexc_file_stream<<endl;
}

//cout<<sum_<<endl;


for(int site=0;site<L1*L2;site++){
for(int sitep=0;sitep<L1*L2;sitep++){
KondoExchange_file_stream<<JK_exc_[site][sitep]<<"  ";
}
KondoExchange_file_stream<<endl;
}






int site_j;
int site_i1, site_i2, site_i;
int site_l1, site_l2, site_l;
for(int site_j1=0;site_j1<L1;site_j1++){
for(int site_j2=0;site_j2<L2;site_j2++){
 site_j = site_j1 + site_j2*(L1);


 //neighbouring-pair-1 (i,l)
 site_i = site_j;
 site_l1 = site_j1;
 site_l2 = (site_j2 + 1 +L2)%L2;  
 site_l = site_l1 + site_l2*L1;

 KondoHopping_file_stream<<site_i<<"  "<<site_l<<"  "<<site_j<<"  "<<KondoHopping_val<<endl; 
 KondoHopping_file_stream<<site_l<<"  "<<site_i<<"  "<<site_j<<"  "<<KondoHopping_val<<endl; 


 //neighbouring-pair-2 (i,l)
 site_i = site_j;
 site_l1 = (site_j1+1+L1)%L1;
 site_l2 = (site_j2);
 site_l = site_l1 + site_l2*L1;

 KondoHopping_file_stream<<site_i<<"  "<<site_l<<"  "<<site_j<<"  "<<KondoHopping_val<<endl;
 KondoHopping_file_stream<<site_l<<"  "<<site_i<<"  "<<site_j<<"  "<<KondoHopping_val<<endl;


 //neighbouring-pair-2 (i,l)
 site_i1 = (site_j1+1+L1)%L1;
 site_i2 = (site_j2);
 site_i = site_i1 + site_i2*L1;

 site_l1 = (site_j1);
 site_l2 = (site_j2+1+L2)%L2;
 site_l = site_l1 + site_l2*L1;

 KondoHopping_file_stream<<site_i<<"  "<<site_l<<"  "<<site_j<<"  "<<KondoHopping_val<<endl;
 KondoHopping_file_stream<<site_l<<"  "<<site_i<<"  "<<site_j<<"  "<<KondoHopping_val<<endl;



}
}






return 0;
}

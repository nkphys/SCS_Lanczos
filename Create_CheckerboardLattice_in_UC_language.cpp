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
double J2_exc=1.2;

int L1=3;
int L2=3;
int L;

int B_=0;
int C_=1;

string Hopp_file_str = "J1file.dat";
ofstream Hopp_file_stream(Hopp_file_str.c_str());


Mat_2_doub Tx; Tx.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Tx[i].resize(2*L1*L2);
}
string Translation_x_file_str = "Translation_x.txt";
ofstream Translation_x_file_stream(Translation_x_file_str.c_str());

Mat_2_doub Ty; Ty.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Ty[i].resize(2*L1*L2);
}
string Translation_y_file_str = "Translation_y.txt";
ofstream Translation_y_file_stream(Translation_y_file_str.c_str());


//Reflection from y axis : x <----> -x
Mat_2_doub Ry; Ry.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Ry[i].resize(2*L1*L2);
}
string Reflection_y_file_str = "Reflection_y.txt";
ofstream Reflection_y_file_stream(Reflection_y_file_str.c_str());

//Reflection from x axis : y <----> -y
Mat_2_doub Rx; Rx.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Rx[i].resize(2*L1*L2);
}
string Reflection_x_file_str = "Reflection_x.txt";
ofstream Reflection_x_file_stream(Reflection_x_file_str.c_str());


//Reflection from x=-y axis : (x,y) <----> (-y,-x)
Mat_2_doub Rpxmy; Rpxmy.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Rpxmy[i].resize(2*L1*L2);
}
string Reflection_pxmy_file_str = "Reflection_pxmy.txt";
ofstream Reflection_pxmy_file_stream(Reflection_pxmy_file_str.c_str());


//Reflection from x=y axis : (x,y) <----> (y,x)
Mat_2_doub Rpxpy; Rpxpy.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Rpxpy[i].resize(2*L1*L2);
}
string Reflection_pxpy_file_str = "Reflection_pxpy.txt";
ofstream Reflection_pxpy_file_stream(Reflection_pxpy_file_str.c_str());


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
Tx[site_new][site]=1.0;


site2_new = (site2 + 1 + L2)%L2; 
site1_new = site1;
s_new = s;
site_new = s + (2*site1_new) + site2_new*2*L1;
Translation_y_file_stream<<site<<"  "<<site_new<<endl;
Ty[site_new][site]=1.0;



if(s==B_){
site1_new = (-site1 + (L1-1))%L1;
site2_new = site2;
}
if(s==C_){
site1_new = (-site1 + L1)%L1;
site2_new = site2;
}
site_new = s + (2*site1_new) + site2_new*2*L1;
Reflection_y_file_stream<<site<<"  "<<site_new<<endl;
Ry[site_new][site]=1.0;


if(s==B_){
site2_new = (-site2 + (L2))%L2;
site1_new = site1;
}
if(s==C_){
site2_new = (-site2 + (L2-1))%L2;
site1_new = site1;
}
site_new = s + (2*site1_new) + site2_new*2*L1;
Reflection_x_file_stream<<site<<"  "<<site_new<<endl;
Rx[site_new][site]=1.0;


//Rpxmy and //Rpxpy
if(L1==L2){
L=L1;

//Rpxmy
if(s==B_){
site2_new = ( (L-2)-site1 + L )%L;
site1_new = ( (L-1)-site2 + L )%L;
s_new = C_;
}
if(s==C_){
site2_new = ( (L-1)-site1 + L )%L;
site1_new = ( (L-2)-site2 + L )%L;
s_new = B_;
}
site_new = s_new + (2*site1_new) + site2_new*2*L1;
Reflection_pxmy_file_stream<<site<<"  "<<site_new<<endl;
Rpxmy[site_new][site]=1.0;


//Rpxpy
if(s==B_){
site2_new = ( site1 + L )%L;
site1_new = ( site2 + L )%L;
s_new = C_;
}
if(s==C_){
site2_new = ( site1 + L )%L;
site1_new = ( site2 + L )%L;
s_new = B_;
}
site_new = s_new + (2*site1_new) + site2_new*2*L1;
Reflection_pxpy_file_stream<<site<<"  "<<site_new<<endl;
Rpxpy[site_new][site]=1.0;


}



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







cout<<"#Checking Lattice Transformation are symmetric or not?-------------"<<endl;

Mat_2_doub Tx_J, Ty_J, Ry_J, Rx_J, Rpxmy_J, Rpxpy_J;
Ty_J.resize(2*L1*L2);Tx_J.resize(2*L1*L2);
Ry_J.resize(2*L1*L2);Rx_J.resize(2*L1*L2);
Rpxmy_J.resize(2*L1*L2);Rpxpy_J.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Tx_J[i].resize(2*L1*L2);Ty_J[i].resize(2*L1*L2);
Ry_J[i].resize(2*L1*L2);Rx_J[i].resize(2*L1*L2);
Rpxmy_J[i].resize(2*L1*L2);Rpxpy_J[i].resize(2*L1*L2);
}


//Tx_J=Tx X J_mat
for(int i=0;i<2*L1*L2;i++){
for(int j=0;j<2*L1*L2;j++){
Tx_J[i][j] =0.0;Ty_J[i][j] =0.0;
Ry_J[i][j] =0.0;Rx_J[i][j] =0.0;
Rpxmy_J[i][j] =0.0;Rpxpy_J[i][j] =0.0;
for(int l=0;l<2*L1*L2;l++){
Tx_J[i][j] += Tx[i][l]*J1_[l][j];
Ty_J[i][j] += Ty[i][l]*J1_[l][j];
Ry_J[i][j] += Ry[i][l]*J1_[l][j];
Rx_J[i][j] += Rx[i][l]*J1_[l][j];
Rpxmy_J[i][j] += Rpxmy[i][l]*J1_[l][j];
Rpxpy_J[i][j] += Rpxpy[i][l]*J1_[l][j];
}
}
}

//Tx_J_TxDag
Mat_2_doub Tx_J_TxDag, Ty_J_TyDag, Ry_J_RyDag, Rx_J_RxDag, Rpxmy_J_RpxmyDag, Rpxpy_J_RpxpyDag;
Ty_J_TyDag.resize(2*L1*L2);Tx_J_TxDag.resize(2*L1*L2);
Ry_J_RyDag.resize(2*L1*L2);Rx_J_RxDag.resize(2*L1*L2);
Rpxmy_J_RpxmyDag.resize(2*L1*L2);Rpxpy_J_RpxpyDag.resize(2*L1*L2);
for(int i=0;i<2*L1*L2;i++){
Tx_J_TxDag[i].resize(2*L1*L2);Ty_J_TyDag[i].resize(2*L1*L2);
Ry_J_RyDag[i].resize(2*L1*L2);Rx_J_RxDag[i].resize(2*L1*L2);
Rpxmy_J_RpxmyDag[i].resize(2*L1*L2);Rpxpy_J_RpxpyDag[i].resize(2*L1*L2);
}
for(int i=0;i<2*L1*L2;i++){
for(int j=0;j<2*L1*L2;j++){
Tx_J_TxDag[i][j] =0.0;Ty_J_TyDag[i][j] =0.0;
Ry_J_RyDag[i][j] =0.0;Rx_J_RxDag[i][j] =0.0;
Rpxmy_J_RpxmyDag[i][j] =0.0;Rpxpy_J_RpxpyDag[i][j] =0.0;
for(int l=0;l<2*L1*L2;l++){
Tx_J_TxDag[i][j] += Tx_J[i][l]*Tx[j][l];
Ty_J_TyDag[i][j] += Ty_J[i][l]*Ty[j][l];
Ry_J_RyDag[i][j] += Ry_J[i][l]*Ry[j][l];
Rx_J_RxDag[i][j] += Rx_J[i][l]*Rx[j][l];
Rpxmy_J_RpxmyDag[i][j] += Rpxmy_J[i][l]*Rpxmy[j][l];
Rpxpy_J_RpxpyDag[i][j] += Rpxpy_J[i][l]*Rpxpy[j][l];
}
}
}

double mat_diff_Tx=0.0;
double mat_diff_Ty=0.0;
double mat_diff_Ry=0.0;
double mat_diff_Rx=0.0;
double mat_diff_Rpxmy=0.0;
double mat_diff_Rpxpy=0.0;
for(int i=0;i<2*L1*L2;i++){
for(int j=0;j<2*L1*L2;j++){
mat_diff_Tx += abs(J1_[i][j] - Tx_J_TxDag[i][j]);
mat_diff_Ty += abs(J1_[i][j] - Ty_J_TyDag[i][j]);
mat_diff_Ry += abs(J1_[i][j] - Ry_J_RyDag[i][j]);
mat_diff_Rx += abs(J1_[i][j] - Rx_J_RxDag[i][j]);
mat_diff_Rpxmy += abs(J1_[i][j] - Rpxmy_J_RpxmyDag[i][j]);
mat_diff_Rpxpy += abs(J1_[i][j] - Rpxpy_J_RpxpyDag[i][j]);
}}


cout<<"|TJTDag - J| For Tx = "<<mat_diff_Tx<<endl;
cout<<"|TJTDag - J| For Ty = "<<mat_diff_Ty<<endl;
cout<<"|TJTDag - J| For Ry = "<<mat_diff_Ry<<endl;
cout<<"|TJTDag - J| For Rx = "<<mat_diff_Rx<<endl;
cout<<"|TJTDag - J| For Rpxmy = "<<mat_diff_Rpxmy<<endl;
cout<<"|TJTDag - J| For Rpxpy = "<<mat_diff_Rpxpy<<endl;


return 0;
}

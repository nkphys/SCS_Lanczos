/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_KondoModel.h"
using namespace std;


void BASIS_KondoModel::Construct_basis(){


    if(Target_Total_Ne>Length){
       cout<<"Target_Total_Ne cannot be larger than No. of sites"<<endl;
        assert(false);
    }

    int TwoTimesTarget_Total_Sz = floor((2.0*Target_Total_Sz)+0.5);

    if(TwoTimesTarget_Total_Sz>(Length+Target_Total_Ne)){
      cout<<"Target_Total_Sz cannot be larger than 0.5*(Length+Target_Total_Ne)"<<endl;
        assert(false);
    }

    if(TwoTimesTarget_Total_Sz<(-1*(Length+Target_Total_Ne))){
      cout<<"Target_Total_Sz cannot be smaller than -0.5*(Length+Target_Total_Ne)"<<endl;
        assert(false);
    }

    Target_Total_SzWithOffset=int(((2.0*Target_Total_Sz + (Length+Target_Total_Ne))*0.5) + 0.5);

    Get_SzSectors_Allowed();

    Dec_LocalizedSpins.resize(SzWithOffsetAllowed.size());
    InverseDec_LocalizedSpins.resize(SzWithOffsetAllowed.size());


    Dec_dnup_Fermions.resize(SzWithOffsetAllowed.size());
    InverseDec_Fermions.resize(SzWithOffsetAllowed.size());

    No_of_SzSets = SzWithOffsetAllowed.size();

    Dmax_dn.resize(No_of_SzSets);
    Dmax_up.resize(No_of_SzSets);

    basis_size=0;
    for(int SzSetNo=0;SzSetNo<SzWithOffsetAllowed.size();SzSetNo++){


    //--------------------------------------------------
    //For Localized Spins--------------------------
    //-------------------------------------------------

    int LocalSpinsSzWithOffset=SzWithOffsetAllowed[SzSetNo].first;
    int BASE_NEW=Length-LocalSpinsSzWithOffset+2;

    //int Dmax_=(BASE_NEW-1)*pow(BASE_NEW,LocalSpinsSzWithOffset);

    char state_str[2*Length];
    for(int i=0;i<Length;i++){
        if(i<LocalSpinsSzWithOffset){
        state_str[i]='1';
        }
        else{
        state_str[i]='0';
        }
    }
    state_str[Length] = '\0';
    int length_state_str = strlen(state_str);
    findPermutations_type2(state_str, 0, length_state_str, Dec_LocalizedSpins[SzSetNo], BASE_NEW);

    int Dmax_temp=Get_max(Dec_LocalizedSpins[SzSetNo]);
    InverseDec_LocalizedSpins[SzSetNo].resize(Dmax_temp+1);
    for(int i=0;i<Dec_LocalizedSpins[SzSetNo].size();i++){
    int d_temp = Dec_LocalizedSpins[SzSetNo][i];
    InverseDec_LocalizedSpins[SzSetNo][d_temp] = i;
    }

    //cout<<"here"<<endl;

    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------


    //-----------------------------------------------------
    //For Fermions-----------------------------------------
    //-----------------------------------------------------

    int FermionsSzWithOffset=SzWithOffsetAllowed[SzSetNo].second;
    assert(SzWithOffsetAllowed[SzSetNo].second<=Target_Total_Ne);

    int N_up = FermionsSzWithOffset;
    int N_dn = Target_Total_Ne-N_up;
    assert(N_dn>=0);

    int BASE_NEW_DN = Length-N_dn+2;
    int BASE_NEW_UP = Length-N_up+2;


    char state_str_dnup[2*Length];
    for(int i=0;i<N_dn;i++){
        state_str_dnup[i]='1';
    }
    for(int i=N_dn;i<(N_up+N_dn);i++){
        state_str_dnup[i]='2';
    }
    for(int i=(N_dn+N_up);i<Length;i++){
        state_str_dnup[i]='0';
    }

    state_str_dnup[Length] = '\0';
    int length_state_str_dnup = strlen(state_str_dnup);
    findPermutations_type3(state_str_dnup, 0, length_state_str_dnup, Dec_dnup_Fermions[SzSetNo], BASE_NEW_DN, BASE_NEW_UP);


    int Dmax_dn_temp=Get_max_first_type(Dec_dnup_Fermions[SzSetNo]);
    int Dmax_up_temp=Get_max_second_type(Dec_dnup_Fermions[SzSetNo]);

    Dmax_dn[SzSetNo]=Dmax_dn_temp;
    Dmax_up[SzSetNo]=Dmax_up_temp;

//    int temp_size = (Dmax_dn_temp+1)*(Dmax_up_temp+1);
//    cout<<temp_size<<endl;
    InverseDec_Fermions[SzSetNo].resize((Dmax_dn_temp+1)*(Dmax_up_temp+1));
    for(int i=0;i<Dec_dnup_Fermions[SzSetNo].size();i++){
    int d_dn_temp = Dec_dnup_Fermions[SzSetNo][i].first;
    int d_up_temp = Dec_dnup_Fermions[SzSetNo][i].second;
    InverseDec_Fermions[SzSetNo][d_dn_temp + (Dmax_dn_temp+1)*d_up_temp] = i;
    }


    //inverse mapping
    //

    basis_size += Dec_LocalizedSpins[SzSetNo].size()*Dec_dnup_Fermions[SzSetNo].size();

}//SzSetNo


    cout<<"Total Basis Size = "<<basis_size<<endl;

}


int BASIS_KondoModel::Concatenate_Spins_and_Fermions(int i_fermion, int i_Lspins, int SzSetNo){

int i_ = i_fermion + i_Lspins*Dec_dnup_Fermions[SzSetNo].size();

return i_;
}


void BASIS_KondoModel::Split_to_Spins_and_Fermions(int &i_fermion, int &i_Lspins, int SzSetNo, int i_local){

i_fermion = i_local%Dec_dnup_Fermions[SzSetNo].size();
i_Lspins = int((1.0*((i_local - i_fermion)/Dec_dnup_Fermions[SzSetNo].size()))+0.5);

}



ulli BASIS_KondoModel::Get_basis_ind(int i_local, int SzSetNo){

ulli m;
m=i_local;

for(int sz_set_no=0;sz_set_no<SzSetNo;sz_set_no++){
m += Dec_LocalizedSpins[sz_set_no].size()*Dec_dnup_Fermions[sz_set_no].size();
}
return m;
}


void BASIS_KondoModel::Get_SzSec_and_localindex(ulli m, int &i_local, int &SzSetNo){

ulli m_offset, m_temp;
m_offset=0;
m_temp=0;
for(int sz_set_no=0;sz_set_no<SzWithOffsetAllowed.size();sz_set_no++){

m_temp += Dec_LocalizedSpins[sz_set_no].size()*Dec_dnup_Fermions[sz_set_no].size();

if(m_temp<=m){
m_offset = m_temp;
SzSetNo = sz_set_no;
}
else{
break;
}
}

i_local = m - m_offset;

}






int BASIS_KondoModel::Get_max(Mat_1_int A_vec){

    int max=0;
    for(int i=0;i<A_vec.size();i++){
        if(A_vec[i]>=max){
            max=A_vec[i];
        }
    }

    return max;
}


int BASIS_KondoModel::Get_max_first_type(Mat_1_intpair A_vec){

    int max=0;
    for(int i=0;i<A_vec.size();i++){
        if(A_vec[i].first>=max){
            max=A_vec[i].first;
        }
    }

    return max;
}


int BASIS_KondoModel::Get_max_second_type(Mat_1_intpair A_vec){

    int max=0;
    for(int i=0;i<A_vec.size();i++){
        if(A_vec[i].second>=max){
            max=A_vec[i].second;
        }
    }

    return max;
}

void BASIS_KondoModel::Get_SzSectors_Allowed(){

    //TwoTimesSzAllowed

    Inverse_SzWithOffsetAllowed.clear();
    SzWithOffsetAllowed.clear();

    int SzFermionsMaxWithOffset = Target_Total_Ne;
    int SzFermionsMinWithOffset = 0;

    int LocalSpinsSzMaxWithOffset = Length;
    int LocalSpinsSzMinWithOffset = 0;


    pair_int Szpair;
    int index=0;

    Inverse_SzWithOffsetAllowed.resize((Length+1)*(Target_Total_Ne+1));
    for(int LocalSpinsSzWithOffset=0;LocalSpinsSzWithOffset<=Length;LocalSpinsSzWithOffset++){
        for(int FermionsSzWithOffset=0;FermionsSzWithOffset<=Target_Total_Ne;FermionsSzWithOffset++){

            if((LocalSpinsSzWithOffset + FermionsSzWithOffset)==Target_Total_SzWithOffset){
                Szpair.first=LocalSpinsSzWithOffset;
                Szpair.second=FermionsSzWithOffset;

                SzWithOffsetAllowed.push_back(Szpair);
                Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset + (Length+1)*(FermionsSzWithOffset)]=index;
                index +=1;
            }
        }
    }



}


void BASIS_KondoModel::clear(){
    //D_up_basis.clear();
    //D_dn_basis.clear();
}

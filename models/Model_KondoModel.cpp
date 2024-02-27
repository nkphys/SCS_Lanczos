/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "Model_KondoModel.h"
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
#define PI 3.14159265

void MODEL_KondoModel::Act_Hamil(BASIS_KondoModel &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){}

void MODEL_KondoModel::Add_diagonal_terms(BASIS_KondoModel &basis){


    EPS_=0.000000001;
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;


    //Remember H[l][m]=<l|H|m>
    ulli m;
    double value;
    int dec_i_Ls;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int bit_i, bit_j;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    //int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    //Mat_1_int bit_dn_array, bit_up_array;
    //Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset;
    int BASE_NEW_DN, BASE_NEW_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------


    double Sz_Localized, sz_fermion;


    //Local Spins SzSz interaction
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

    value=0;
    dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

    LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
    base_Ls=Length-LocalSpinsSzWithOffset+2;
    from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    for(int site_i=0;site_i<Length;site_i++){
    for(int site_j=0;site_j<Length;site_j++){

    if(site_j>site_i){

    bit_i = bit_val_at_site_n_array(site_i, Length, n_array);
    bit_j = bit_val_at_site_n_array(site_j, Length, n_array);

    //0-->dn, 1-->up
    value += J_LSpins_mat[site_i][site_j]*((1.0*bit_i)-0.5)*((1.0*bit_j)-0.5);

    }
    }}

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){
        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        if(value!=0){
            Hamil.value.push_back(value*one);
            Hamil.rows.push_back(m);
            Hamil.columns.push_back(m);
        }
    }


    }
    }


    //Kondo coupling SzXsz
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){


        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_NEW_DN = Length-N_dn+2;
        BASE_NEW_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_NEW_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_NEW_UP, n_up_array);

        //Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        //Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        value=0.0;
        for(int site_i=0;site_i<Length;site_i++){
            for(int site_j=0;site_j<Length;site_j++){
           //sz_i Sz_j

           if(abs(KondoExchange_mat[site_i][site_j])>EPS_){
            bit_j = bit_val_at_site_n_array(site_j, Length, n_array);
            Sz_Localized = ((1.0*bit_j)-0.5);

            bit_up_fermion = bit_val_at_site_n_array(site_i, Length, n_up_array);
            bit_dn_fermion = bit_val_at_site_n_array(site_i, Length, n_dn_array);
            sz_fermion = (0.5*(bit_up_fermion - bit_dn_fermion));

            value += KondoExchange_mat[site_i][site_j]*Sz_Localized*sz_fermion;
                }
        }}
            if(value!=0){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m);
                Hamil.columns.push_back(m);
            }




    }
    }
    }



}


void MODEL_KondoModel::Add_LocalSpin_couplings(BASIS_KondoModel &basis){


    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    int bits_in_bw_up2, bits_in_bw_dn2;
    double FM_SIGN;


    //Kondo coupling Sminus_i X Splus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);


        //Kondo
        for(int site_i=0;site_i<Length;site_i++){
            for(int site_j=site_i+1;site_j<Length;site_j++){

           if(abs(J_LSpins_mat[site_i][site_j])>EPS_){


            //Sminus_i Splus_j
            if( bit_array[site_j]==0  &&
                bit_array[site_i]==1
                    ){

                 //FermionsSzWithOffset_new = FermionsSzWithOffset -1;
                 //LocalSpinsSzWithOffset_new=LocalSpinsSzWithOffset+1;

                 //SzSetNo_new = basis.Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset_new + (Length+1)*(FermionsSzWithOffset_new)];

                 bit_array_new=bit_array;

                 bit_array_new[site_j]=1;
                 bit_array_new[site_i]=0;

                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls);

                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins_new, SzSetNo);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

                 value = 0.5*J_LSpins_mat[site_i][site_j];

                 if(value!=0){
                     assert(m_new>=m);
                     Hamil.value.push_back(value*one);
                     Hamil.rows.push_back(m);
                     Hamil.columns.push_back(m_new);
                 }

            }




            }
        }
    }

    }
    }
    }

}

void MODEL_KondoModel::Add_FermionHopping(BASIS_KondoModel &basis){

    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    ulli m, m_new;
    int i_fermion_new;
    double value;
    Mat_1_int n_dn_array, n_up_array;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    int LocalSpinsSzWithOffset, base_Ls;
    int bit_i, bit_j;
    int i_local, i_local_new;
    int FermionsSzWithOffset;
    int BASE_NEW_DN, BASE_NEW_UP;
    int bits_in_bw;
    double FM_SIGN;

    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    //for(int SzSetNo=0;SzSetNo<=0;SzSetNo++){
    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

    d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
    d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;


    LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
    FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
    N_up = FermionsSzWithOffset;
    N_dn = Target_Total_Ne-N_up;

    BASE_NEW_DN = Length-N_dn+2;
    BASE_NEW_UP = Length-N_up+2;

    from_deci_type2_to_n_array(d_dn_, BASE_NEW_DN, n_dn_array);
    from_deci_type2_to_n_array(d_up_, BASE_NEW_UP, n_up_array);

    Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
    Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);

    //dn hopping
    for(int site_i=0;site_i<Length;site_i++){
        for(int site_j=0;site_j<Length;site_j++){
            if(abs(Hopping_mat[site_i][site_j])>0.000000001){
            //c_{i,dn}^{dagger}c_{j,dn}

            if((bit_dn_array[site_j]==1 && bit_dn_array[site_i]==0) &&
                bit_up_array[site_i]==0){

             bit_dn_array_new = bit_dn_array;
             bit_dn_array_new[site_j]=0;
             bit_dn_array_new[site_i]=1;

            d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_NEW_DN);

            //[SzSetNo][d_dn_temp + (Dmax_dn_temp+1)*d_up_temp]
            i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_new_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_];

            for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){
            i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
            i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);

            m = basis.Get_basis_ind(i_local, SzSetNo);
            m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

            bits_in_bw = Count_bits_in_bw(site_i, site_j, bit_dn_array);


            FM_SIGN = pow(-1.0,(1.0*bits_in_bw));

            value = -1.0*FM_SIGN*Hopping_mat[site_i][site_j];
            if(value!=0){
                assert(m_new<=m);
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m_new);
                Hamil.columns.push_back(m);
            }
            }

            }}
        }}

    //up hopping
    for(int site_i=0;site_i<Length;site_i++){
        for(int site_j=0;site_j<Length;site_j++){
            if(abs(Hopping_mat[site_i][site_j])>0.000000001){
            //c_{i,up}^{dagger}c_{j,up}

            if((bit_up_array[site_j]==1 && bit_up_array[site_i]==0) &&
                bit_dn_array[site_i]==0){

             bit_up_array_new = bit_up_array;
             bit_up_array_new[site_j]=0;
             bit_up_array_new[site_i]=1;

            d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_NEW_UP);

            //[SzSetNo][d_dn_temp + (Dmax_dn_temp+1)*d_up_temp]
            i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_new_];

            for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

            i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
            i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);

            m = basis.Get_basis_ind(i_local, SzSetNo);
            m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

            bits_in_bw = Count_bits_in_bw(site_i, site_j, bit_up_array);


            FM_SIGN = pow(-1.0,(1.0*bits_in_bw));

            value = -1.0*FM_SIGN*Hopping_mat[site_i][site_j];
            if(value!=0){
                assert(m_new<=m);
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m_new);
                Hamil.columns.push_back(m);
            }
            }

            }
            }
        }}


    }
    }


}

void MODEL_KondoModel::Add_Kondocouplings(BASIS_KondoModel &basis){


    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    int bits_in_bw_up2, bits_in_bw_dn2;
    double FM_SIGN;


    //Kondo coupling Splus X sminus
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo


        for(int site_i=0;site_i<Length;site_i++){
            for(int site_j=0;site_j<Length;site_j++){

           if(abs(KondoExchange_mat[site_i][site_j])>EPS_){


            //sminus_i Splus_j
            if( bit_array[site_j]==0  &&
                bit_up_array[site_i]==1
                    ){
                assert(bit_dn_array[site_i]==0);

                 FermionsSzWithOffset_new = FermionsSzWithOffset -1;
                 LocalSpinsSzWithOffset_new=LocalSpinsSzWithOffset+1;

                 SzSetNo_new = basis.Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset_new + (Length+1)*(FermionsSzWithOffset_new)];

                 bit_array_new=bit_array;
                 bit_dn_array_new=bit_dn_array;
                 bit_up_array_new=bit_up_array;

                 bit_array_new[site_j]=1;
                 bit_dn_array_new[site_i]=1;
                 bit_up_array_new[site_i]=0;

                 BASE_NEW_DN = Length-(N_dn+1)+2;
                 BASE_NEW_UP = Length-(N_up-1)+2;

                 d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_NEW_UP);
                 d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_NEW_DN);

                 i_fermion_new = basis.InverseDec_Fermions[SzSetNo_new][d_dn_new_ + (basis.Dmax_dn[SzSetNo_new]+1)*d_up_new_];

                 base_Ls_new=Length-LocalSpinsSzWithOffset_new+2;
                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls_new);

                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo_new][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins_new, SzSetNo_new);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo_new);


                 //basis: up ...down
                 //bits_in_bw_up = Count_bits_in_bw(site_i, "closed", Length, "open", bit_up_array) + 1;
                 //bits_in_bw_dn = Count_bits_in_bw(0, "closed", site_i, "open",bit_dn_array);

                 //basis: down ...up
                 bits_in_bw_dn = Count_bits_in_bw(site_i, "closed", Length, "open", bit_dn_array);
                 bits_in_bw_up = Count_bits_in_bw(0, "closed", site_i, "open", bit_up_array);

                 FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                 value = 0.5*KondoExchange_mat[site_i][site_j]*FM_SIGN;

                 if(value!=0){
                     assert(m_new>=m);
                     Hamil.value.push_back(value*one);
                     Hamil.rows.push_back(m);
                     Hamil.columns.push_back(m_new);
                 }

            }




            }
        }
    }

    }
    }
    }


    //Kondo hopping  +0.5 X c(up,i)_dag X c(up,l) X Sz(j)
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);

        for(int index=0;index<KondoHoppings.size();index++){
            int site_i=KondoHoppings_sites[index].first_ ;
            int site_l=KondoHoppings_sites[index].second_;
            int site_j=KondoHoppings_sites[index].third_;

            bit_j = bit_val_at_site_n_array(site_j, Length, n_array);
            Sz_Localized = ((1.0*bit_j)-0.5);

            if(site_l>=site_i){
            if((bit_up_array[site_l]==1 && bit_up_array[site_i]==0) &&
                bit_dn_array[site_i]==0){

             bit_up_array_new = bit_up_array;
             bit_up_array_new[site_l]=0;
             bit_up_array_new[site_i]=1;

            d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_UP);
            i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_new_];
            i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);
            m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

            bits_in_bw = Count_bits_in_bw(site_i, site_l, bit_up_array);
            FM_SIGN = pow(-1.0,(1.0*bits_in_bw));

            value = 0.5*FM_SIGN*KondoHoppings[index]*Sz_Localized;
            if(value!=0){
                assert(m_new<=m);
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m_new);
                Hamil.columns.push_back(m);
            }
            }
            }
        }
    }
    }
    }



    //Kondo hopping  -0.5 X c(dn,i)_dag X c(dn,l) X Sz(j)
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);
        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);

        for(int index=0;index<KondoHoppings.size();index++){
            int site_i=KondoHoppings_sites[index].first_ ;
            int site_l=KondoHoppings_sites[index].second_;
            int site_j=KondoHoppings_sites[index].third_;

            bit_j = bit_val_at_site_n_array(site_j, Length, n_array);
            Sz_Localized = ((1.0*bit_j)-0.5);

            if(site_l>=site_i){
            if((bit_dn_array[site_l]==1 && bit_dn_array[site_i]==0) &&
                bit_up_array[site_i]==0){

             bit_dn_array_new = bit_dn_array;
             bit_dn_array_new[site_l]=0;
             bit_dn_array_new[site_i]=1;

            d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_DN);
            i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_new_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_];
            i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);
            m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

            bits_in_bw = Count_bits_in_bw(site_i, site_l, bit_dn_array);

            FM_SIGN = pow(-1.0,(1.0*bits_in_bw));

            value = -0.5*FM_SIGN*KondoHoppings[index]*Sz_Localized;
            if(value!=0){
                assert(m_new<=m);
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m_new);
                Hamil.columns.push_back(m);
            }
            }
        }
        }

    }
    }
    }




    //Kondo hopping  c(dn,i)_dag X c(up,l) X Splus(j)
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);



        for(int index=0;index<KondoHoppings.size();index++){
            int site_i=KondoHoppings_sites[index].first_ ;
            int site_l=KondoHoppings_sites[index].second_;
            int site_j=KondoHoppings_sites[index].third_;


            if( (bit_up_array[site_l]==1 && (bit_dn_array[site_i]==0 && bit_up_array[site_i]==0)) &&
                bit_array[site_j]==0){


                FermionsSzWithOffset_new = FermionsSzWithOffset -1;
                LocalSpinsSzWithOffset_new=LocalSpinsSzWithOffset+1;

                SzSetNo_new = basis.Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset_new + (Length+1)*(FermionsSzWithOffset_new)];

                bit_array_new=bit_array;
                bit_dn_array_new=bit_dn_array;
                bit_up_array_new=bit_up_array;

                bit_array_new[site_j]=1;
                bit_dn_array_new[site_i]=1;
                bit_up_array_new[site_l]=0;

                BASE_NEW_DN = Length-(N_dn+1)+2;
                BASE_NEW_UP = Length-(N_up-1)+2;

                d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_NEW_UP);
                d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_NEW_DN);

                i_fermion_new = basis.InverseDec_Fermions[SzSetNo_new][d_dn_new_ + (basis.Dmax_dn[SzSetNo_new]+1)*d_up_new_];

                base_Ls=Length-LocalSpinsSzWithOffset_new+2;
                dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls);

                i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo_new][dec_i_Ls_new];

                i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins_new, SzSetNo_new);
                m_new = basis.Get_basis_ind(i_local_new, SzSetNo_new);

                bits_in_bw_up = Count_bits_in_bw(0, "closed", site_l, "open", bit_up_array);
                bits_in_bw_dn = Count_bits_in_bw(site_i, "closed", Length, "open", bit_dn_array);

                FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                value = 0.5*KondoHoppings[index]*FM_SIGN;

                if(value!=0){
                    assert(m_new>=m);
                    Hamil.value.push_back(value*one);
                    Hamil.rows.push_back(m);
                    Hamil.columns.push_back(m_new);
                }

            }

        }

    }
    }
    }






}


void MODEL_KondoModel::Initializer_opr_LocalizedSpin_SzSz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){

    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;

    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();



    //Remember H[l][m]=<l|H|m>
    ulli m;
    double value;
    int dec_i_Ls;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int bit_i, bit_j;
    int i_local;


    //Sz_i Sz_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

    value=0;
    dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

    LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
    base_Ls=Length-LocalSpinsSzWithOffset+2;
    from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    bit_i = bit_val_at_site_n_array(site_i, Length, n_array);
    bit_j = bit_val_at_site_n_array(site_j, Length, n_array);

    //0-->dn, 1-->up
    value = ((1.0*bit_i)-0.5)*((1.0*bit_j)-0.5);


    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){
        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        if(value!=0){
            OPR.value.push_back(value*one);
            OPR.rows.push_back(m);
            OPR.columns.push_back(m);
        }
    }

    }
    }


}



void MODEL_KondoModel::Initializer_opr_LocalizedSpin_SpSm_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){

    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;

    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();



    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    Mat_1_int bit_array, bit_array_new;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    Mat_1_int bit_dn_array, bit_up_array;
    int FermionsSzWithOffset;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int i_local_new, i_Lspins_new;




    //Kondo coupling Splus_i X Sminus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);
        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

            //Splus_i X Sminus_j
            if( (bit_array[site_j]==1  &&
                bit_array[site_i]==0)
                    ||
                    ((site_j==site_i)  && bit_array[site_j]==1 )
                    ){

                 bit_array_new=bit_array;

                 if(site_j!=site_i){
                 bit_array_new[site_j]=0;
                 bit_array_new[site_i]=1;
                    }
                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls);
                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins_new, SzSetNo);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

                 value = 1.0;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }
    }
    }
    }


}




void MODEL_KondoModel::Initializer_opr_LocalizedSpin_SmSp_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){

    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;

    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();



    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    Mat_1_int bit_array, bit_array_new;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    Mat_1_int bit_dn_array, bit_up_array;
    int FermionsSzWithOffset;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int i_local_new, i_Lspins_new;




    //Kondo coupling Sminus_i X Splus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);
        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

            //Sminus_i X Splus_j
            if( (bit_array[site_i]==1  &&
                bit_array[site_j]==0)
                    ||
                    ((site_j==site_i)  && bit_array[site_j]==0 )
                    ){

                 bit_array_new=bit_array;

                 if(site_j!=site_i){
                 bit_array_new[site_j]=1;
                 bit_array_new[site_i]=0;
                    }
                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls);
                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins_new, SzSetNo);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

                 value = 1.0;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }
    }
    }
    }


}



void MODEL_KondoModel::Initializer_opr_Fermions_LocalOprs(string opr_type, BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i){


    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();

    //Remember H[l][m]=<l|H|m>
    ulli m;
    double value;
    int dec_i_Ls;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    //int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    //Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------


    double n_site_i, sz_fermion_site_j;


   if(opr_type=="n_local"){
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){

       // dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        //LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
       // base_Ls=Length-LocalSpinsSzWithOffset+2;
        //from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


           //n(site_i)

            bit_up_fermion = bit_val_at_site_n_array(site_i, Length, n_up_array);
            bit_dn_fermion = bit_val_at_site_n_array(site_i, Length, n_dn_array);
            n_site_i = (1.0*(bit_up_fermion + bit_dn_fermion));

            value = n_site_i;

            if(abs(value)>EPS_){
             for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){
                 i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
                 m = basis.Get_basis_ind(i_local, SzSetNo);

                OPR.value.push_back(value*one);
                OPR.rows.push_back(m);
                OPR.columns.push_back(m);
            }
    }
    }
    }
        }


}



void MODEL_KondoModel::Initializer_opr_Fermions_szsz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){


    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();

    //Remember H[l][m]=<l|H|m>
    ulli m;
    double value;
    int dec_i_Ls;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    //int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    //Mat_1_int bit_dn_array, bit_up_array;
    //Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------


    double sz_fermion_site_i, sz_fermion_site_j;


    //Kondo coupling szXsz
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){

       // dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        //LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
       // base_Ls=Length-LocalSpinsSzWithOffset+2;
        //from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        //Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        //Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);




           //sz_i sz_j

            bit_up_fermion = bit_val_at_site_n_array(site_i, Length, n_up_array);
            bit_dn_fermion = bit_val_at_site_n_array(site_i, Length, n_dn_array);
            sz_fermion_site_i = (0.5*(bit_up_fermion - bit_dn_fermion));

            bit_up_fermion = bit_val_at_site_n_array(site_j, Length, n_up_array);
            bit_dn_fermion = bit_val_at_site_n_array(site_j, Length, n_dn_array);
            sz_fermion_site_j = (0.5*(bit_up_fermion - bit_dn_fermion));

            value = sz_fermion_site_i*sz_fermion_site_j;

            if(value!=0){
             for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){
                 i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
                 m = basis.Get_basis_ind(i_local, SzSetNo);

                OPR.value.push_back(value*one);
                OPR.rows.push_back(m);
                OPR.columns.push_back(m);
            }
    }
    }
    }



}


void MODEL_KondoModel::Initializer_opr_Fermions_spsm_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){

    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();

    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    double FM_SIGN;


    //Kondo coupling splus_i X sminus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

        //splus_i X sminus_j
            if( (bit_up_array[site_j]==1  &&
                bit_dn_array[site_i]==1)
                    ||

                ((site_i==site_j) && bit_up_array[site_j]==1)
                    ){

                 bit_dn_array_new=bit_dn_array;
                 bit_up_array_new=bit_up_array;

                 if(site_i!=site_j){
                 assert(bit_dn_array[site_j]==0);
                 assert(bit_up_array[site_i]==0);

                 bit_dn_array_new[site_i]=0;
                 bit_up_array_new[site_i]=1;

                 bit_up_array_new[site_j]=0;
                 bit_dn_array_new[site_j]=1;
                }
                 else{
                     assert(bit_dn_array[site_j]==0);

                     //bit_up_array_new[site_i]=0;
                    // bit_dn_array_new[site_i]=1;
                 }


                 d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_UP);
                 d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_DN);

                 i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_new_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_new_];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

                 bits_in_bw_up = Count_bits_in_bw(0, "closed", site_i, "open", bit_up_array) + Count_bits_in_bw(0, "closed", site_j, "open", bit_up_array);
                 bits_in_bw_dn = Count_bits_in_bw(0, "closed", site_i, "open", bit_dn_array) + Count_bits_in_bw(0, "closed", site_j, "open", bit_dn_array);

                 FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                 value = FM_SIGN;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }

    }
    }
    }




}



void MODEL_KondoModel::Initializer_opr_Fermions_smsp_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){

    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();

    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    double FM_SIGN;


    //Kondo coupling splus_i X sminus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

        //sminus_i X splus_j
            if( (bit_dn_array[site_j]==1  &&
                bit_up_array[site_i]==1)
                    ||

                ((site_i==site_j) && bit_dn_array[site_j]==1)
                    ){

                 bit_dn_array_new=bit_dn_array;
                 bit_up_array_new=bit_up_array;

                 if(site_i!=site_j){
                 assert(bit_up_array[site_j]==0);
                 assert(bit_dn_array[site_i]==0);

                 bit_dn_array_new[site_i]=1;
                 bit_up_array_new[site_i]=0;

                 bit_up_array_new[site_j]=1;
                 bit_dn_array_new[site_j]=0;
                }
                 else{
                     assert(bit_up_array[site_j]==0);

                     //bit_up_array_new[site_i]=0;
                    // bit_dn_array_new[site_i]=1;
                 }


                 d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_UP);
                 d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_DN);

                 i_fermion_new = basis.InverseDec_Fermions[SzSetNo][d_dn_new_ + (basis.Dmax_dn[SzSetNo]+1)*d_up_new_];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins, SzSetNo);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo);

                 //bits_in_bw_up = Count_bits_in_bw(0, site_i, bit_up_array) + Count_bits_in_bw(0, site_j, bit_up_array);
                // bits_in_bw_dn = Count_bits_in_bw(0, site_i, bit_dn_array) + Count_bits_in_bw(0, site_j, bit_dn_array);

                 bits_in_bw_up = Count_bits_in_bw(0, "closed", site_i, "open", bit_up_array) + Count_bits_in_bw(0, "closed", site_j, "open", bit_up_array);
                 bits_in_bw_dn = Count_bits_in_bw(0, "closed", site_i, "open", bit_dn_array) + Count_bits_in_bw(0, "closed", site_j, "open", bit_dn_array);

                 FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                 value = FM_SIGN;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }

    }
    }
    }




}






void MODEL_KondoModel::Initializer_opr_Fermionsz_LocalizedSz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){


    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();

    ulli m;
    double value;
    int dec_i_Ls;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int bit_i, bit_j;
    int i_local;


    //-----------------
    int d_dn_, d_up_;
    //int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    //Mat_1_int bit_dn_array, bit_up_array;
    //Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset;
    int BASE_NEW_DN, BASE_NEW_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------


    double Sz_Localized, sz_fermion;


    //Kondo coupling sz_iXSz_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){


        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_NEW_DN = Length-N_dn+2;
        BASE_NEW_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_NEW_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_NEW_UP, n_up_array);

        //Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        //Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        value=0.0;

           //sz_i Sz_j
            bit_j = bit_val_at_site_n_array(site_j, Length, n_array);
            Sz_Localized = ((1.0*bit_j)-0.5);

            bit_up_fermion = bit_val_at_site_n_array(site_i, Length, n_up_array);
            bit_dn_fermion = bit_val_at_site_n_array(site_i, Length, n_dn_array);
            sz_fermion = (0.5*(bit_up_fermion - bit_dn_fermion));

            value = Sz_Localized*sz_fermion;

            if(value!=0){
                OPR.value.push_back(value*one);
                OPR.rows.push_back(m);
                OPR.columns.push_back(m);
            }

    }
    }
    }


}



void MODEL_KondoModel::Initializer_opr_Fermionsminus_LocalizedSplus_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){


    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();


    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    double FM_SIGN;


    //Kondo coupling sminus_i Splus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

            //sminus_i Splus_j
            if( bit_array[site_j]==0  &&
                bit_up_array[site_i]==1
                    ){
                assert(bit_dn_array[site_i]==0);

                 FermionsSzWithOffset_new = FermionsSzWithOffset -1;
                 LocalSpinsSzWithOffset_new=LocalSpinsSzWithOffset+1;

                 SzSetNo_new = basis.Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset_new + (Length+1)*(FermionsSzWithOffset_new)];

                 bit_array_new=bit_array;
                 bit_dn_array_new=bit_dn_array;
                 bit_up_array_new=bit_up_array;

                 bit_array_new[site_j]=1;
                 bit_dn_array_new[site_i]=1;
                 bit_up_array_new[site_i]=0;

                 BASE_NEW_DN = Length-(N_dn+1)+2;
                 BASE_NEW_UP = Length-(N_up-1)+2;

                 d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_NEW_UP);
                 d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_NEW_DN);

                 i_fermion_new = basis.InverseDec_Fermions[SzSetNo_new][d_dn_new_ + (basis.Dmax_dn[SzSetNo_new]+1)*d_up_new_];

                 base_Ls_new=Length-LocalSpinsSzWithOffset_new+2;
                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls_new);

                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo_new][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins_new, SzSetNo_new);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo_new);


                 //basis: up ...down
                 //bits_in_bw_up = Count_bits_in_bw(site_i, "closed", Length, "open", bit_up_array) + 1;
                 //bits_in_bw_dn = Count_bits_in_bw(0, "closed", site_i, "open",bit_dn_array);

                 //basis: down ...up
                 bits_in_bw_dn = Count_bits_in_bw(site_i, "closed", Length, "open", bit_dn_array);
                 bits_in_bw_up = Count_bits_in_bw(0, "closed", site_i, "open", bit_up_array);

                 FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                 value = FM_SIGN;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }

    }
    }
    }


}




void MODEL_KondoModel::Initializer_opr_Fermionsplus_LocalizedSminus_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j){


    OPR.nrows = basis.basis_size;
    OPR.ncols = OPR.nrows;
    OPR.columns.clear();
    OPR.rows.clear();
    OPR.value.clear();


    //Remember H[l][m]=<l|H|m>
    ulli m, m_new;
    double value;
    int dec_i_Ls, dec_i_Ls_new;
    Mat_1_int n_array;
    int LocalSpinsSzWithOffset, base_Ls;
    int LocalSpinsSzWithOffset_new, base_Ls_new;
    int bit_i, bit_j;
    double Sz_Localized;
    Mat_1_int bit_array, bit_array_new;
    int i_local;

    int bits_in_bw;

    //-----------------
    int d_dn_, d_up_;
    int d_dn_new_, d_up_new_;
    int N_up, N_dn;
    //ulli m, m_new;
    //int i_fermion_new;
    Mat_1_int n_dn_array, n_up_array;
    int bit_up_fermion, bit_dn_fermion;
    Mat_1_int bit_dn_array, bit_up_array;
    Mat_1_int bit_dn_array_new, bit_up_array_new;
    //int i_local, i_local_new;
    int FermionsSzWithOffset, FermionsSzWithOffset_new;
    int BASE_NEW_DN, BASE_NEW_UP;
    int BASE_DN, BASE_UP;
    //int bits_in_bw;
    //double FM_SIGN;
    //-------------------

    int SzSetNo_new;
    int i_fermion_new;
    int i_local_new, i_Lspins_new;

    int bits_in_bw_up, bits_in_bw_dn;
    double FM_SIGN;


    //Kondo coupling splus_i Sminus_j
    for(int SzSetNo=0;SzSetNo<basis.No_of_SzSets;SzSetNo++){
    for(int i_Lspins=0;i_Lspins<basis.Dec_LocalizedSpins[SzSetNo].size();i_Lspins++){

        dec_i_Ls=basis.Dec_LocalizedSpins[SzSetNo][i_Lspins];

        LocalSpinsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].first;
        base_Ls=Length-LocalSpinsSzWithOffset+2;
        from_deci_type2_to_n_array(dec_i_Ls, base_Ls, n_array);

        Convert_n_array_to_bit_array(Length, n_array, bit_array);

    for(int i_fermion=0;i_fermion<basis.Dec_dnup_Fermions[SzSetNo].size();i_fermion++){

        i_local = basis.Concatenate_Spins_and_Fermions(i_fermion, i_Lspins, SzSetNo);
        m = basis.Get_basis_ind(i_local, SzSetNo);

        d_dn_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].first;
        d_up_= basis.Dec_dnup_Fermions[SzSetNo][i_fermion].second;

        FermionsSzWithOffset=basis.SzWithOffsetAllowed[SzSetNo].second;
        N_up = FermionsSzWithOffset;
        N_dn = Target_Total_Ne-N_up;

        BASE_DN = Length-N_dn+2;
        BASE_UP = Length-N_up+2;

        from_deci_type2_to_n_array(d_dn_, BASE_DN, n_dn_array);
        from_deci_type2_to_n_array(d_up_, BASE_UP, n_up_array);

        Convert_n_array_to_bit_array(Length, n_dn_array, bit_dn_array);
        Convert_n_array_to_bit_array(Length, n_up_array, bit_up_array);


        //local Kondo

            //splus_i Sminus_j
            if( bit_array[site_j]==1  &&
                bit_dn_array[site_i]==1
                    ){
                assert(bit_up_array[site_i]==0);

                 FermionsSzWithOffset_new = FermionsSzWithOffset +1;
                 LocalSpinsSzWithOffset_new=LocalSpinsSzWithOffset -1;

                 SzSetNo_new = basis.Inverse_SzWithOffsetAllowed[LocalSpinsSzWithOffset_new + (Length+1)*(FermionsSzWithOffset_new)];

                 bit_array_new=bit_array;
                 bit_dn_array_new=bit_dn_array;
                 bit_up_array_new=bit_up_array;

                 bit_array_new[site_j]=0;
                 bit_dn_array_new[site_i]=0;
                 bit_up_array_new[site_i]=1;

                 BASE_NEW_DN = Length-(N_dn-1)+2;
                 BASE_NEW_UP = Length-(N_up+1)+2;

                 d_up_new_ = FromBitArray_toDeci_type2(bit_up_array_new, BASE_NEW_UP);
                 d_dn_new_ = FromBitArray_toDeci_type2(bit_dn_array_new, BASE_NEW_DN);

                 i_fermion_new = basis.InverseDec_Fermions[SzSetNo_new][d_dn_new_ + (basis.Dmax_dn[SzSetNo_new]+1)*d_up_new_];

                 base_Ls_new=Length-LocalSpinsSzWithOffset_new+2;
                 dec_i_Ls_new = FromBitArray_toDeci_type2(bit_array_new, base_Ls_new);

                 i_Lspins_new = basis.InverseDec_LocalizedSpins[SzSetNo_new][dec_i_Ls_new];

                 i_local_new = basis.Concatenate_Spins_and_Fermions(i_fermion_new, i_Lspins_new, SzSetNo_new);
                 m_new = basis.Get_basis_ind(i_local_new, SzSetNo_new);


                 //basis: up ...down
                 //bits_in_bw_up = Count_bits_in_bw(site_i, "closed", Length, "open", bit_up_array);
                 //bits_in_bw_dn = Count_bits_in_bw(0, "closed", site_i, "open",bit_dn_array);

                 //basis: down ...up
                 bits_in_bw_dn = Count_bits_in_bw(site_i, "closed", Length, "open", bit_dn_array) + 1;
                 bits_in_bw_up = Count_bits_in_bw(0, "closed", site_i, "open", bit_up_array);

                 FM_SIGN = pow(-1.0,(1.0*(bits_in_bw_up + bits_in_bw_dn)));

                 value = FM_SIGN;

                 if(value!=0){
                     //assert(m_new>=m);
                     OPR.value.push_back(value*one);
                     OPR.rows.push_back(m);
                     OPR.columns.push_back(m_new);
                 }

            }

    }
    }
    }


}


void MODEL_KondoModel::Read_parameters(BASIS_KondoModel &basis, string filename){


    //    bool read_basis, write_basis;
    //    string read_basis_file, write_basis_file;

    string filepath = filename;

    string read_basis_, Read_Basis_= "Read_Basis = ";
    string write_basis_, Write_Basis_= "Write_Basis = ";
    string Read_Basis_File_= "Read_Basis_File = ";
    string Write_Basis_File_= "Write_Basis_File = ";

    string length, Length_str = "Length = ";

    string targettotalsz, TargetTotalSz = "TotalSzTarget = ";
    string targettotalne, TargetTotalNe = "Target_Total_Ne = ";

    string hmag, Hmag = "H_mag = ";

    string local_kondo, Local_Kondo = "LocalKondo = ";

    string KondoExchange_filepath, KondoExchange_file_ = "KondoExchange_file = ";
    string KondoHopping_filepath, KondoHopping_file_ = "KondoHopping_file = ";

    string LongRangeJ_LS_filepath ,LongRangeJ_LS_file_ = "LongRangeJ_LocalSpins_file = ";
    string LongRange_Hopping_filepath ,LongRange_Hopping_file_ = "LongRange_Hopping_file = ";



    int offset;
    string line;
    ifstream inputfile(filepath.c_str());

    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(LongRange_Hopping_file_, 0)) != string::npos) {
                LongRange_Hopping_filepath = line.substr (offset+LongRange_Hopping_file_.length());  }

            if ((offset = line.find(LongRangeJ_LS_file_, 0)) != string::npos) {
                LongRangeJ_LS_filepath = line.substr (offset+LongRangeJ_LS_file_.length());  }

            if ((offset = line.find(KondoExchange_file_, 0)) != string::npos) {
                KondoExchange_filepath = line.substr (offset+KondoExchange_file_.length());  }

            if ((offset = line.find(KondoHopping_file_, 0)) != string::npos) {
                KondoHopping_filepath = line.substr (offset+KondoHopping_file_.length());  }


            if ((offset = line.find(Read_Basis_, 0)) != string::npos) {
                read_basis_ = line.substr (offset+Read_Basis_.length());				}

            if ((offset = line.find(Write_Basis_, 0)) != string::npos) {
                write_basis_ = line.substr (offset+Write_Basis_.length());				}

            if ((offset = line.find(Read_Basis_File_, 0)) != string::npos) {
                basis.read_basis_file = line.substr (offset+Read_Basis_File_.length());				}

            if ((offset = line.find(Write_Basis_File_, 0)) != string::npos) {
                basis.write_basis_file = line.substr (offset+Write_Basis_File_.length());				}


            if ((offset = line.find(Length_str, 0)) != string::npos) {
                length = line.substr (offset + Length_str.length());		}


            if ((offset = line.find(TargetTotalSz, 0)) != string::npos) {
                targettotalsz = line.substr (offset + TargetTotalSz.length());		}


            if ((offset = line.find(TargetTotalNe, 0)) != string::npos) {
                targettotalne = line.substr (offset + TargetTotalNe.length());		}


            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(Local_Kondo, 0)) != string::npos) {
                local_kondo = line.substr (offset + Local_Kondo.length());		}



        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}



    if(read_basis_ == "true"){
        basis.read_basis =true;
    }
    else{
        basis.read_basis =false;
    }

    if(write_basis_ == "true"){
        basis.write_basis =true;
    }
    else{
        basis.write_basis =false;
    }


    basis.Length=atoi(length.c_str());
    basis.Target_Total_Ne=atoi(targettotalne.c_str());
    basis.Target_Total_Sz = (atof(targettotalsz.c_str()));

    LocalKondo = atof(local_kondo.c_str());

    Length=basis.Length;
    Target_Total_Ne=basis.Target_Total_Ne;
    Target_Total_Sz=basis.Target_Total_Sz;


    J_LSpins_mat.resize(basis.Length);
    for(int i =0;i<basis.Length;i++){
        J_LSpins_mat[i].resize(basis.Length);
        for(int j =0;j<basis.Length;j++){
            J_LSpins_mat[i][j]=0.0;
        }
    }


    Hopping_mat.resize(basis.Length);
    for(int i =0;i<basis.Length;i++){
        Hopping_mat[i].resize(basis.Length);
        for(int j =0;j<basis.Length;j++){
            Hopping_mat[i][j]=0.0;
        }
    }

    KondoExchange_mat.resize(basis.Length);
    for(int i =0;i<basis.Length;i++){
        KondoExchange_mat[i].resize(basis.Length);
        for(int j =0;j<basis.Length;j++){
            KondoExchange_mat[i][j]=0.0;
        }
    }

    ifstream inputJ1(LongRangeJ_LS_filepath.c_str());

    string line_temp;
    int i1, i2;
    double value_temp;
/*
    while(getline(inputJ1, line_temp)){
        stringstream line_temp_ss(line_temp);
        line_temp_ss >> i1 >> i2 >> value_temp;
        J_LSpins_mat[i1][i2]=value_temp;
    }
*/

    for(int i1=0;i1<Length;i1++){
        for(int i2=0;i2<Length;i2++){
            inputJ1>>J_LSpins_mat[i1][i2];
        }
    }


    ifstream inputHopping(LongRange_Hopping_filepath.c_str());
    for(int i1=0;i1<Length;i1++){
        for(int i2=0;i2<Length;i2++){
            inputHopping>>Hopping_mat[i1][i2];
        }
    }



    ifstream inputKondoExchange(KondoExchange_filepath.c_str());
    for(int i1=0;i1<Length;i1++){
        for(int i2=0;i2<Length;i2++){
            inputKondoExchange>>KondoExchange_mat[i1][i2];
        }
    }

    ifstream inputKondoHopping(KondoHopping_filepath.c_str());
    int site_i, site_l,site_j;
    double val_;
    triad_int temp_int_triad;
    KondoHoppings_sites.clear();
    KondoHoppings.clear();
    while(getline(inputKondoHopping,line_temp)){
    stringstream line_stream;
    line_stream<<line_temp;
    line_stream>>site_i>>site_l>>site_j>>val_;
    temp_int_triad.first_=site_i;
    temp_int_triad.second_=site_l;
    temp_int_triad.third_=site_j;
    KondoHoppings_sites.push_back(temp_int_triad);
    KondoHoppings.push_back(val_);
    }


}


/*
void MODEL_1_orb_tJ::Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis){

//      int T_no_oprs=2;
//    int orb;
//    int spin;



//   //  0 n
//   //  1 Sz

//    One_point_oprts.resize(T_no_oprs);
//    for(int i=0;i<T_no_oprs;i++){
//        One_point_oprts[i].resize(basis.Length);
//    }

}


void MODEL_1_orb_tJ::Initialize_two_point_to_calculate(BASIS_1_orb_tJ &basis){
    two_point_obs.resize(1);
    two_point_obs[0]="SpSm";


    Two_point_oprts.resize(1);
    int T_no_oprs=1;



    for(int i=0;i<T_no_oprs;i++){
        Two_point_oprts[i].resize(basis.Length);
        for(int j=0;j<basis.Length;j++){
            Two_point_oprts[i][j].resize(basis.Length);
        }
    }




    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        if(two_point_obs[opr_no]=="SpSm"){

            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;
            double value;


            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){

                        m=i ;

                        value=0;

                        value+=0.25*( ( bit_value(basis.D_up_basis[i], site) -
                                        bit_value(basis.D_dn_basis[j], site) )*
                                      ( bit_value(basis.D_up_basis[i], site2) -
                                        bit_value(basis.D_dn_basis[j], site2) )
                                      );

                        if(value!=0){
                            Two_point_oprts[opr_no][site][site2].value.push_back(value*one);
                            Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                            Two_point_oprts[opr_no][site][site2].columns.push_back(m);
                        }
                    }
                }
            }
        }
    }
}

void MODEL_1_orb_tJ::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR,
                                                                  int site, int site2, BASIS_1_orb_tJ &basis){

    OPR.nrows = basis.D_up_basis.size();
    OPR.ncols = OPR.nrows;


    //Remember OPR[l][m]=<l|OPR|m>
    int m,j,nup_temp,ndn_temp;
    double value;


    if(opr_type=="SzSz"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            j=i;
            m=i;
            value=0;
            value+=0.25*( ( bit_value(basis.D_up_basis[i], site) -
                            bit_value(basis.D_dn_basis[j], site) )*
                          ( bit_value(basis.D_up_basis[i], site2) -
                            bit_value(basis.D_dn_basis[j], site2) )
                          );
            if(value!=0){
                OPR.value.push_back(value);
                OPR.rows.push_back(m);
                OPR.columns.push_back(m);
            }

        }
    }


    if(opr_type=="SpSm"){
        //Remember OPR[l][m]=<l|OPR|m>
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;
        for (int i=0;i<basis.D_up_basis.size();i++){
            j=i;
            m=i;

            //Sp_site[site]*Sm_site[site2]:
            //there have to be ony up electron at site2
            //there have to be only down electron at site

            if(((bit_value(basis.D_dn_basis[j], site)==1)
                &&
                (bit_value(basis.D_up_basis[i], site)==0)
                )
                    &&
                    ((bit_value(basis.D_up_basis[i], site2)==1)
                     &&
                     (bit_value(basis.D_dn_basis[j], site2)==0)
                     ))
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2, site2)
                              + pow(2, site) );
                D_dn = (int) (basis.D_dn_basis[j] + pow(2, site2)
                              - pow(2, site) );

                nup_temp = __builtin_popcount(D_up);
                ndn_temp = __builtin_popcount(D_dn);
                m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                l= site;
                lp= site2;

                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);

                //assert(m_new<m);
                OPR.value.push_back(sign_FM);
                OPR.rows.push_back(m_new);
                OPR.columns.push_back(m);
            }

            if((site==site2)){
                if(
                        ((bit_value(basis.D_up_basis[i], site2)==1)
                         &&
                         (bit_value(basis.D_dn_basis[j], site2)==0)
                         )
                        )
                {
                    OPR.value.push_back(1.0);
                    OPR.rows.push_back(m);
                    OPR.columns.push_back(m);
                }
            }
        }
    }

    if(opr_type=="SmSp"){
        //Remember OPR[l][m]=<l|OPR|m>
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;
        for (int i=0;i<basis.D_up_basis.size();i++){
            j=i;
            m=i;

//            if site !=site2
//                Sm[site]*Sp[site2]=Sp[site2]*Sm[site]

            //there have to be ony up electron at site2
            //there have to be only down electron at site

            if(((bit_value(basis.D_dn_basis[j], site2)==1)
                &&
                (bit_value(basis.D_up_basis[i], site2)==0)
                )
                    &&
                    ((bit_value(basis.D_up_basis[i], site)==1)
                     &&
                     (bit_value(basis.D_dn_basis[j], site)==0)
                     ))
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2, site)
                              + pow(2, site2) );
                D_dn = (int) (basis.D_dn_basis[j] + pow(2, site)
                              - pow(2, site2) );

                nup_temp = __builtin_popcount(D_up);
                ndn_temp = __builtin_popcount(D_dn);
                m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                l= site2;
                lp= site;

                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);

                //assert(m_new<m);
                OPR.value.push_back(sign_FM);
                OPR.rows.push_back(m_new);
                OPR.columns.push_back(m);
            }

            if((site==site2)){
                if(
                        ((bit_value(basis.D_up_basis[i], site2)==0)
                         &&
                         (bit_value(basis.D_dn_basis[j], site2)==1)
                         )
                        )
                {
                    OPR.value.push_back(1.0);
                    OPR.rows.push_back(m);
                    OPR.columns.push_back(m);
                }
            }
        }
    }


}

void MODEL_1_orb_tJ::Initialize_Opr_for_Dynamics(BASIS_1_orb_tJ &basis){



}

*/
#endif

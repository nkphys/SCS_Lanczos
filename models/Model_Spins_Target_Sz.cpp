/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "Model_Spins_Target_Sz.h"
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
#define PI 3.14159265



void MODEL_Spins_Target_Sz::Act_Hamil(BASIS_Spins_Target_Sz &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    cout<<"NOT WORKING AT PRESENT"<<endl;

}

void MODEL_Spins_Target_Sz::Add_diagonal_terms(BASIS_Spins_Target_Sz &basis){



    //    //Remember H[l][m]=<l|H|m>

        double_type value;
      for (int i=basis.D_min;i<basis.D_max + 1;i++){

            value=zero;

    //        //magnetic Field in Z-direction
    //        for(int site=0;site<basis.Length;site++){
    //            value+=one*(H_field[2])*
    //                    ( ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );
    //        }


            //(Sz_local)^2 exchange
            for(int site_i=0;site_i<basis.Length;site_i++){
                if(abs(Dz_anisotropy[site_i])>0.0000001){
                    value+=  Dz_anisotropy[site_i]*
                            ( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                              ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                              );
                }
            }

            if(abs(value)>0.000001){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(i);
                Hamil.columns.push_back(i);
            }

       }


}


void MODEL_Spins_Target_Sz::Add_non_diagonal_terms(BASIS_Spins_Target_Sz &basis){


}


void MODEL_Spins_Target_Sz::Add_connections_strictly2point(BASIS_Spins_Target_Sz &basis){

    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;

    Mat_1_doub J1_vals, J2_vals, J3_vals;
    Mat_1_tetra_int J2_sites,J3_sites;
    Mat_1_intpair J1_sites;
    pair_int temp_pair_int;
    tetra_int temp_tetra_int;

    //TERM J1: J1[k][l] x Svec[k].Svec[l] |m>
    for(int site_k=0;site_k<basis.Length;site_k++){
        for(int site_l=0;site_l<basis.Length;site_l++){
            if(J1_mat[site_k][site_l]!=0.0){
                J1_vals.push_back(one*J1_mat[site_k][site_l]);
                temp_pair_int.first=site_k;
                temp_pair_int.second=site_l;
                J1_sites.push_back(temp_pair_int);
            }
        }
    }


    //TERM J2 and J3
    for(int site_i=0;site_i<basis.Length;site_i++){
        for(int site_j=0;site_j<basis.Length;site_j++){
            for(int site_k=0;site_k<basis.Length;site_k++){
                for(int site_l=0;site_l<basis.Length;site_l++){

                    if(J2_mat[site_i][site_j][site_k][site_l]!=0.0){
                        J2_vals.push_back(one*J2_mat[site_i][site_j][site_k][site_l]);
                        temp_tetra_int.first=site_i;
                        temp_tetra_int.second=site_j;
                        temp_tetra_int.third=site_k;
                        temp_tetra_int.fourth=site_l;
                        J2_sites.push_back(temp_tetra_int);
                    }

                    if(J3_mat[site_i][site_j][site_k][site_l]!=0.0){
                        J3_vals.push_back(one*J3_mat[site_i][site_j][site_k][site_l]);
                        temp_tetra_int.first=site_i;
                        temp_tetra_int.second=site_j;
                        temp_tetra_int.third=site_k;
                        temp_tetra_int.fourth=site_l;
                        J3_sites.push_back(temp_tetra_int);
                    }

                }
            }
        }
    }


    assert(J3_sites.size()==0);
    assert(J2_sites.size()==0);



    Hamiltonian_1_COO Hamil_private;
    int N_threads=1;
#ifdef _OPENMP
    N_threads=no_of_proc;
#endif

    Hamil_private.resize(N_threads);


#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id;
    ulli dec_new, dec_new_temp;
    ulli dec_m;
    ulli dec_max;
    ulli dec_part0_, dec_part1_;
    int m_new;
    int n_index;
    Mat_1_int state_vec;
    double_type value;
    int val_site_k, val_site_l;
    int val_site_k_new, val_site_l_new;


#ifdef _OPENMP
#pragma omp for
#endif
    for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){

        thread_id=0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif

        //dec_m = basis.D_basis[m];
        dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];


        //Sz(i)Sz(j)----------------------------------------------------------------
        value=0.0;
        for(int site_k=0;site_k<basis.Length;site_k++){
            for(int site_l=0;site_l<basis.Length;site_l++){

                if(J1_mat[site_k][site_l]!=0.0){
                    value += J1_mat[site_k][site_l]*((1.0*value_at_pos(dec_m, site_k, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                            ((1.0*value_at_pos(dec_m, site_l, basis.BASE)) - (0.5*basis.TwoTimesSpin));
                }

        }}

        // Hamil.value.push_back(value*one);
        // Hamil.rows.push_back(m);
        // Hamil.columns.push_back(m);
        Hamil_private[thread_id].value.push_back(value*one);
        Hamil_private[thread_id].rows.push_back(m);
        Hamil_private[thread_id].columns.push_back(m);

        //0.5*S+(k)S-(l)---------------------------------------------------
        value=0.0;

        for(int site_k=0;site_k<basis.Length;site_k++){
            for(int site_l=0;site_l<basis.Length;site_l++){
                val_site_k = value_at_pos(dec_m, site_k, basis.BASE);
                val_site_l = value_at_pos(dec_m, site_l, basis.BASE);
                if(J1_mat[site_k][site_l]!=0.0){
                assert(site_k!=site_l);
                bool allowed = ((val_site_l != 0) && (val_site_k != (basis.BASE - 1)));
                    if(allowed){
                    val_site_k_new = val_site_k + 1;
                    val_site_l_new = val_site_l - 1;

                    dec_new_temp = Updated_decimal_with_value_at_pos(dec_m, site_k, basis.BASE, val_site_k_new);
                    dec_new = Updated_decimal_with_value_at_pos(dec_new_temp, site_l, basis.BASE, val_site_l_new);

                    value = one*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                       ((val_site_k - (0.5*basis.TwoTimesSpin))*
                                        (val_site_k_new - (0.5*basis.TwoTimesSpin))  ) );
                    value = J1_mat[site_k][site_l]*value*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                               ((val_site_l - (0.5*basis.TwoTimesSpin))*
                                                (val_site_l_new - (0.5*basis.TwoTimesSpin))  ) );




                    //m_new = Find_int_in_intarray(dec_new, basis.D_basis);
                    //m_new = Find_int_in_intarray_using_bisection(dec_new, basis.D_basis);

                    //m_new = Find_int_in_intarray_using_multisections(dec_new, basis.D_basis, MultisectionSearch_int);


                    dec_part1_ = dec_new/pow(basis.BASE, basis.Partition_Length[0]);
                    dec_part0_ = dec_new - (dec_part1_*pow(basis.BASE, basis.Partition_Length[0]));
                    m_new = basis.Dec_to_Index_part0_[dec_part0_] + basis.Dec_to_Index_part1_[dec_part1_];

                    /*
                    fromDeci_to_Vecint(state_vec, basis.BASE, dec_new , basis.Length);
                    quicksort(state_vec, 0, state_vec.size() -1);
                    fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                    n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                    m_new = Find_int_in_part_of_intarray(dec_new, basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);
                    */

                    assert(m_new<m);
                    if(m_new<m){
                        Hamil_private[thread_id].value.push_back(value*one);
                        Hamil_private[thread_id].rows.push_back(m_new);
                        Hamil_private[thread_id].columns.push_back(m);
                        // Hamil.value.push_back(value*one);
                        // Hamil.rows.push_back(m);
                        // Hamil.columns.push_back(m_new);
                    }

                    }


                    }
        }}


        if(m%1000==0){
            cout<<m<<" basis done in Hamil construction by thread "<< thread_id <<endl;
        }

    }
#ifdef _OPENMP
        }
#endif


        // cout<<"HERE"<<endl;
        for(int thread=0;thread<N_threads;thread++){
            Hamil.value.insert(Hamil.value.end(),Hamil_private[thread].value.begin(), Hamil_private[thread].value.end() );
            Hamil.rows.insert(Hamil.rows.end(),Hamil_private[thread].rows.begin(), Hamil_private[thread].rows.end() );
            Hamil.columns.insert(Hamil.columns.end(),Hamil_private[thread].columns.begin(), Hamil_private[thread].columns.end() );
        }



}


void MODEL_Spins_Target_Sz::Add_arbitraryconnections_from_files(BASIS_Spins_Target_Sz &basis){

    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;


    double EPS_=0.000001;
    int thread_id;
    ulli dec_m;
    ulli dec_max;
    int m_new;
    int n_index;
    Mat_1_int state_vec;

    Mat_1_ullint dec_vec_out_given_m;
    Mat_1_doub val_vec_out_given_m;
    Mat_1_ullint dec_vec_out_final;
    Mat_1_doub val_vec_out_final;

    Mat_1_ullint dec_vec_out;
    Mat_1_doub val_vec_out;
    Mat_1_string oprs_list;
    Mat_1_int oprs_site;
    //stringstream connection_stream;
    string temp_opr_str;


    //abort();


     for (int m=0;m<basis.D_basis.size();m++){

         dec_m = basis.D_basis[m];
         dec_vec_out_given_m.clear();
         val_vec_out_given_m.clear();

         for(int FileNo=0;FileNo<ConnectionFiles.size();FileNo++){
             for(int connection_no=0;connection_no<Connections[FileNo].size();connection_no++){
               // cout<<"here 2"<<endl;
       // connection_stream.str("");

                stringstream connection_stream;
                connection_stream<<Connections[FileNo][connection_no];

                int n_oprs;
                oprs_list.clear();
                oprs_site.clear();


                double_type connection_val;
                connection_stream>>n_oprs;
                oprs_list.resize(n_oprs);
                oprs_site.resize(n_oprs);

                for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                    connection_stream>>temp_opr_str;
                    oprs_list[opr_no]=temp_opr_str;
                }
                for(int opr_no=(n_oprs-1);opr_no>=0;opr_no--){
                    connection_stream>>oprs_site[opr_no];
                }
                connection_stream>>connection_val;

                dec_vec_out.clear();
                val_vec_out.clear();
                dec_vec_out.push_back(dec_m);
                val_vec_out.push_back(connection_val);


                Act_LocalOprString_by_Recursion(oprs_site, oprs_list, dec_vec_out, val_vec_out, basis, 0);


                for(int index_temp=0;index_temp<dec_vec_out.size();index_temp++){
                    dec_vec_out_given_m.push_back(dec_vec_out[index_temp]);
                    val_vec_out_given_m.push_back(val_vec_out[index_temp]);
                }



               // cout<<"here 1"<<endl;
             }
         }



         Remove_repetitions(dec_vec_out_given_m, val_vec_out_given_m, dec_vec_out_final, val_vec_out_final);

         for(int j=0;j<dec_vec_out_final.size();j++){

             if(abs(val_vec_out_final[j])>EPS_){
             fromDeci_to_Vecint(state_vec, basis.BASE, dec_vec_out_final[j] , basis.Length);
             quicksort(state_vec, 0, state_vec.size() -1);
             fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
             n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
             m_new = Find_int_in_part_of_intarray(dec_vec_out_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


             if(m_new<=m){
                 Hamil.value.push_back(val_vec_out_final[j]*one);
                 Hamil.rows.push_back(m_new);
                 Hamil.columns.push_back(m);
             }
             }

         }




     }


     state_vec.clear();
    dec_vec_out.clear();
     val_vec_out.clear();
   // connection_stream.str(std::string());


}


void MODEL_Spins_Target_Sz::Get_LocalOPR(Matrix<double_type> &OPR_, string opr_type){

    if(opr_type=="Sz"){
        OPR_=SzLocal;
    }
    else if(opr_type=="Sx"){
        OPR_=SxLocal;
    }
    else if(opr_type=="SyTimesIota"){
        OPR_=SyLocalTimesIota;
    }
    else if(opr_type=="Sp"){
        OPR_=SplusLocal;
    }
    else if(opr_type=="Sm"){
        OPR_=SminusLocal;
    }
    else if (opr_type=="Sx2"){
        OPR_=Sx2Local;
    }
    else if (opr_type=="Sy2"){
        OPR_=Sy2Local;
    }
    else if (opr_type=="Sz2"){
        OPR_=Sz2Local;
    }
    else if (opr_type=="Qxz"){
        OPR_=QxzLocal;
    }
    else if (opr_type=="QyzI"){
        OPR_=QyzLocalTimesIota;
    }
    else if (opr_type=="QxyI"){
        OPR_=QxyLocalTimesIota;
    }
    else{
        cout<<"OPR TYPE "<<opr_type<<" not present"<<endl;
    }

}

void MODEL_Spins_Target_Sz::Act_LocalOprString_by_Recursion(Mat_1_int oprs_site, Mat_1_string oprs_list, Mat_1_ullint & dec_vec_out, Mat_1_doub &val_vec_out,BASIS_Spins_Target_Sz & basis, int opr_no){

    double EPS_=0.0000001;


    if(opr_no<oprs_list.size()){

        Matrix<double_type> OPR_;
        Get_LocalOPR(OPR_, oprs_list[opr_no]);

    //act opr on dec_vec and update dec_vec_out and val_vec_out
        Mat_1_ullint dec_vec_out_temp_final;
        Mat_1_doub val_vec_out_temp_final;
        Mat_1_ullint dec_vec_out_temp;
        Mat_1_doub val_vec_out_temp;

        ulli dec_old, dec_new_temp;
        int val_site_old;
        for(int dec_index=0;dec_index<dec_vec_out.size();dec_index++){
            dec_old = dec_vec_out[dec_index];
            val_site_old = value_at_pos(dec_old, oprs_site[opr_no], basis.BASE);

            for(int val_col=0;val_col<basis.BASE;val_col++){
                if(abs(OPR_(val_col,val_site_old))>EPS_){
                dec_new_temp = Updated_decimal_with_value_at_pos(dec_old, oprs_site[opr_no], basis.BASE, val_col);
                dec_vec_out_temp.push_back(dec_new_temp);
                val_vec_out_temp.push_back(OPR_(val_col,val_site_old)*val_vec_out[dec_index]);
                }
            }
        }

        Remove_repetitions(dec_vec_out_temp, val_vec_out_temp, dec_vec_out_temp_final, val_vec_out_temp_final);
        dec_vec_out=dec_vec_out_temp_final;
        val_vec_out=val_vec_out_temp_final;

    }

    if((opr_no+1)<oprs_list.size()){
        Act_LocalOprString_by_Recursion(oprs_site, oprs_list, dec_vec_out, val_vec_out, basis, opr_no+1);
    }



}


void MODEL_Spins_Target_Sz::Add_connections_new(BASIS_Spins_Target_Sz &basis){



    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;

    Mat_1_doub J1_vals, J2_vals, J3_vals;
    Mat_1_tetra_int J2_sites,J3_sites;
    Mat_1_intpair J1_sites;
    pair_int temp_pair_int;
    tetra_int temp_tetra_int;

    //TERM J1: J1[k][l] x Svec[k].Svec[l] |m>
    for(int site_k=0;site_k<basis.Length;site_k++){
        for(int site_l=0;site_l<basis.Length;site_l++){
            if(J1_mat[site_k][site_l]!=0.0){
                J1_vals.push_back(one*J1_mat[site_k][site_l]);
                temp_pair_int.first=site_k;
                temp_pair_int.second=site_l;
                J1_sites.push_back(temp_pair_int);
            }
        }
    }

    //TERM J2 and J3
    for(int site_i=0;site_i<basis.Length;site_i++){
        for(int site_j=0;site_j<basis.Length;site_j++){
            for(int site_k=0;site_k<basis.Length;site_k++){
                for(int site_l=0;site_l<basis.Length;site_l++){

                    if(J2_mat[site_i][site_j][site_k][site_l]!=0.0){
                        J2_vals.push_back(one*J2_mat[site_i][site_j][site_k][site_l]);
                        temp_tetra_int.first=site_i;
                        temp_tetra_int.second=site_j;
                        temp_tetra_int.third=site_k;
                        temp_tetra_int.fourth=site_l;
                        J2_sites.push_back(temp_tetra_int);
                    }

                    if(J3_mat[site_i][site_j][site_k][site_l]!=0.0){
                        J3_vals.push_back(one*J3_mat[site_i][site_j][site_k][site_l]);
                        temp_tetra_int.first=site_i;
                        temp_tetra_int.second=site_j;
                        temp_tetra_int.third=site_k;
                        temp_tetra_int.fourth=site_l;
                        J3_sites.push_back(temp_tetra_int);
                    }

                }
            }
        }
    }





    Hamiltonian_1_COO Hamil_private;
    int N_threads=1;
#ifdef _OPENMP
    N_threads=no_of_proc;
#endif

    Hamil_private.resize(N_threads);


#ifdef _OPENMP
#pragma omp parallel
    {
#endif
        int thread_id;
        ulli dec_m;
        Mat_1_ullint m_connected_final;
        Mat_1_doub coeffs_final;
        Mat_1_ullint m_connected_temp2;
        Mat_1_doub coeffs_temp2;
        Mat_1_ullint m_connected_temp3;
        Mat_1_doub coeffs_temp3;
        ulli mi,mj;
        ulli dec_max;
        int m_new;
        int n_index;
        Mat_1_int state_vec;
        Mat_1_ullint m_connected;
        Mat_1_doub coeffs;
        Mat_1_ullint m_connected_temp;
        Mat_1_doub coeffs_temp;
        //abort();

#ifdef _OPENMP
#pragma omp for
#endif
        for (int m=0;m<basis.D_basis.size();m++){

            thread_id=0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif

            dec_m = basis.D_basis[m];
            m_connected.clear();
            coeffs.clear();


            //TERM J1: J1[k][l] x Svec[k].Svec[l] |m>
            for(int site_k=0;site_k<basis.Length;site_k++){
                for(int site_l=0;site_l<basis.Length;site_l++){

                    if(J1_mat[site_k][site_l]!=0.0){
                        Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                        value_multiply_vector(one*J1_mat[site_k][site_l], coeffs_temp);
                        m_connected.insert(m_connected.end(), m_connected_temp.begin(), m_connected_temp.end());
                        coeffs.insert(coeffs.end(), coeffs_temp.begin(), coeffs_temp.end());
                    }
                }
            }



            //TERM J2: J2 [i][j][k][l] X (Svec[i].Svec[j]) X (Svec[k].Svec[l]) |m>

            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){
                    for(int site_k=0;site_k<basis.Length;site_k++){
                        for(int site_l=0;site_l<basis.Length;site_l++){
                            if(J2_mat[site_i][site_j][site_k][site_l]!=0.0){

                                Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                                value_multiply_vector(one*J2_mat[site_i][site_j][site_k][site_l], coeffs_temp);
                                for(int i_=0;i_<m_connected_temp.size();i_++){
                                    mi = m_connected_temp[i_];
                                    Act_SiSj(site_i, site_j, mi, m_connected_temp2, coeffs_temp2, basis);

                                    value_multiply_vector(coeffs_temp[i_], coeffs_temp2);
                                    m_connected.insert(m_connected.end(), m_connected_temp2.begin(), m_connected_temp2.end());
                                    coeffs.insert(coeffs.end(), coeffs_temp2.begin(), coeffs_temp2.end());

                                }

                            }
                        }
                    }
                }
            }



            //TERM J3: J3 [i][j][k][l] X (Svec[i].Svec[j]) X (Svec[k].Svec[l]) X (Svec[k].Svec[l])|m>
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){
                    for(int site_k=0;site_k<basis.Length;site_k++){
                        for(int site_l=0;site_l<basis.Length;site_l++){
                            if(J3_mat[site_i][site_j][site_k][site_l]!=0.0){

                                Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                                value_multiply_vector(one*J3_mat[site_i][site_j][site_k][site_l], coeffs_temp);
                                for(int i_=0;i_<m_connected_temp.size();i_++){
                                    mi = m_connected_temp[i_];
                                    Act_SiSj(site_k, site_l, mi, m_connected_temp2, coeffs_temp2, basis);
                                    value_multiply_vector(coeffs_temp[i_], coeffs_temp2);

                                    for(int j_=0;j_<m_connected_temp2.size();j_++){
                                        mj = m_connected_temp2[j_];
                                        Act_SiSj(site_i, site_j, mj, m_connected_temp3, coeffs_temp3, basis);
                                        value_multiply_vector(coeffs_temp2[j_], coeffs_temp3);

                                        m_connected.insert(m_connected.end(), m_connected_temp3.begin(), m_connected_temp3.end());
                                        coeffs.insert(coeffs.end(), coeffs_temp3.begin(), coeffs_temp3.end());

                                    }
                                }
                            }
                        }
                    }
                }
            }



            //remove_repetitions
            Remove_repetitions(m_connected, coeffs, m_connected_final, coeffs_final);


            for(int j=0;j<m_connected_final.size();j++){

                //40sec for L=6
                //m_new = Find_int_in_intarray(m_connected_final[j], basis.D_basis);


                //25sec for L=6
                fromDeci_to_Vecint(state_vec, basis.BASE,m_connected_final[j] , basis.Length);
                quicksort(state_vec, 0, state_vec.size() -1);
                fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                m_new = Find_int_in_part_of_intarray(m_connected_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


                if(m_new<=m){
                    Hamil_private[thread_id].value.push_back(coeffs_final[j]*one);
                    Hamil_private[thread_id].rows.push_back(m_new);
                    Hamil_private[thread_id].columns.push_back(m);
                }
            }

            if(m%500==0){
                cout<<m<<" basis done in Hamil construction by thread "<< thread_id <<endl;
            }


        } // m i.e. basis index, columns of H


#ifdef _OPENMP
    }
#endif


    // cout<<"HERE"<<endl;
    for(int thread=0;thread<N_threads;thread++){
        Hamil.value.insert(Hamil.value.end(),Hamil_private[thread].value.begin(), Hamil_private[thread].value.end() );
        Hamil.rows.insert(Hamil.rows.end(),Hamil_private[thread].rows.begin(), Hamil_private[thread].rows.end() );
        Hamil.columns.insert(Hamil.columns.end(),Hamil_private[thread].columns.begin(), Hamil_private[thread].columns.end() );
    }

    //assert(false);

}

void MODEL_Spins_Target_Sz::Add_connections(BASIS_Spins_Target_Sz &basis){

    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;

    Mat_1_ullint m_connected;
    Mat_1_doub coeffs;
    Mat_1_ullint m_connected_final;
    Mat_1_doub coeffs_final;

    Mat_1_ullint m_connected_temp;
    Mat_1_doub coeffs_temp;
    Mat_1_ullint m_connected_temp2;
    Mat_1_doub coeffs_temp2;
    Mat_1_ullint m_connected_temp3;
    Mat_1_doub coeffs_temp3;
    ulli mi,mj;
    ulli dec_max, dec_m;
    int m_new;
    int n_index;
    Mat_1_int state_vec;


    for (int m=0;m<basis.D_basis.size();m++){

        dec_m = basis.D_basis[m];
        m_connected.clear();
        coeffs.clear();


        //TERM J1: J1[k][l] x Svec[k].Svec[l] |m>
        for(int site_k=0;site_k<basis.Length;site_k++){
            for(int site_l=0;site_l<basis.Length;site_l++){

                if(J1_mat[site_k][site_l]!=0.0){
                    Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                    value_multiply_vector(one*J1_mat[site_k][site_l], coeffs_temp);
                    m_connected.insert(m_connected.end(), m_connected_temp.begin(), m_connected_temp.end());
                    coeffs.insert(coeffs.end(), coeffs_temp.begin(), coeffs_temp.end());
                }
            }
        }



        //TERM J2: J2 [i][j][k][l] X (Svec[i].Svec[j]) X (Svec[k].Svec[l]) |m>

        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                for(int site_k=0;site_k<basis.Length;site_k++){
                    for(int site_l=0;site_l<basis.Length;site_l++){
                        if(J2_mat[site_i][site_j][site_k][site_l]!=0.0){

                            Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                            value_multiply_vector(one*J2_mat[site_i][site_j][site_k][site_l], coeffs_temp);
                            for(int i_=0;i_<m_connected_temp.size();i_++){
                                mi = m_connected_temp[i_];
                                Act_SiSj(site_i, site_j, mi, m_connected_temp2, coeffs_temp2, basis);

                                value_multiply_vector(coeffs_temp[i_], coeffs_temp2);
                                m_connected.insert(m_connected.end(), m_connected_temp2.begin(), m_connected_temp2.end());
                                coeffs.insert(coeffs.end(), coeffs_temp2.begin(), coeffs_temp2.end());

                            }

                        }
                    }
                }
            }
        }



        //TERM J3: J3 [i][j][k][l] X (Svec[i].Svec[j]) X (Svec[k].Svec[l]) X (Svec[k].Svec[l])|m>
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                for(int site_k=0;site_k<basis.Length;site_k++){
                    for(int site_l=0;site_l<basis.Length;site_l++){
                        if(J3_mat[site_i][site_j][site_k][site_l]!=0.0){

                            Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
                            value_multiply_vector(one*J3_mat[site_i][site_j][site_k][site_l], coeffs_temp);
                            for(int i_=0;i_<m_connected_temp.size();i_++){
                                mi = m_connected_temp[i_];
                                Act_SiSj(site_k, site_l, mi, m_connected_temp2, coeffs_temp2, basis);
                                value_multiply_vector(coeffs_temp[i_], coeffs_temp2);

                                for(int j_=0;j_<m_connected_temp2.size();j_++){
                                    mj = m_connected_temp2[j_];
                                    Act_SiSj(site_i, site_j, mj, m_connected_temp3, coeffs_temp3, basis);
                                    value_multiply_vector(coeffs_temp2[j_], coeffs_temp3);

                                    m_connected.insert(m_connected.end(), m_connected_temp3.begin(), m_connected_temp3.end());
                                    coeffs.insert(coeffs.end(), coeffs_temp3.begin(), coeffs_temp3.end());

                                }
                            }
                        }
                    }
                }
            }
        }



        //remove_repetitions
        Remove_repetitions(m_connected, coeffs, m_connected_final, coeffs_final);


        for(int j=0;j<m_connected_final.size();j++){

            //40sec for L=6
            //m_new = Find_int_in_intarray(m_connected_final[j], basis.D_basis);


            //25sec for L=6
            fromDeci_to_Vecint(state_vec, basis.BASE,m_connected_final[j] , basis.Length);
            quicksort(state_vec, 0, state_vec.size() -1);
            fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
            n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
            m_new = Find_int_in_part_of_intarray(m_connected_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


            if(m_new<=m){
                Hamil.value.push_back(coeffs_final[j]*one);
                Hamil.rows.push_back(m_new);
                Hamil.columns.push_back(m);
            }
        }

        if(m%500==0){
            cout<<m<<" basis done in Hamil construction"<<endl;
        }


    } // m i.e. basis index, columns of H

    //assert(false);

}


void MODEL_Spins_Target_Sz::Act_translational_opr(BASIS_Spins_Target_Sz &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    assert(Vec_in.size()==basis.D_basis.size());
    Vec_out.resize(Vec_in.size());
    ulli dec_m, dec_new, dec_new_temp, dec_max;
    Mat_1_int state_vec;

    int val_site;
    int n_index, m_new;
    int site, site_new;
    int Length_SL = (int)( basis.Length/2);

    cout <<"This routine works assuming a S-L chain"<<endl;
    for (int m=0;m<basis.D_basis.size();m++){
        dec_m = basis.D_basis[m];

        dec_new=dec_m;
        for(int l=0;l<Length_SL;l++){
            int l_new = l+1;
            l_new = l_new%Length_SL;

            //translating S
             site=2*l;
             site_new=2*l_new;
             val_site = value_at_pos(dec_m, site, basis.BASE);
             dec_new_temp = Updated_decimal_with_value_at_pos(dec_new, site_new, basis.BASE, val_site);
             dec_new = dec_new_temp;


            //translating L
             site=(2*l)+1;
             site_new=(2*l_new)+1;
             val_site = value_at_pos(dec_m, site, basis.BASE);
             dec_new_temp = Updated_decimal_with_value_at_pos(dec_new, site_new, basis.BASE, val_site);
             dec_new = dec_new_temp;

        }


        fromDeci_to_Vecint(state_vec, basis.BASE, dec_new , basis.Length);
        quicksort(state_vec, 0, state_vec.size() -1);
        fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
        n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
        m_new = Find_int_in_part_of_intarray(dec_new, basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


        Vec_out[m_new] = Vec_in[m];


    }


}




double_type MODEL_Spins_Target_Sz::Get_SzSz(int site1, int site2, Mat_1_doub &Vec_, BASIS_Spins_Target_Sz &basis){


    ulli dec_m;
    double_type value=0.0;
   // for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){
        //dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];

        for (int m=0;m<basis.D_basis.size();m++){
        dec_m = basis.D_basis[m];


        //Sz(i)Sz(j)----------------------------------------------------------------
        value +=conjugate(Vec_[m])*Vec_[m]*((1.0*value_at_pos(dec_m, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                ((1.0*value_at_pos(dec_m, site2, basis.BASE)) - (0.5*basis.TwoTimesSpin));

    }

    return value;

}

double_type MODEL_Spins_Target_Sz::Get_Sz(int site1, Mat_1_doub &Vec_, BASIS_Spins_Target_Sz &basis){


    ulli dec_m;
    double_type value=0.0;
    //for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){
    for (int m=0;m<basis.D_basis.size();m++){
        dec_m = basis.D_basis[m];

        //Using Li tables
        //dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];


        //Sz(i)Sz(j)----------------------------------------------------------------
        value +=conjugate(Vec_[m])*Vec_[m]*((1.0*value_at_pos(dec_m, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin));

    }

    return value;

}

void MODEL_Spins_Target_Sz::Act_SiSj(int &site_i, int &site_j, ulli &m, Mat_1_ullint &m_out_array, Mat_1_doub &Coeff_out_Array, BASIS_Spins_Target_Sz &basis){


    m_out_array.clear();
    Coeff_out_Array.clear();
    double_type value;
    int val_site_i, val_site_j;
    int val_site_i_new, val_site_j_new;
    ulli dec_old, dec_new, dec_new_temp;
    bool allowed;


    dec_old = m;
    val_site_i = value_at_pos(dec_old, site_i, basis.BASE);
    val_site_j = value_at_pos(dec_old, site_j, basis.BASE);


    //Sz(i)Sz(j)----------------------------------------------------------------
    value = ((1.0*value_at_pos(dec_old, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
            ( (1.0*value_at_pos(dec_old, site_j, basis.BASE)) - (0.5*basis.TwoTimesSpin));

    m_out_array.push_back(dec_old);
    Coeff_out_Array.push_back(value);
    //----------------------------------------------------------------------------


    //0.5*S+(i)S-(j)---------------------------------------------------------------
    //if i notequal to j : then site_i cannot be in Spin, and site_j cannot be in -Spin
    if(site_i!=site_j){
        allowed = ((val_site_j != 0) && (val_site_i != (basis.BASE - 1)));
        if(allowed)
        {

            val_site_i_new = val_site_i + 1;
            val_site_j_new = val_site_j - 1;

            dec_new_temp = Updated_decimal_with_value_at_pos(dec_old, site_i, basis.BASE, val_site_i_new);
            dec_new = Updated_decimal_with_value_at_pos(dec_new_temp, site_j, basis.BASE, val_site_j_new);

            value = one*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                              ((val_site_i - (0.5*basis.TwoTimesSpin))*
                               (val_site_i_new - (0.5*basis.TwoTimesSpin))  ) );
            value = value*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                    ((val_site_j - (0.5*basis.TwoTimesSpin))*
                                     (val_site_j_new - (0.5*basis.TwoTimesSpin))  ) );

            m_out_array.push_back(dec_new);
            Coeff_out_Array.push_back(value);

        }
    }
    else{ //i.e. i==j
        //if i equal to j : then j(or i) cannot be in -Spin
        allowed = (val_site_j != 0);
        val_site_j_new = val_site_j - 1;

        value = one*0.5*( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                          ((val_site_j - (0.5*basis.TwoTimesSpin))*
                           (val_site_j_new - (0.5*basis.TwoTimesSpin))  ) );

        m_out_array.push_back(dec_old);
        Coeff_out_Array.push_back(value);

    }
    //------------------------------------------------------------------------------------


    //0.5*S-(i)S+(j)--------------------------------------------------------------------------
    //if i notequal to j : then site_i cannot be in -Spin, and site_j cannot be in Spin
    if(site_i!=site_j){
        allowed = ((val_site_i != 0) && (val_site_j != (basis.BASE - 1)));
        if(allowed)
        {

            val_site_j_new = val_site_j + 1;
            val_site_i_new = val_site_i - 1;

            dec_new_temp = Updated_decimal_with_value_at_pos(dec_old, site_j, basis.BASE, val_site_j_new);
            dec_new = Updated_decimal_with_value_at_pos(dec_new_temp, site_i, basis.BASE, val_site_i_new);

            value = one*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                              ((val_site_i - (0.5*basis.TwoTimesSpin))*
                               (val_site_i_new - (0.5*basis.TwoTimesSpin))  ) );
            value = value*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                    ((val_site_j - (0.5*basis.TwoTimesSpin))*
                                     (val_site_j_new - (0.5*basis.TwoTimesSpin))  ) );

            m_out_array.push_back(dec_new);
            Coeff_out_Array.push_back(value);

        }
    }
    else{ //i.e. i==j
        //if i equal to j : then j(or i) cannot be in Spin
        allowed = (val_site_j != (basis.BASE - 1));
        val_site_j_new = val_site_j + 1;

        value = one*0.5*( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                          ((val_site_j - (0.5*basis.TwoTimesSpin))*
                           (val_site_j_new - (0.5*basis.TwoTimesSpin))  ) );

        m_out_array.push_back(dec_old);
        Coeff_out_Array.push_back(value);

    }
    //------------------------------------------------------------------------------------------


}


void MODEL_Spins_Target_Sz::Initialize_four_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site_i, int site_j, int site_k, int site_l, BASIS_Spins_Target_Sz &basis){


    OPR_.nrows = basis.basis_size;
    OPR_.ncols = OPR_.nrows;

    ulli dec_m;
    Mat_1_ullint m_connected_final;
    Mat_1_doub coeffs_final;
    Mat_1_ullint m_connected_temp2;
    Mat_1_doub coeffs_temp2;
    Mat_1_ullint m_connected_temp3;
    Mat_1_doub coeffs_temp3;
    ulli mi,mj;
    ulli dec_max;
    int m_new;
    int n_index;
    Mat_1_int state_vec;
    Mat_1_ullint m_connected;
    Mat_1_doub coeffs;
    Mat_1_ullint m_connected_temp;
    Mat_1_doub coeffs_temp;


    for (int m=0;m<basis.D_basis.size();m++){

        dec_m = basis.D_basis[m];
        m_connected.clear();
        coeffs.clear();


        if(opr_type=="J2_type"){
            Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
            value_multiply_vector(one*J2_mat[site_i][site_j][site_k][site_l], coeffs_temp);
            for(int i_=0;i_<m_connected_temp.size();i_++){
                mi = m_connected_temp[i_];
                Act_SiSj(site_i, site_j, mi, m_connected_temp2, coeffs_temp2, basis);

                value_multiply_vector(coeffs_temp[i_], coeffs_temp2);
                m_connected.insert(m_connected.end(), m_connected_temp2.begin(), m_connected_temp2.end());
                coeffs.insert(coeffs.end(), coeffs_temp2.begin(), coeffs_temp2.end());
            }
        }



        if(opr_type=="J3_type"){
            Act_SiSj(site_k, site_l, dec_m, m_connected_temp, coeffs_temp, basis);
            value_multiply_vector(one*J3_mat[site_i][site_j][site_k][site_l], coeffs_temp);
            for(int i_=0;i_<m_connected_temp.size();i_++){
                mi = m_connected_temp[i_];
                Act_SiSj(site_k, site_l, mi, m_connected_temp2, coeffs_temp2, basis);
                value_multiply_vector(coeffs_temp[i_], coeffs_temp2);

                for(int j_=0;j_<m_connected_temp2.size();j_++){
                    mj = m_connected_temp2[j_];
                    Act_SiSj(site_i, site_j, mj, m_connected_temp3, coeffs_temp3, basis);
                    value_multiply_vector(coeffs_temp2[j_], coeffs_temp3);

                    m_connected.insert(m_connected.end(), m_connected_temp3.begin(), m_connected_temp3.end());
                    coeffs.insert(coeffs.end(), coeffs_temp3.begin(), coeffs_temp3.end());

                }
            }
        }



        //remove_repetitions
        Remove_repetitions(m_connected, coeffs, m_connected_final, coeffs_final);


        for(int j=0;j<m_connected_final.size();j++){

            fromDeci_to_Vecint(state_vec, basis.BASE,m_connected_final[j] , basis.Length);
            quicksort(state_vec, 0, state_vec.size() -1);
            fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
            n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
            m_new = Find_int_in_part_of_intarray(m_connected_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


            if(m_new<=m){
                OPR_.value.push_back(coeffs_final[j]*one);
                OPR_.rows.push_back(m_new);
                OPR_.columns.push_back(m);
            }
        }
    }







}


void MODEL_Spins_Target_Sz::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site1, int site2, BASIS_Spins_Target_Sz &basis){



    OPR_.nrows = basis.basis_size;
    OPR_.ncols = OPR_.nrows;

    if(opr_type=="Svec.Svec"){

        //-----------------------

        Mat_1_ullint m_connected;
        Mat_1_doub coeffs;
        Mat_1_ullint m_connected_final;
        Mat_1_doub coeffs_final;

        Mat_1_ullint m_connected_temp;
        Mat_1_doub coeffs_temp;

        ulli dec_m, dec_max;
        ulli dec_part1_, dec_part0_;
        int m_new;
        int n_index;
        Mat_1_int state_vec;
        int site_i, site_j;


      //  for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){

        for (int m=0;m<basis.D_basis.size();m++){
            dec_m = basis.D_basis[m];
           // dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];

            m_connected.clear();
            coeffs.clear();

            //TERM  Svec[i].Svec[j] |m>
            site_i=site1;
            site_j=site2;


            Act_SiSj(site_i, site_j, dec_m, m_connected_temp, coeffs_temp, basis);

            value_multiply_vector(one, coeffs_temp);
            m_connected.insert(m_connected.end(), m_connected_temp.begin(), m_connected_temp.end());
            coeffs.insert(coeffs.end(), coeffs_temp.begin(), coeffs_temp.end());


            //remove_repetitions
            Remove_repetitions(m_connected, coeffs, m_connected_final, coeffs_final);


            for(int j=0;j<m_connected_final.size();j++){

                // dec_part1_ = m_connected_final[j]/pow(basis.BASE, basis.Partition_Length[0]);
                // dec_part0_ = m_connected_final[j] - (dec_part1_*pow(basis.BASE, basis.Partition_Length[0]));
                // m_new = basis.Dec_to_Index_part0_[dec_part0_] + basis.Dec_to_Index_part1_[dec_part1_];

                //m_new = Find_int_in_intarray_using_bisection(m_connected_final[j], basis.D_basis);
                //or
                // m_new = Find_int_in_intarray_using_multisections(m_connected_final[j], basis.D_basis, MultisectionSearch_int);

                //or
                // m_new = Find_int_in_intarray(m_connected_final[j], basis.D_basis);

                //or

                 fromDeci_to_Vecint(state_vec, basis.BASE,m_connected_final[j] , basis.Length);
                 quicksort(state_vec, 0, state_vec.size() -1);
                 fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                 n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                 m_new = Find_int_in_part_of_intarray(m_connected_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                OPR_.value.push_back(coeffs_final[j]*one);
                OPR_.rows.push_back(m_new);
                OPR_.columns.push_back(m);

            }

        } // m i.e. basis index, columns of H


        //----------------------


    }

}

void MODEL_Spins_Target_Sz::Read_parameters(BASIS_Spins_Target_Sz &basis, string filename){


    //    bool read_basis, write_basis;
    //    string read_basis_file, write_basis_file;

    string filepath = filename;

    string read_basis_, Read_Basis_= "Read_Basis = ";
    string write_basis_, Write_Basis_= "Write_Basis = ";
    string Read_Basis_File_= "Read_Basis_File = ";
    string Write_Basis_File_= "Write_Basis_File = ";


    string multisectionsearchint_, MultiSectionSearchInt_ = "MultiSectionBasisSearchInt = ";
    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string twotimesspin, TwoTimesSpin = "TwoTimesSpin = ";

    string twotimestotalsztarget, TwoTimesTotalSzTarget = "TwoTimesTotalSzTarget = ";

    string hmag, Hmag = "H_mag = ";
    string d_anisotropy_, D_Anisotropy_ = "Dz_anisotropy = ";

    string LongRangeJ1file_ = "LongRangeJ1_file = ";

    string LongRangeJ2file_ = "LongRangeJ2_file = ";

    string LongRangeJ3file_ = "LongRangeJ3_file = ";

    string genericconnectionsfiles_, GenericConnectionsFiles_ = "GenericConnectionsFiles = ";

    string fourpointobservablessitesfile_ ,FourPointObservablesSitesFile_ = "FourPointObservablesSites file = ";

    //string no_of_onepoint_obs_, No_Of_Onepoint_Obs_ = "No_of_onepoint_obs = ";

    int offset;
    string line;
    ifstream inputfile(filepath.c_str());

    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(GenericConnectionsFiles_, 0)) != string::npos) {
               genericconnectionsfiles_ = line.substr (offset+GenericConnectionsFiles_.length());  }

            // if ((offset = line.find(No_Of_Onepoint_Obs_, 0)) != string::npos) {
            //     no_of_onepoint_obs_ = line.substr (offset+No_Of_Onepoint_Obs_.length());  }

            if ((offset = line.find(FourPointObservablesSitesFile_, 0)) != string::npos) {
                fourpointobservablessitesfile_ = line.substr (offset+FourPointObservablesSitesFile_.length());  }

            if ((offset = line.find(LongRangeJ1file_, 0)) != string::npos) {
                LongRangeJ1_filepath = line.substr (offset+LongRangeJ1file_.length());  }

            if ((offset = line.find(LongRangeJ2file_, 0)) != string::npos) {
                LongRangeJ2_filepath = line.substr (offset+LongRangeJ2file_.length());  }

            if ((offset = line.find(LongRangeJ3file_, 0)) != string::npos) {
                LongRangeJ3_filepath = line.substr (offset+LongRangeJ3file_.length());  }

            if ((offset = line.find(PBC_, 0)) != string::npos) {
                pbc_ = line.substr (offset+PBC_.length());				}

            if ((offset = line.find(Read_Basis_, 0)) != string::npos) {
                read_basis_ = line.substr (offset+Read_Basis_.length());				}

            if ((offset = line.find(Write_Basis_, 0)) != string::npos) {
                write_basis_ = line.substr (offset+Write_Basis_.length());				}

            if ((offset = line.find(Read_Basis_File_, 0)) != string::npos) {
                basis.read_basis_file = line.substr (offset+Read_Basis_File_.length());				}

            if ((offset = line.find(Write_Basis_File_, 0)) != string::npos) {
                basis.write_basis_file = line.substr (offset+Write_Basis_File_.length());				}


            if ((offset = line.find(MultiSectionSearchInt_, 0)) != string::npos) {
                multisectionsearchint_ = line.substr (offset + MultiSectionSearchInt_.length());		}

            if ((offset = line.find(Length, 0)) != string::npos) {
                length = line.substr (offset + Length.length());		}

            if ((offset = line.find(TwoTimesTotalSzTarget, 0)) != string::npos) {
                twotimestotalsztarget = line.substr (offset + TwoTimesTotalSzTarget.length());		}

            if ((offset = line.find(TwoTimesSpin, 0)) != string::npos) {
                twotimesspin = line.substr (offset + TwoTimesSpin.length());		}

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(D_Anisotropy_, 0)) != string::npos) {
                d_anisotropy_ = line.substr (offset + D_Anisotropy_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    if(pbc_ == "true"){
        PBC =true;
    }
    else{
        PBC=false;
    }

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
    basis.Target_Total_Sz = 0.5*(atof(twotimestotalsztarget.c_str()));
    basis.TwoTimesSpin=atoi(twotimesspin.c_str());

    MultisectionSearch_int = atoi(multisectionsearchint_.c_str());

    //No_of_onepoint_obs=atoi(no_of_onepoint_obs_.c_str());


    stringstream genericconnectionsfiles_stream;
    genericconnectionsfiles_stream<<genericconnectionsfiles_;
    genericconnectionsfiles_stream>>N_ConnectionsFiles;
    ConnectionFiles.resize(N_ConnectionsFiles);
    string filename_temp;
    for(int i=0;i<N_ConnectionsFiles;i++){
        genericconnectionsfiles_stream>>filename_temp;
        ConnectionFiles[i]=filename_temp;
    }


    Connections.resize(N_ConnectionsFiles);
    for(int FileNo=0;FileNo<N_ConnectionsFiles;FileNo++){
    string line_connection;
    ifstream inputfileConnection(ConnectionFiles[FileNo].c_str());
    Connections[FileNo].clear();
    while(getline(inputfileConnection,line_connection)){
    Connections[FileNo].push_back(line_connection);
    }
    }


    double temp_h;
    stringstream hmag_stream;
    hmag_stream<<hmag;
    H_field.clear();
    H_field.resize(3); //Hx, Hy, Hz
    for(int i=0;i<3;i++){
        hmag_stream >> temp_h;
        H_field[i]=temp_h;
    }


    double temp_anis;
    stringstream anis_stream;
    anis_stream<<d_anisotropy_;
    int terms_anis, site_temp;
    double value_temp_;
    anis_stream>>terms_anis;
    Dz_anisotropy.clear();
    Dz_anisotropy.resize(basis.Length); //Dz
    for(int i=0;i<terms_anis;i++){
        anis_stream >> site_temp>>value_temp_;
        Dz_anisotropy[site_temp]=value_temp_;
    }



    J1_mat.resize(basis.Length);
    J2_mat.resize(basis.Length);
    J3_mat.resize(basis.Length);
    for(int i =0;i<basis.Length;i++){
        J1_mat[i].resize(basis.Length);
        J2_mat[i].resize(basis.Length);
        J3_mat[i].resize(basis.Length);
        for(int j =0;j<basis.Length;j++){
            J1_mat[i][j]=0.0;
            J2_mat[i][j].resize(basis.Length);
            J3_mat[i][j].resize(basis.Length);
            for(int k =0;k<basis.Length;k++){
                J2_mat[i][j][k].resize(basis.Length);
                J3_mat[i][j][k].resize(basis.Length);
                for(int l =0;l<basis.Length;l++){
                    J2_mat[i][j][k][l]=0.0;
                    J3_mat[i][j][k][l]=0.0;
                }
            }
        }
    }


    ifstream inputJ1(LongRangeJ1_filepath.c_str());
    ifstream inputJ2(LongRangeJ2_filepath.c_str());
    ifstream inputJ3(LongRangeJ3_filepath.c_str());

    string line_temp;
    int i1, i2, i3, i4;
    double value_temp;

    while(getline(inputJ1, line_temp)){
        stringstream line_temp_ss(line_temp);
        line_temp_ss >> i1 >> i2 >> value_temp;
        cout<<i1<<" "<<i2<<" "<<value_temp<<endl;
        J1_mat[i1][i2]=value_temp;
    }

    while(getline(inputJ2, line_temp)){
        stringstream line_temp_ss(line_temp);
        line_temp_ss >> i1 >> i2 >> i3>>i4>> value_temp;
        J2_mat[i1][i2][i3][i4]=value_temp;
    }

    while(getline(inputJ3, line_temp)){
        stringstream line_temp_ss(line_temp);
        line_temp_ss >> i1 >> i2 >> i3>>i4>> value_temp;
        J3_mat[i1][i2][i3][i4]=value_temp;
    }





    /*
    One_point_strs.resize(No_of_onepoint_obs);
    string One_point_oprt_path = "One_point_oprt_path_";
    string one_point_oprt_path_;


    ifstream inputfile2(filepath.c_str());

    int i_temp;

    if(inputfile2.is_open())
    {
        while(!inputfile2.eof())
        {
            getline(inputfile2,line);

            for(int i=0;i<No_of_onepoint_obs;i++){

                stringstream ss_temp;
                string str_temp;

                i_temp = i+1;
                ss_temp << i_temp;
                ss_temp >> str_temp;
                one_point_oprt_path_ = One_point_oprt_path + str_temp + " = ";

                if ((offset = line.find(one_point_oprt_path_, 0)) != string::npos) {
                    One_point_strs[i] = line.substr (offset+one_point_oprt_path_.length());  }

            }


        }
        inputfile2.close();
    }
    else
    {cout<<"Unable to open input file 2nd time while in the Model class."<<endl;}



    ifstream input4site(fourpointobservablessitesfile_.c_str());
    Four_point_sites_set.clear();
    tetra_int Tetra_int_temp;
    double double_temp;
    while (getline(input4site, line))
    {
        stringstream ss_temp;
        ss_temp<<line;
        ss_temp>>double_temp;
        Tetra_int_temp.first=int(double_temp);
        ss_temp>>double_temp;
        Tetra_int_temp.second=int(double_temp);
        ss_temp>>double_temp;
        Tetra_int_temp.third=int(double_temp);
        ss_temp>>double_temp;
        Tetra_int_temp.fourth=int(double_temp);


        Four_point_sites_set.push_back(Tetra_int_temp);
        // process pair (a,b)
    }

    */


}


void MODEL_Spins_Target_Sz::CreateLocalOprs_in_LocalHilbertBasis(BASIS_Spins_Target_Sz &basis){

    assert(basis.TwoTimesSpin>0);
    int HS; //Hilbert Space Size
    HS=basis.TwoTimesSpin+1;
    double Spin=0.5*basis.TwoTimesSpin;
    //basis convention = [-S, -S+1,..., S-1, S]

    //Dipole Oprs
    SxLocal.resize(HS,HS);SyLocal.resize(HS,HS);SzLocal.resize(HS,HS);
    SplusLocal.resize(HS,HS);SminusLocal.resize(HS,HS);SyLocalTimesIota.resize(HS,HS);
    int mp_temp;
    for(int m=0;m<HS;m++){
        SzLocal(m,m) = -Spin + 1.0*m;
        mp_temp=m+1;

        if(mp_temp<HS){
        SplusLocal(mp_temp,m) = sqrt( Spin*(Spin+1.0) - ( (-Spin+mp_temp*1.0)*(-Spin+m*1.0) ) );
        SminusLocal(m,mp_temp) = SplusLocal(mp_temp,m);
        }
    }

    for(int m=0;m<HS;m++){
        for(int mp=0;mp<HS;mp++){
            SxLocal(mp,m) = 0.5*(SplusLocal(mp,m) + SminusLocal(mp,m));
            SyLocalTimesIota(mp,m) = 0.5*(SplusLocal(mp,m) - SminusLocal(mp,m));
        }
    }

#ifdef USE_COMPLEX
    for(int m=0;m<HS;m++){
        for(int mp=0;mp<HS;mp++){
        SyLocal(mp,m) = -0.5*iota_comp*(SplusLocal(mp,m) - SminusLocal(mp,m));
        }
    }
#endif



    //Quadropolar Oprs
    if(basis.TwoTimesSpin>=2){
    QxxLocal.resize(HS,HS);QxyLocal.resize(HS,HS);QxzLocal.resize(HS,HS);
    QyyLocal.resize(HS,HS);QyzLocal.resize(HS,HS);QzzLocal.resize(HS,HS);
    QxyLocalTimesIota.resize(HS,HS);QyzLocalTimesIota.resize(HS,HS);
    Sx2Local.resize(HS,HS);Sy2Local.resize(HS,HS);Sz2Local.resize(HS,HS);



    for(int m=0;m<HS;m++){
        for(int mp=0;mp<HS;mp++){
        for(int k=0;k<HS;k++){
                Sx2Local(mp,m) += SxLocal(mp,k)*SxLocal(k,m);
                Sy2Local(mp,m) += -1.0*SyLocalTimesIota(mp,k)*SyLocalTimesIota(k,m);
                Sz2Local(mp,m) += SzLocal(mp,k)*SzLocal(k,m);
                QxzLocal(mp,m) += SxLocal(mp,k)*SzLocal(k,m) + SzLocal(mp,k)*SxLocal(k,m);
                QyzLocalTimesIota(mp,m) += SyLocalTimesIota(mp,k)*SzLocal(k,m) +
                                           SzLocal(mp,k)*SyLocalTimesIota(k,m);
                QxyLocalTimesIota(mp,m) += SxLocal(mp,k)*SyLocalTimesIota(k,m) +
                                           SyLocalTimesIota(mp,k)*SxLocal(k,m);
            }
    }}


    }

}



void MODEL_Spins_Target_Sz::Read_parameters_for_dynamics(string filename){


    string DYN_STRING, Dyn_opr_string_  = "Opr_for_Dynamics = ";
    int No_of_oprts;


    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(Dyn_opr_string_, 0)) != string::npos) {
                DYN_STRING = line.substr (offset + Dyn_opr_string_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    stringstream DYN_SSTREAM(DYN_STRING);

    DYN_SSTREAM >> No_of_oprts;
    DYN_SSTREAM >> Dyn_opr_string;
    Dyn_opr_int.resize(No_of_oprts);
    Dyn_opr_coeffs.resize(No_of_oprts);

    for(int opr_no=0;opr_no<No_of_oprts;opr_no++){
        DYN_SSTREAM >>  Dyn_opr_coeffs[opr_no];
        DYN_SSTREAM >> Dyn_opr_int[opr_no];
    }


}


void MODEL_Spins_Target_Sz::Direct_product_of_localJ_Ansatz(Mat_1_pair_realInt JJzBasis_state, Mat_1_real &Coefficients, Mat_1_int &m_basis, BASIS_Spins_Target_Sz &basis){

    cout<<"THIS ROUTINE IS STRICTLY USED ONLY FOR S-L MODEL"<<endl;
    Coefficients.clear();
    m_basis.clear();
    Mat_2_int SzLzBasis_state;

    ulli dec_, dec_max;
    int m_, n_index;
    Mat_2_real Coefficients_local;
    Mat_2_pair_realint SzLzBasis_state_local;
    Coefficients_local.resize(JJzBasis_state.size());
    SzLzBasis_state_local.resize(JJzBasis_state.size());

    assert((basis.Length/2) == JJzBasis_state.size());
    for(int local_state=0;local_state<JJzBasis_state.size();local_state++){

        if( JJzBasis_state[local_state].first==0 ){  //|J=0,Jz=0>
            assert(JJzBasis_state[local_state].second==0);
            Coefficients_local[local_state].resize(3);
            SzLzBasis_state_local[local_state].resize(3);
            Coefficients_local[local_state][0]=1.0/sqrt(3.0);
            Coefficients_local[local_state][1]=-1.0/sqrt(3.0);
            Coefficients_local[local_state][2]=-1.0/sqrt(3.0);
            SzLzBasis_state_local[local_state][0].first=0;SzLzBasis_state_local[local_state][0].second=0;
            SzLzBasis_state_local[local_state][1].first=-1;SzLzBasis_state_local[local_state][1].second=1;
            SzLzBasis_state_local[local_state][2].first=1;SzLzBasis_state_local[local_state][2].second=-1;
        }
        if( JJzBasis_state[local_state].first==1 && JJzBasis_state[local_state].second==0){  //|J=1, Jz=0>
            Coefficients_local[local_state].resize(2);
            SzLzBasis_state_local[local_state].resize(2);
            Coefficients_local[local_state][0]=1.0/sqrt(2.0);
            Coefficients_local[local_state][1]=-1.0/sqrt(2.0);
            SzLzBasis_state_local[local_state][0].first=-1;SzLzBasis_state_local[local_state][0].second=1;
            SzLzBasis_state_local[local_state][1].first=1;SzLzBasis_state_local[local_state][1].second=-1;
        }
        if( JJzBasis_state[local_state].first==1 && JJzBasis_state[local_state].second!=0){  //|J=1, Jz=+-1>
            assert(abs(JJzBasis_state[local_state].second)==1);
            Coefficients_local[local_state].resize(2);
            SzLzBasis_state_local[local_state].resize(2);
            Coefficients_local[local_state][0]=1.0/sqrt(2.0);
            Coefficients_local[local_state][1]=-1.0/sqrt(2.0);
            SzLzBasis_state_local[local_state][0].first=0;SzLzBasis_state_local[local_state][0].second=JJzBasis_state[local_state].second;
            SzLzBasis_state_local[local_state][1].first=JJzBasis_state[local_state].second;SzLzBasis_state_local[local_state][1].second=0;
        }
        if( JJzBasis_state[local_state].first==2 && JJzBasis_state[local_state].second==0){  //|J=2,Jz=0>
            Coefficients_local[local_state].resize(3);
            SzLzBasis_state_local[local_state].resize(3);
            Coefficients_local[local_state][0]=sqrt(2.0/3.0);
            Coefficients_local[local_state][1]=1.0/sqrt(6.0);
            Coefficients_local[local_state][2]=1.0/sqrt(6.0);
            SzLzBasis_state_local[local_state][0].first=0;SzLzBasis_state_local[local_state][0].second=0;
            SzLzBasis_state_local[local_state][1].first=-1;SzLzBasis_state_local[local_state][1].second=1;
            SzLzBasis_state_local[local_state][2].first=1;SzLzBasis_state_local[local_state][2].second=-1;
        }
        if( JJzBasis_state[local_state].first==2 && abs(JJzBasis_state[local_state].second)==1){  //|J=2, Jz=+-1>
            Coefficients_local[local_state].resize(2);
            SzLzBasis_state_local[local_state].resize(2);
            Coefficients_local[local_state][0]=1.0/sqrt(2.0);
            Coefficients_local[local_state][1]=1.0/sqrt(2.0);
            SzLzBasis_state_local[local_state][0].first=0;SzLzBasis_state_local[local_state][0].second=JJzBasis_state[local_state].second;
            SzLzBasis_state_local[local_state][1].first=JJzBasis_state[local_state].second;SzLzBasis_state_local[local_state][1].second=0;
        }
        if( JJzBasis_state[local_state].first==2 && abs(JJzBasis_state[local_state].second)==2){  //|J=2, Jz=+-2>
            Coefficients_local[local_state].resize(1);
            SzLzBasis_state_local[local_state].resize(1);
            Coefficients_local[local_state][0]=1.0;
            SzLzBasis_state_local[local_state][0].first=JJzBasis_state[local_state].second/abs(JJzBasis_state[local_state].second);
            SzLzBasis_state_local[local_state][0].second=JJzBasis_state[local_state].second/abs(JJzBasis_state[local_state].second);
        }

    }



    if(SzLzBasis_state_local.size()==2){
        for(int b0_=0;b0_<Coefficients_local[0].size();b0_++){

            for(int b1_=0;b1_<Coefficients_local[1].size();b1_++){

                Mat_1_int temp_szlz_state;
                double val_temp;
                val_temp = Coefficients_local[0][b0_]*Coefficients_local[1][b1_];
                Coefficients.push_back(val_temp);

                temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].first); //sz
                temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].second); //lz
                temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].first);
                temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].second);

                for(int i_=0;i_<temp_szlz_state.size();i_++){
                    temp_szlz_state[i_] = temp_szlz_state[i_] + basis.SPIN;
                }

                fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_, basis.Length);
                quicksort(temp_szlz_state, 0, temp_szlz_state.size() -1);
                fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_max, basis.Length);
                n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                m_ = Find_int_in_part_of_intarray(dec_, basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                m_basis.push_back(m_);

            }
        }
    }


    if(SzLzBasis_state_local.size()==4){
        for(int b0_=0;b0_<Coefficients_local[0].size();b0_++){

            for(int b1_=0;b1_<Coefficients_local[1].size();b1_++){

                for(int b2_=0;b2_<Coefficients_local[2].size();b2_++){

                    for(int b3_=0;b3_<Coefficients_local[3].size();b3_++){
                        Mat_1_int temp_szlz_state;
                        double val_temp;
                        val_temp = Coefficients_local[0][b0_]*Coefficients_local[1][b1_]*Coefficients_local[2][b2_]*Coefficients_local[3][b3_];
                        Coefficients.push_back(val_temp);

                        temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].first); //sz
                        temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].second); //lz
                        temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].first);
                        temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].second);
                        temp_szlz_state.push_back(SzLzBasis_state_local[2][b2_].first);
                        temp_szlz_state.push_back(SzLzBasis_state_local[2][b2_].second);
                        temp_szlz_state.push_back(SzLzBasis_state_local[3][b3_].first);
                        temp_szlz_state.push_back(SzLzBasis_state_local[3][b3_].second);

                        for(int i_=0;i_<temp_szlz_state.size();i_++){
                            temp_szlz_state[i_] = temp_szlz_state[i_] + basis.SPIN;
                        }

                        fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_, basis.Length);
                        quicksort(temp_szlz_state, 0, temp_szlz_state.size() -1);
                        fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_max, basis.Length);
                        n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                        m_ = Find_int_in_part_of_intarray(dec_, basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                        m_basis.push_back(m_);


                    }
                }
            }
        }
    }


    if(SzLzBasis_state_local.size()==6){
        for(int b0_=0;b0_<Coefficients_local[0].size();b0_++){

            for(int b1_=0;b1_<Coefficients_local[1].size();b1_++){

                for(int b2_=0;b2_<Coefficients_local[2].size();b2_++){

                    for(int b3_=0;b3_<Coefficients_local[3].size();b3_++){

                        for(int b4_=0;b4_<Coefficients_local[4].size();b4_++){

                            for(int b5_=0;b5_<Coefficients_local[5].size();b5_++){
                                Mat_1_int temp_szlz_state;
                                double val_temp;
                                val_temp = Coefficients_local[0][b0_]*Coefficients_local[1][b1_]*Coefficients_local[2][b2_]*Coefficients_local[3][b3_]*Coefficients_local[4][b4_]*Coefficients_local[5][b5_];
                                Coefficients.push_back(val_temp);

                                temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].first); //sz
                                temp_szlz_state.push_back(SzLzBasis_state_local[0][b0_].second); //lz
                                temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].first);
                                temp_szlz_state.push_back(SzLzBasis_state_local[1][b1_].second);
                                temp_szlz_state.push_back(SzLzBasis_state_local[2][b2_].first);
                                temp_szlz_state.push_back(SzLzBasis_state_local[2][b2_].second);
                                temp_szlz_state.push_back(SzLzBasis_state_local[3][b3_].first);
                                temp_szlz_state.push_back(SzLzBasis_state_local[3][b3_].second);
                                temp_szlz_state.push_back(SzLzBasis_state_local[4][b4_].first);
                                temp_szlz_state.push_back(SzLzBasis_state_local[4][b4_].second);
                                temp_szlz_state.push_back(SzLzBasis_state_local[5][b5_].first);
                                temp_szlz_state.push_back(SzLzBasis_state_local[5][b5_].second);

                                for(int i_=0;i_<temp_szlz_state.size();i_++){
                                    temp_szlz_state[i_] = temp_szlz_state[i_] + basis.SPIN;
                                }

                                fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_, basis.Length);
                                quicksort(temp_szlz_state, 0, temp_szlz_state.size() -1);
                                fromVecint_to_Deci(temp_szlz_state, basis.BASE, dec_max, basis.Length);
                                n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                                m_ = Find_int_in_part_of_intarray(dec_, basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                                m_basis.push_back(m_);


                            }
                        }
                    }
                }
            }
        }
    }



}

void MODEL_Spins_Target_Sz::Create_JJz_Trial_States(Mat_2_pair_realint &JJzBasis_states, Mat_2_real &Coefficients, Mat_2_int &m_basis, BASIS_Spins_Target_Sz &basis){

    cout<<"THIS ROUTINE IS STRICTLY USED ONLY FOR S-L MODEL"<<endl;

    JJzBasis_states.clear();
    Coefficients.clear();
    m_basis.clear();

    int SL_sites;
    SL_sites=basis.Length/2;

    Mat_1_pair_realInt JJz_basis_temp;
    pair_real_int temp_pair_int;

    temp_pair_int.first=2;
    for(int temp_i=-2;temp_i<=2;temp_i++){
        temp_pair_int.second=temp_i;
        JJz_basis_temp.push_back(temp_pair_int);
    }

    temp_pair_int.first=1;
    for(int temp_i=-1;temp_i<=1;temp_i++){
        temp_pair_int.second=temp_i;
        JJz_basis_temp.push_back(temp_pair_int);
    }

    temp_pair_int.first=0;
    temp_pair_int.second=0;
    JJz_basis_temp.push_back(temp_pair_int);







    if(SL_sites==2){

        for(int J0_=0;J0_<JJz_basis_temp.size();J0_++){ //site0
            for(int J1_=0;J1_<JJz_basis_temp.size();J1_++){ //site1

                Mat_1_pair_realInt JJz_basis_state_var;
                Mat_1_real Coefficients_var;
                Mat_1_int m_basis_var;
                int total_Jz;

                JJz_basis_state_var.push_back(JJz_basis_temp[J0_]); //site0;
                JJz_basis_state_var.push_back(JJz_basis_temp[J1_]); //site1;

                total_Jz =JJz_basis_temp[J0_].second + JJz_basis_temp[J1_].second;
                total_Jz = int(total_Jz + (basis.SPIN*basis.Length));

                if(total_Jz == basis.Target_Total_Value){
                    Direct_product_of_localJ_Ansatz(JJz_basis_state_var, Coefficients_var, m_basis_var, basis);
                    JJzBasis_states.push_back(JJz_basis_state_var);
                    Coefficients.push_back(Coefficients_var);
                    m_basis.push_back(m_basis_var);
                }

            }
        }
    }


    if(SL_sites==4){

        for(int J0_=0;J0_<JJz_basis_temp.size();J0_++){ //site0
            for(int J1_=0;J1_<JJz_basis_temp.size();J1_++){ //site1

                for(int J2_=0;J2_<JJz_basis_temp.size();J2_++){ //site2
                    for(int J3_=0;J3_<JJz_basis_temp.size();J3_++){ //site3

                        Mat_1_pair_realInt JJz_basis_state_var;
                        Mat_1_real Coefficients_var;
                        Mat_1_int m_basis_var;
                        int total_Jz;

                        JJz_basis_state_var.push_back(JJz_basis_temp[J0_]); //site0;
                        JJz_basis_state_var.push_back(JJz_basis_temp[J1_]); //site1;
                        JJz_basis_state_var.push_back(JJz_basis_temp[J2_]); //site2;
                        JJz_basis_state_var.push_back(JJz_basis_temp[J3_]); //site3;

                        total_Jz =JJz_basis_temp[J0_].second + JJz_basis_temp[J1_].second + JJz_basis_temp[J2_].second + JJz_basis_temp[J3_].second;
                        total_Jz = int(total_Jz + (basis.SPIN*basis.Length));

                        if(total_Jz == basis.Target_Total_Value){
                            Direct_product_of_localJ_Ansatz(JJz_basis_state_var, Coefficients_var, m_basis_var, basis);

                            JJzBasis_states.push_back(JJz_basis_state_var);
                            Coefficients.push_back(Coefficients_var);
                            m_basis.push_back(m_basis_var);
                        }

                    }
                }
            }
        }
    }

    if(SL_sites==6){

        for(int J0_=0;J0_<JJz_basis_temp.size();J0_++){ //site0
            for(int J1_=0;J1_<JJz_basis_temp.size();J1_++){ //site1

                for(int J2_=0;J2_<JJz_basis_temp.size();J2_++){ //site2
                    for(int J3_=0;J3_<JJz_basis_temp.size();J3_++){ //site3

                        for(int J4_=0;J4_<JJz_basis_temp.size();J4_++){ //site4
                            for(int J5_=0;J5_<JJz_basis_temp.size();J5_++){ //site5

                                Mat_1_pair_realInt JJz_basis_state_var;
                                Mat_1_real Coefficients_var;
                                Mat_1_int m_basis_var;
                                int total_Jz;

                                JJz_basis_state_var.push_back(JJz_basis_temp[J0_]); //site0;
                                JJz_basis_state_var.push_back(JJz_basis_temp[J1_]); //site1;
                                JJz_basis_state_var.push_back(JJz_basis_temp[J2_]); //site2;
                                JJz_basis_state_var.push_back(JJz_basis_temp[J3_]); //site3;
                                JJz_basis_state_var.push_back(JJz_basis_temp[J4_]); //site4;
                                JJz_basis_state_var.push_back(JJz_basis_temp[J5_]); //site5;

                                total_Jz =JJz_basis_temp[J0_].second + JJz_basis_temp[J1_].second + JJz_basis_temp[J2_].second + JJz_basis_temp[J3_].second + JJz_basis_temp[J4_].second + JJz_basis_temp[J5_].second;
                                total_Jz = int(total_Jz + (basis.SPIN*basis.Length));

                                if(total_Jz == basis.Target_Total_Value){
                                    Direct_product_of_localJ_Ansatz(JJz_basis_state_var, Coefficients_var, m_basis_var, basis);

                                    JJzBasis_states.push_back(JJz_basis_state_var);
                                    Coefficients.push_back(Coefficients_var);
                                    m_basis.push_back(m_basis_var);
                                }

                            }
                        }
                    }
                }
            }
        }
    }

}


void MODEL_Spins_Target_Sz::Overlap_of_JJzBasis_with_State(Mat_1_doub &Vec_ , Mat_1_doub &Overlaps_ , Mat_1_int &sorted_indices,  Mat_2_pair_realint &JJzBasis_states, Mat_2_real &Coefficients, Mat_2_int &m_basis, BASIS_Spins_Target_Sz &basis){


    Mat_1_doub Overlaps_old;

    Overlaps_old.resize(JJzBasis_states.size());
    Overlaps_.resize(JJzBasis_states.size());

    for(int ansatz_basis=0;ansatz_basis<JJzBasis_states.size();ansatz_basis++){

        double_type temp_doub =0.0;
        for(int i_=0;i_<m_basis[ansatz_basis].size();i_++){
            temp_doub += conjugate(Coefficients[ansatz_basis][i_])*Vec_[m_basis[ansatz_basis][i_]];
        }

        Overlaps_old[ansatz_basis]=temp_doub;
    }

    Sort_vector_in_decreasing_order_in_file(Overlaps_old, Overlaps_ , sorted_indices);


}

void MODEL_Spins_Target_Sz::Initialize_Opr_for_Dynamics(BASIS_Spins_Target_Sz &basis){

    Dyn_opr.value.clear();
    Dyn_opr.rows.clear();
    Dyn_opr.columns.clear();
    Dyn_opr.nrows = basis.MainIndex_to_Dec_part0_.size();//basis.D_basis.size();
    Dyn_opr.ncols = basis.MainIndex_to_Dec_part0_.size();


    double_type value_;
    int dec_;

    if(Dyn_opr_string=="Sz"){
        for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){
            dec_ = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];

            //dec_ = basis.D_basis[m];
            value_ =0.0;
            for(int opr_no=0;opr_no<Dyn_opr_int.size();opr_no++){
                value_ += Dyn_opr_coeffs[opr_no]*(
                            ((1.0*value_at_pos(dec_, Dyn_opr_int[opr_no], basis.BASE)) - (0.5*basis.TwoTimesSpin))
                            );

            }
            value_ = value_*(1.0/(sqrt(1.0*Dyn_opr_int.size())));
            Dyn_opr.value.push_back(value_);
            Dyn_opr.rows.push_back(m);
            Dyn_opr.columns.push_back(m);

        }
    }


    //#ifdef USE_COMPLEX
    //            value2=exp(iota_comp*(1.0*(site+1))*Dyn_Momentum*PI)*sqrt(1.0/(basis.Length));
    //#endif
    //#ifndef USE_COMPLEX
    //            cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
    //#endif

    //    for(int site_=0;site_<basis.Length;site_++){
    //        vector< int >().swap( Oprs_local[site_].columns );
    //        vector< int >().swap( Oprs_local[site_].rows );
    //        vector< double_type >().swap( Oprs_local[site_].value );
    //    }



}


void MODEL_Spins_Target_Sz::Initialize_State_for_Dynamics(BASIS_Spins_Target_Sz &basis_new, BASIS_Spins_Target_Sz &basis, Mat_1_doub &Vec_, Mat_1_doub &Vec_Dyn){



    double_type value_;
    int dec_, dec_new, dec_max, m_new;
    Vec_Dyn.resize(basis_new.D_basis.size());
    int site_;
    bool allowed;
    double val_site_, val_site_new;
    int n_index;
    Mat_1_int state_vec;

    if(Dyn_opr_string=="Splus"){
        for (int m=0;m<basis.D_basis.size();m++){
            dec_ = basis.D_basis[m];


            for(int opr_no=0;opr_no<Dyn_opr_int.size();opr_no++){
                site_=Dyn_opr_int[opr_no];
                val_site_ = value_at_pos(dec_, site_, basis.BASE);
                allowed=(val_site_ != (basis.BASE - 1));

                if(allowed){
                    val_site_new = val_site_ + 1;
                    dec_new = Updated_decimal_with_value_at_pos(dec_, site_, basis.BASE, val_site_new);
                    value_ = Vec_[m]*Dyn_opr_coeffs[opr_no]*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                                                  ((val_site_ - (0.5*basis.TwoTimesSpin))*
                                                                   (val_site_new - (0.5*basis.TwoTimesSpin))  ) );


                    fromDeci_to_Vecint(state_vec, basis.BASE, dec_new , basis.Length);
                    quicksort(state_vec, 0, state_vec.size() -1);
                    fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                    n_index = Find_int_in_intarray(dec_max, basis_new.Partitions_Dec);
                    m_new = Find_int_in_part_of_intarray(dec_new, basis_new.D_basis, basis_new.Partitions_pos[n_index].first, basis_new.Partitions_pos[n_index].second);

                    assert(m_new<= basis_new.D_basis.size());
                    Vec_Dyn[m_new] +=value_;
                }

            }
        }
    }

    if(Dyn_opr_string=="Sminus"){
        for (int m=0;m<basis.D_basis.size();m++){
            dec_ = basis.D_basis[m];


            for(int opr_no=0;opr_no<Dyn_opr_int.size();opr_no++){
                site_=Dyn_opr_int[opr_no];
                val_site_ = value_at_pos(dec_, site_, basis.BASE);
                allowed=(val_site_ != 0);

                if(allowed){
                    val_site_new = val_site_ - 1;
                    dec_new = Updated_decimal_with_value_at_pos(dec_, site_, basis.BASE, val_site_new);
                    value_ = Vec_[m]*Dyn_opr_coeffs[opr_no]*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                                                  ((val_site_ - (0.5*basis.TwoTimesSpin))*
                                                                   (val_site_new - (0.5*basis.TwoTimesSpin))  ) );


                    fromDeci_to_Vecint(state_vec, basis.BASE, dec_new , basis.Length);
                    quicksort(state_vec, 0, state_vec.size() -1);
                    fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                    n_index = Find_int_in_intarray(dec_max, basis_new.Partitions_Dec);
                    m_new = Find_int_in_part_of_intarray(dec_new, basis_new.D_basis, basis_new.Partitions_pos[n_index].first, basis_new.Partitions_pos[n_index].second);

                    assert(m_new<= basis_new.D_basis.size());
                    Vec_Dyn[m_new] +=value_;
                }

            }
        }
    }

}

void MODEL_Spins_Target_Sz::Initialize_one_point_to_calculate_from_file(BASIS_Spins_Target_Sz &basis){

    One_point_oprts_onsite.resize(No_of_onepoint_obs);

    for(int i=0;i<No_of_onepoint_obs;i++){
        Read_matrix_from_file(One_point_strs[i],
                              One_point_oprts_onsite[i],basis.BASE, basis.BASE);
    }


    int T_no_oprs=No_of_onepoint_obs;

    One_point_oprts.resize(T_no_oprs);
    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }


    one_point_obs=One_point_strs;
    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.basis_size;
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }


        //Remember OPR[l][m]=<l|OPR|m>
        int j;
        double_type value;


        for(int site=0;site<basis.Length;site++){

            One_point_oprts[opr_no][site].value.clear();
            One_point_oprts[opr_no][site].rows.clear();
            One_point_oprts[opr_no][site].columns.clear();

            for(int col_local=0;col_local<basis.BASE;col_local++){
                for(int row_local=0;row_local<basis.BASE;row_local++){

                    if(One_point_oprts_onsite[opr_no][row_local][col_local] != zero){


                        if(row_local !=col_local){
                            for (int i=0;i<basis.basis_size;i++){
                                if(value_at_pos(i, site, basis.BASE) == col_local){
                                    j = Updated_decimal_with_value_at_pos(i, site, basis.BASE, row_local);
                                    value = One_point_oprts_onsite[opr_no][row_local][col_local];

                                    One_point_oprts[opr_no][site].value.push_back(value);
                                    One_point_oprts[opr_no][site].rows.push_back(j);
                                    One_point_oprts[opr_no][site].columns.push_back(i);
                                }
                            }
                        }
                        else{
                            for (int i=0;i<basis.basis_size;i++){
                                if(value_at_pos(i, site, basis.BASE) == col_local){
                                    value = One_point_oprts_onsite[opr_no][row_local][col_local];
                                    One_point_oprts[opr_no][site].value.push_back(value);
                                    One_point_oprts[opr_no][site].rows.push_back(i);
                                    One_point_oprts[opr_no][site].columns.push_back(i);
                                }
                            }
                        }
                    }
                }
            }
        }
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

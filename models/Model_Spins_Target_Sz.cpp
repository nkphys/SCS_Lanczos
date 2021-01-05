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


void MODEL_Spins_Target_Sz::Add_diagonal_terms(BASIS_Spins_Target_Sz &basis){



    //    //Remember H[l][m]=<l|H|m>

    //    double_type value;
    //    for (int i=basis.D_min;i<basis.D_max + 1;i++){

    //        value=zero;

    //        //magnetic Field in Z-direction
    //        for(int site=0;site<basis.Length;site++){
    //            value+=one*(H_field[2])*
    //                    ( ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );
    //        }


    //        //(Sz_local)^2 exchange
    //        for(int site_i=0;site_i<basis.Length;site_i++){
    //            if(D_anisotropy[2]!=zero){
    //                value+=  D_anisotropy[2]*
    //                        ( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
    //                          ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))
    //                          );
    //            }
    //        }

    //        if(value!=zero){
    //            Hamil.value.push_back(value*one);
    //            Hamil.rows.push_back(i);
    //            Hamil.columns.push_back(i);
    //        }

    //    }


}


void MODEL_Spins_Target_Sz::Add_non_diagonal_terms(BASIS_Spins_Target_Sz &basis){


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
                        J3_vals.push_back(one*J2_mat[site_i][site_j][site_k][site_l]);
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


    //0.5*S+(i)S-(j)------------------------------------------------------------------
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
        int m_new;
        int n_index;
        Mat_1_int state_vec;
        int site_i, site_j;


        for (int m=0;m<basis.D_basis.size();m++){
            dec_m = basis.D_basis[m];
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

                // m_new = Find_int_in_intarray(m_connected_final[j], basis.D_basis);
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


    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string twotimesspin, TwoTimesSpin = "TwoTimesSpin = ";

    string twotimestotalsztarget, TwoTimesTotalSzTarget = "TwoTimesTotalSzTarget = ";

    string hmag, Hmag = "H_mag = ";
    string d_anisotropy_, D_Anisotropy_ = "D_anisotropy = ";

    string LongRangeJ1file_ = "LongRangeJ1_file = ";

    string LongRangeJ2file_ = "LongRangeJ2_file = ";

    string LongRangeJ3file_ = "LongRangeJ3_file = ";


    string fourpointobservablessitesfile_ ,FourPointObservablesSitesFile_ = "FourPointObservablesSites file = ";

    string no_of_onepoint_obs_, No_Of_Onepoint_Obs_ = "No_of_onepoint_obs = ";

    int offset;
    string line;
    ifstream inputfile(filepath.c_str());

    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(No_Of_Onepoint_Obs_, 0)) != string::npos) {
                no_of_onepoint_obs_ = line.substr (offset+No_Of_Onepoint_Obs_.length());  }

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

    No_of_onepoint_obs=atoi(no_of_onepoint_obs_.c_str());



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
    D_anisotropy.clear();
    D_anisotropy.resize(3); //Dx, Dy, Dz
    for(int i=0;i<3;i++){
        anis_stream >> temp_anis;
        D_anisotropy[i]=temp_anis;
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



void MODEL_Spins_Target_Sz::Initialize_Opr_for_Dynamics(BASIS_Spins_Target_Sz &basis){

    Dyn_opr.value.clear();
    Dyn_opr.rows.clear();
    Dyn_opr.columns.clear();
    Dyn_opr.nrows = basis.D_basis.size();
    Dyn_opr.ncols = basis.D_basis.size();


    double_type value_;
    int dec_;

    if(Dyn_opr_string=="Sz"){
    for (int m=0;m<basis.D_basis.size();m++){
        dec_ = basis.D_basis[m];
        value_ =0.0;
        for(int opr_no=0;opr_no<Dyn_opr_int.size();opr_no++){
            value_ += Dyn_opr_coeffs[opr_no]*(
                        ((1.0*value_at_pos(dec_, Dyn_opr_int[opr_no], basis.BASE)) - (0.5*basis.TwoTimesSpin))
                        );

        }
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

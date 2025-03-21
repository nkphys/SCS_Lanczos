/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "Model_Spins.h"
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
#define PI 3.14159265

/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ]
    site----->[012....................(L-1)]


2)  similarly for "down spin" basis

3)  For total
    m=basis.D_up_basis.size();
*/


void MODEL_Spins::Act_Hamil(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert (Vec_in.size() == basis.D_max + 1 - basis.D_min);
    Vec_out.clear();
    Vec_out.resize(basis.D_max + 1 - basis.D_min);
    Act_diagonal_terms(basis, Vec_in, Vec_out);
    cout<<"Diagonal done"<<endl;
   // Act_non_diagonal_terms(basis, Vec_in, Vec_out);
   // cout<<"Non diagonal done"<<endl;
    Act_connections(basis, Vec_in, Vec_out);
    cout<<"Connections done"<<endl;

}


void MODEL_Spins::Act_diagonal_terms(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert (Vec_out.size() == basis.D_max + 1 - basis.D_min);
    assert (Vec_in.size() == Vec_out.size());


    //Remember H[l][i]=<l|H|i>

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
   for (int i=basis.D_min;i<basis.D_max + 1;i++){
            double_type value;
            value=zero;
            
        //magnetic Field in Z-direction
        for(int site=0;site<basis.Length;site++){
            value+=one*(Hz_field[site])*
                    ( ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );
        }



        //SzSz exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                if(Jzz_Exchange_mat[site_i][site_j]!=zero){
                    assert(Jzz_Exchange_mat[site_i][site_j] == Jzz_Exchange_mat[site_j][site_i]);
                    value+=  Jzz_Exchange_mat[site_i][site_j]*
                            (0.5*( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                                   ( (1.0*value_at_pos(i, site_j, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                                   ));

                }
            }
        }

        //(Sz_local)^2 exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            if(D_anisotropy[2]!=zero){
                value+=  D_anisotropy[2]*
                        ( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                          ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                          );
            }
        }




            if(abs(value)>0){
                Vec_out[i] +=value*one*Vec_in[i];
            }
        
    }

}




void MODEL_Spins::Add_arbitraryconnections_from_files_new(BASIS_Spins &basis){

    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;



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
    double EPS_=0.000001;
    //Mat_1_int state_vec;

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


#ifdef _OPENMP
#pragma omp for
#endif
    for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){

        thread_id=0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif

        //dec_m = basis.D_basis[m];
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

                // fromDeci_to_Vecint(state_vec, basis.BASE, dec_vec_out_final[j] , basis.Length);
                // quicksort(state_vec, 0, state_vec.size() -1);
                // fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                // n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                // m_new = Find_int_in_part_of_intarray(dec_vec_out_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


                if(dec_vec_out_final[j]<=dec_m){
                    Hamil_private[thread_id].value.push_back(val_vec_out_final[j]*one);
                    Hamil_private[thread_id].rows.push_back(dec_vec_out_final[j]);
                    Hamil_private[thread_id].columns.push_back(dec_m);
                }
            }

        }


        if(dec_m%10000==0){
            cout<<"dec_m = "<<dec_m<<" done "<<endl;
        }

    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        Hamil.value.insert(Hamil.value.end(),Hamil_private[thread].value.begin(), Hamil_private[thread].value.end() );
        Hamil.rows.insert(Hamil.rows.end(),Hamil_private[thread].rows.begin(), Hamil_private[thread].rows.end() );
        Hamil.columns.insert(Hamil.columns.end(),Hamil_private[thread].columns.begin(), Hamil_private[thread].columns.end() );
    }

    //state_vec.clear();
    //dec_vec_out.clear();
    //val_vec_out.clear();
    // connection_stream.str(std::string());


}


void MODEL_Spins::Add_arbitraryconnections_from_files(BASIS_Spins &basis){

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


    //for (int m=0;m<basis.D_basis.size();m++){
    for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){


        //dec_m = basis.D_basis[m];
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

                // fromDeci_to_Vecint(state_vec, basis.BASE, dec_vec_out_final[j] , basis.Length);
                // quicksort(state_vec, 0, state_vec.size() -1);
                // fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                // n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                // m_new = Find_int_in_part_of_intarray(dec_vec_out_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);


                if(dec_vec_out_final[j]<=dec_m){
                    Hamil.value.push_back(val_vec_out_final[j]*one);
                    Hamil.rows.push_back(dec_vec_out_final[j]);
                    Hamil.columns.push_back(dec_m);
                }
            }

        }


        if(dec_m%5000==0){
            cout<<"dec_m = "<<dec_m<<" done "<<endl;
        }

    }


    state_vec.clear();
    dec_vec_out.clear();
    val_vec_out.clear();
    // connection_stream.str(std::string());


}


void MODEL_Spins::Get_LocalOPR(Matrix<double_type> &OPR_, string opr_type){

    if(opr_type=="Id"){
        OPR_=Idendity;
    }
    else if(opr_type=="Sz"){
        OPR_=SzLocal;
    }
    else if(opr_type=="Sx"){
        OPR_=SxLocal;
    }
    else if(opr_type=="SyI"){
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
    else if (opr_type=="TauZ"){
        OPR_ = TauZLocal;
    }
    else{
        cout<<"OPR TYPE "<<opr_type<<" not present"<<endl;
    }

}

void MODEL_Spins::Act_LocalOprString_by_Recursion(Mat_1_int oprs_site, Mat_1_string oprs_list, Mat_1_ullint & dec_vec_out, Mat_1_doub &val_vec_out,BASIS_Spins & basis, int opr_no){

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



void MODEL_Spins::Act_connections(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert (Vec_out.size() == basis.D_max + 1 - basis.D_min);
    assert (Vec_in.size() == Vec_out.size());
    


    double_type phase;
    double_type value;
    int m,j;
    unsigned long long int i_new, i_new_temp;
    unsigned long long int m_new;
    int Val_site, Val_site_p;
    int Val_site_new, Val_site_p_new;
    double_type Factor;
    complex<double> iota_(0.0,1.0);


    for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){

        value=zero;

        for(int site_p=0;site_p<basis.Length ;site_p++){
            Val_site_p = value_at_pos(i, site_p, basis.BASE);

            for(int site=site_p+1;site<basis.Length ;site++){
                Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);


                if((Jpm_Exchange_mat[site_p][site]!=zero)){
                    //SpSm-connections
                    //Sp_site_p*Sm_site  Exchange coupling:
                    //Jpm[site_p][site] x Sp_site_p x Sm_site
                    // + (Jpm[site_p][site])* x Sm_site_p x Sp_site

                    //site cannot be in -Spin, and
                    //site_p cannot be in Spin

                    if(
                            (Val_site != 0)
                            &&
                            (Val_site_p != (basis.BASE - 1))
                            )
                    {

                        Val_site_p_new = Val_site_p + 1;
                        Val_site_new = Val_site - 1;


                        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);
                        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site_p, basis.BASE, Val_site_p_new);


                        Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                       ((Val_site - (0.5*basis.TwoTimesSpin))*
                                        (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );
                        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                              ((Val_site_p - (0.5*basis.TwoTimesSpin))*
                                               (Val_site_p_new - (0.5*basis.TwoTimesSpin))  ) );

                        assert(i_new<i);

                        Hamil.value.push_back((0.5*(Jpm_Exchange_mat[site_p][site]))*Factor*one);
                        Hamil.rows.push_back(i_new);
                        Hamil.columns.push_back(i);

                    } // if SpSm possible

                }

            }//site_p

        } // site


     //Hx field
        for(int site=0;site<basis.Length ;site++){
            Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);

         //Sminus*0.5*Hx_field
            if(Val_site != 0){
                Val_site_new = Val_site - 1;
                i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);

                Factor = Hx_field[site]*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                               ((Val_site - (0.5*basis.TwoTimesSpin))*
                                (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );

                assert(i_new<i);

                Hamil.value.push_back(Factor*one);
                Hamil.rows.push_back(i_new);
                Hamil.columns.push_back(i);

            }

        }



    } // "i" i.e up_decimals
    
    
}


void MODEL_Spins::Update_Hamiltonian_Params(BASIS_Spins &basis, double Hx_factor,double Hz_factor, double Jpm_factor, double Jzz_factor){





    for(int site=0;site<basis.Length;site++){
        Hz_field[site] = Hz_field_const[site]*Hz_factor;
    }

    for(int site_p=0;site_p<basis.Length ;site_p++){
        for(int site=0;site<basis.Length ;site++){
        Jpm_Exchange_mat[site_p][site] = Jpm_Exchange_mat_const[site_p][site]*Jpm_factor;
        }}


    for(int site=0;site<basis.Length;site++){
        Hx_field[site] = Hx_field_const[site]*Hx_factor;
    }

    for(int site_p=0;site_p<basis.Length ;site_p++){
        for(int site=0;site<basis.Length ;site++){
        Jzz_Exchange_mat[site_p][site] = Jzz_Exchange_mat_const[site_p][site]*Jzz_factor;
        }}



}


void MODEL_Spins::Update_Hamiltonian_Params_with_multicolors(BASIS_Spins &basis, Mat_1_real& A_factor, Mat_1_real& B_factor){


    for(int site=0;site<basis.Length;site++){
        Hz_field[site] = Hz_field_const[site]*B_factor[site];
    }

    for(int site_p=0;site_p<basis.Length ;site_p++){
        for(int site=0;site<basis.Length ;site++){
        Jpm_Exchange_mat[site_p][site] = Jpm_Exchange_mat_const[site_p][site]*0.0;
        }}


    for(int site=0;site<basis.Length;site++){
        Hx_field[site] = Hx_field_const[site]*A_factor[site];
    }

    for(int site_p=0;site_p<basis.Length ;site_p++){
        for(int site=0;site<basis.Length ;site++){
        assert(B_factor[site]>=0.0);
        assert(B_factor[site_p]>=0.0);
        Jzz_Exchange_mat[site_p][site] = Jzz_Exchange_mat_const[site_p][site]*
                                         sqrt(B_factor[site]*B_factor[site_p]);
        }}



}


void MODEL_Spins::Get_BipartiteEntanglement(BASIS_Spins &basis, Mat_1_int &Sys1_, Mat_1_int &Sys2_, Mat_1_doub &Vec_, double & VonNuemannEntropy){


double VonNuemannEntropy_temp;
BASIS_Spins basis_sys1;
BASIS_Spins basis_sys2;

Matrix<double_type> Psi_LB;

/*
    int Length;
    int TwoTimesSpin;
    double SPIN;
    int BASE;
    int basis_size;
*/

basis_sys1.Length=Sys1_.size();
basis_sys1.TwoTimesSpin=basis.TwoTimesSpin;
basis_sys1.Construct_basis();

basis_sys2.Length=Sys2_.size();
basis_sys2.TwoTimesSpin=basis.TwoTimesSpin;
basis_sys2.Construct_basis();

assert( (basis_sys1.Length+basis_sys2.Length) == basis.Length);
//also check mutually exclusive systems


Psi_LB.resize(basis_sys1.basis_size, basis_sys2.basis_size);

Mat_1_int State_Full, State_Sys1_, State_Sys2_;

int dec1_, dec2_;
for(int i=basis.D_min;i<basis.D_max+1;i++){
fromDeci_to_Vecint(State_Full, basis.BASE, i, basis.Length);
State_Sys1_.clear();
State_Sys2_.clear();
for(int i1_=0;i1_<basis_sys1.Length;i1_++){
State_Sys1_.push_back(State_Full[Sys1_[i1_]]);
}
for(int i2_=0;i2_<basis_sys2.Length;i2_++){
State_Sys2_.push_back(State_Full[Sys2_[i2_]]);
}
fromVecint_to_Deci(State_Sys1_, basis_sys1.BASE , dec1_, basis_sys1.Length);
fromVecint_to_Deci(State_Sys2_, basis_sys2.BASE , dec2_, basis_sys2.Length);
Psi_LB(dec1_,dec2_) = Vec_[i];
}



//Perform SVD
int r_;
r_=min(basis_sys1.basis_size,basis_sys2.basis_size);
Matrix<double_type> VT_; //mxm
Matrix<double_type> U_;  //nxn
vector<double> Sigma_; //atmost non-zero min(n,m) values

Perform_SVD(Psi_LB,VT_,U_,Sigma_);

VonNuemannEntropy_temp=0.0;
for(int n=0;n<r_;n++){
VonNuemannEntropy_temp += -1.0*Sigma_[n]*Sigma_[n]*(log2(Sigma_[n]*Sigma_[n]));
}


VonNuemannEntropy=VonNuemannEntropy_temp;

}


void MODEL_Spins::Add_diagonal_terms(BASIS_Spins &basis){


    Hamil.value.clear();
    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;

    //Remember H[l][m]=<l|H|m>

    double_type value;
    for (int i=basis.D_min;i<basis.D_max + 1;i++){

        value=zero;

        //magnetic Field in Z-direction
        for(int site=0;site<basis.Length;site++){
            value+=one*(Hz_field[site])*
                    ( ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );
        }



        //SzSz exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                if(Jzz_Exchange_mat[site_i][site_j]!=zero){
                    assert(Jzz_Exchange_mat[site_i][site_j] == Jzz_Exchange_mat[site_j][site_i]);
                    value+=  Jzz_Exchange_mat[site_i][site_j]*
                            (0.5*( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                                   ( (1.0*value_at_pos(i, site_j, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                                   ));

                }
            }
        }

        //(Sz_local)^2 exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            if(D_anisotropy[2]!=zero){
                value+=  D_anisotropy[2]*
                        ( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                          ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                          );
            }
        }

        if(value!=zero){
            Hamil.value.push_back(value*one);
            Hamil.rows.push_back(i);
            Hamil.columns.push_back(i);
        }

    }


}


void MODEL_Spins::Add_non_diagonal_terms(BASIS_Spins &basis){


}

void MODEL_Spins::Add_connections(BASIS_Spins &basis){

    double_type phase;
    double_type value;
    int m,j;
    unsigned long long int i_new, i_new_temp;
    unsigned long long int m_new;
    int Val_site, Val_site_p;
    int Val_site_new, Val_site_p_new;
    double_type Factor;
    complex<double> iota_(0.0,1.0);

    for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){

        value=zero;

        for(int site_p=0;site_p<basis.Length ;site_p++){
            Val_site_p = value_at_pos(i, site_p, basis.BASE);

            for(int site=site_p+1;site<basis.Length ;site++){
                Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);


                if((Jpm_Exchange_mat[site_p][site]!=zero)){
                    //SpSm-connections
                    //Sp_site_p*Sm_site  Exchange coupling:
                    //Jpm[site_p][site] x Sp_site_p x Sm_site
                    // + (Jpm[site_p][site])* x Sm_site_p x Sp_site

                    //site cannot be in -Spin, and
                    //site_p cannot be in Spin

                    if(
                            (Val_site != 0)
                            &&
                            (Val_site_p != (basis.BASE - 1))
                            )
                    {

                        Val_site_p_new = Val_site_p + 1;
                        Val_site_new = Val_site - 1;


                        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);
                        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site_p, basis.BASE, Val_site_p_new);


                        Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                       ((Val_site - (0.5*basis.TwoTimesSpin))*
                                        (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );
                        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                              ((Val_site_p - (0.5*basis.TwoTimesSpin))*
                                               (Val_site_p_new - (0.5*basis.TwoTimesSpin))  ) );

                        assert(i_new<i);

                        Hamil.value.push_back((0.5*(Jpm_Exchange_mat[site_p][site]))*Factor*one);
                        Hamil.rows.push_back(i_new);
                        Hamil.columns.push_back(i);

                    } // if SpSm possible

                }

            }//site_p

        } // site


        //Hx field
           for(int site=0;site<basis.Length ;site++){
               Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);

            //Sminus*0.5*Hx_field
               if(Val_site != 0){
                   Val_site_new = Val_site - 1;
                   i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);

                   Factor = Hx_field[site]*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                  ((Val_site - (0.5*basis.TwoTimesSpin))*
                                   (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );
                   assert(i_new<i);
                   Hamil.value.push_back(Factor*one);
                   Hamil.rows.push_back(i_new);
                   Hamil.columns.push_back(i);
               }
           }


    } // "i" i.e up_decimals

}



void MODEL_Spins::MeasureEnergy(BASIS_Spins &basis, Mat_1_doub &Vec, double_type &Energy_){

    Energy_=0;

//Connections
    double_type value;
    unsigned long long int i_new, i_new_temp;
    int Val_site, Val_site_p;
    int Val_site_new, Val_site_p_new;
    double_type Factor;

    for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){

        value=zero;

        for(int site_p=0;site_p<basis.Length ;site_p++){
            Val_site_p = value_at_pos(i, site_p, basis.BASE);

            for(int site=site_p+1;site<basis.Length ;site++){
                Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);


                if((Jpm_Exchange_mat[site_p][site]!=zero)){
                    //SpSm-connections
                    //Sp_site_p*Sm_site  Exchange coupling:
                    //Jpm[site_p][site] x Sp_site_p x Sm_site
                    // + (Jpm[site_p][site])* x Sm_site_p x Sp_site

                    //site cannot be in -Spin, and
                    //site_p cannot be in Spin

                    if(
                            (Val_site != 0)
                            &&
                            (Val_site_p != (basis.BASE - 1))
                            )
                    {

                        Val_site_p_new = Val_site_p + 1;
                        Val_site_new = Val_site - 1;


                        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);
                        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site_p, basis.BASE, Val_site_p_new);


                        Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                       ((Val_site - (0.5*basis.TwoTimesSpin))*
                                        (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );
                        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                              ((Val_site_p - (0.5*basis.TwoTimesSpin))*
                                               (Val_site_p_new - (0.5*basis.TwoTimesSpin))  ) );

                        assert(i_new<i);

                        Energy_ += (0.5*(Jpm_Exchange_mat[site_p][site]))*Factor*conjugate(Vec[i_new])*Vec[i];
                        Energy_ += conjugate((0.5*(Jpm_Exchange_mat[site_p][site]))*Factor*conjugate(Vec[i_new])*Vec[i]);
                    } // if SpSm possible

                }

            }//site_p

        } // site


        //Hx field
           for(int site=0;site<basis.Length ;site++){
               Val_site = value_at_pos<unsigned long long int>(i, site, basis.BASE);

            //Sminus*0.5*Hx_field
               if(Val_site != 0){
                   Val_site_new = Val_site - 1;
                   i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);

                   Factor = Hx_field[site]*0.5*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                  ((Val_site - (0.5*basis.TwoTimesSpin))*
                                   (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );
                   assert(i_new<i);

                   Energy_ +=Factor*one*conjugate(Vec[i_new])*Vec[i];
                   Energy_ +=conjugate(Factor*one*conjugate(Vec[i_new])*Vec[i]);
               }
           }


    } // "i" i.e up_decimals



    //Diagonal terms
    for (int i=basis.D_min;i<basis.D_max + 1;i++){

        //magnetic Field in Z-direction
        for(int site=0;site<basis.Length;site++){
            Energy_ +=one*(Hz_field[site])*
                    ( ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin)) )*
                      conjugate(Vec[i])*Vec[i];
        }



        //SzSz exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                if(Jzz_Exchange_mat[site_i][site_j]!=zero){
                    assert(Jzz_Exchange_mat[site_i][site_j] == Jzz_Exchange_mat[site_j][site_i]);
                    Energy_ +=  Jzz_Exchange_mat[site_i][site_j]*
                            (0.5*( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                                   ( (1.0*value_at_pos(i, site_j, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                                   ))*conjugate(Vec[i])*Vec[i];

                }
            }
        }

        //(Sz_local)^2 exchange
        for(int site_i=0;site_i<basis.Length;site_i++){
            if(D_anisotropy[2]!=zero){
                Energy_ +=  D_anisotropy[2]*
                        ( ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                          ( (1.0*value_at_pos(i, site_i, basis.BASE)) - (0.5*basis.TwoTimesSpin))
                          )*conjugate(Vec[i])*Vec[i];
            }
        }
    }


}




void MODEL_Spins::MeasureTwoPointOprs(BASIS_Spins &basis, Mat_1_doub &Vec, string tag){

assert(Vec.size()==basis.basis_size);


double_type CHSH_std_01;
CHSH_std_01=0.0;

double_type CHSH_std_12;
CHSH_std_12=0.0;

//SzSz
double_type temp_obs;
double_type Sz2_total;
Sz2_total=0;
cout<<" SzSz["<<tag<<"][site1][site2]--------------------"<<endl;
for(int site1=0;site1<basis.Length;site1++){
for(int site2=0;site2<basis.Length;site2++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
     temp_obs+=Vec[i]*conjugate(Vec[i])*
     ( (1.0*value_at_pos(i, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
     ( (1.0*value_at_pos(i, site2, basis.BASE)) - (0.5*basis.TwoTimesSpin));
 }
cout<<temp_obs<<"  ";
Sz2_total +=temp_obs;

if(site1==0 &&site2==1){
    CHSH_std_01 = CHSH_std_01 - temp_obs*4.0;
}
if(site1==1 &&site2==2){
    CHSH_std_12 = CHSH_std_12 - temp_obs*4.0;
}
    }
cout<<endl;
}

cout<<"Sz2_Total ["<<tag<<"] : "<<Sz2_total<<endl;
cout<<"---------------------------"<<endl<<endl;



cout<<" SzSm["<<tag<<"][site1][site2]--------------------"<<endl;
int Val_site1, Val_site2;
int Val_site1_new, Val_site2_new;
int i_new, i_new_temp;
double Factor;
for(int site1=0;site1<basis.Length;site1++){
for(int site2=0;site2<basis.Length;site2++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
Val_site1 = value_at_pos(i, site1, basis.BASE);
Val_site2 = value_at_pos(i, site2, basis.BASE);
if(Val_site2 != 0){

        Val_site2_new = Val_site2 - 1;
        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site2, basis.BASE, Val_site2_new);

        Factor = ( (1.0*value_at_pos(i, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin));
        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site2 - (0.5*basis.TwoTimesSpin))*
                        (Val_site2_new - (0.5*basis.TwoTimesSpin))  ) );

        temp_obs+=Vec[i]*conjugate(Vec[i_new])*Factor;

}

 }
cout<<temp_obs<<"  ";
 }
cout<<endl;
}



cout<<" SzSp["<<tag<<"][site1][site2]--------------------"<<endl;
for(int site1=0;site1<basis.Length;site1++){
for(int site2=0;site2<basis.Length;site2++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
Val_site1 = value_at_pos(i, site1, basis.BASE);
Val_site2 = value_at_pos(i, site2, basis.BASE);
if(Val_site2 != (basis.BASE - 1)){

        Val_site2_new = Val_site2 + 1;
        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site2, basis.BASE, Val_site2_new);

        Factor = ( (1.0*value_at_pos(i, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin));
        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site2 - (0.5*basis.TwoTimesSpin))*
                        (Val_site2_new - (0.5*basis.TwoTimesSpin))  ) );

        temp_obs+=Vec[i]*conjugate(Vec[i_new])*Factor;

}

 }
cout<<temp_obs<<"  ";
 }
cout<<endl;
}



bool allowed;
cout<<" SpSm["<<tag<<"][site1][site2]--------------------"<<endl;
for(int site1=0;site1<basis.Length;site1++){
for(int site2=0;site2<basis.Length;site2++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
Val_site1 = value_at_pos(i, site1, basis.BASE);
Val_site2 = value_at_pos(i, site2, basis.BASE);

if(site1!=site2){
allowed = ((Val_site1 != (basis.BASE - 1))) && (Val_site2 != 0);
}
else{
//allowed = ((Val_site1 != (basis.BASE - 1))) && (Val_site2 != 0);
assert(site1==site2);
allowed = (Val_site2 != 0);
}

    if(allowed){

        Val_site2_new = Val_site2 - 1;
        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site2, basis.BASE, Val_site2_new);

        if(site1==site2){
        Val_site1=Val_site2_new;
        }

        Val_site1_new = Val_site1 + 1;
        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site1, basis.BASE, Val_site1_new);

        Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site1 - (0.5*basis.TwoTimesSpin))*
                        (Val_site1_new - (0.5*basis.TwoTimesSpin))  ) );
        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site2 - (0.5*basis.TwoTimesSpin))*
                        (Val_site2_new - (0.5*basis.TwoTimesSpin))  ) );

        temp_obs+=Vec[i]*conjugate(Vec[i_new])*Factor;

}
}
cout<<temp_obs<<"  ";
if(site1==0 &&site2==1){
    CHSH_std_01 = CHSH_std_01 - temp_obs*1.0;
}
if(site2==0 &&site1==1){
    CHSH_std_01 = CHSH_std_01 - temp_obs*1.0;
}
if(site1==1 &&site2==2){
    CHSH_std_12 = CHSH_std_12 - temp_obs*1.0;
}
if(site2==1 &&site1==2){
    CHSH_std_12 = CHSH_std_12 - temp_obs*1.0;
}

}
cout<<endl;
}


cout<<" SpSp["<<tag<<"][site1][site2]--------------------"<<endl;
for(int site1=0;site1<basis.Length;site1++){
for(int site2=0;site2<basis.Length;site2++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
Val_site1 = value_at_pos(i, site1, basis.BASE);
Val_site2 = value_at_pos(i, site2, basis.BASE);
if(Val_site1 != (basis.BASE - 1)){
    if(Val_site2 != (basis.BASE - 1)){

        if(site1!=site2){
        Val_site1_new = Val_site1 + 1;
        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site1, basis.BASE, Val_site1_new);

        Val_site2_new = Val_site2 + 1;
        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site2, basis.BASE, Val_site2_new);

        Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site1 - (0.5*basis.TwoTimesSpin))*
                        (Val_site1_new - (0.5*basis.TwoTimesSpin))  ) );
        Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                       ((Val_site2 - (0.5*basis.TwoTimesSpin))*
                        (Val_site2_new - (0.5*basis.TwoTimesSpin))  ) );

        temp_obs+=Vec[i]*conjugate(Vec[i_new])*Factor;

        }
        else{
            assert(site1==site2);
            Val_site1_new = Val_site1 + 1;
            i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site1, basis.BASE, Val_site1_new);
            Val_site1_new = Val_site1_new + 1;
           if(Val_site1_new<=(basis.BASE - 1)){
            i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site1, basis.BASE, Val_site1_new);

            Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                           ((Val_site1 - (0.5*basis.TwoTimesSpin))*
                            (Val_site1_new - (0.5*basis.TwoTimesSpin))  ) );
            Factor = Factor*sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                           ((Val_site2 - (0.5*basis.TwoTimesSpin))*
                            (Val_site2_new - (0.5*basis.TwoTimesSpin))  ) );

            temp_obs+=Vec[i]*conjugate(Vec[i_new])*Factor;
           }
        }



}}

 }
cout<<temp_obs<<"  ";
if(site1==0 &&site2==1){
    CHSH_std_01 = CHSH_std_01 - temp_obs*1.0 - conjugate(temp_obs)*1.0;
}
if(site1==1 &&site2==2){
    CHSH_std_12 = CHSH_std_12 - temp_obs*1.0 - conjugate(temp_obs)*1.0;
}
}
cout<<endl;
}


CHSH_std_01= CHSH_std_01*sqrt(2.0);
cout<< "CHSH_std_01 : "<<CHSH_std_01<<endl;

CHSH_std_12= CHSH_std_12*sqrt(2.0);
cout<< "CHSH_std_12 : "<<CHSH_std_12<<endl;

//HERE SxSz


}

void MODEL_Spins::MeasureLocalOprs(BASIS_Spins &basis, Mat_1_doub &Vec, string tag){

assert(Vec.size()==basis.basis_size);


//Sz
double_type temp_obs;
cout<<" Local Sz ("<<tag<<"): site  value--------------------"<<endl;
for(int site=0;site<basis.Length;site++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
     temp_obs+=Vec[i]*conjugate(Vec[i])*
     ( (1.0*value_at_pos(i, site, basis.BASE)) - (0.5*basis.TwoTimesSpin));
 }

cout<<"Sz local ["<<tag<<"] : "<<site<<"  "<<temp_obs<<endl;

}



//Splus
double Factor;
int Val_site, Val_site_new, i_new;
cout<<" Local Splus : site  value --------------------"<<endl;
for(int site=0;site<basis.Length;site++){
temp_obs=0.0;
for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
Val_site = value_at_pos(i, site, basis.BASE);

if(Val_site != (basis.BASE - 1)){
    Val_site_new = Val_site + 1;
    i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, Val_site_new);

    Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                   ((Val_site - (0.5*basis.TwoTimesSpin))*
                    (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );

    temp_obs+=Factor*Vec[i]*conjugate(Vec[i_new]);
}

 }

cout<<"Splus local ["<<tag<<"] : "<<site<<"  "<<temp_obs<<endl;

}




}


void MODEL_Spins::Act_Operator(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub &Vec_out, string opr_str, int opr_site){

    Vec_out.resize(Vec_in.size());

    if(opr_str=="Sz"){
    for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){
           Vec_out[i]=Vec_in[i]*
             ( (1.0*value_at_pos(i, opr_site, basis.BASE)) - (0.5*basis.TwoTimesSpin));
    }
    }

}

void MODEL_Spins::Read_parameters(BASIS_Spins &basis, string filename){


    string filepath = filename;
    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string twotimesspin, TwoTimesSpin = "TwoTimesSpin = ";

    string hmag, Hmag = "H_mag = ";
    string d_anisotropy_, D_Anisotropy_ = "D_anisotropy = ";

    string LongRangeExchangeZZfile_ = "LongRangeExchangeZZ file = ";

    string LongRangeExchangePMfile_ = "LongRangeExchangePM file = ";

    string MagFieldVectorFile_ = "MagFieldVectorFile = ";

    string fourpointobservablessitesfile_ ,FourPointObservablesSitesFile_ = "FourPointObservablesSites file = ";

    string no_of_onepoint_obs_, No_Of_Onepoint_Obs_ = "No_of_onepoint_obs = ";

    string genericconnectionsfiles_, GenericConnectionsFiles_ = "GenericConnectionsFiles = ";



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

            if ((offset = line.find(No_Of_Onepoint_Obs_, 0)) != string::npos) {
                no_of_onepoint_obs_ = line.substr (offset+No_Of_Onepoint_Obs_.length());  }

            if ((offset = line.find(FourPointObservablesSitesFile_, 0)) != string::npos) {
                fourpointobservablessitesfile_ = line.substr (offset+FourPointObservablesSitesFile_.length());  }

            if ((offset = line.find(LongRangeExchangeZZfile_, 0)) != string::npos) {
                LongRangeExchangeZZfilepath = line.substr (offset+LongRangeExchangeZZfile_.length());  }

            if ((offset = line.find(LongRangeExchangePMfile_, 0)) != string::npos) {
                LongRangeExchangePMfilepath = line.substr (offset+LongRangeExchangePMfile_.length());  }


            if ((offset = line.find(MagFieldVectorFile_, 0)) != string::npos) {
                MagFieldVectorfilepath = line.substr (offset+MagFieldVectorFile_.length());  }

            if ((offset = line.find(PBC_, 0)) != string::npos) {
                pbc_ = line.substr (offset+PBC_.length());				}

            if ((offset = line.find(Length, 0)) != string::npos) {
                length = line.substr (offset + Length.length());		}

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

    basis.Length=atoi(length.c_str());
    basis.TwoTimesSpin=atoi(twotimesspin.c_str());


    No_of_onepoint_obs=atoi(no_of_onepoint_obs_.c_str());



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


    MagFieldVectorfilepath = MagFieldVectorfilepath;// + Extenstion_to_FilePaths;

    //cout<<"Reading Magnetic field :"<<endl;
   /* ifstream inMagfile(MagFieldVectorfilepath.c_str());
    double temp_hx, temp_hy , temp_hz;
    int temp_site;
    Hx_field_const.clear();Hy_field_const.clear();Hz_field_const.clear();
    Hx_field_const.resize(basis.Length);Hy_field_const.resize(basis.Length);Hz_field_const.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        inMagfile >> temp_site >> temp_hx >> temp_hy >> temp_hz ;
        Hx_field_const[i]=temp_hx;Hy_field_const[i]=temp_hy;Hz_field_const[i]=temp_hz;
        assert(i==temp_site);
        //cout<<i<<"  "<<Hx_field_const[i]<<"  "<<Hy_field_const[i]<<"  "<<Hz_field_const[i]<<endl;
    }*/
    //cout<<endl;

    /*
    Hx_field = Hx_field_const;
    Hy_field = Hy_field_const;
    Hz_field = Hz_field_const;

    double temp_anis;
    stringstream anis_stream;
    anis_stream<<d_anisotropy_;
    D_anisotropy_const.clear();
    D_anisotropy_const.resize(3); //Dx, Dy, Dz
    for(int i=0;i<3;i++){
        anis_stream >> temp_anis;
        D_anisotropy_const[i]=temp_anis;
    }
    D_anisotropy=D_anisotropy_const;

    LongRangeExchangePMfilepath = LongRangeExchangePMfilepath;// + Extenstion_to_FilePaths;
    LongRangeExchangeZZfilepath = LongRangeExchangeZZfilepath;// + Extenstion_to_FilePaths;

    Read_matrix_from_file(LongRangeExchangePMfilepath,
                          Jpm_Exchange_mat_const ,basis.Length,basis.Length);
    Read_matrix_from_file(LongRangeExchangeZZfilepath,
                          Jzz_Exchange_mat_const ,basis.Length,basis.Length);
    Jpm_Exchange_mat = Jpm_Exchange_mat_const;
    Jzz_Exchange_mat = Jzz_Exchange_mat_const;


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




void MODEL_Spins::CreateLocalOprs_in_LocalHilbertBasis(BASIS_Spins &basis){

    assert(basis.TwoTimesSpin>0);
    int HS; //Hilbert Space Size
    HS=basis.TwoTimesSpin+1;
    double Spin=0.5*basis.TwoTimesSpin;
    //basis convention = [-S, -S+1,..., S-1, S]

    Idendity.resize(HS,HS);

    //Dipole Oprs
    SxLocal.resize(HS,HS);SyLocal.resize(HS,HS);SzLocal.resize(HS,HS);
    SplusLocal.resize(HS,HS);SminusLocal.resize(HS,HS);SyLocalTimesIota.resize(HS,HS);
    int mp_temp;
    for(int m=0;m<HS;m++){
        Idendity(m,m)=1.0;
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
        TauZLocal.resize(HS,HS);


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


        for(int m=0;m<HS;m++){
            for(int mp=0;mp<HS;mp++){
                TauZLocal(mp,m) = 0.5*(Sx2Local(mp,m) - Sy2Local(mp,m));
            }}


    }

}






void MODEL_Spins::Read_parameters_for_dynamics(string filename){

    string dyn_momentum_, Dyn_Momentum_ = "k = ";
    string dyn_momentum_resolved_, Dyn_Momentum_Resolved_ = "Momentum_resolved = ";
    string Dyn_opr_string_  = "Opr_for_Dynamics = ";


    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(Dyn_Momentum_Resolved_, 0)) != string::npos) {
                dyn_momentum_resolved_ = line.substr (offset + Dyn_Momentum_Resolved_.length());		}

            if ((offset = line.find(Dyn_Momentum_, 0)) != string::npos) {
                dyn_momentum_ = line.substr (offset + Dyn_Momentum_.length());		}

            if ((offset = line.find(Dyn_opr_string_, 0)) != string::npos) {
                Dyn_opr_string = line.substr (offset + Dyn_opr_string_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    Dyn_Momentum=atof(dyn_momentum_.c_str());

    if(dyn_momentum_resolved_=="true"){
        Dyn_Momentum_Resolved=true;
    }
    else{
        Dyn_Momentum_Resolved=false;
    }

}



void MODEL_Spins::Initialize_Opr_for_Dynamics(BASIS_Spins &basis, int site_){

    vector< int >().swap( Dyn_opr.columns );
    vector< int >().swap( Dyn_opr.rows );
    vector< double_type >().swap( Dyn_opr.value );
    Dyn_opr.ncols = Hamil.ncols;
    Dyn_opr.nrows = Hamil.nrows;

    double_type value;
    int Val_site, Val_site_new;
    double Factor;
    unsigned long long int i_new;


    if(Dyn_opr_string=="Sz"){
        for (unsigned long long int i=basis.D_min;i<basis.D_max + 1;i++){
            value=zero;
            //Sz[site]
            value=one*
                    ( ( (1.0*value_at_pos<unsigned long long int>(i, site_, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );


            if(value!=zero){
                Dyn_opr.value.push_back(value*one);
                Dyn_opr.rows.push_back(i);
                Dyn_opr.columns.push_back(i);
            }
        }
    }
    else if(Dyn_opr_string=="Sm"){

        for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){

            value=zero;
            Val_site = value_at_pos<unsigned long long int>(i, site_, basis.BASE);

            //Sm[site]:
            //site cannot be in -Spin

            if(Val_site != 0)
            {
                Val_site_new = Val_site - 1;
                i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site_, basis.BASE, Val_site_new);

                Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                               ((Val_site - (0.5*basis.TwoTimesSpin))*
                                (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );

                assert(i_new<i);

                Dyn_opr.value.push_back(Factor*one);
                Dyn_opr.rows.push_back(i_new);
                Dyn_opr.columns.push_back(i);

            } // if Sm possible
        } // "i" i.e up_decimals

    }
    else{
        cout<<"Dyn Opr: only Sz or Sm are allowed"<<endl;
    }

}




void MODEL_Spins::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site1, int site2, BASIS_Spins &basis){



    OPR_.nrows = basis.basis_size;
    OPR_.ncols = OPR_.nrows;

    Mat_1_int state_vec;

    double EPS_=0.000001;
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





    if(opr_type=="Svec.Svec"){
        Mat_1_string connection_type;
        connection_type.push_back("Sz Sz");
        connection_type.push_back("Sp Sm");
        connection_type.push_back("Sm Sp");
        Mat_1_doub connection_val;
        connection_val.push_back(1.0);
        connection_val.push_back(0.5);
        connection_val.push_back(0.5);

    for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){


        //dec_m = basis.D_basis[m];
        dec_vec_out_given_m.clear();
        val_vec_out_given_m.clear();


                // cout<<"here 2"<<endl;
                // connection_stream.str("");

        for(int connection_type_ind=0;connection_type_ind<connection_type.size();connection_type_ind++){
                string connections="2 "+ connection_type[connection_type_ind] + " "+to_string(site1) +" "+to_string(site2)+" "+to_string(connection_val[connection_type_ind]);
                stringstream connection_stream;
                connection_stream<<connections;
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

        Remove_repetitions(dec_vec_out_given_m, val_vec_out_given_m, dec_vec_out_final, val_vec_out_final);

        for(int j=0;j<dec_vec_out_final.size();j++){

            if(abs(val_vec_out_final[j])>EPS_){

                // fromDeci_to_Vecint(state_vec, basis.BASE, dec_vec_out_final[j] , basis.Length);
                // quicksort(state_vec, 0, state_vec.size() -1);
                // fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                // n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                // m_new = Find_int_in_part_of_intarray(dec_vec_out_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                    OPR_.value.push_back(val_vec_out_final[j]*one);
                    OPR_.rows.push_back(dec_vec_out_final[j]);
                    OPR_.columns.push_back(dec_m);

            }

        }

    }

    }


    state_vec.clear();
    dec_vec_out.clear();
    val_vec_out.clear();

    if(opr_type=="TauZTauZ"){
        Mat_1_string connection_type;
        connection_type.push_back("TauZ TauZ");
        Mat_1_doub connection_val;
        connection_val.push_back(1.0);
        for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){


            //dec_m = basis.D_basis[m];
            dec_vec_out_given_m.clear();
            val_vec_out_given_m.clear();


            // cout<<"here 2"<<endl;
            // connection_stream.str("");

            for(int connection_type_ind=0;connection_type_ind<connection_type.size();connection_type_ind++){
                string connections="2 "+ connection_type[connection_type_ind] + " "+to_string(site1) +" "+to_string(site2)+" "+to_string(connection_val[connection_type_ind]);
                stringstream connection_stream;
                connection_stream<<connections;
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

            Remove_repetitions(dec_vec_out_given_m, val_vec_out_given_m, dec_vec_out_final, val_vec_out_final);

            for(int j=0;j<dec_vec_out_final.size();j++){

                if(abs(val_vec_out_final[j])>EPS_){

                        // fromDeci_to_Vecint(state_vec, basis.BASE, dec_vec_out_final[j] , basis.Length);
                        // quicksort(state_vec, 0, state_vec.size() -1);
                        // fromVecint_to_Deci(state_vec, basis.BASE, dec_max, basis.Length);
                        // n_index = Find_int_in_intarray(dec_max, basis.Partitions_Dec);
                        // m_new = Find_int_in_part_of_intarray(dec_vec_out_final[j], basis.D_basis, basis.Partitions_pos[n_index].first, basis.Partitions_pos[n_index].second);

                    OPR_.value.push_back(val_vec_out_final[j]*one);
                    OPR_.rows.push_back(dec_vec_out_final[j]);
                    OPR_.columns.push_back(dec_m);

                }

            }

        }

    }

    state_vec.clear();
    dec_vec_out.clear();
    val_vec_out.clear();




}



double_type MODEL_Spins::Get_SzSz(int site1, int site2, Mat_1_doub &Vec_, BASIS_Spins &basis){



    double_type value=0.0;
    // for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){
    //dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];

    for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){


        //Sz(i)Sz(j)----------------------------------------------------------------
        value +=conjugate(Vec_[dec_m])*Vec_[dec_m]*((1.0*value_at_pos(dec_m, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin))*
                 ((1.0*value_at_pos(dec_m, site2, basis.BASE)) - (0.5*basis.TwoTimesSpin));

    }

    return value;

}

double_type MODEL_Spins::Get_Sz(int site1, Mat_1_doub &Vec_, BASIS_Spins &basis){


    ulli dec_m;
    double_type value=0.0;
    //for (int m=0;m<basis.MainIndex_to_Dec_part0_.size();m++){
    for (unsigned long long int dec_m=basis.D_min;dec_m<basis.D_max+1;dec_m++){


        //Using Li tables
        //dec_m = basis.MainIndex_to_Dec_part0_[m] +  pow(basis.BASE, basis.Partition_Length[0])*basis.MainIndex_to_Dec_part1_[m];


        //Sz(i)Sz(j)----------------------------------------------------------------
        value +=conjugate(Vec_[dec_m])*Vec_[dec_m]*((1.0*value_at_pos(dec_m, site1, basis.BASE)) - (0.5*basis.TwoTimesSpin));

    }

    return value;

}


void MODEL_Spins::Initialize_Opr_for_Dynamics(BASIS_Spins &basis){



    assert(Dyn_Momentum_Resolved);

    vector< int >().swap( Dyn_opr.columns );
    vector< int >().swap( Dyn_opr.rows );
    vector< double_type >().swap( Dyn_opr.value );
    Dyn_opr.ncols = Hamil.ncols;
    Dyn_opr.nrows = Hamil.nrows;


    Hamiltonian_1_COO Oprs_local;
    Oprs_local.resize(basis.Length);
    for(int site=0;site<basis.Length;site++){
        Oprs_local[site].nrows = basis.basis_size ;
        Oprs_local[site].ncols = Oprs_local[site].nrows;
    }






    double_type value;
    int Val_site, Val_site_new;
    double Factor;
    unsigned long long int i_new;


    if(Dyn_opr_string=="Sz"){

        for(int site_=0;site_<basis.Length;site_++){
            for (unsigned long long int i=basis.D_min;i<basis.D_max + 1;i++){
                value=zero;
                //Sz[site]
                value=one*
                        ( ( (1.0*value_at_pos<unsigned long long int>(i, site_, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );

                if(value!=zero){
                    Oprs_local[site_].value.push_back(value*one);
                    Oprs_local[site_].rows.push_back(i);
                    Oprs_local[site_].columns.push_back(i);
                }
            }
        }
    }
    else if(Dyn_opr_string=="Sm"){

        for(int site_=0;site_<basis.Length;site_++){
            for (unsigned long long int i=basis.D_min;i<basis.D_max+1;i++){

                value=zero;
                Val_site = value_at_pos<unsigned long long int>(i, site_, basis.BASE);

                //Sm[site]:
                //site cannot be in -Spin

                if(Val_site != 0)
                {
                    Val_site_new = Val_site - 1;
                    i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site_, basis.BASE, Val_site_new);

                    Factor = sqrt( (1.0*basis.SPIN*(1.0+basis.SPIN))  -
                                   ((Val_site - (0.5*basis.TwoTimesSpin))*
                                    (Val_site_new - (0.5*basis.TwoTimesSpin))  ) );

                    assert(i_new<i);

                    Oprs_local[site_].value.push_back(Factor*one);
                    Oprs_local[site_].rows.push_back(i_new);
                    Oprs_local[site_].columns.push_back(i);

                } // if Sm possible
            } // "i" i.e up_decimals
        }

    }
    else{
        cout<<"Dyn Opr: only Sz or Sm are allowed"<<endl;
    }


    //In Momentum space--------For PBC use cosine, for OBC use sin----------



    Matrix_COO temp;
    temp=Oprs_local[0];
    double_type value1, value2;
    for(int site=0;site<basis.Length-1;site++){
        if(PBC==false){
            value2=one*sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
        }
        else{
#ifdef USE_COMPLEX
            value2=exp(iota_comp*(1.0*(site+1))*Dyn_Momentum*PI)*sqrt(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
            cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif
        }
        if(site==0){
            if(PBC==false){
                value1=one*sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
            }
            else{
#ifdef USE_COMPLEX
                value1=exp(iota_comp*(1.0*site)*Dyn_Momentum*PI)*sqrt(1.0/(basis.Length));
#endif
            }
            Sum(temp, Oprs_local[site+1], temp, value1, value2);}
        else{
            Sum(temp, Oprs_local[site+1], temp, 1.0, value2);
        }

    }

    Dyn_opr=temp;

    //----------------------------------------------------------------------

    vector< int >().swap( temp.columns );
    vector< int >().swap( temp.rows );
    vector< double_type >().swap( temp.value );

    for(int site_=0;site_<basis.Length;site_++){
        vector< int >().swap( Oprs_local[site_].columns );
        vector< int >().swap( Oprs_local[site_].rows );
        vector< double_type >().swap( Oprs_local[site_].value );
    }



}


void MODEL_Spins::Initialize_one_point_to_calculate_from_file(BASIS_Spins &basis){

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
        unsigned long long int j;
        double_type value;


        for(int site=0;site<basis.Length;site++){

            One_point_oprts[opr_no][site].value.clear();
            One_point_oprts[opr_no][site].rows.clear();
            One_point_oprts[opr_no][site].columns.clear();

            for(int col_local=0;col_local<basis.BASE;col_local++){
                for(int row_local=0;row_local<basis.BASE;row_local++){

                    if(One_point_oprts_onsite[opr_no][row_local][col_local] != zero){


                        if(row_local !=col_local){
                            for (unsigned long long int i=0;i<basis.basis_size;i++){
                                if(value_at_pos<unsigned long long int>(i, site, basis.BASE) == col_local){
                                    j = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.BASE, row_local);
                                    value = One_point_oprts_onsite[opr_no][row_local][col_local];

                                    One_point_oprts[opr_no][site].value.push_back(value);
                                    One_point_oprts[opr_no][site].rows.push_back(j);
                                    One_point_oprts[opr_no][site].columns.push_back(i);
                                }
                            }
                        }
                        else{
                            for (int i=0;i<basis.basis_size;i++){
                                if(value_at_pos<unsigned long long int>(i, site, basis.BASE) == col_local){
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

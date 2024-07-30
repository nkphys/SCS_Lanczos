/*
This class includes the Model for which Lanczos is being done
*/
#include "Model_Moire_Kspace.h"
#include <stdlib.h>
using namespace std;
#define PI 3.14159265
#ifdef _OPENMP
#include <omp.h>
#endif

/*convention for basis:

1)  for "up-spin" basis
                     [_______________________  _  ] [_______________________  _  ]
    site(K1,K2)----->[012....................(L-1)] [012....................(L-1)]
    orbital    ----->[.........orb - "0"..........] [.........orb - "1"..........]

2)  similarly for "down spin" basis

*/


void MODEL_Moire_Kspace::Act_Hamil(BASIS_Moire_Kspace &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

}

void MODEL_Moire_Kspace::Add_Dispersion_term(BASIS_Moire_Kspace &basis){

    Hamil.columns.clear();
    Hamil.rows.clear();
    Hamil.value.clear();

    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;



    int no_of_proc;
    no_of_proc=1;

#ifdef _OPENMP
    int temp_int;
    temp_int = basis.D_up_basis.size();
    no_of_proc=min(temp_int, NProcessors_);
    omp_set_num_threads(no_of_proc);
    cout<<"Dispersion term adding: "<<no_of_proc<<" processors"<<endl;
#endif




//#ifdef _OPENMP
//#pragma omp parallel
//    {
//#endif
//#ifdef _OPENMP
//#pragma omp for
//#endif
//    for (int m=0;m<basis.D_up_basis.size();m++){
//        cout <<"m = "<<m<<endl;
//    }

//#ifdef _OPENMP
//    }
//#endif



#ifdef _OPENMP
#pragma omp parallel
    {
        Matrix_COO Hamil_parallel;
        Hamil_parallel.nrows = basis.basis_size;
        Hamil_parallel.ncols = Hamil_parallel.nrows;
#endif



#ifdef _OPENMP
#pragma omp for nowait
#endif
    for (int m=0;m<basis.D_up_basis.size();m++){
       // cout<<"m = "<< m<<endl;
         assert(m< basis.D_up_basis.size());

        int SPIN_UP=0;
        int SPIN_DN=1;

    for(int k1_1=0;k1_1<Length1;k1_1++){
    for(int k1_2=0;k1_2<Length2;k1_2++){
        int k1_ =  k1_1 + k1_2*Length1;

        for(int spin=0;spin<2;spin++){
            for(int orb=0;orb<n_orb;orb++){

        bool electron_in_k1_spin_orb;
         if(spin==SPIN_UP){
           electron_in_k1_spin_orb=(bit_value(basis.D_up_basis[m],orb*basis.Length + k1_)==1);
          }
         else{
             assert(spin==SPIN_DN);
             electron_in_k1_spin_orb=(bit_value(basis.D_dn_basis[m],orb*basis.Length + k1_)==1);
         }


         if(electron_in_k1_spin_orb){

#ifdef _OPENMP
        Hamil_parallel.value.push_back(Ek_dispersion[spin][orb][k1_]*one);
        Hamil_parallel.rows.push_back(m);
        Hamil_parallel.columns.push_back(m);
#endif
#ifndef _OPENMP
             Hamil.value.push_back(Ek_dispersion[spin][orb][k1_]*one);
             Hamil.rows.push_back(m);
             Hamil.columns.push_back(m);
#endif
         }



            }
        }

    }}


    }
    cout<<"done"<<endl;
#ifdef _OPENMP
#pragma omp critical
    {
Hamil.value.insert(Hamil.value.end(), Hamil_parallel.value.begin(), Hamil_parallel.value.end());
Hamil.rows.insert(Hamil.rows.end(), Hamil_parallel.rows.begin(), Hamil_parallel.rows.end());
Hamil.columns.insert(Hamil.columns.end(), Hamil_parallel.columns.begin(), Hamil_parallel.columns.end());
  }
}
#endif




}


void MODEL_Moire_Kspace::Add_Interaction_terms(BASIS_Moire_Kspace &basis){



#ifdef _OPENMP
#pragma omp parallel
    {

        Matrix_COO Hamil_parallel;
        Hamil_parallel.nrows = basis.basis_size;
        Hamil_parallel.ncols = Hamil_parallel.nrows;
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif

    for (int m=0;m<basis.D_up_basis.size();m++){
        //cout<< "m = "<<m<<endl;
        assert(m < basis.D_up_basis.size());
        int SPIN_UP=0;
        int SPIN_DN=1;
        double EnergyEps=0.00000001;
        complex<double> value;
        ulli D_up_0,D_dn_0;
        ulli D_up_1,D_dn_1;
        ulli D_up_2,D_dn_2;
        ulli D_up_3,D_dn_3;
        ulli D_up_4,D_dn_4;
        int i_new,j_new_min, j_new_max;
        int m_new;
        double sign_FM ;
        int sign_pow_up, sign_pow_dn ;

     D_up_0=basis.D_up_basis[m];
     D_dn_0=basis.D_dn_basis[m];

     value=0;
    for(int k1_1=0;k1_1<Length1;k1_1++){
    for(int k1_2=0;k1_2<Length2;k1_2++){
        int k1_ =  k1_1 + k1_2*Length1;

        for(int k2_1=0;k2_1<Length1;k2_1++){
        for(int k2_2=0;k2_2<Length2;k2_2++){
          int k2_ =  k2_1 + k2_2*Length1;

            for(int q_1=0;q_1<Length1;q_1++){
            for(int q_2=0;q_2<Length2;q_2++){
            int q_ =  q_1 + q_2*Length1;

            int kkq_dof = k1_ + Length*k2_  + Length*Length*q_;

            int k2_1_plus_q_1 = (k2_1 + q_1)%Length1;
            int k2_2_plus_q_2 = (k2_2 + q_2)%Length2;
            int k1_1_minus_q_1 = (k1_1 - q_1 + Length1)%Length1;
            int k1_2_minus_q_2 = (k1_2 - q_2 + Length2)%Length2;

            int k1_minus_q_ = k1_1_minus_q_1 + k1_2_minus_q_2*Length1;
            int k2_plus_q_ = k2_1_plus_q_1 + k2_2_plus_q_2*Length1;


          for(int n1_orb=0;n1_orb<n_orb;n1_orb++){
              for(int n2_orb=0;n2_orb<n_orb;n2_orb++){
                  for(int n3_orb=0;n3_orb<n_orb;n3_orb++){
                      for(int n4_orb=0;n4_orb<n_orb;n4_orb++){
                int orbs_dof = n1_orb + n_orb*n2_orb + n_orb*n_orb*n3_orb + n_orb*n_orb*n_orb*n4_orb;

                          for(int tau=0;tau<2;tau++){
                              for(int tau_p=0;tau_p<2;tau_p++){
                int spins_dof = tau + 2*tau_p;


//                if(spins_dof==1 && kkq_dof==107){
//                    cout<<"just a check"<<endl;
//                }

              //  int n4_k1_tau = (n4_orb*basis.Length + k1_) + tau*(n_orb*Length);
              //  int n3_k2_tau_p = (n3_orb*basis.Length + k2_) + tau_p*(n_orb*Length);
              //  int n2_k2_plus_q_tau_p = (n2_orb*basis.Length + k2_plus_q_) + tau_p*(n_orb*Length);
              //  int n1_k1_minus_q_tau = (n1_orb*basis.Length + k1_minus_q_) + tau*(n_orb*Length);


              if(abs(W_int[spins_dof][orbs_dof][kkq_dof])>EnergyEps){
                sign_FM=1.0;
                sign_pow_up=0;
                sign_pow_dn=0;
              //n1_k1_tau ne n3_k2_tau_p
             // assert(n4_k1_tau!=n3_k2_tau_p);
             // assert(n2_k2_plus_q_tau_p!=n1_k1_minus_q_tau);



              //electron in n1,k1,tau
              bool electron_in_n4_k1_tau;
               if(tau==SPIN_UP){
                electron_in_n4_k1_tau=(bit_value(basis.D_up_basis[m],n4_orb*basis.Length + k1_)==1);
               }
               else{
                assert(tau==SPIN_DN);
                electron_in_n4_k1_tau=(bit_value(basis.D_dn_basis[m],n4_orb*basis.Length + k1_)==1);
               }

               //acting a_{n1 k1 tau} on [|d_up(m)> X |d_dn(m)>]
               if(electron_in_n4_k1_tau){

               if(tau==SPIN_UP){
                D_up_1 = (ulli) (D_up_0 - pow(2,n4_orb*basis.Length + k1_));
                D_dn_1 = D_dn_0;
                if( (n4_orb*basis.Length + k1_)!=0 ){
                sign_pow_up += one_bits_in_bw(0,n4_orb*basis.Length + k1_ , D_up_0)
                               + bit_value(D_up_0,0);
                }
               }
               else{
               assert(tau==SPIN_DN);
               D_up_1 = D_up_0;
               D_dn_1 = (ulli) (D_dn_0 - pow(2,n4_orb*basis.Length + k1_));
               if( (n4_orb*basis.Length + k1_)!=0 ){
               sign_pow_dn += one_bits_in_bw(0,n4_orb*basis.Length + k1_ , D_dn_0)
                              + bit_value(D_dn_0,0);
               }
               sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_0)
                             + bit_value(D_up_0,0)
                             + bit_value(D_up_0,Length*n_orb -1);
               }


               //electron in n3,k2,tau_p
               bool electron_in_n3_k2_tau_p;
                if(tau_p==SPIN_UP){
                 electron_in_n3_k2_tau_p=(bit_value(D_up_1,n3_orb*basis.Length + k2_)==1);
                }
                else{
                 assert(tau_p==SPIN_DN);
                 electron_in_n3_k2_tau_p=(bit_value(D_dn_1,n3_orb*basis.Length + k2_)==1);
                }


                if(electron_in_n3_k2_tau_p){
                    if(tau_p==SPIN_UP){
                     D_up_2 = (ulli) (D_up_1 - pow(2,n3_orb*basis.Length + k2_));
                     D_dn_2 = D_dn_1;
                     if( (n3_orb*basis.Length + k2_)!=0 ){
                     sign_pow_up += one_bits_in_bw(0,n3_orb*basis.Length + k2_ , D_up_1)
                                   + bit_value(D_up_1,0);
                     }
                    }
                    else{
                    assert(tau_p==SPIN_DN);
                    D_up_2 = D_up_1;
                    D_dn_2 = (ulli) (D_dn_1 - pow(2,n3_orb*basis.Length + k2_));
                    if( (n3_orb*basis.Length + k2_)!=0 ){
                    sign_pow_dn += one_bits_in_bw(0,n3_orb*basis.Length + k2_ , D_dn_1)
                                   + bit_value(D_dn_1,0);
                    }
                    sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_1)
                                  + bit_value(D_up_1,0)
                                  + bit_value(D_up_1,Length*n_orb -1);
                    }



                    //hole in n2,k2_plus_q,tau_p
                    bool hole_in_n2_k2_plus_q_tau_p;
                     if(tau_p==SPIN_UP){
                      hole_in_n2_k2_plus_q_tau_p=(bit_value(D_up_2,n2_orb*basis.Length + k2_plus_q_)==0);
                     }
                     else{
                      assert(tau_p==SPIN_DN);
                      hole_in_n2_k2_plus_q_tau_p=(bit_value(D_dn_2,n2_orb*basis.Length + k2_plus_q_)==0);
                     }


                     if(hole_in_n2_k2_plus_q_tau_p){
                         if(tau_p==SPIN_UP){
                          D_up_3 = (ulli) (D_up_2 + pow(2,n2_orb*basis.Length + k2_plus_q_));
                          D_dn_3 = D_dn_2;
                          if( (n2_orb*basis.Length + k2_plus_q_)!=0 ){
                          sign_pow_up += one_bits_in_bw(0,n2_orb*basis.Length + k2_plus_q_ , D_up_2)
                                        + bit_value(D_up_2,0);
                          }
                         }
                         else{
                         assert(tau_p==SPIN_DN);
                         D_up_3 = D_up_2;
                         D_dn_3 = (ulli) (D_dn_2 + pow(2,n2_orb*basis.Length + k2_plus_q_));

                         if( (n2_orb*basis.Length + k2_plus_q_)!=0 ){
                         sign_pow_dn += one_bits_in_bw(0,n2_orb*basis.Length + k2_plus_q_ , D_dn_2)
                                        + bit_value(D_dn_2,0);
                         }
                         sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_2)
                                       + bit_value(D_up_2,0)
                                       + bit_value(D_up_2,Length*n_orb -1);
                         }


                         //hole in n1,k1_minus_q,tau
                         bool hole_in_n1_k1_minus_q_tau;
                          if(tau==SPIN_UP){
                           hole_in_n1_k1_minus_q_tau=(bit_value(D_up_3,n1_orb*basis.Length + k1_minus_q_)==0);
                          }
                          else{
                           assert(tau==SPIN_DN);
                           hole_in_n1_k1_minus_q_tau=(bit_value(D_dn_3,n1_orb*basis.Length + k1_minus_q_)==0);
                          }

                          if(hole_in_n1_k1_minus_q_tau){
                           if(tau==SPIN_UP){
                             D_up_4 = (ulli) (D_up_3 + pow(2,n1_orb*basis.Length + k1_minus_q_));
                             D_dn_4 = D_dn_3;
                             if( (n1_orb*basis.Length + k1_minus_q_)!=0 ){
                             sign_pow_up += one_bits_in_bw(0,n1_orb*basis.Length + k1_minus_q_ , D_up_3)
                                           + bit_value(D_up_3,0);
                             }
                             }
                            else{
                            assert(tau==SPIN_DN);
                            D_up_4 = D_up_3;
                            D_dn_4 = (ulli) (D_dn_3 + pow(2,n1_orb*basis.Length + k1_minus_q_));
                            if( (n1_orb*basis.Length + k1_minus_q_)!=0 ){
                            sign_pow_dn += one_bits_in_bw(0,n1_orb*basis.Length + k1_minus_q_ , D_dn_3)
                                           + bit_value(D_dn_3,0);
                            }
                            sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_3)
                                          + bit_value(D_up_3,0)
                                          + bit_value(D_up_3,Length*n_orb -1);
                             }


                i_new = Find_int_in_intarray(D_up_4,basis.D_up_basis_range);
                j_new_min = basis.D_up_basis_range_min[i_new];
                j_new_max = basis.D_up_basis_range_max[i_new];
                m_new = Find_int_in_part_of_intarray(D_dn_4, basis.D_dn_basis, j_new_min, j_new_max);
                assert(D_up_4==basis.D_up_basis[m_new]);


                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

//                if(sign_FM<0){
//                 cout<<"something wrong"<<endl;
//                }

                //append to Hamil(m_new,m)
                if(m>=m_new){
                 value = sign_FM*((W_int[spins_dof][orbs_dof][kkq_dof]))*one;
                    if(m==m_new){
                    value = complex<double>(value.real(),0.0);
                    }

#ifdef _OPENMP
                    Hamil_parallel.value.push_back(value);
                    Hamil_parallel.rows.push_back(m_new);
                    Hamil_parallel.columns.push_back(m);
#endif
#ifndef _OPENMP
        Hamil.value.push_back(value);
        Hamil.rows.push_back(m_new);
        Hamil.columns.push_back(m);
#endif


                    }


                }//hole_in_n1_k1_minus_q_tau
                }//hole_in_n2_k2_plus_q_tau_p
                }//electron_in_n3_k2_tau_p
               }//electron_in_n4_k1_tau





              }//W_int ne 0



                              }
                          }

                      }
                  }
              }
          }

            }}
    }}
}}



     } // "m" i.e decimal for up and dn both


#ifdef _OPENMP
#pragma omp critical
    {
Hamil.value.insert(Hamil.value.end(), Hamil_parallel.value.begin(), Hamil_parallel.value.end());
Hamil.rows.insert(Hamil.rows.end(), Hamil_parallel.rows.begin(), Hamil_parallel.rows.end());
Hamil.columns.insert(Hamil.columns.end(), Hamil_parallel.columns.begin(), Hamil_parallel.columns.end());
    }
    }
#endif


}





void MODEL_Moire_Kspace::Initialize_Sk_opr(BASIS_Moire_Kspace &basis, Matrix_COO &OPR_, int q_1, int q_2, string alpha, string beta){



    Matrix<complex<double>> Pauli_z, Pauli_plus, Pauli_minus;
    Pauli_z.resize(2,2);
    Pauli_plus.resize(2,2);
    Pauli_minus.resize(2,2);


    Pauli_z(0,0)=1.0;Pauli_z(1,1)=-1.0;
    Pauli_plus(0,1)=2.0;
    Pauli_minus(1,0)=2.0;

    Matrix<complex<double>> Pauli_alpha, Pauli_beta;
    if(alpha=="z"){
        Pauli_alpha=Pauli_z;
    }
    if(alpha=="plus"){
        Pauli_alpha=Pauli_plus;
    }
    if(alpha=="minus"){
        Pauli_alpha=Pauli_minus;
    }

    if(beta=="z"){
        Pauli_beta=Pauli_z;
    }
    if(beta=="plus"){
        Pauli_beta=Pauli_plus;
    }
    if(beta=="minus"){
        Pauli_beta=Pauli_minus;
    }



    assert(n_orb==1);



    OPR_.columns.clear();
    OPR_.rows.clear();
    OPR_.value.clear();

    OPR_.nrows = basis.basis_size;
    OPR_.ncols = OPR_.nrows;



#ifdef _OPENMP
#pragma omp parallel
    {

        Matrix_COO OPR_parallel;
        OPR_parallel.nrows = basis.basis_size;
        OPR_parallel.ncols = OPR_parallel.nrows;
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
    for (int m=0;m<basis.D_up_basis.size();m++){

        int SPIN_UP=0;
        int SPIN_DN=1;
        double EnergyEps=0.00000001;
        complex<double> value;
        ulli D_up_0,D_dn_0;
        ulli D_up_1,D_dn_1;
        ulli D_up_2,D_dn_2;
        ulli D_up_3,D_dn_3;
        ulli D_up_4,D_dn_4;
        int i_new,j_new_min, j_new_max;
        int m_new;
        double sign_FM;
        int sign_pow_up, sign_pow_dn;
        int q_ =  q_1 + q_2*Length1;

        int n1_orb=0;
        int n2_orb=0;
        int n3_orb=0;
        int n4_orb=0;

     D_up_0=basis.D_up_basis[m];
     D_dn_0=basis.D_dn_basis[m];

     value=0;
    for(int k1_1=0;k1_1<Length1;k1_1++){
    for(int k1_2=0;k1_2<Length2;k1_2++){
        int k1_ =  k1_1 + k1_2*Length1;

        for(int k2_1=0;k2_1<Length1;k2_1++){
        for(int k2_2=0;k2_2<Length2;k2_2++){
          int k2_ =  k2_1 + k2_2*Length1;

            int k2_1_plus_q_1 = (k2_1 + q_1)%Length1;
            int k2_2_plus_q_2 = (k2_2 + q_2)%Length2;
            int k1_1_minus_q_1 = (k1_1 - q_1 + Length1)%Length1;
            int k1_2_minus_q_2 = (k1_2 - q_2 + Length2)%Length2;

            int k1_minus_q_ = k1_1_minus_q_1 + k1_2_minus_q_2*Length1;
            int k2_plus_q_ = k2_1_plus_q_1 + k2_2_plus_q_2*Length1;


            for(int s_=0;s_<2;s_++){
                for(int s_p_=0;s_p_<2;s_p_++){
                for(int tau=0;tau<2;tau++){
                    for(int tau_p=0;tau_p<2;tau_p++){
               // int spins_dof = tau + 2*tau_p;


//                if(spins_dof==1 && kkq_dof==107){
//                    cout<<"just a check"<<endl;
//                }

              //  int n4_k1_tau = (n4_orb*basis.Length + k1_) + tau*(n_orb*Length);
              //  int n3_k2_tau_p = (n3_orb*basis.Length + k2_) + tau_p*(n_orb*Length);
              //  int n2_k2_plus_q_tau_p = (n2_orb*basis.Length + k2_plus_q_) + tau_p*(n_orb*Length);
              //  int n1_k1_minus_q_tau = (n1_orb*basis.Length + k1_minus_q_) + tau*(n_orb*Length);


              if(abs(Pauli_alpha(tau,tau_p)*Pauli_beta(s_,s_p_))>EnergyEps){
                sign_FM=1.0;
                sign_pow_up=0;
                sign_pow_dn=0;
              //n1_k1_tau ne n3_k2_tau_p
             // assert(n4_k1_tau!=n3_k2_tau_p);
             // assert(n2_k2_plus_q_tau_p!=n1_k1_minus_q_tau);



              //electron in n1,k1,tau
              bool electron_in_n4_k1_m_q_sp;
               if(s_p_==SPIN_UP){
                electron_in_n4_k1_m_q_sp=(bit_value(basis.D_up_basis[m],n4_orb*basis.Length + k1_minus_q_)==1);
               }
               else{
                assert(s_p_==SPIN_DN);
                electron_in_n4_k1_m_q_sp=(bit_value(basis.D_dn_basis[m],n4_orb*basis.Length + k1_minus_q_)==1);
               }

               //acting a_{n1 k1 tau} on [|d_up(m)> X |d_dn(m)>]
               if(electron_in_n4_k1_m_q_sp){

               if(s_p_==SPIN_UP){
                D_up_1 = (ulli) (D_up_0 - pow(2,n4_orb*basis.Length + k1_minus_q_));
                D_dn_1 = D_dn_0;
                if( (n4_orb*basis.Length + k1_minus_q_)!=0 ){
                sign_pow_up += one_bits_in_bw(0,n4_orb*basis.Length + k1_minus_q_ , D_up_0)
                               + bit_value(D_up_0,0);
                }
               }
               else{
               assert(s_p_==SPIN_DN);
               D_up_1 = D_up_0;
               D_dn_1 = (ulli) (D_dn_0 - pow(2,n4_orb*basis.Length + k1_minus_q_));
               if( (n4_orb*basis.Length + k1_minus_q_)!=0 ){
               sign_pow_dn += one_bits_in_bw(0,n4_orb*basis.Length + k1_minus_q_ , D_dn_0)
                              + bit_value(D_dn_0,0);
               }
               sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_0)
                             + bit_value(D_up_0,0)
                             + bit_value(D_up_0,Length*n_orb -1);
               }


               //hole in n3,k1,s_
               bool hole_in_n3_k1_s;
                if(s_==SPIN_UP){
                 hole_in_n3_k1_s =(bit_value(D_up_1, n3_orb*basis.Length + k1_)==0);
                }
                else{
                 assert(s_==SPIN_DN);
                 hole_in_n3_k1_s=(bit_value(D_dn_1,n3_orb*basis.Length + k1_)==0);
                }


                if(hole_in_n3_k1_s){
                    if(s_==SPIN_UP){
                     D_up_2 = (ulli) (D_up_1 + pow(2,n3_orb*basis.Length + k1_));
                     D_dn_2 = D_dn_1;
                     if( (n3_orb*basis.Length + k1_)!=0 ){
                     sign_pow_up += one_bits_in_bw(0,n3_orb*basis.Length + k1_ , D_up_1)
                                   + bit_value(D_up_1,0);
                     }
                    }
                    else{
                    assert(s_==SPIN_DN);
                    D_up_2 = D_up_1;
                    D_dn_2 = (ulli) (D_dn_1 + pow(2,n3_orb*basis.Length + k1_));
                    if( (n3_orb*basis.Length + k1_)!=0 ){
                    sign_pow_dn += one_bits_in_bw(0,n3_orb*basis.Length + k1_ , D_dn_1)
                                   + bit_value(D_dn_1,0);
                    }
                    sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_1)
                                  + bit_value(D_up_1,0)
                                  + bit_value(D_up_1,Length*n_orb -1);
                    }



                    //electron in n2,k2_plus_q,tau_p
                    bool electron_in_n2_k2_plus_q_tau_p;
                     if(tau_p==SPIN_UP){
                      electron_in_n2_k2_plus_q_tau_p=(bit_value(D_up_2,n2_orb*basis.Length + k2_plus_q_)==1);
                     }
                     else{
                      assert(tau_p==SPIN_DN);
                      electron_in_n2_k2_plus_q_tau_p=(bit_value(D_dn_2,n2_orb*basis.Length + k2_plus_q_)==1);
                     }


                     if(electron_in_n2_k2_plus_q_tau_p){
                         if(tau_p==SPIN_UP){
                          D_up_3 = (ulli) (D_up_2 - pow(2,n2_orb*basis.Length + k2_plus_q_));
                          D_dn_3 = D_dn_2;
                          if( (n2_orb*basis.Length + k2_plus_q_)!=0 ){
                          sign_pow_up += one_bits_in_bw(0,n2_orb*basis.Length + k2_plus_q_ , D_up_2)
                                        + bit_value(D_up_2,0);
                          }
                         }
                         else{
                         assert(tau_p==SPIN_DN);
                         D_up_3 = D_up_2;
                         D_dn_3 = (ulli) (D_dn_2 - pow(2,n2_orb*basis.Length + k2_plus_q_));

                         if( (n2_orb*basis.Length + k2_plus_q_)!=0 ){
                         sign_pow_dn += one_bits_in_bw(0,n2_orb*basis.Length + k2_plus_q_ , D_dn_2)
                                        + bit_value(D_dn_2,0);
                         }
                         sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_2)
                                       + bit_value(D_up_2,0)
                                       + bit_value(D_up_2,Length*n_orb -1);
                         }


                         //hole in n1,k2,tau
                         bool hole_in_n1_k2_tau;
                          if(tau==SPIN_UP){
                           hole_in_n1_k2_tau=(bit_value(D_up_3,n1_orb*basis.Length + k2_)==0);
                          }
                          else{
                           assert(tau==SPIN_DN);
                           hole_in_n1_k2_tau=(bit_value(D_dn_3,n1_orb*basis.Length + k2_)==0);
                          }

                          if(hole_in_n1_k2_tau){
                           if(tau==SPIN_UP){
                             D_up_4 = (ulli) (D_up_3 + pow(2,n1_orb*basis.Length + k2_));
                             D_dn_4 = D_dn_3;
                             if( (n1_orb*basis.Length + k2_)!=0 ){
                             sign_pow_up += one_bits_in_bw(0,n1_orb*basis.Length + k2_ , D_up_3)
                                           + bit_value(D_up_3,0);
                             }
                             }
                            else{
                            assert(tau==SPIN_DN);
                            D_up_4 = D_up_3;
                            D_dn_4 = (ulli) (D_dn_3 + pow(2,n1_orb*basis.Length + k2_));
                            if( (n1_orb*basis.Length + k2_)!=0 ){
                            sign_pow_dn += one_bits_in_bw(0,n1_orb*basis.Length + k2_ , D_dn_3)
                                           + bit_value(D_dn_3,0);
                            }
                            sign_pow_dn += one_bits_in_bw(0, n_orb*Length-1, D_up_3)
                                          + bit_value(D_up_3,0)
                                          + bit_value(D_up_3,Length*n_orb -1);
                             }


                i_new = Find_int_in_intarray(D_up_4,basis.D_up_basis_range);
                j_new_min = basis.D_up_basis_range_min[i_new];
                j_new_max = basis.D_up_basis_range_max[i_new];
                m_new = Find_int_in_part_of_intarray(D_dn_4, basis.D_dn_basis, j_new_min, j_new_max);
                assert(D_up_4==basis.D_up_basis[m_new]);


                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

//                if(sign_FM<0){
//                 cout<<"something wrong"<<endl;
//                }

                //append to Hamil(m_new,m)
                //if(m>=m_new){
                 value = 0.25*sign_FM*(Pauli_alpha(tau,tau_p)*Pauli_beta(s_,s_p_))*one;
//                    if(m==m_new){
//                    value = complex<double>(value.real(),0.0);
//                    }


#ifdef _OPENMP
                    OPR_parallel.value.push_back(value);
                    OPR_parallel.rows.push_back(m_new);
                    OPR_parallel.columns.push_back(m);
#endif
#ifndef _OPENMP
                    OPR_.value.push_back(value);
                    OPR_.rows.push_back(m_new);
                    OPR_.columns.push_back(m);
#endif


                    //}


                }//hole_in_n1_k1_minus_q_tau
                }//hole_in_n2_k2_plus_q_tau_p
                }//electron_in_n3_k2_tau_p
               }//electron_in_n4_k1_tau





              }//W_int ne 0



                              }
                          }
                }}


    }}
}}



     } // "m" i.e decimal for up and dn both


#ifdef _OPENMP
#pragma omp critical
    {
OPR_.value.insert(OPR_.value.end(), OPR_parallel.value.begin(), OPR_parallel.value.end());
OPR_.rows.insert(OPR_.rows.end(), OPR_parallel.rows.begin(), OPR_parallel.rows.end());
OPR_.columns.insert(OPR_.columns.end(), OPR_parallel.columns.begin(), OPR_parallel.columns.end());
  }
}
#endif


}

void MODEL_Moire_Kspace::Read_Dispersion(BASIS_Moire_Kspace &basis){



    Ek_dispersion.resize(2);
    for(int spin=0;spin<2;spin++){
      Ek_dispersion[spin].resize(n_orb);
      for(int orb_no=0;orb_no<n_orb;orb_no++){
      Ek_dispersion[spin][orb_no].resize(Length);
      }
    }

    bool UseSquareLatticeDispersion=false;
    double t1_hop=1.0;
    if(UseSquareLatticeDispersion){
        for(int spin=0;spin<2;spin++){
          for(int orb_no=0;orb_no<n_orb;orb_no++){
              for(int k1=0;k1<Length1;k1++){
                  for(int k2=0;k2<Length2;k2++){
                Ek_dispersion[spin][orb_no][k1+k2*Length1] = -2.0*t1_hop*
                                            ( cos(2.0*k1*PI/Length1) + cos(2.0*k2*PI/Length2) );
                  }
              }
          }
        }
    }
    else{ //Reading
    for(int orb=0;orb<n_orb;orb++){
    ifstream file_in(DispersionFilepaths[orb].c_str());
    string line_temp;
    getline(file_in,line_temp);
    int k_temp;
    double E_k_temp;
    while(getline(file_in,line_temp)){
     stringstream line_stream(line_temp);
     line_stream >> k_temp >> E_k_temp;
     for(int spin=0;spin<2;spin++){
     Ek_dispersion[spin][orb][k_temp]=E_k_temp;
    }
    }
    }//orb
    }






}

void MODEL_Moire_Kspace::Read_Interations(BASIS_Moire_Kspace &basis){

    assert(n_orb==1);

    W_int.resize(4); //spins
    for(int spins_dof=0;spins_dof<4;spins_dof++){
    W_int[spins_dof].resize(n_orb*n_orb*n_orb*n_orb);
    for(int orbs_dof=0;orbs_dof<W_int[spins_dof].size();orbs_dof++){
    W_int[spins_dof][orbs_dof].resize(Length*Length*Length);
    }
    }


    bool OnsiteHubbard=false;


    //For onsite Hubbard interaction i.e U\sum_{i}n_{i up}n_{i dn}
    if(OnsiteHubbard){
    double U0_=10.0;
    assert(n_orb==1);
    for(int spin=0;spin<2;spin++){
        for(int spin_p=0;spin_p<2;spin_p++){
        int spin_dof = spin +  2*spin_p;
        for(int kkq=0;kkq<Length*Length*Length;kkq++){
            if(spin==0 && spin_p==1){
         W_int[spin_dof][0][kkq]=U0_/Length;
            }
            else{
           W_int[spin_dof][0][kkq]=0;
            }
        }
        }
    }
    }
    else{
    ifstream InteractionFileStream(InteractionFilepath.c_str());
    string line_temp;
    getline(InteractionFileStream, line_temp);
    int m1, m2 ,m3, m4, spin, spin_p, k1, k2, q;
    double WReal, WImag;

    complex<double> W_temp;
    while(getline(InteractionFileStream, line_temp)){
    stringstream line_temp_stream(line_temp);
    line_temp_stream>>m1>>m2>>m3>>m4>>spin>>spin_p>>k1>>k2>>q>>WReal>>WImag;
    W_temp = complex<double> (WReal, WImag);

    int kkq = k1 + Length*k2  + Length*Length*q;

    if( (m2==m3)  && (m1==m4)){
        if( (m2==(spin_p+2))  && (m1==(spin+2))){

        W_int[spin +  2*spin_p][0][kkq] = W_temp;

//            if(spin==0 && spin_p==0){
//        W_int[spin +  2*spin_p][0][kkq] = W_temp;
//        W_int[(spin+1) +  2*(spin_p+1)][0][kkq] = W_temp;
//            }
//            if(spin==0 && spin_p==1){
//        W_int[spin +  2*spin_p][0][kkq] = W_temp;
//        W_int[spin_p +  2*(spin)][0][kkq] = W_temp;
//            }
        }
    }

    }

    }



}

void MODEL_Moire_Kspace::Read_parameters(BASIS_Moire_Kspace &basis, string filename){


    string filepath = filename;

    string length1_, Length1_ = "Length1 = ";
    string length2_, Length2_ = "Length2 = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string norb, Norb = "Norbitals = ";
    string k1_target, K1_target = "K1_target = ";
    string k2_target, K2_target = "K2_target = ";

    string dispersionfile, DispersionFile = "DispersionFiles = ";
    string interactionfile, InteractionFile = "InteractionFile = ";

    string processors_, Processors_ = "Processors = ";

    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Length1_, 0)) != string::npos) {
                length1_ = line.substr (offset + Length1_.length());		}

            if ((offset = line.find(Length2_, 0)) != string::npos) {
                length2_ = line.substr (offset + Length2_.length());		}


            if ((offset = line.find(K1_target, 0)) != string::npos) {
                k1_target = line.substr (offset + K1_target.length());		}

            if ((offset = line.find(K2_target, 0)) != string::npos) {
                k2_target = line.substr (offset + K2_target.length());		}

            if ((offset = line.find(Ndn, 0)) != string::npos) {
                ndn = line.substr (offset + Ndn.length());		}

            if ((offset = line.find(Nup, 0)) != string::npos) {
                nup= line.substr (offset + Nup.length());		}

            if ((offset = line.find(Norb, 0)) != string::npos) {
                norb = line.substr (offset + Norb.length());		}

            if ((offset = line.find(DispersionFile, 0)) != string::npos) {
                dispersionfile = line.substr (offset + DispersionFile.length());		}


            if ((offset = line.find(InteractionFile, 0)) != string::npos) {
                interactionfile = line.substr (offset + InteractionFile.length());		}

            if ((offset = line.find(Processors_ , 0)) != string::npos) {
                processors_ = line.substr (offset+Processors_ .length());	}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}



    basis.Length1=atoi(length1_.c_str());
    basis.Length2=atoi(length2_.c_str());
    basis.K1_target=atoi(k1_target.c_str());
    basis.K2_target=atoi(k2_target.c_str());
    basis.Ndn=atoi(ndn.c_str());
    basis.Nup=atoi(nup.c_str());
    basis.n_orb=atoi(norb.c_str());


    n_orb=basis.n_orb;
    Length1 = basis.Length1;
    Length2 = basis.Length2;
    Length = Length1*Length2;



    stringstream DispersionFile_stream;
    DispersionFile_stream<<dispersionfile;
    DispersionFilepaths.resize(n_orb);
    for(int orb=0;orb<n_orb;orb++){
    string temp_str;
    DispersionFile_stream>>temp_str;
    DispersionFilepaths[orb]=temp_str;
    }


    InteractionFilepath = interactionfile;

    NProcessors_=atoi(processors_.c_str());

}




#ifndef Model_1_orb_Hubbard_GC_Functions
#define Model_1_orb_Hubbard_GC_Functions

#include "Model_1_orb_Hubbard_GC.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
#define PI 3.14159265

namespace {

double_type Read_FourPoint_Coefficient(const string& token){

#ifdef USE_COMPLEX
    if(token.size()>2 && token[0]=='(' && token[token.size()-1]==')'){
        string coeff_string = token.substr(1, token.size()-2);
        size_t comma_pos = coeff_string.find(',');
        if(comma_pos != string::npos){
            double real_part = atof(coeff_string.substr(0, comma_pos).c_str());
            double imag_part = atof(coeff_string.substr(comma_pos + 1).c_str());
            return complex<double>(real_part, imag_part);
        }
    }

    return complex<double>(atof(token.c_str()), 0.0);
#endif

#ifndef USE_COMPLEX
    if(token.size()>2 && token[0]=='(' && token[token.size()-1]==')'){
        string coeff_string = token.substr(1, token.size()-2);
        size_t comma_pos = coeff_string.find(',');
        if(comma_pos != string::npos){
            return atof(coeff_string.substr(0, comma_pos).c_str());
        }
    }

    return atof(token.c_str());
#endif
}

double GetDoubleTypeMagnitude(const double_type& value_){
    return abs(value_);
}

double GetDoubleTypeSignedPart(const double_type& value_){
#ifdef USE_COMPLEX
    return value_.real();
#endif

#ifndef USE_COMPLEX
    return value_;
#endif
}

string DoubleTypeToString(const double_type& value_){
    stringstream value_stream;
    value_stream<<setprecision(6)<<value_;
    return value_stream.str();
}

}

/*convention for basis:

    1)  for "up-spin" basis
                  [_______________________  _  ]
        site----->[012....................(L-1)]


    2)  similarly for "down-spin" basis

    */


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Act_Hamil(BASIS_1_orb_Hubbard_GC &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert(Vec_in.size() == basis.D_up_basis.size());

    Vec_out.clear();
    Vec_out.resize(basis.D_up_basis.size());
    for(int i=0;i<basis.D_up_basis.size();i++){
        Vec_out[i] = zero;
    }

    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Mat_1_doub> Vec_out_private;
    Vec_out_private.resize(N_threads);
    for(int thread=0;thread<N_threads;thread++){
        Vec_out_private[thread].resize(basis.D_up_basis.size());
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out_private[thread][i] = zero;
        }
    }

    
#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    int alpha_, alpha_p_;
    int j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int max_up, max_dn, min_dn;
    int l,lp;
    double value;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){

        j=i;
        value=0;

        //intra-orbital coulomb repulsion:
        value += U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);

        //Long range density-density interaction: n_{i,s} n_{j,s'}
        for(int spin1=0;spin1<2;spin1++){
            for(int site1=0;site1<basis.Length;site1++){
                int alpha1 = basis.Length*spin1 + site1;
                int n_alpha1;
                if(spin1==0){
                    n_alpha1 = bit_value(basis.D_up_basis[i], site1);
                }
                else{
                    n_alpha1 = bit_value(basis.D_dn_basis[j], site1);
                }

                for(int spin2=0;spin2<2;spin2++){
                    for(int site2=0;site2<basis.Length;site2++){
                        int alpha2 = basis.Length*spin2 + site2;
                        int n_alpha2;
                        if(spin2==0){
                            n_alpha2 = bit_value(basis.D_up_basis[i], site2);
                        }
                        else{
                            n_alpha2 = bit_value(basis.D_dn_basis[j], site2);
                        }
                        value += DenDenInt_mat_LongRange[alpha1][alpha2]*(n_alpha1*n_alpha2);
                    }
                }
            }
        }

        //Crystal Field Splitting (CFE):
        for(int site=0;site<basis.Length;site++){
            value += (CFS[site])*
                    ( ( bit_value(basis.D_up_basis[i], site) +
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }

        //magnetic Field * [Sz]
        for(int site=0;site<basis.Length;site++){
            value += 0.5*(H_field[site])*
                    ( ( bit_value(basis.D_up_basis[i],site) -
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }

        Vec_out_private[thread_id][i] += (value*one)*Vec_in[i];


        for(int sigma=0;sigma<2 ;sigma++){
            for(int site=0;site<basis.Length ;site++){
                alpha_ = basis.Length*sigma + site;

                for(int sigma_p=0;sigma_p<2 ;sigma_p++){
                    for(int site_p=0;site_p<basis.Length ;site_p++){
                        alpha_p_ = basis.Length*sigma_p + site_p;


                        if(Hopping_mat_LongRange[alpha_p_][alpha_]!=zero){

                            //HOPPING COEFFICIENT IN FRONT OF
                            // (?) X c_{site_p,sigma_p}^{\dagger} c_{site,sigma}



                            if(sigma==0 && sigma_p==0){

                                //---------------Hopping: up to up electrons-------------------//
                                //there have to be one up electron on site
                                //there have to be no up electron on site_p
                                if(
                                        ( (bit_value(basis.D_up_basis[i],site)==1)
                                          &&
                                          (bit_value(basis.D_up_basis[i],site_p)==0)
                                          )
                                        &&
                                        (site_p<site)
                                        )
                                {


                                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                                  - pow(2, site) );
                                    D_dn = basis.D_dn_basis[j];

                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                    sign_FM = pow(-1.0, sign_pow_up);


                                    assert(m_new<i);
                                    double_type value_ = sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one;
                                    Vec_out_private[thread_id][m_new] += value_*Vec_in[i];
                                    Vec_out_private[thread_id][i] += conjugate(value_)*Vec_in[m_new];


                                } // if up-up hopping possible

                            }


                            if(sigma==1 && sigma_p==1){

                                //---------------Hopping: dn to dn electrons-------------------//
                                //there have to be one dn electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (
                                            (bit_value(basis.D_dn_basis[j], site)==1)
                                            &&
                                            (bit_value(basis.D_dn_basis[j], site_p)==0)
                                            )
                                        &&
                                        (site_p<site)
                                        )
                                {

                                    D_up = basis.D_up_basis[i];
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                                  - pow(2, site) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn);

                                    assert(m_new<i);
                                    double_type value_ = 1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one;
                                    Vec_out_private[thread_id][m_new] += value_*Vec_in[i];
                                    Vec_out_private[thread_id][i] += conjugate(value_)*Vec_in[m_new];


                                } // if dn-dn hopping possible

                            }


                            if(sigma==0 && sigma_p==1){

                                //---------------Hopping: up to dn state-------------------//
                                //there have to be one up electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (bit_value(basis.D_up_basis[i], site)==1)
                                        &&
                                        (bit_value(basis.D_dn_basis[j], site_p)==0)

                                        )
                                {

                                    D_up = (int) (basis.D_up_basis[i] - pow(2, site)   );
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    //-----
                                    max_up = basis.Length -1;
                                    max_dn = basis.Length -1;
                                    min_dn = 0;


                                    sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                    if(l != max_up){
                                        sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                    }
                                    sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[j]);
                                    if(lp != min_dn){
                                        sign_pow_dn += bit_value(basis.D_dn_basis[j],min_dn);
                                    }

                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                                    //-----

                                    assert(m_new<i);
                                    double_type value_ = 1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one;
                                    Vec_out_private[thread_id][m_new] += value_*Vec_in[i];
                                    Vec_out_private[thread_id][i] += conjugate(value_)*Vec_in[m_new];


                                } // if up ---to---> dn hopping possible

                            }


                        }//if hopping matrix element is non-zero

                    }//site_p
                }//sigma_p

            } // site
        }//sigma
    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out[i] += Vec_out_private[thread][i];
        }
    }

}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_diagonal_terms_old(){

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;

    //Remember H[l][m]=<l|H|m>
    int m,j;
    double value;


    for (int i=0;i<basis.D_up_basis.size();i++){

        m=i;
        j=i;

        value=0;
        //intra-orbital coulomb repulsion:
        value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);


        //Long range density-density interaction: n_{i,s} n_{j,s'}
        for(int spin1=0;spin1<2;spin1++){
            for(int site1=0;site1<basis.Length;site1++){
                int alpha1 = basis.Length*spin1 + site1;
                int n_alpha1;
                if(spin1==0){
                    n_alpha1 = bit_value(basis.D_up_basis[i], site1);
                }
                else{
                    n_alpha1 = bit_value(basis.D_dn_basis[j], site1);
                }

                for(int spin2=0;spin2<2;spin2++){
                    for(int site2=0;site2<basis.Length;site2++){
                        int alpha2 = basis.Length*spin2 + site2;
                        int n_alpha2;
                        if(spin2==0){
                            n_alpha2 = bit_value(basis.D_up_basis[i], site2);
                        }
                        else{
                            n_alpha2 = bit_value(basis.D_dn_basis[j], site2);
                        }
                        value+=DenDenInt_mat_LongRange[alpha1][alpha2]*(n_alpha1*n_alpha2);
                    }
                }
            }
        }

        //Crystal Field Splitting (CFE):
        for(int site=0;site<basis.Length;site++){
            value+=(CFS[site])*
                    ( ( bit_value(basis.D_up_basis[i], site) +
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }

        //magnetic Field * [Sz]
        for(int site=0;site<basis.Length;site++){
            value+=0.5*(H_field[site])*
                    ( ( bit_value(basis.D_up_basis[i],site) -
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }


        if(value!=0){
            Hamil.value.push_back(value*one);
            Hamil.rows.push_back(m);
            Hamil.columns.push_back(m);
        }

    }

}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_diagonal_terms(){

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;

    //Remember H[l][m]=<l|H|m>
    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Matrix_COO> Hamil_private;
    Hamil_private.resize(N_threads);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int m,j;
    double value;
    int thread_id = 0;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){

        m=i;
        j=i;

        value=0;
        //intra-orbital coulomb repulsion:
        value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);


        //Long range density-density interaction: n_{i,s} n_{j,s'}
        for(int spin1=0;spin1<2;spin1++){
            for(int site1=0;site1<basis.Length;site1++){
                int alpha1 = basis.Length*spin1 + site1;
                int n_alpha1;
                if(spin1==0){
                    n_alpha1 = bit_value(basis.D_up_basis[i], site1);
                }
                else{
                    n_alpha1 = bit_value(basis.D_dn_basis[j], site1);
                }

                for(int spin2=0;spin2<2;spin2++){
                    for(int site2=0;site2<basis.Length;site2++){
                        int alpha2 = basis.Length*spin2 + site2;
                        int n_alpha2;
                        if(spin2==0){
                            n_alpha2 = bit_value(basis.D_up_basis[i], site2);
                        }
                        else{
                            n_alpha2 = bit_value(basis.D_dn_basis[j], site2);
                        }
                        value+=DenDenInt_mat_LongRange[alpha1][alpha2]*(n_alpha1*n_alpha2);
                    }
                }
            }
        }

        //Crystal Field Splitting (CFE):
        for(int site=0;site<basis.Length;site++){
            value+=(CFS[site])*
                    ( ( bit_value(basis.D_up_basis[i], site) +
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }

        //magnetic Field * [Sz]
        for(int site=0;site<basis.Length;site++){
            value+=0.5*(H_field[site])*
                    ( ( bit_value(basis.D_up_basis[i],site) -
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }


        if(value!=0){
            Hamil_private[thread_id].value.push_back(value*one);
            Hamil_private[thread_id].rows.push_back(m);
            Hamil_private[thread_id].columns.push_back(m);
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

}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_non_diagonal_terms(){
    //NOTHING FOR THIS MODEL
}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_connections_old(){


    double_type value;
    int alpha_, alpha_p_;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int max_up, max_dn, min_up, min_dn;
    int l,lp;
    complex<double> iota_ (0.0,1.0);




    for (int i=0;i<basis.D_up_basis.size();i++){
        //cout<<i<<" done"<<endl;
        m=i;
        j=i;


        for(int sigma=0;sigma<2 ;sigma++){
            for(int site=0;site<basis.Length ;site++){
                alpha_ = basis.Length*sigma + site;

                for(int sigma_p=0;sigma_p<2 ;sigma_p++){
                    for(int site_p=0;site_p<basis.Length ;site_p++){
                        alpha_p_ = basis.Length*sigma_p + site_p;


                        if(Hopping_mat_LongRange[alpha_p_][alpha_]!=zero){

                            //HOPPING COEFFICIENT IN FRONT OF
                            // (?) X c_{site_p,sigma_p}^{\dagger} c_{site,sigma}



                            if(sigma==0 && sigma_p==0){

                                //---------------Hopping: up to up electrons-------------------//
                                //there have to be one up electron on site
                                //there have to be no up electron on site_p
                                if(
                                        ( (bit_value(basis.D_up_basis[i],site)==1)
                                          &&
                                          (bit_value(basis.D_up_basis[i],site_p)==0)
                                          )
                                        &&
                                        (site_p<site)
                                        )
                                {


                                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                                  - pow(2, site) );
                                    D_dn = basis.D_dn_basis[j];

                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                    sign_FM = pow(-1.0, sign_pow_up);


                                    assert(m_new<m);
                                    Hamil.value.push_back(sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));


                                } // if up-up hopping possible

                            }


                            if(sigma==1 && sigma_p==1){

                                //---------------Hopping: dn to dn electrons-------------------//
                                //there have to be one dn electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (
                                            (bit_value(basis.D_dn_basis[j], site)==1)
                                            &&
                                            (bit_value(basis.D_dn_basis[j], site_p)==0)
                                            )
                                        &&
                                        (site_p<site)
                                        )
                                {

                                    D_up = basis.D_up_basis[i];
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                                  - pow(2, site) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn);

                                    assert(m_new<m);
                                    Hamil.value.push_back(1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));


                                } // if dn-dn hopping possible

                            }


                            if(sigma==0 && sigma_p==1){

                                //---------------Hopping: up to dn state-------------------//
                                //there have to be one up electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (bit_value(basis.D_up_basis[i], site)==1)
                                        &&
                                        (bit_value(basis.D_dn_basis[j], site_p)==0)

                                        )
                                {

                                    D_up = (int) (basis.D_up_basis[i] - pow(2, site)   );
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    //                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                    //                                sign_FM = pow(-1.0, sign_pow_dn);

                                    //-----
                                    max_up = basis.Length -1;
                                    min_up = 0;
                                    max_dn = basis.Length -1;
                                    min_dn = 0;


                                    sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                    if(l != max_up){
                                        sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                    }
                                    sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                    if(lp != min_dn){
                                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                    }

                                    //try this as well
                                    /*
                                            sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                            sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                            */

                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                                    //-----

                                    assert(m_new<m);
                                    Hamil.value.push_back(1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));


                                } // if up ---to---> dn hopping possible

                            }


                        }//if hopping matrix element is non-zero

                    }//site_p
                }//sigma_p

            } // site
        }//sigma

    } // "i" i.e up_decimals

}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_connections(){


   
    complex<double> iota_ (0.0,1.0);



    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Matrix_COO> Hamil_private;
    Hamil_private.resize(N_threads);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    double_type value;
    int alpha_, alpha_p_;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int max_up, max_dn, min_up, min_dn;
    int l,lp;
#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        //cout<<i<<" done"<<endl;
        m=i;
        j=i;


        for(int sigma=0;sigma<2 ;sigma++){
            for(int site=0;site<basis.Length ;site++){
                alpha_ = basis.Length*sigma + site;

                for(int sigma_p=0;sigma_p<2 ;sigma_p++){
                    for(int site_p=0;site_p<basis.Length ;site_p++){
                        alpha_p_ = basis.Length*sigma_p + site_p;


                        if(Hopping_mat_LongRange[alpha_p_][alpha_]!=zero){

                            //HOPPING COEFFICIENT IN FRONT OF
                            // (?) X c_{site_p,sigma_p}^{\dagger} c_{site,sigma}


                            //cout<<i<<"  "<<sigma<<"  "<<site<<"  "<<sigma_p<<"  "<<site_p<<"  "<<Hopping_mat_LongRange[alpha_p_][alpha_]<<endl;

                            if(sigma==0 && sigma_p==0){

                                //---------------Hopping: up to up electrons-------------------//
                                //there have to be one up electron on site
                                //there have to be no up electron on site_p
                                if(
                                        ( (bit_value(basis.D_up_basis[i],site)==1)
                                          &&
                                          (bit_value(basis.D_up_basis[i],site_p)==0)
                                          )
                                        &&
                                        (site_p<site)
                                        )
                                {


                                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                                  - pow(2, site) );
                                    D_dn = basis.D_dn_basis[j];

                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                    sign_FM = pow(-1.0, sign_pow_up);


                                    assert(m_new<m);
                                    Hamil_private[thread_id].value.push_back(sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil_private[thread_id].rows.push_back((m_new));
                                    Hamil_private[thread_id].columns.push_back((m));


                                } // if up-up hopping possible

                            }


                            if(sigma==1 && sigma_p==1){

                                //---------------Hopping: dn to dn electrons-------------------//
                                //there have to be one dn electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (
                                            (bit_value(basis.D_dn_basis[j], site)==1)
                                            &&
                                            (bit_value(basis.D_dn_basis[j], site_p)==0)
                                            )
                                        &&
                                        (site_p<site)
                                        )
                                {

                                    D_up = basis.D_up_basis[i];
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                                  - pow(2, site) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn);

                                    assert(m_new<m);
                                    Hamil_private[thread_id].value.push_back(1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil_private[thread_id].rows.push_back((m_new));
                                    Hamil_private[thread_id].columns.push_back((m));


                                } // if dn-dn hopping possible

                            }


                            if(sigma==0 && sigma_p==1){

                                //---------------Hopping: up to dn state-------------------//
                                //there have to be one up electron on site
                                //there have to be no dn electron on site_p
                                if(
                                        (bit_value(basis.D_up_basis[i], site)==1)
                                        &&
                                        (bit_value(basis.D_dn_basis[j], site_p)==0)

                                        )
                                {

                                    D_up = (int) (basis.D_up_basis[i] - pow(2, site)   );
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p) );


                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;


                                    l= site;
                                    lp= site_p;

                                    //                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                    //                                sign_FM = pow(-1.0, sign_pow_dn);

                                    //-----
                                    max_up = basis.Length -1;
                                    min_up = 0;
                                    max_dn = basis.Length -1;
                                    min_dn = 0;


                                    sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                    if(l != max_up){
                                        sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                    }
                                    sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                    if(lp != min_dn){
                                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                    }

                                    //try this as well
                                    /*
                                            sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                            sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                            */

                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                                    //-----

                                    assert(m_new<m);
                                    Hamil_private[thread_id].value.push_back(1.0*sign_FM*(Hopping_mat_LongRange[alpha_p_][alpha_])*one);
                                    Hamil_private[thread_id].rows.push_back((m_new));
                                    Hamil_private[thread_id].columns.push_back((m));


                                } // if up ---to---> dn hopping possible

                            }


                        }//if hopping matrix element is non-zero

                    }//site_p
                }//sigma_p

            } // site
        }//sigma

    } // "i" i.e up_decimals

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        Hamil.value.insert(Hamil.value.end(),Hamil_private[thread].value.begin(), Hamil_private[thread].value.end() );
        Hamil.rows.insert(Hamil.rows.end(),Hamil_private[thread].rows.begin(), Hamil_private[thread].rows.end() );
        Hamil.columns.insert(Hamil.columns.end(),Hamil_private[thread].columns.begin(), Hamil_private[thread].columns.end() );
    }

}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Read_parameters(string filename){


    string filepath = filename;


    double temp_val;
    string length, Length = "Length = ";
    string ntotal, Ntotal = "N_Total = ";

    string ucoul, Ucoul = "U = ";
    string hmag_site_resolved, Hmag_site_resolved = "H_mag_site_resolved = ";
    string cfs_site_resolved, CFS_SITE_RESOLVED = "CFS_SITE_RESOLVED = ";

    string LongRangeHoppingfile_ = "LongRangeHopping_file = ";
    string LongRangeDenDenIntfile_ = "LongRangeDenDenInt_file = ";
    string BipartiteEntanglement_ = "Calculate_Bipartite_Entanglement = ";
    string Sys1_sites_ = "Sys1 = ";
    string Sys2_sites_ = "Sys2 = ";

    string FourPointObsSet_file_ = "FourPointObsSet_file = ";

    string bipartite_entanglement_;
    string sys1_sites_string_;
    string sys2_sites_string_;

    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Length, 0)) != string::npos) {
                length = line.substr (offset + Length.length());		}

            if ((offset = line.find(Ntotal, 0)) != string::npos) {
                ntotal = line.substr (offset + Ntotal.length());		}

            if ((offset = line.find(Ucoul, 0)) != string::npos) {
                ucoul= line.substr (offset + Ucoul.length());		}

            if ((offset = line.find(Hmag_site_resolved, 0)) != string::npos) {
                hmag_site_resolved = line.substr (offset + Hmag_site_resolved.length());		}

            if ((offset = line.find(CFS_SITE_RESOLVED, 0)) != string::npos) {
                cfs_site_resolved = line.substr (offset+CFS_SITE_RESOLVED.length());				}

            if ((offset = line.find(LongRangeHoppingfile_, 0)) != string::npos) {
                LongRangeHoppingfilepath = line.substr (offset+LongRangeHoppingfile_.length());  }

            if ((offset = line.find(LongRangeDenDenIntfile_, 0)) != string::npos) {
                LongRangeDenDenIntfilepath = line.substr (offset+LongRangeDenDenIntfile_.length());  }

            if ((offset = line.find(BipartiteEntanglement_, 0)) != string::npos) {
                bipartite_entanglement_ = line.substr (offset+BipartiteEntanglement_.length());  }

            if ((offset = line.find(Sys1_sites_, 0)) != string::npos) {
                sys1_sites_string_ = line.substr (offset+Sys1_sites_.length());  }

            if ((offset = line.find(Sys2_sites_, 0)) != string::npos) {
                sys2_sites_string_ = line.substr (offset+Sys2_sites_.length());  }

            if ((offset = line.find(FourPointObsSet_file_, 0)) != string::npos) {
                FourPointObsSet_filepath = line.substr (offset+FourPointObsSet_file_.length());  }

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}



    basis.Length=atoi(length.c_str());
    cout<<"Length = "<<basis.Length<<endl;

    basis.N_total=atoi(ntotal.c_str());
    cout<<"N_total = "<<basis.N_total<<endl;

    U=atof(ucoul.c_str());
    cout<<"U = "<<U<<endl;

    Calculate_Bipartite_Entanglement = false;
    Sys1_Ent_sites.clear();
    Sys2_Ent_sites.clear();

    if(bipartite_entanglement_=="true"){
        Calculate_Bipartite_Entanglement = true;
    }

    if(sys1_sites_string_.size()>0){
        stringstream sys1_stream(sys1_sites_string_);
        int site_temp;
        while(sys1_stream >> site_temp){
            Sys1_Ent_sites.push_back(site_temp);
        }
    }
    if(sys2_sites_string_.size()>0){
        stringstream sys2_stream(sys2_sites_string_);
        int site_temp;
        while(sys2_stream >> site_temp){
            Sys2_Ent_sites.push_back(site_temp);
        }
    }

    //if(Sys1_Ent_sites.size()>0 && Sys2_Ent_sites.size()>0){
    //    Calculate_Bipartite_Entanglement = true;
   // }

    if(Calculate_Bipartite_Entanglement){
        cout<<"Bipartite entanglement calculation enabled"<<endl;
        cout<<"Sys1 = ";
        for(int i=0;i<Sys1_Ent_sites.size();i++){
            cout<<Sys1_Ent_sites[i]<<" ";
        }
        cout<<endl;
        cout<<"Sys2 = ";
        for(int i=0;i<Sys2_Ent_sites.size();i++){
            cout<<Sys2_Ent_sites[i]<<" ";
        }
        cout<<endl;
    }

    //double h;
    string temp_string_h;
    stringstream h_site_resolved_stream;
    h_site_resolved_stream<<hmag_site_resolved;
    h_site_resolved_stream>>temp_string_h;

    if(temp_string_h == "true"){
        H_MAG_SITE_RESOLVED_bool = true;
    }
    else{
        H_MAG_SITE_RESOLVED_bool = false;
    }


    H_field.clear();
    H_field.resize(basis.Length);
     for(int i=0;i<basis.Length;i++){
        if(H_MAG_SITE_RESOLVED_bool==true){

            h_site_resolved_stream >> temp_val;
            H_field[i]=temp_val;
        }
        else{
            H_field[i]=0.0;
        }
    }


    cout<<"H_mag_site_resolved = "<<endl;
    for(int i=0;i<basis.Length;i++){
        cout<<H_field[i]<<" ";
    }
    cout<<endl;


    string temp_string;
    stringstream cfs_site_resolved_stream;
    cfs_site_resolved_stream<<cfs_site_resolved;
    cfs_site_resolved_stream>>temp_string;

    if(temp_string == "true"){
        CFS_SITE_RESOLVED_bool = true;
    }
    else{
        CFS_SITE_RESOLVED_bool = false;
    }


    CFS.clear();
    CFS.resize(basis.Length);


    for(int i=0;i<basis.Length;i++){
        if(CFS_SITE_RESOLVED_bool==true){

            cfs_site_resolved_stream >> temp_val;
            CFS[i]=temp_val;
        }
        else{
            CFS[i]=0.0;
        }

    }

    cout<<"CFS read from inputfile"<<endl;

    cout<<"Reading hopping matrix from : "<<LongRangeHoppingfilepath<<endl;

    Read_matrix_from_file(LongRangeHoppingfilepath, Hopping_mat_LongRange , 2*basis.Length, 2*basis.Length);

    cout<<"Hopping matrix read from file : "<<endl;
    Print_Matrix(Hopping_mat_LongRange);




    cout<<"Reading DenDenInt matrix from : "<<LongRangeDenDenIntfilepath<<endl;

    Read_matrix_from_file(LongRangeDenDenIntfilepath, DenDenInt_mat_LongRange , 2*basis.Length, 2*basis.Length);


    cout<<"DenDenInt matrix read from file : "<<endl;
    Print_Matrix(DenDenInt_mat_LongRange);




    ifstream FourPointInfile(FourPointObsSet_filepath.c_str());
    fourpointSitesSet.clear();
    fourpointSpinsSet.clear();
    fourpointValuesSet.clear();
    if(FourPointInfile.is_open())
    {
        string line;
        while(getline(FourPointInfile, line))
        {
            stringstream line_stream(line);
            string first_token;
            line_stream >> first_token;

            if(first_token.size()==0 || first_token[0]=='#'){
                continue;
            }

            int no_of_oprs = atoi(first_token.c_str());
            Mat_1_tetra_int sites_set_temp;
            Mat_1_tetra_int spins_set_temp;
            Mat_1_doub values_set_temp;

            for(int opr_no=0;opr_no<no_of_oprs;opr_no++){
                string coeff_token;
                int site1, spin1, site2, spin2, site3, spin3, site4, spin4;
                tetra_int sites_temp;
                tetra_int spins_temp;

                line_stream >> coeff_token;
                line_stream >> site1 >> spin1 >> site2 >> spin2 >> site3 >> spin3 >> site4 >> spin4;

                if(line_stream.fail()){
                    cout<<"Malformed four point observable entry: "<<line<<endl;
                    break;
                }

                sites_temp.first = site1;
                sites_temp.second = site2;
                sites_temp.third = site3;
                sites_temp.fourth = site4;

                spins_temp.first = spin1;
                spins_temp.second = spin2;
                spins_temp.third = spin3;
                spins_temp.fourth = spin4;

                sites_set_temp.push_back(sites_temp);
                spins_set_temp.push_back(spins_temp);
                values_set_temp.push_back(Read_FourPoint_Coefficient(coeff_token));
            }

            if(sites_set_temp.size()>0){
                fourpointSitesSet.push_back(sites_set_temp);
                fourpointSpinsSet.push_back(spins_set_temp);
                fourpointValuesSet.push_back(values_set_temp);
            }
        }
        FourPointInfile.close();
    }
    else    {cout<<"Unable to open input file for four point obs set."<<endl;}

    cout<<"No. of four point observable sets read = "<<fourpointSitesSet.size()<<endl;


    //THINK ABOUT IT LATER :)
    /*for(int site=0;site<basis.Length ;site++){
        for(int site_p=0;site_p<basis.Length ;site_p++){
            if(site_p>=site){
                Hopping_mat[site_p][site]=zero;
            }

        }}
        */

    cout<<"PARAMETERS READ"<<endl;


}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Read_parameters_for_dynamics(string filename){

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


//NOT WORKING PROPERLY AT PRESENT, THINK ABOUT IT LATER
template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Create_Lattice_Graph(string output_filename){

    ofstream graph_file(output_filename.c_str());

    if(!graph_file.is_open()){
        cout<<"Unable to open lattice graph output file."<<endl;
        return;
    }

    const double canvas_x = 900.0;
    const double canvas_y = 900.0;
    const double center_x = canvas_x*0.5;
    const double center_y = canvas_y*0.5;
    const double radius = 0.36*((canvas_x<canvas_y)?canvas_x:canvas_y);
    const double node_radius = 18.0;

    Mat_1_real site_x, site_y;
    site_x.resize(basis.Length);
    site_y.resize(basis.Length);

    if(basis.Length==1){
        site_x[0]=center_x;
        site_y[0]=center_y;
    }
    else{
        for(int site=0;site<basis.Length;site++){
            double angle = (2.0*PI*site)/(1.0*basis.Length);
            site_x[site] = center_x + radius*cos(angle);
            site_y[site] = center_y + radius*sin(angle);
        }
    }

    double max_strength = 0.0;
    for(int site1=0;site1<basis.Length;site1++){
        for(int site2=site1+1;site2<basis.Length;site2++){
            for(int spin_from=0;spin_from<2;spin_from++){
                for(int spin_to=0;spin_to<2;spin_to++){
                    int alpha_12 = basis.Length*spin_from + site1;
                    int alpha_21 = basis.Length*spin_to + site2;

                    double strength1 = GetDoubleTypeMagnitude(Hopping_mat_LongRange[alpha_21][alpha_12]);
                    double strength2 = GetDoubleTypeMagnitude(Hopping_mat_LongRange[alpha_12][alpha_21]);

                    if(strength1>max_strength){max_strength=strength1;}
                    if(strength2>max_strength){max_strength=strength2;}
                }
            }
        }
    }

    if(max_strength==0.0){
        max_strength=1.0;
    }

    graph_file<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<canvas_x
              <<"\" height=\""<<canvas_y<<"\" viewBox=\"0 0 "<<canvas_x<<" "<<canvas_y<<"\">"<<endl;
    graph_file<<"<rect x=\"0\" y=\"0\" width=\""<<canvas_x<<"\" height=\""<<canvas_y
              <<"\" fill=\"white\"/>"<<endl;
    graph_file<<"<text x=\""<<(0.5*canvas_x)<<"\" y=\"40\" text-anchor=\"middle\" font-size=\"24\" font-family=\"Arial\">"
              <<"Lattice graph from Hopping_mat_LongRange</text>"<<endl;

    for(int site1=0;site1<basis.Length;site1++){
        for(int site2=site1+1;site2<basis.Length;site2++){
            double_type uu_12 = Hopping_mat_LongRange[site2][site1];
            double_type uu_21 = Hopping_mat_LongRange[site1][site2];
            double_type dd_12 = Hopping_mat_LongRange[basis.Length + site2][basis.Length + site1];
            double_type dd_21 = Hopping_mat_LongRange[basis.Length + site1][basis.Length + site2];
            double_type ud_12 = Hopping_mat_LongRange[basis.Length + site2][site1];
            double_type ud_21 = Hopping_mat_LongRange[site1][basis.Length + site2];
            double_type du_12 = Hopping_mat_LongRange[site2][basis.Length + site1];
            double_type du_21 = Hopping_mat_LongRange[basis.Length + site1][site2];

            double edge_strength = 0.0;
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(uu_12));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(uu_21));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(dd_12));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(dd_21));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(ud_12));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(ud_21));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(du_12));
            edge_strength = max(edge_strength, GetDoubleTypeMagnitude(du_21));

            if(edge_strength==0.0){
                continue;
            }

            double_type dominant_hopping = uu_12;
            if(GetDoubleTypeMagnitude(uu_21) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = uu_21;
            }
            if(GetDoubleTypeMagnitude(dd_12) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = dd_12;
            }
            if(GetDoubleTypeMagnitude(dd_21) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = dd_21;
            }
            if(GetDoubleTypeMagnitude(ud_12) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = ud_12;
            }
            if(GetDoubleTypeMagnitude(ud_21) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = ud_21;
            }
            if(GetDoubleTypeMagnitude(du_12) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = du_12;
            }
            if(GetDoubleTypeMagnitude(du_21) > GetDoubleTypeMagnitude(dominant_hopping)){
                dominant_hopping = du_21;
            }

            double signed_value = GetDoubleTypeSignedPart(dominant_hopping);

            bool spin_conserving = ((GetDoubleTypeMagnitude(uu_12)!=0.0) || (GetDoubleTypeMagnitude(uu_21)!=0.0)
                                    || (GetDoubleTypeMagnitude(dd_12)!=0.0) || (GetDoubleTypeMagnitude(dd_21)!=0.0));
            bool spin_flip = ((GetDoubleTypeMagnitude(ud_12)!=0.0) || (GetDoubleTypeMagnitude(ud_21)!=0.0)
                              || (GetDoubleTypeMagnitude(du_12)!=0.0) || (GetDoubleTypeMagnitude(du_21)!=0.0));

            string edge_color = "#6c757d";
            string dash_style = "none";
            if(signed_value > 0.0){
                edge_color = "#1d4ed8";
            }
            if(signed_value < 0.0){
                edge_color = "#dc2626";
            }

            if(spin_conserving && spin_flip){
                dash_style = "8,5";
            }
            else if(spin_flip){
                dash_style = "8,5";
            }

            double normalized_strength = edge_strength/max_strength;
            double line_width = 1.0 + 10.0*normalized_strength;

            graph_file<<"<line x1=\""<<site_x[site1]<<"\" y1=\""<<site_y[site1]
                      <<"\" x2=\""<<site_x[site2]<<"\" y2=\""<<site_y[site2]
                      <<"\" stroke=\""<<edge_color<<"\" stroke-width=\""<<line_width<<"\" ";
            if(dash_style!="none"){
                graph_file<<"stroke-dasharray=\""<<dash_style<<"\" ";
            }
            graph_file<<">"<<endl;

            graph_file<<"<title>sites "<<site1<<" and "<<site2
                      <<"; uu_21="<<DoubleTypeToString(uu_12)
                      <<"; uu_12="<<DoubleTypeToString(uu_21)
                      <<"; dd_21="<<DoubleTypeToString(dd_12)
                      <<"; dd_12="<<DoubleTypeToString(dd_21)
                      <<"; ud_21="<<DoubleTypeToString(ud_12)
                      <<"; ud_12="<<DoubleTypeToString(ud_21)
                      <<"; du_21="<<DoubleTypeToString(du_12)
                      <<"; du_12="<<DoubleTypeToString(du_21)
                      <<"; dominant="<<DoubleTypeToString(dominant_hopping)
                      <<"</title>"<<endl;
            graph_file<<"</line>"<<endl;
        }
    }

    for(int site=0;site<basis.Length;site++){
        graph_file<<"<circle cx=\""<<site_x[site]<<"\" cy=\""<<site_y[site]
                  <<"\" r=\""<<node_radius<<"\" fill=\"#cfe8ff\" stroke=\"#1f4e79\" stroke-width=\"2\"/>"<<endl;
        graph_file<<"<text x=\""<<site_x[site]<<"\" y=\""<<(site_y[site] + 5.0)
                  <<"\" text-anchor=\"middle\" font-size=\"16\" font-family=\"Arial\" fill=\"#102a43\">"
                  <<site<<"</text>"<<endl;
    }

    graph_file<<"<rect x=\"25\" y=\""<<(canvas_y-125.0)<<"\" width=\"250\" height=\"90\" fill=\"#fbfbfb\" stroke=\"#bbbbbb\"/>"<<endl;
    graph_file<<"<text x=\"40\" y=\""<<(canvas_y-95.0)<<"\" font-size=\"16\" font-family=\"Arial\">Legend</text>"<<endl;
    graph_file<<"<line x1=\"45\" y1=\""<<(canvas_y-70.0)<<"\" x2=\"105\" y2=\""<<(canvas_y-70.0)
              <<"\" stroke=\"#1d4ed8\" stroke-width=\"3\"/>"<<endl;
    graph_file<<"<text x=\"115\" y=\""<<(canvas_y-64.0)<<"\" font-size=\"14\" font-family=\"Arial\">positive hopping</text>"<<endl;
    graph_file<<"<line x1=\"45\" y1=\""<<(canvas_y-45.0)<<"\" x2=\"105\" y2=\""<<(canvas_y-45.0)
              <<"\" stroke=\"#dc2626\" stroke-width=\"3\"/>"<<endl;
    graph_file<<"<text x=\"115\" y=\""<<(canvas_y-39.0)<<"\" font-size=\"14\" font-family=\"Arial\">negative hopping</text>"<<endl;
    graph_file<<"<text x=\"40\" y=\""<<(canvas_y-18.0)<<"\" font-size=\"14\" font-family=\"Arial\">dashed lines indicate spin-flip terms; width tracks magnitude</text>"<<endl;

    graph_file<<"</svg>"<<endl;
    graph_file.close();

    cout<<"Lattice graph written to "<<output_filename<<endl;

}






template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_BipartiteEntanglement(Mat_1_int &Sys1_, Mat_1_int &Sys2_,
                                                                    Mat_1_doub &Vec_, double & VonNuemannEntropy){

    assert((int)Vec_.size()==(int)basis.D_up_basis.size());
    assert((int)Sys1_.size() + (int)Sys2_.size() == basis.Length);

    Mat_1_int part_check;
    part_check.resize(basis.Length);
    for(int site=0;site<basis.Length;site++){
        part_check[site]=0;
    }
    for(int s1=0;s1<Sys1_.size();s1++){
        assert(Sys1_[s1]>=0 && Sys1_[s1]<basis.Length);
        part_check[Sys1_[s1]] += 1;
    }
    for(int s2=0;s2<Sys2_.size();s2++){
        assert(Sys2_[s2]>=0 && Sys2_[s2]<basis.Length);
        part_check[Sys2_[s2]] += 1;
    }
    for(int site=0;site<basis.Length;site++){
        assert(part_check[site]==1);
    }

    int len1 = Sys1_.size();
    int len2 = Sys2_.size();

    int dim_up_1 = (int)pow(2, len1);
    int dim_dn_1 = (int)pow(2, len1);
    int dim1 = dim_up_1*dim_dn_1;

    int dim_up_2 = (int)pow(2, len2);
    int dim_dn_2 = (int)pow(2, len2);
    int dim2 = dim_up_2*dim_dn_2;

    Matrix<double_type> Psi_LB;
    Psi_LB.resize(dim1, dim2);
    for(int i1=0;i1<dim1;i1++){
        for(int i2=0;i2<dim2;i2++){
            Psi_LB(i1,i2)=zero;
        }
    }

    for(int m=0;m<basis.D_up_basis.size();m++){
        int D_up = basis.D_up_basis[m];
        int D_dn = basis.D_dn_basis[m];

        int up1_dec = 0;
        int dn1_dec = 0;
        int up2_dec = 0;
        int dn2_dec = 0;

        for(int x=0;x<len1;x++){
            int site = Sys1_[x];
            up1_dec += bit_value(D_up, site)*(int)pow(2,x);
            dn1_dec += bit_value(D_dn, site)*(int)pow(2,x);
        }

        for(int x=0;x<len2;x++){
            int site = Sys2_[x];
            up2_dec += bit_value(D_up, site)*(int)pow(2,x);
            dn2_dec += bit_value(D_dn, site)*(int)pow(2,x);
        }

        int row = up1_dec + dim_up_1*dn1_dec;
        int col = up2_dec + dim_up_2*dn2_dec;

        Psi_LB(row,col) += Vec_[m];
    }

    int r_ = min(dim1,dim2);
    Matrix<double_type> VT_;
    Matrix<double_type> U_;
    vector<double> Sigma_;
    Perform_SVD(Psi_LB,VT_,U_,Sigma_);

    double VonNuemannEntropy_temp=0.0;
    for(int n=0;n<r_ && n<Sigma_.size();n++){
        double prob_n = Sigma_[n]*Sigma_[n];
        if(prob_n>1.0e-16){
            VonNuemannEntropy_temp += -1.0*prob_n*(log2(prob_n));
        }
    }

    VonNuemannEntropy = VonNuemannEntropy_temp;
}




template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_two_point_observables(Mat_1_doub &Vec_){

    Matrix_COO OPR_;
    Mat_1_doub VecL_, VecR_;
    double_type value_;

    Mat_2_doub AMat0; //[sigma][sigma_p]
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
    }


    int TOTAL_NO_OBS=5;
    Mat_3_doub AMatL, AMatR;
    AMatL.resize(TOTAL_NO_OBS);
    AMatR.resize(TOTAL_NO_OBS);

    two_point_obs.resize(TOTAL_NO_OBS);

    two_point_obs[0]="<n[i].n[j]>";
    AMatL[0]=AMat0;AMatR[0]=AMat0;
    AMatL[0][0][0]=one;AMatL[0][1][1]=one;
    AMatR[0][0][0]=one;AMatR[0][1][1]=one;

    two_point_obs[1]="<Sz[i].Sz[j]>";
    AMatL[1]=AMat0;AMatR[1]=AMat0;
    AMatL[1][0][0]=one*(0.5);AMatL[1][1][1]=one*(-0.5);
    AMatR[1][0][0]=one*(0.5);AMatR[1][1][1]=one*(-0.5);

    two_point_obs[2]="<Splus[i].Sminus[j]>";
    AMatL[2]=AMat0;AMatR[2]=AMat0;
    AMatL[2][1][0]=one; //S-
    AMatR[2][1][0]=one; //S-

    two_point_obs[3]="<Sminus[i].Splus[j]>";
    AMatL[3]=AMat0;AMatR[3]=AMat0;
    AMatL[3][0][1]=one; //S+
    AMatR[3][0][1]=one; //S+

    two_point_obs[4]="<nup[i].ndn[j]>";
    AMatL[4]=AMat0;AMatR[4]=AMat0;
    AMatL[4][0][0]=one;
    AMatR[4][1][1]=one;



    double_type sum_;
    for(int obs_no=0;obs_no<TOTAL_NO_OBS;obs_no++){
        cout<<"--------------"<<two_point_obs[obs_no]<<"-------------------"<<endl;


        sum_=zero;
        for(int siteL=0;siteL<basis.Length;siteL++){
            Get_CdaggerC_type_Opr(AMatL[obs_no], OPR_, siteL);
            Matrix_COO_vector_multiplication("cx", OPR_, Vec_, VecL_);

            vector< int >().swap( OPR_.columns );
            vector< int >().swap( OPR_.rows );
            vector< double_type >().swap( OPR_.value );


            for(int siteR=0;siteR<basis.Length;siteR++){
                Get_CdaggerC_type_Opr(AMatR[obs_no], OPR_, siteR);
                Matrix_COO_vector_multiplication("cx", OPR_, Vec_, VecR_);

                value_ = dot_product(VecR_,VecL_);
                sum_ += value_;

                cout<<value_<<"  ";

                vector< int >().swap( OPR_.columns );
                vector< int >().swap( OPR_.rows );
                vector< double_type >().swap( OPR_.value );

            }
            cout<<endl;
        }

        cout<<"-------------------------------------------------------"<<endl;
        cout<<"sum = "<<sum_<<endl<<endl;
        cout<<endl;


    }



}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_one_point_observables(Mat_1_doub &Vec_){


    Matrix_COO OPR_;
    Mat_1_doub Vec_temp_;
    double_type value_;

    Mat_2_doub AMat0; //[sigma][sigma_p]
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
    }


    int TOTAL_NO_OBS=3;
    Mat_3_doub AMat;
    AMat.resize(TOTAL_NO_OBS);

    one_point_obs.resize(TOTAL_NO_OBS);
    one_point_obs[0]="n_up";
    AMat[0]=AMat0;
    AMat[0][0][0]=one;

    one_point_obs[1]="n_dn";
    AMat[1]=AMat0;
    AMat[1][1][1]=one;

    one_point_obs[2]="S_plus";
    AMat[2]=AMat0;
    AMat[2][0][1]=one;


    for(int obs_no=0;obs_no<one_point_obs.size();obs_no++){

        for(int site=0;site<basis.Length;site++){

            Get_CdaggerC_type_Opr(AMat[obs_no], OPR_, site);
            Matrix_COO_vector_multiplication("cx", OPR_, Vec_, Vec_temp_);

            value_ = dot_product(Vec_temp_,Vec_);
#ifdef USE_COMPLEX
            cout<<one_point_obs[obs_no]<<"["<<site<<"] = "<<value_.real() << "  "<<value_.imag()<<endl;
#endif
#ifndef USE_COMPLEX
            cout<<one_point_obs[obs_no]<<"["<<site<<"] = "<<value_<<endl;
#endif

        }

        cout<<endl;

        vector< int >().swap( OPR_.columns );
        vector< int >().swap( OPR_.rows );
        vector< double_type >().swap( OPR_.value );
    }

}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_four_point_observables(Mat_1_doub &Vec_){

    Matrix_COO OPR1_, OPR2_;
    Mat_1_doub Vec_temp_, Vec_temp2_,Vec_final_;
    double_type value_, value_sum_;
    double_type value_sum_quantum;

    double_type Value1_class, Value2_class;

    double_type Total_Value_Sum=zero;
    double_type Total_Value_Sum_quantum=zero;

    Mat_2_doub AMat0;
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
        for(int ind2=0;ind2<2;ind2++){
            AMat0[ind][ind2]=zero;
        }
    }

    cout<<"-------------- <cdag c cdag c> and <cdagc><cdagc>  and <cdag c cdag c> - <cdagc><cdagc>-------------------"<<endl;

    for(int set_no=0;set_no<fourpointSitesSet.size();set_no++){
        assert(fourpointSitesSet[set_no].size()==fourpointSpinsSet[set_no].size());
        assert(fourpointSitesSet[set_no].size()==fourpointValuesSet[set_no].size());

        value_sum_=zero;

        cout<<"Set = "<<set_no<<endl;
        for(int term_no=0;term_no<fourpointSitesSet[set_no].size();term_no++){
            Mat_2_doub AMat1, AMat2;
            tetra_int sites_;
            tetra_int spins_;

            AMat1=AMat0;
            AMat2=AMat0;

            sites_ = fourpointSitesSet[set_no][term_no];
            spins_ = fourpointSpinsSet[set_no][term_no];

            AMat1[spins_.first][spins_.second]=one;
            AMat2[spins_.third][spins_.fourth]=one;

            Get_CdaggerC_type_Opr(AMat2, OPR2_, sites_.third, sites_.fourth);
            Matrix_COO_vector_multiplication("cx", OPR2_, Vec_, Vec_temp_);


            Get_CdaggerC_type_Opr(AMat1, OPR1_, sites_.first, sites_.second);
            Matrix_COO_vector_multiplication("cx", OPR1_, Vec_temp_, Vec_final_);

            //For classical <cdag1 c2>
            Matrix_COO_vector_multiplication("cx", OPR1_, Vec_, Vec_temp2_);


            value_ = fourpointValuesSet[set_no][term_no]*dot_product(Vec_final_, Vec_);


            Value2_class = dot_product(Vec_temp_, Vec_);
            Value1_class = dot_product(Vec_temp2_, Vec_);
            

            value_sum_ += value_;

            value_sum_quantum +=  value_ - (fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class);

            cout<<"term="<<term_no<<"  coeff="<<fourpointValuesSet[set_no][term_no]
                <<"  sites=("<<sites_.first<<","<<sites_.second<<","<<sites_.third<<","<<sites_.fourth<<")"
                <<"  spins=("<<spins_.first<<","<<spins_.second<<","<<spins_.third<<","<<spins_.fourth<<")"
                <<"  value="<<value_
                <<"  value_classical="<<fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class
                <<"  value_quantum="<<value_ - (fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class)
                <<endl;

            vector< int >().swap( OPR1_.columns );
            vector< int >().swap( OPR1_.rows );
            vector< double_type >().swap( OPR1_.value );

            vector< int >().swap( OPR2_.columns );
            vector< int >().swap( OPR2_.rows );
            vector< double_type >().swap( OPR2_.value );

            vector< double_type >().swap( Vec_temp_ );
            vector< double_type >().swap( Vec_temp2_ );
            vector< double_type >().swap( Vec_final_ );
        }

        cout<<"Total for set "<<set_no<<" = "<<value_sum_<<endl;
        cout<<"Total quantum for set "<<set_no<<" = "<<value_sum_quantum<<endl;
        cout<<endl;
        Total_Value_Sum += value_sum_;
        Total_Value_Sum_quantum += value_sum_quantum;
    }

    cout<<"------------------------------------------------"<<endl;
    cout<<"Total Value Sum = "<<Total_Value_Sum<<endl;
    cout<<"Total Value Sum Quantum = "<<Total_Value_Sum_quantum<<endl;

}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_two_point_observables_acting(Mat_1_doub &Vec_){

    Mat_1_doub VecL_, VecR_;
    double_type value_;

    Mat_2_doub AMat0; //[sigma][sigma_p]
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
    }


    int TOTAL_NO_OBS=5;
    Mat_3_doub AMatL, AMatR;
    AMatL.resize(TOTAL_NO_OBS);
    AMatR.resize(TOTAL_NO_OBS);

    two_point_obs.resize(TOTAL_NO_OBS);

    two_point_obs[0]="<n[i].n[j]>";
    AMatL[0]=AMat0;AMatR[0]=AMat0;
    AMatL[0][0][0]=one;AMatL[0][1][1]=one;
    AMatR[0][0][0]=one;AMatR[0][1][1]=one;

    two_point_obs[1]="<Sz[i].Sz[j]>";
    AMatL[1]=AMat0;AMatR[1]=AMat0;
    AMatL[1][0][0]=one*(0.5);AMatL[1][1][1]=one*(-0.5);
    AMatR[1][0][0]=one*(0.5);AMatR[1][1][1]=one*(-0.5);

    two_point_obs[2]="<Splus[i].Sminus[j]>";
    AMatL[2]=AMat0;AMatR[2]=AMat0;
    AMatL[2][1][0]=one; //S-
    AMatR[2][1][0]=one; //S-

    two_point_obs[3]="<Sminus[i].Splus[j]>";
    AMatL[3]=AMat0;AMatR[3]=AMat0;
    AMatL[3][0][1]=one; //S+
    AMatR[3][0][1]=one; //S+

    two_point_obs[4]="<nup[i].ndn[j]>";
    AMatL[4]=AMat0;AMatR[4]=AMat0;
    AMatL[4][0][0]=one;
    AMatR[4][1][1]=one;



    double_type sum_;
    for(int obs_no=0;obs_no<TOTAL_NO_OBS;obs_no++){
        cout<<"--------------"<<two_point_obs[obs_no]<<"-------------------"<<endl;


        sum_=zero;
        for(int siteL=0;siteL<basis.Length;siteL++){
            Get_CdaggerC_type_Opr(AMatL[obs_no], Vec_, VecL_, siteL);


            for(int siteR=0;siteR<basis.Length;siteR++){
                Get_CdaggerC_type_Opr(AMatR[obs_no], Vec_, VecR_, siteR);

                value_ = dot_product(VecR_,VecL_);
                sum_ += value_;

                cout<<value_<<"  ";

            }
            cout<<endl;
        }

        cout<<"-------------------------------------------------------"<<endl;
        cout<<"sum = "<<sum_<<endl<<endl;
        cout<<endl;


    }



}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_one_point_observables_acting(Mat_1_doub &Vec_){


    Mat_1_doub Vec_temp_;
    double_type value_;

    Mat_2_doub AMat0; //[sigma][sigma_p]
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
    }


    int TOTAL_NO_OBS=3;
    Mat_3_doub AMat;
    AMat.resize(TOTAL_NO_OBS);

    one_point_obs.resize(TOTAL_NO_OBS);
    one_point_obs[0]="n_up";
    AMat[0]=AMat0;
    AMat[0][0][0]=one;

    one_point_obs[1]="n_dn";
    AMat[1]=AMat0;
    AMat[1][1][1]=one;

    one_point_obs[2]="S_plus";
    AMat[2]=AMat0;
    AMat[2][0][1]=one;


    for(int obs_no=0;obs_no<one_point_obs.size();obs_no++){

        for(int site=0;site<basis.Length;site++){

            Get_CdaggerC_type_Opr(AMat[obs_no], Vec_, Vec_temp_, site);

            value_ = dot_product(Vec_temp_,Vec_);
#ifdef USE_COMPLEX
            cout<<one_point_obs[obs_no]<<"["<<site<<"] = "<<value_.real() << "  "<<value_.imag()<<endl;
#endif
#ifndef USE_COMPLEX
            cout<<one_point_obs[obs_no]<<"["<<site<<"] = "<<value_<<endl;
#endif

        }

        cout<<endl;
    }

}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Calculate_four_point_observables_acting(Mat_1_doub &Vec_){

    Mat_1_doub Vec_temp_, Vec_temp2_,Vec_final_;
    double_type value_, value_sum_;
    double_type value_sum_quantum;

    double_type Value1_class, Value2_class;

    double_type Total_Value_Sum=zero;
    double_type Total_Value_Sum_quantum=zero;

    Mat_2_doub AMat0;
    AMat0.resize(2);
    for(int ind=0;ind<2;ind++){
        AMat0[ind].resize(2);
        for(int ind2=0;ind2<2;ind2++){
            AMat0[ind][ind2]=zero;
        }
    }

    cout<<"-------------- <cdag c cdag c> and <cdagc><cdagc>  and <cdag c cdag c> - <cdagc><cdagc>-------------------"<<endl;

    for(int set_no=0;set_no<fourpointSitesSet.size();set_no++){
        assert(fourpointSitesSet[set_no].size()==fourpointSpinsSet[set_no].size());
        assert(fourpointSitesSet[set_no].size()==fourpointValuesSet[set_no].size());

        value_sum_=zero;

        cout<<"Set = "<<set_no<<endl;
        for(int term_no=0;term_no<fourpointSitesSet[set_no].size();term_no++){
            Mat_2_doub AMat1, AMat2;
            tetra_int sites_;
            tetra_int spins_;

            AMat1=AMat0;
            AMat2=AMat0;

            sites_ = fourpointSitesSet[set_no][term_no];
            spins_ = fourpointSpinsSet[set_no][term_no];

            AMat1[spins_.first][spins_.second]=one;
            AMat2[spins_.third][spins_.fourth]=one;

            Get_CdaggerC_type_Opr(AMat2, Vec_, Vec_temp_, sites_.third, sites_.fourth);


            Get_CdaggerC_type_Opr(AMat1, Vec_temp_, Vec_final_, sites_.first, sites_.second);

            //For classical <cdag1 c2>
            Get_CdaggerC_type_Opr(AMat1, Vec_, Vec_temp2_, sites_.first, sites_.second);


            value_ = fourpointValuesSet[set_no][term_no]*dot_product(Vec_final_, Vec_);


            Value2_class = dot_product(Vec_temp_, Vec_);
            Value1_class = dot_product(Vec_temp2_, Vec_);


            value_sum_ += value_;

            value_sum_quantum +=  value_ - (fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class);

            cout<<"term="<<term_no<<"  coeff="<<fourpointValuesSet[set_no][term_no]
                <<"  sites=("<<sites_.first<<","<<sites_.second<<","<<sites_.third<<","<<sites_.fourth<<")"
                <<"  spins=("<<spins_.first<<","<<spins_.second<<","<<spins_.third<<","<<spins_.fourth<<")"
                <<"  value="<<value_
                <<"  value_classical="<<fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class
                <<"  value_quantum="<<value_ - (fourpointValuesSet[set_no][term_no]*Value1_class*Value2_class)
                <<endl;

            vector< double_type >().swap( Vec_temp_ );
            vector< double_type >().swap( Vec_temp2_ );
            vector< double_type >().swap( Vec_final_ );
        }

        cout<<"Total for set "<<set_no<<" = "<<value_sum_<<endl;
        cout<<"Total quantum for set "<<set_no<<" = "<<value_sum_quantum<<endl;
        cout<<endl;
        Total_Value_Sum += value_sum_;
        Total_Value_Sum_quantum += value_sum_quantum;
    }

    cout<<"------------------------------------------------"<<endl;
    cout<<"Total Value Sum = "<<Total_Value_Sum<<endl;
    cout<<"Total Value Sum Quantum = "<<Total_Value_Sum_quantum<<endl;

}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Mat_1_doub &Vec_in, Mat_1_doub &Vec_out, int site){


    assert(Vec_in.size()==basis.D_up_basis.size());
    Vec_out.clear();
    Vec_out.resize(basis.D_up_basis.size());
    for(int i=0;i<basis.D_up_basis.size();i++){
        Vec_out[i]=zero;
    }
    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Mat_1_doub> Vec_out_private;
    Vec_out_private.resize(N_threads);
    for(int thread=0;thread<N_threads;thread++){
        Vec_out_private[thread].resize(basis.D_up_basis.size());
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out_private[thread][i]=zero;
        }
    }

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    int i_new, j_new, m_new;
    int l,lp, max_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;
    int value_;
    double_type value_diagonal;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ((sigma)) == ((sigma_p)) ){
                        assert(sigma == sigma_p);
                        if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i], site);}
                        if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i], site);}
                        if(value_ != 0){
                            value_diagonal += AMat[(sigma)][(sigma_p)]*one;
                        }
                    }
                }
            }
        }

        Vec_out_private[thread_id][i] += value_diagonal*Vec_in[i];

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ((sigma)) != ((sigma_p)) ){
                        if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i], site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }
                        if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }

                        if(check)
                        {
                            D_up = (int) (basis.D_up_basis[i]
                                          + ((1-sigma)*pow(2,site))
                                          - ((1-sigma_p)*pow(2,site)) );
                            D_dn = (int) (basis.D_dn_basis[i]
                                          + (sigma*pow(2,site))
                                          - (sigma_p*pow(2,site)) );

                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;

                            l=site;
                            lp=site;
                            max_up = basis.Length -1;
                            min_dn = 0;

                            if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                if(l != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                if(lp != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                                if(lp != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                                if(l != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                sign_pow_dn=0;
                                sign_pow_up = one_bits_in_bw(l ,lp,basis.D_up_basis[i]);
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                sign_pow_up=0;
                                sign_pow_dn = one_bits_in_bw(l ,lp,basis.D_dn_basis[i]);
                            }

                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                            Vec_out_private[thread_id][m_new] += (sign_FM*AMat[(sigma)][(sigma_p)]*one)*Vec_in[i];
                        }
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out[i] += Vec_out_private[thread][i];
        }
    }


}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Mat_1_doub &Vec_in, Mat_1_doub &Vec_out, int site, int site_p){

    //sum_{sigma,sigma_p} AMat[sigma][sigma_p]
    //c_{sigma}^{dagger,site}*c_{sigma_p,site_p}
    assert(Vec_in.size()==basis.D_up_basis.size());
    Vec_out.clear();
    Vec_out.resize(basis.D_up_basis.size());
    for(int i=0;i<basis.D_up_basis.size();i++){
        Vec_out[i]=zero;
    }
    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Mat_1_doub> Vec_out_private;
    Vec_out_private.resize(N_threads);
    for(int thread=0;thread<N_threads;thread++){
        Vec_out_private[thread].resize(basis.D_up_basis.size());
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out_private[thread][i]=zero;
        }
    }

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    int i_new, j_new, m_new;
    int l,lp, max_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;
    int value_;
    double_type value_diagonal;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;
        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if( AMat[(sigma)][(sigma_p)] != zero){
                    if( ( ((sigma)) == ((sigma_p)) ) &&
                            (site == site_p)
                            ){
                        assert(sigma == sigma_p);
                        if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i],site);}
                        if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i],site);}
                        if(value_ != 0){
                            value_diagonal += AMat[(sigma)][(sigma_p)]*one;
                        }
                    }
                }
            }
        }

        Vec_out_private[thread_id][i] += value_diagonal*Vec_in[i];

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ( ((sigma)) != ((sigma_p)) ) ||
                            (site != site_p)
                            ){
                        if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }
                        if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i], site)==0));
                        }

                        if(check)
                        {
                            D_up = (int) (basis.D_up_basis[i]
                                          + ((1-sigma)*pow(2, site))
                                          - ((1-sigma_p)*pow(2, site_p)) );
                            D_dn = (int) (basis.D_dn_basis[i]
                                          + (sigma*pow(2, site))
                                          - (sigma_p*pow(2, site_p)) );

                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;

                            l=site;
                            lp=site_p;
                            max_up = basis.Length -1;
                            min_dn = 0;

                            if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                if(l != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                if(lp != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                                if(lp != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                                if(l != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                sign_pow_dn=0;
                                sign_pow_up = one_bits_in_bw(l ,lp,basis.D_up_basis[i]);
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                sign_pow_up=0;
                                sign_pow_dn = one_bits_in_bw(l ,lp,basis.D_dn_basis[i]);
                            }

                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                            Vec_out_private[thread_id][m_new] += (sign_FM*AMat[(sigma)][(sigma_p)]*one)*Vec_in[i];
                        }
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        for(int i=0;i<basis.D_up_basis.size();i++){
            Vec_out[i] += Vec_out_private[thread][i];
        }
    }



}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site){


    OPR.value.clear();
    OPR.rows.clear();
    OPR.columns.clear();
    OPR.nrows = basis.D_up_basis.size();
    OPR.ncols = OPR.nrows;
    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Matrix_COO> OPR_private;
    OPR_private.resize(N_threads);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    int i_new, j_new, m_new;
    int l,lp, max_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;
    int value_;
    double_type value_diagonal;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ((sigma)) == ((sigma_p)) ){
                        assert(sigma == sigma_p);
                        if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i], site);}
                        if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i], site);}
                        if(value_ != 0){
                            value_diagonal += AMat[(sigma)][(sigma_p)]*one;
                        }
                    }
                }
            }
        }

        OPR_private[thread_id].value.push_back(value_diagonal*one);
        OPR_private[thread_id].rows.push_back(i);
        OPR_private[thread_id].columns.push_back(i);

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ((sigma)) != ((sigma_p)) ){
                        if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i], site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }
                        if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }

                        if(check)
                        {
                            D_up = (int) (basis.D_up_basis[i]
                                          + ((1-sigma)*pow(2,site))
                                          - ((1-sigma_p)*pow(2,site)) );
                            D_dn = (int) (basis.D_dn_basis[i]
                                          + (sigma*pow(2,site))
                                          - (sigma_p*pow(2,site)) );

                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;

                            l=site;
                            lp=site;
                            max_up = basis.Length -1;
                            min_dn = 0;

                            if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                if(l != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                if(lp != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                                if(lp != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                                if(l != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                sign_pow_dn=0;
                                sign_pow_up = one_bits_in_bw(l ,lp,basis.D_up_basis[i]);
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                sign_pow_up=0;
                                sign_pow_dn = one_bits_in_bw(l ,lp,basis.D_dn_basis[i]);
                            }

                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                            OPR_private[thread_id].value.push_back(sign_FM*AMat[(sigma)][(sigma_p)]*one);
                            OPR_private[thread_id].rows.push_back(m_new);
                            OPR_private[thread_id].columns.push_back(i);
                        }
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        OPR.value.insert(OPR.value.end(),OPR_private[thread].value.begin(), OPR_private[thread].value.end() );
        OPR.rows.insert(OPR.rows.end(),OPR_private[thread].rows.begin(), OPR_private[thread].rows.end() );
        OPR.columns.insert(OPR.columns.end(),OPR_private[thread].columns.begin(), OPR_private[thread].columns.end() );
    }


}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site, int site_p){

    //sum_{sigma,sigma_p} AMat[sigma][sigma_p]
    //c_{sigma}^{dagger,site}*c_{sigma_p,site_p}
    OPR.value.clear();
    OPR.rows.clear();
    OPR.columns.clear();
    OPR.nrows = basis.D_up_basis.size();
    OPR.ncols = OPR.nrows;
    int N_threads = 1;
#ifdef _OPENMP
    N_threads = omp_get_max_threads();
#endif
    vector<Matrix_COO> OPR_private;
    OPR_private.resize(N_threads);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
    int thread_id = 0;
    int i_new, j_new, m_new;
    int l,lp, max_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;
    int value_;
    double_type value_diagonal;

#ifdef _OPENMP
    thread_id = omp_get_thread_num();

#pragma omp for
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;
        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if( AMat[(sigma)][(sigma_p)] != zero){
                    if( ( ((sigma)) == ((sigma_p)) ) &&
                            (site == site_p)
                            ){
                        assert(sigma == sigma_p);
                        if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i],site);}
                        if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i],site);}
                        if(value_ != 0){
                            value_diagonal += AMat[(sigma)][(sigma_p)]*one;
                        }
                    }
                }
            }
        }

        OPR_private[thread_id].value.push_back(value_diagonal*one);
        OPR_private[thread_id].rows.push_back(i);
        OPR_private[thread_id].columns.push_back(i);

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                if(AMat[(sigma)][(sigma_p)] != zero){
                    if( ( ((sigma)) != ((sigma_p)) ) ||
                            (site != site_p)
                            ){
                        if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i],site)==0));
                        }
                        if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                            check=(bit_value(basis.D_dn_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_up_basis[i],site)==0));
                        }
                        if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                            check=(bit_value(basis.D_up_basis[i],site_p)==1);
                            check = (check &&
                                     (bit_value(basis.D_dn_basis[i], site)==0));
                        }

                        if(check)
                        {
                            D_up = (int) (basis.D_up_basis[i]
                                          + ((1-sigma)*pow(2, site))
                                          - ((1-sigma_p)*pow(2, site_p)) );
                            D_dn = (int) (basis.D_dn_basis[i]
                                          + (sigma*pow(2, site))
                                          - (sigma_p*pow(2, site_p)) );

                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;

                            l=site;
                            lp=site_p;
                            max_up = basis.Length -1;
                            min_dn = 0;

                            if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                                if(l != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                                if(lp != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                                if(lp != max_up){
                                    sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                                }
                                sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                                if(l != min_dn){
                                    sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                                }
                            }
                            if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                sign_pow_dn=0;
                                sign_pow_up = one_bits_in_bw(l ,lp,basis.D_up_basis[i]);
                            }
                            if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                sign_pow_up=0;
                                sign_pow_dn = one_bits_in_bw(l ,lp,basis.D_dn_basis[i]);
                            }

                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);
                            OPR_private[thread_id].value.push_back(sign_FM*AMat[(sigma)][(sigma_p)]*one);
                            OPR_private[thread_id].rows.push_back(m_new);
                            OPR_private[thread_id].columns.push_back(i);
                        }
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    }
#endif

    for(int thread=0;thread<N_threads;thread++){
        OPR.value.insert(OPR.value.end(),OPR_private[thread].value.begin(), OPR_private[thread].value.end() );
        OPR.rows.insert(OPR.rows.end(),OPR_private[thread].rows.begin(), OPR_private[thread].rows.end() );
        OPR.columns.insert(OPR.columns.end(),OPR_private[thread].columns.begin(), OPR_private[thread].columns.end() );
    }



}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Initialize_Opr_for_Dynamics(){

}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_c_on_GS(Mat_1_doub &EigVec_, BASIS_1_orb_Hubbard_GC & basis_Nm1,
                                                     Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


    int site_val ;
    int orb_val ;
    int spin_val ;
    int D_dn_new, D_up_new;
    int i_new, i_up,i_dn;
    int max_up, max_dn, min_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    int l;
    double sign_FM;
    double_type value_in, value;

    State_c_on_GS.clear();
    State_c_on_GS.resize(basis_Nm1.D_up_basis.size());

    assert(TRIO_VEC.size() == values.size());

    for(int n=0;n<TRIO_VEC.size();n++){
        site_val = TRIO_VEC[n].site_;
        orb_val = 0;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        for (int i=0;i<basis.D_up_basis.size();i++){


            if(spin_val==0){
                //For c_up|GS>
                if(bit_value(basis.D_up_basis[i],site_val)==1){
                    l = site_val;

                    D_up_new = (int) (basis.D_up_basis[i] - pow(2,site_val) );
                    D_dn_new = basis.D_dn_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Nm1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Nm1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_up = basis.Length -1;
                    min_up =0 ;
                    sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                    if(l != min_up){
                        sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                    }

                    sign_FM = pow(-1.0, sign_pow_up);

                    value = sign_FM*EigVec_[i]*value_in;

                    State_c_on_GS[i_new] += value;


                }
            }

            if(spin_val==1){
                //For c_dn|GS>
                if(bit_value(basis.D_dn_basis[i], site_val)==1){
                    l = site_val;

                    D_dn_new = (int) (basis.D_dn_basis[i] - pow(2, site_val) );
                    D_up_new = basis.D_up_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Nm1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Nm1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_dn = basis.Length -1;
                    min_dn=0;
                    sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[i]) ;
                    if(l != min_dn){
                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                    }
                    sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                    sign_FM = pow(-1.0, sign_pow_dn);

                    value = sign_FM*EigVec_[i]*value_in;

                    State_c_on_GS[i_new] += value;


                }
            }



        }


    }



}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_cdagger_on_GS(Mat_1_doub & EigVec_, BASIS_1_orb_Hubbard_GC & basis_Np1,
                                                           Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


    int site_val ;
    int orb_val ;
    int spin_val ;
    int D_dn_new, D_up_new;
    int i_new, i_up,i_dn;
    int max_up,max_dn, min_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    int l;
    double sign_FM;
    double_type value_in, value;

    State_cdagger_on_GS.clear();
    State_cdagger_on_GS.resize(basis_Np1.D_up_basis.size());

    assert(TRIO_VEC.size() == values.size());

    for(int n=0;n<TRIO_VEC.size();n++){
        site_val = TRIO_VEC[n].site_;
        orb_val = 0;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        for (int i=0;i<basis.D_up_basis.size();i++){


            if(spin_val==0){
                //For c_dagger_up|GS>
                if(bit_value(basis.D_up_basis[i], site_val)==0){
                    l =  site_val;

                    D_up_new = (int) (basis.D_up_basis[i] + pow(2, site_val) );
                    D_dn_new = basis.D_dn_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Np1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Np1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_up = basis.Length -1;
                    min_up=0;
                    sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                    if(l != min_up){
                        sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                    }

                    sign_FM = pow(-1.0, sign_pow_up);

#ifdef USE_COMPLEX
                    value = sign_FM*EigVec_[i]*conj(value_in);
#endif
#ifndef USE_COMPLEX
                    value = sign_FM*EigVec_[i]*(value_in);
#endif


                    State_cdagger_on_GS[i_new] += value;


                }
            }

            if(spin_val==1){
                //For c_dn|GS>
                if(bit_value(basis.D_dn_basis[i], site_val)==0){
                    l = site_val;

                    D_dn_new = (int) (basis.D_dn_basis[i] + pow(2,site_val) );
                    D_up_new = basis.D_up_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Np1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Np1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_dn = basis.Length -1;
                    min_dn=0;
                    sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[i]) ;
                    if(l != min_dn){
                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                    }
                    sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                    sign_FM = pow(-1.0, sign_pow_dn);


#ifdef USE_COMPLEX
                    value = sign_FM*EigVec_[i]*conj(value_in);
#endif
#ifndef USE_COMPLEX
                    value = sign_FM*EigVec_[i]*(value_in);
#endif

                    State_cdagger_on_GS[i_new] += value;


                }
            }



        }


    }



}



#endif


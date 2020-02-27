#include "Model_3_orb_Hubbard_chain_GC.h"
#include <stdlib.h>
using namespace std;
#define PI 3.14159265
/*convention for basis:

    1)  for "up-spin" basis
                  [_______________________  _  ] [_______________________  _  ]
        site----->[012....................(L-1)] [012....................(L-1)]
        orbital-->[.........orb - "0"..........] [.........orb - "1"..........]

    2)  similarly for "down-spin" basis

    4) orb_0 = d_xz, orb_1 = d_yz, orb_2 = d_xy
    */
template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Add_diagonal_terms(){


    //----------------MACRO OPRS INITIALIZATION-------------------//
    macro_obs.clear();
    Macro_oprts.clear();


    macro_obs.push_back("Hopping");macro_obs.push_back("CFS");
    macro_obs.push_back("InterOrbRepulsion");
    macro_obs.push_back("IntraOrbRepulsion");
    macro_obs.push_back("Hunds_zz");
    macro_obs.push_back("Hunds_pm_mp");
    macro_obs.push_back("PairHopping");
    assert(macro_obs.size()==(num_PairHopping+1));

    Macro_oprts.resize(num_PairHopping+1);
    for(int type_=num_Hopping;type_<=num_PairHopping;type_++){
        Macro_oprts[type_].nrows = basis.D_up_basis.size();
        Macro_oprts[type_].ncols = Macro_oprts[type_].nrows;
    }

    //---------------------------------------------------------------------//

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;

    int gamma, gamma_p;
    int D_up, D_dn;
    int nup_temp, ndn_temp;
    int m_new, i_new, j_new;
    int l, lp;
    int sign_pow_up;
    double sign_FM;


    //Remember H[l][m]=<l|H|m>
    int m,j;
    double value1, value2, value3, value4, value5, value;
    for (int i=0;i<basis.D_up_basis.size();i++){

        m=i;
        j=i;

        value=0;
        //intra-orbital coulomb repulsion:
        value1=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);

        value2=0;
        //inter-orbital coulomb repulsion:
        for(int gamma=0;gamma<3;gamma++){
            for(int gamma_p=gamma+1;gamma_p<3;gamma_p++){
                for(int site=0;site<basis.Length;site++){
                    value2+=(U_p - (J_H*0.5))*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                              ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site) +
                                bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site) )
                              );
                }
            }
        }

        value3=0;
        //SzSz Hunds coupling:
        for(int gamma=0;gamma<3;gamma++){
            for(int gamma_p=gamma+1;gamma_p<3;gamma_p++){
                for(int site=0;site<basis.Length;site++){
                    value3+=0.25*(-J_H*2.0)*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                              ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site) -
                                bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site) )
                              );
                }
            }
        }


        value4=0;
        //Crystal Field Splitting (CFE):
        for(int gamma=0;gamma<3;gamma++){
            for(int site=0;site<basis.Length;site++){
                value4+=(CFS_SITE_RESOLVED_[site][gamma])*
                        ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                            bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                          );
            }
        }



        value5=0;
        //magnetic Field * [Sz]
        for(int gamma=0;gamma<3;gamma++){
            for(int site=0;site<basis.Length;site++){
                value5+=0.5*(H_field[site])*
                        ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                            bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                          );
            }
        }


        value=value1+value2+value3+value4+value5;
        if(value!=0){
            Hamil.value.push_back(value*one);
            Hamil.rows.push_back(m);
            Hamil.columns.push_back(m);
        }

        //Macro_oprs
        if(value1!=0){
            Macro_oprts[num_IntraOrbRepulsion].value.push_back(value1*one);
            Macro_oprts[num_IntraOrbRepulsion].rows.push_back(m);
            Macro_oprts[num_IntraOrbRepulsion].columns.push_back(m);
        }
        if(value2!=0){
            Macro_oprts[num_InterOrbRepulsion].value.push_back(value2*one);
            Macro_oprts[num_InterOrbRepulsion].rows.push_back(m);
            Macro_oprts[num_InterOrbRepulsion].columns.push_back(m);
        }
        if(value3!=0){
            Macro_oprts[num_Hunds_zz].value.push_back(value3*one);
            Macro_oprts[num_Hunds_zz].rows.push_back(m);
            Macro_oprts[num_Hunds_zz].columns.push_back(m);

        }
        if(value4!=0){
            Macro_oprts[num_CFS].value.push_back(value4*one);
            Macro_oprts[num_CFS].rows.push_back(m);
            Macro_oprts[num_CFS].columns.push_back(m);
        }

    }

}


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Add_non_diagonal_terms(){


    double value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM ;
    int sign_pow_up, sign_pow_dn ;
    int l,lp;
    int _gamma, _gamma_p;
    int nup_temp, ndn_temp;
    for (int i=0;i<basis.D_up_basis.size();i++){

        m=i;
        j=i;

        value=0;


        for(int gamma=0;gamma<3;gamma++){
            for(int gamma_p=gamma+1;gamma_p<3;gamma_p++){
                for(int site=0;site<basis.Length;site++){

                    //Sp_site_gamma*Sm_site_gamma_p  Hunds coupling:
                    //there have to be ony up electron in gamma_p, site
                    //there have to be only down electron in gamma, site
                    if(
                            ((bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)
                             &&
                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0)
                             )
                            &&
                            ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site)==0)
                             )
                            )
                    {

                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );
                        D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site)
                                      - pow(2,gamma*basis.Length + site) );



                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,basis.D_dn_basis);

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            //m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //      basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];
                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }
                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);


                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                        assert(m_new<m);

                        Hamil.value.push_back(sign_FM*(0.5*(-J_H*2.0))*one);
                        Hamil.rows.push_back(m_new);
                        Hamil.columns.push_back(m);

                        Macro_oprts[num_Hunds_pm_mp].value.push_back(sign_FM*(0.5*(-J_H*2.0))*one);
                        Macro_oprts[num_Hunds_pm_mp].rows.push_back(m_new);
                        Macro_oprts[num_Hunds_pm_mp].columns.push_back(m);



                    } // if SpSm possible



                    //Pair hopping: P_i\gamma = c_i\gamma\dn c_i\gamma\dn
                    //there have to be pair present at i,gamma_p
                    //there have to nothing at i,gamma

                    if(
                            ((bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==0)
                             &&
                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0)
                             )
                            &&
                            ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site)==1)
                             )

                            )
                    {

                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );
                        D_dn = (int) (basis.D_dn_basis[j] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );



                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,basis.D_dn_basis);


                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            //m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //      basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);


                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);



                        assert(m_new<m);
                        Hamil.value.push_back(sign_FM*J_H*one);
                        Hamil.rows.push_back((m_new));
                        Hamil.columns.push_back((m));

                        Macro_oprts[num_PairHopping].value.push_back(sign_FM*J_H*one);
                        Macro_oprts[num_PairHopping].rows.push_back((m_new));
                        Macro_oprts[num_PairHopping].columns.push_back((m));



                    } //Pair-Hopping


                } // site
            } //gamma_p
        } //gamma



#ifdef USE_COMPLEX

        //magnetic Field * [-Lz]
        //---------------------------------------
        _gamma=1;_gamma_p=0;
        for(int site=0;site<basis.Length;site++){
            //for up spin
            if(
                    (bit_value(basis.D_up_basis[i],_gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],_gamma_p*basis.Length + site)==0)
                    )
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2,_gamma*basis.Length + site) + pow(2,_gamma_p*basis.Length + site) );
                D_dn = (int) (basis.D_dn_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=_gamma*basis.Length + site;
                lp=_gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                assert(m_new<m);
                Hamil.value.push_back((1.0*H_field[site])*sign_FM*iota_comp);
                Hamil.rows.push_back((m_new));
                Hamil.columns.push_back((m));

            }

            //for down spin
            if(
                    (bit_value(basis.D_dn_basis[i],_gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_dn_basis[i],_gamma_p*basis.Length + site)==0)
                    )
            {

                D_dn = (int) (basis.D_dn_basis[i] - pow(2,_gamma*basis.Length + site) + pow(2,_gamma_p*basis.Length + site) );
                D_up = (int) (basis.D_up_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=_gamma*basis.Length + site;
                lp=_gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                assert(m_new<m);
                Hamil.value.push_back((1.0*H_field[site])*sign_FM*iota_comp);
                Hamil.rows.push_back((m_new));
                Hamil.columns.push_back((m));

            }

        }
        //---------------------------------------
#endif


    } // "i" i.e up_decimals

}



template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Add_connections(){

    //REMEMBER TO PUT IT FALSE
    bool make_blocks=false;

    double_type value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int nup_temp, ndn_temp;
    int N1_temp, N2_temp;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int l,lp;
    bool Boundary_Connected;
    complex<double> iota_(0.0,1.0);
    double_type boundary_value;


    if(PBC==false && TBC==false){
        //it is OBC
        Boundary_Connected=false;
        boundary_value=zero;
    }
    else if(PBC==true && TBC==false){
        //it is PBC
        Boundary_Connected=true;
        boundary_value=one;
    }
    else if(PBC==false && TBC==true){
        //it is TBC
        Boundary_Connected=true;
#ifdef USE_COMPLEX
        boundary_value=exp(iota_*((2.0*PI*basis.Boundary_phase_exponent)/(1.0*basis.No_of_sys_copies)));
#endif
#ifndef USE_COMPLEX
        cout<<"To use TBC, please use DUSE_COMPLEX"<<endl;
        assert(1==0);
#endif
    }





    for (int i=0;i<basis.D_up_basis.size();i++){
        //cout<<i<<" done"<<endl;
        m=i;
        j=i;

        for(int site=0;site<basis.Length ;site++){

            for(int gamma=0;gamma<3;gamma++){

                for(int site_p=0;site_p<basis.Length ;site_p++){
                    int neigh =site-site_p;

                    for(int gamma_p=gamma;gamma_p>=0;gamma_p--){

                        if( ((abs(neigh)==1) && (gamma_p<gamma))
                                ||
                                ((neigh==1) && (gamma_p==gamma))
                                ||
                                (
                                    (
                                        ((abs(neigh)==(basis.Length - 1) ) && (gamma_p<gamma))
                                        ||
                                        ((neigh==(basis.Length - 1) ) && (gamma_p==gamma))
                                        )
                                    &&
                                    Boundary_Connected
                                    )

                                )
                        {


                            //IS BOUNDARY CONNECTED??
                            if((
                                        ((abs(neigh)==(basis.Length - 1) ) && (gamma_p<gamma))
                                        ||
                                        ((neigh==(basis.Length - 1) ) && (gamma_p==gamma))
                                        )==true){
                                //Yes, It is boundary
                                value=boundary_value;

                            }
                            else{
                                //No, it is not boundary
                                value=one;
                            }

                            //if(site == 1 && make_blocks==true){
                            //  value=0.0*one;
                            //}



                            //---------------Hopping for up electrons-------------------//
                            //there have to be one up electron in gamma, site
                            //there have to be no up electron in gamma_p, site_p
                            if(
                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                    &&
                                    (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site_p)==0)
                                    )
                            {



                                D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site_p)
                                              - pow(2,gamma*basis.Length + site) );
                                D_dn = basis.D_dn_basis[j];

                                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,basis.D_dn_basis);

                                if(basis.Restricted==true){
                                    nup_temp = __builtin_popcount(D_up);
                                    ndn_temp = __builtin_popcount(D_dn);
                                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    N1_temp =0;
                                    N2_temp =0;
                                    for(int orb_temp=0;orb_temp<3;orb_temp++){
                                        N1_temp += bit_value(D_up,orb_temp*basis.Length + site) +
                                                bit_value(D_dn,orb_temp*basis.Length + site);
                                        N2_temp += bit_value(D_up,orb_temp*basis.Length + site_p) +
                                                bit_value(D_dn,orb_temp*basis.Length + site_p);

                                    }






                                }
                                else{
                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                }

                                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                l=gamma*basis.Length + site;
                                lp=gamma_p*basis.Length + site_p;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                                sign_FM = pow(-1.0, sign_pow_up);


                                if(!basis.Restricted){
                                    if( (Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));

                                        Macro_oprts[num_Hopping].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Macro_oprts[num_Hopping].rows.push_back((m_new));
                                        Macro_oprts[num_Hopping].columns.push_back((m));

                                    }
                                }
                                else{
                                    if( (Hopping_mat_NN[gamma_p][gamma])!=0 &&
                                            Is_int_in_array(N1_temp, basis.Local_occupations_allowed)
                                            &&
                                            Is_int_in_array(N2_temp, basis.Local_occupations_allowed)
                                            ){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));

                                        Macro_oprts[num_Hopping].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Macro_oprts[num_Hopping].rows.push_back((m_new));
                                        Macro_oprts[num_Hopping].columns.push_back((m));


                                    }
                                }

                            } // if up hopping possible


                            //---------------Hopping for dn electrons-------------------//
                            //there have to be one dn electron in gamma, site
                            //there have to be no dn electron in gamma_p, site_p
                            if(
                                    (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)
                                    &&
                                    (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site_p)==0)
                                    )
                            {

                                D_up = basis.D_up_basis[i];
                                D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site_p)
                                              - pow(2,gamma*basis.Length + site) );


                                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,basis.D_dn_basis);
                                if(basis.Restricted==true){
                                    nup_temp = __builtin_popcount(D_up);
                                    ndn_temp = __builtin_popcount(D_dn);
                                    //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                    //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    N1_temp =0;
                                    N2_temp =0;
                                    for(int orb_temp=0;orb_temp<3;orb_temp++){
                                        N1_temp += bit_value(D_up,orb_temp*basis.Length + site) +
                                                bit_value(D_dn,orb_temp*basis.Length + site);
                                        N2_temp += bit_value(D_up,orb_temp*basis.Length + site_p) +
                                                bit_value(D_dn,orb_temp*basis.Length + site_p);

                                    }

                                }
                                else{
                                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                }

                                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                l=gamma*basis.Length + site;
                                lp=gamma_p*basis.Length + site_p;

                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                sign_FM = pow(-1.0, sign_pow_dn);


                                if(!basis.Restricted){
                                    if( (Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));

                                        Macro_oprts[num_Hopping].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Macro_oprts[num_Hopping].rows.push_back((m_new));
                                        Macro_oprts[num_Hopping].columns.push_back((m));


                                    }
                                }
                                else{
                                    if( (Hopping_mat_NN[gamma_p][gamma])!=0 &&
                                            Is_int_in_array(N1_temp, basis.Local_occupations_allowed)
                                            &&
                                            Is_int_in_array(N2_temp, basis.Local_occupations_allowed)
                                            ){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));

                                        Macro_oprts[num_Hopping].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*value);
                                        Macro_oprts[num_Hopping].rows.push_back((m_new));
                                        Macro_oprts[num_Hopping].columns.push_back((m));


                                    }
                                }


                            } // if dn hopping possible


                        }//nearest neighbour
                    } //gamma_p

                }//site_p

            } //gamma


        } // site

    } // "i" i.e up_decimals

}

template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Add_Spin_Orbit_Coupling(){

#ifdef USE_COMPLEX


    double value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int nup_temp, ndn_temp;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int l,lp;
    int max_up, max_dn;
    int min_up, min_dn;
    int gamma ,gamma_p;

    Mat_1_intpair orbs;
    Mat_1_Complex_doub SO_val;






    for (int i=0;i<basis.D_up_basis.size();i++){
        //cout<<i<<" done"<<endl;
        m=i;
        j=i;

        value=0;


        for(int site=0;site<basis.Length ;site++){


            orbs.resize(4);
            orbs[0].first = 1;orbs[0].second = 2;
            orbs[1].first = 0;orbs[1].second = 2;
            orbs[2].first = 2;orbs[2].second = 1;
            orbs[3].first = 2;orbs[3].second = 0;

            SO_val.resize(4);
            SO_val[0].real(-1.0);SO_val[0].imag(0.0);
            SO_val[1].real(0.0);SO_val[1].imag(-1.0);
            SO_val[2].real(1.0);SO_val[2].imag(0.0);
            SO_val[3].real(0.0);SO_val[3].imag(1.0);




            for(int term_no=0;term_no<orbs.size();term_no++){
                gamma=orbs[term_no].first;
                gamma_p=orbs[term_no].second;


                if(
                        (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                        &&
                        (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                        )
                {

                    D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) );
                    D_dn = (int) (basis.D_dn_basis[i] + pow(2,gamma_p*basis.Length + site) );
                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    l=gamma*basis.Length + site;
                    lp=gamma_p*basis.Length + site;
                    max_up = 2*basis.Length + basis.Length -1;
                    min_up = 0;
                    max_dn = 2*basis.Length + basis.Length -1;
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



                    if((lambda_SOC)!=0){
                        assert(m_new<m);
                        Hamil.value.push_back(sign_FM*SO_val[term_no]*(lambda_SOC*0.5));
                        Hamil.rows.push_back((m_new));
                        Hamil.columns.push_back((m));
                    }



                }

            }


            orbs.resize(1);
            orbs[0].first = 1;orbs[0].second = 0;

            SO_val.resize(1);
            SO_val[0].real(0.0);SO_val[0].imag(-1.0);


            for(int term_no=0;term_no<orbs.size();term_no++){
                gamma=orbs[term_no].first;
                gamma_p=orbs[term_no].second;


                if(
                        (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                        &&
                        (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                        )
                {

                    D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                    D_dn = (int) (basis.D_dn_basis[i] );

                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    l=gamma*basis.Length + site;
                    lp=gamma_p*basis.Length + site;


                    sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);


                    //try this as well
                    /*
                                sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                */

                    sign_FM = pow(-1.0, sign_pow_up);



                    if((lambda_SOC)!=0){
                        assert(m_new<m);
                        Hamil.value.push_back(sign_FM*SO_val[term_no]*(lambda_SOC*0.5));
                        Hamil.rows.push_back((m_new));
                        Hamil.columns.push_back((m));
                    }



                }

            }



            orbs.resize(1);
            orbs[0].first = 1;orbs[0].second = 0;

            SO_val.resize(1);
            SO_val[0].real(0.0);SO_val[0].imag(1.0);


            for(int term_no=0;term_no<orbs.size();term_no++){
                gamma=orbs[term_no].first;
                gamma_p=orbs[term_no].second;


                if(
                        (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                        &&
                        (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                        )
                {

                    D_up = (int) (basis.D_up_basis[i]);
                    D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );

                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    l=gamma*basis.Length + site;
                    lp=gamma_p*basis.Length + site;


                    sign_pow_dn = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);


                    //try this as well
                    /*
                                sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                */

                    sign_FM = pow(-1.0, sign_pow_dn);



                    if((lambda_SOC)!=0){
                        assert(m_new<m);
                        Hamil.value.push_back(sign_FM*SO_val[term_no]*(lambda_SOC*0.5));
                        Hamil.rows.push_back((m_new));
                        Hamil.columns.push_back((m));
                    }



                }

            }







        } // site

    } // "i" i.e up_decimals


#endif
}

template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Read_parameters(string filename){


    string filepath = filename;


    double temp_val;
    string pbc_,PBC_ ="PBC = ";
    string nsc_, NSC_ = "No_of_system_copies = ";
    string m_, M_ = "Boundary_phase_exponent = ";
    string length, Length = "Length = ";
    string ntotal, Ntotal = "N_Total = ";

    string lambda_soc_, Lambda_SOC_ = "lambda_SOC = ";

    string ucoul, Ucoul = "U = ";
    string jhund, Jhund = "JHund = ";
    string upcoul, Upcoul = "Uprime = ";
    string hmag, Hmag = "H_mag = ";
    string cfs_, CFS_ = "CFS = ";
    string cfs_site_resolved, CFS_SITE_RESOLVED = "CFS_SITE_RESOLVED = ";
    string hopp0_, Hopp0_ = "Hopping_mat[0][orb] = ";
    string hopp1_, Hopp1_ = "Hopping_mat[1][orb] = ";
    string hopp2_, Hopp2_ = "Hopping_mat[2][orb] = ";

    string restriction_on_local_occupations_, Restriction_On_Local_Occupations_ = "Restriction_on_local_occupations = ";



    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(PBC_, 0)) != string::npos) {
                pbc_ = line.substr (offset+PBC_.length());				}

            if ((offset = line.find(NSC_, 0)) != string::npos) {
                nsc_ = line.substr (offset + NSC_.length());		}

            if ((offset = line.find(M_, 0)) != string::npos) {
                m_ = line.substr (offset + M_.length());		}

            if ((offset = line.find(Length, 0)) != string::npos) {
                length = line.substr (offset + Length.length());		}

            if ((offset = line.find(Ntotal, 0)) != string::npos) {
                ntotal = line.substr (offset + Ntotal.length());		}

            if ((offset = line.find(Lambda_SOC_, 0)) != string::npos) {
                lambda_soc_= line.substr (offset + Lambda_SOC_.length());		}


            if ((offset = line.find(Ucoul, 0)) != string::npos) {
                ucoul= line.substr (offset + Ucoul.length());		}

            if ((offset = line.find(Jhund, 0)) != string::npos) {
                jhund = line.substr (offset + Jhund.length());		}

            if ((offset = line.find(Upcoul, 0)) != string::npos) {
                upcoul = line.substr (offset + Upcoul.length());		}

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(CFS_, 0)) != string::npos) {
                cfs_ = line.substr (offset+CFS_.length());				}

            if ((offset = line.find(CFS_SITE_RESOLVED, 0)) != string::npos) {
                cfs_site_resolved = line.substr (offset+CFS_SITE_RESOLVED.length());				}

            if ((offset = line.find(Hopp0_, 0)) != string::npos) {
                hopp0_ = line.substr (offset+Hopp0_.length());				}

            if ((offset = line.find(Hopp1_, 0)) != string::npos) {
                hopp1_ = line.substr (offset+Hopp1_.length());				}

            if ((offset = line.find(Hopp2_, 0)) != string::npos) {
                hopp2_ = line.substr (offset+Hopp2_.length());				}


            if ((offset = line.find(Restriction_On_Local_Occupations_, 0)) != string::npos) {
                restriction_on_local_occupations_ = line.substr (offset+Restriction_On_Local_Occupations_.length());				}


        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}

    //cout<<"read-1"<<endl;

    if(pbc_ == "true"){
        PBC =true;
        TBC=false;
    }
    else if(pbc_== "false"){
        PBC=false;
        TBC=false;
    }
    else if(pbc_ == "TBC_IsUsed"){
        PBC=false;
        TBC=true;
        basis.No_of_sys_copies=atoi(nsc_.c_str());
        basis.Boundary_phase_exponent=atof(m_.c_str());
    }
    else {
        cout <<"For this model only one of the following options are allowed:"<<endl;
        cout <<"PBC = true"<<endl;
        cout <<"PBC = false"<<endl;
        cout <<"PBC = TBC_IsUsed"<<endl;
        assert(1==0);
    }

    //   cout<<"here 1"<<endl;

    basis.Length=atoi(length.c_str());
    basis.N_total=atoi(ntotal.c_str());

    //cout<<"here 2"<<endl;

    U=atof(ucoul.c_str());
    J_H=atof(jhund.c_str());
    U_p=atof(upcoul.c_str());
    lambda_SOC=atof(lambda_soc_.c_str());

    //   cout<<"here 3"<<endl;

    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }

    //  cout<<"here 4"<<endl;

    stringstream cfs_stream;
    cfs_stream<<cfs_;

    CFS.clear();
    CFS.resize(3);


    for(int n=0;n<3;n++){
        cfs_stream >> temp_val;
        CFS[n]=temp_val;
    }

    //cout<<"here 5"<<endl;

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


    CFS_SITE_RESOLVED_.clear();
    CFS_SITE_RESOLVED_.resize(basis.Length);


    for(int i=0;i<basis.Length;i++){
        CFS_SITE_RESOLVED_[i].resize(3);
        for(int orb=0;orb<3;orb++){
            if(CFS_SITE_RESOLVED_bool==true){
                cfs_site_resolved_stream >> temp_val;
                CFS_SITE_RESOLVED_[i][orb]=temp_val;
            }
            else{
                CFS_SITE_RESOLVED_[i][orb]=CFS[orb];
            }
        }
    }

    //cout<<"here 6"<<endl;


    Hopping_mat_NN.clear();
    Hopping_mat_NN.resize(3);
    for (int i=0;i<3;i++){
        Hopping_mat_NN[i].resize(3);
    }

    //Hopping_mat_NN[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}

    //cout<<"here 7"<<endl;

    stringstream hopp0_stream(hopp0_);

    for(int n=0;n<3;n++){
        hopp0_stream >> Hopping_mat_NN[0][n];
    }
    stringstream hopp1_stream(hopp1_);

    for(int n=0;n<3;n++){
        hopp1_stream >> Hopping_mat_NN[1][n];
    }

    stringstream hopp2_stream(hopp2_);

    for(int n=0;n<3;n++){
        hopp2_stream >> Hopping_mat_NN[2][n];
    }

    cout<<"READING PARAMETERS"<<endl;

    // cout<<"here 8"<<endl;

    stringstream rolo_stream(restriction_on_local_occupations_);
    int temp_int;
    rolo_stream>>temp_int;

    if(temp_int != 0){
        basis.Local_occupations_allowed.resize(temp_int);
        for(int n=0;n<temp_int;n++){
            rolo_stream >> basis.Local_occupations_allowed[n];
        }
    }


}


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Read_parameters_for_dynamics(string filename){

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


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Initialize_one_point_to_calculate(){
    one_point_obs.resize(6);
    one_point_obs[0]="n_0_up";
    one_point_obs[1]="n_1_up";
    one_point_obs[2]="n_2_up";
    one_point_obs[3]="n_0_dn";
    one_point_obs[4]="n_1_dn";
    one_point_obs[5]="n_2_dn";

    One_point_oprts.resize(6);

    int nup_temp, ndn_temp;
    int T_no_oprs=6;
    int orb;
    int spin;



    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }




    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        if(one_point_obs[opr_no]=="n_0_up" || one_point_obs[opr_no]=="n_0_dn"){
            orb=0;
        }
        if(one_point_obs[opr_no]=="n_1_up" || one_point_obs[opr_no]=="n_1_dn"){
            orb=1;
        }
        if(one_point_obs[opr_no]=="n_2_up" || one_point_obs[opr_no]=="n_2_dn"){
            orb=2;
        }


        if(one_point_obs[opr_no]=="n_0_up" || one_point_obs[opr_no]=="n_1_up" || one_point_obs[opr_no]=="n_2_up"){
            spin=0;
        }
        else{
            spin=1;
        }


        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }


        //Remember OPR[l][m]=<l|OPR|m>
        int m,j;
        double value;


        for(int site=0;site<basis.Length;site++){

            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                //n_orb_spin[site]:
                if(spin==0){
                    value=bit_value(basis.D_up_basis[i],orb*basis.Length + site);
                }
                else{
                    value=bit_value(basis.D_dn_basis[j],orb*basis.Length + site);
                }





                if(value!=0){
                    One_point_oprts[opr_no][site].value.push_back(value*one);
                    One_point_oprts[opr_no][site].rows.push_back(m);
                    One_point_oprts[opr_no][site].columns.push_back(m);
                }

            }

        }


    }


#ifdef USE_COMPLEX //HERE
    one_point_obs.resize(14);
    one_point_obs[6]="n_3by2_3by2";
    one_point_obs[7]="n_3by2_m3by2";
    one_point_obs[8]="n_3by2_1by2";
    one_point_obs[9]="n_3by2_m1by2";
    one_point_obs[10]="n_1by2_1by2";
    one_point_obs[11]="n_1by2_m1by2";
    one_point_obs[12]="Delta_3by2_1by2_1by2_1by2";
    one_point_obs[13]="Delta_3by2_m1by2_1by2_m1by2";

    One_point_oprts.resize(14);

    for(int i=6;i<14;i++){
        One_point_oprts[i].resize(basis.Length);
    }

    cout<<"To construct operators n_{jm} eq-4 of PhysRevB.96.155111 is used"<<endl;

    int gamma, gamma_p;
    int D_up, D_dn;
    int l, lp;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    double value;
    int m,j,m_new,i_new, j_new;
    int max_up, max_dn;
    int min_up, min_dn;
    complex<double> iota(0,1);

    //n_{3by2_3by2}
    for(int site=0;site<basis.Length;site++){
        One_point_oprts[6][site].nrows = basis.D_up_basis.size();
        One_point_oprts[6][site].ncols = One_point_oprts[6][site].nrows;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            gamma=0;gamma_p=1;
            if(
                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {


                D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_dn = (int) (basis.D_dn_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[6][site].value.push_back((-0.5)*sign_FM*iota);
                One_point_oprts[6][site].rows.push_back(m_new);
                One_point_oprts[6][site].columns.push_back(m);

                One_point_oprts[6][site].value.push_back((0.5)*sign_FM*iota);
                One_point_oprts[6][site].rows.push_back(m);
                One_point_oprts[6][site].columns.push_back(m_new);

            }




            for(gamma=0;gamma<2;gamma++){
                value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);

                if(value!=0){
                    One_point_oprts[6][site].value.push_back(0.5*value*one);
                    One_point_oprts[6][site].rows.push_back(m);
                    One_point_oprts[6][site].columns.push_back(m);
                }
            }




        }

    }

    //Print_Matrix_COO(One_point_oprts[6][0]);

    //n_{3by2_m3by2}
    for(int site=0;site<basis.Length;site++){
        One_point_oprts[7][site].nrows = basis.D_up_basis.size();
        One_point_oprts[7][site].ncols = One_point_oprts[7][site].nrows;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            gamma=0;gamma_p=1;
            if(
                    (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {


                D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_up = (int) (basis.D_up_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_dn = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);

                sign_FM = pow(-1.0, sign_pow_dn);


                One_point_oprts[7][site].value.push_back((0.5)*sign_FM*iota);
                One_point_oprts[7][site].rows.push_back(m_new);
                One_point_oprts[7][site].columns.push_back(m);

                One_point_oprts[7][site].value.push_back((-0.5)*sign_FM*iota);
                One_point_oprts[7][site].rows.push_back(m);
                One_point_oprts[7][site].columns.push_back(m_new);

            }




            for(gamma=0;gamma<2;gamma++){
                value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);

                if(value!=0){
                    One_point_oprts[7][site].value.push_back(0.5*value*one);
                    One_point_oprts[7][site].rows.push_back(m);
                    One_point_oprts[7][site].columns.push_back(m);
                }
            }

        }

    }

    //n_{3by2_pm(sby2)}
    int s_val;
    double factor;
    complex<double> factor_comp;
    for(int opr_no=8;opr_no<=9;opr_no++){
        if(opr_no==8){
            s_val=-1;
        }
        else{
            s_val=1;
        }

        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;

            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                for(gamma=0;gamma<3;gamma++){

                    if(s_val==1){
                        if(gamma==2){
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);
                        }
                    }
                    else{
                        if(gamma==2){
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);
                        }
                    }

                    if(gamma==2){
                        factor=4.0/6.0;
                    }
                    else{
                        factor=1.0/6.0;
                    }

                    if(value!=0){
                        One_point_oprts[opr_no][site].value.push_back(factor*value*one);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }




                gamma=0;gamma_p=1;
                if(
                        (
                            (s_val==1)
                            &&
                            (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                            )
                        ||
                        (
                            (s_val==-1)
                            &&
                            (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                            )

                        )
                {


                    if(s_val==1){
                        D_dn = (int) (basis.D_dn_basis[i]);
                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                    }
                    else{
                        D_up = (int) (basis.D_up_basis[i]);
                        D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                    }


                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        //m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //      basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    l=gamma*basis.Length + site;
                    lp=gamma_p*basis.Length + site;

                    if(s_val==1){
                        sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);
                        sign_FM = pow(-1.0, sign_pow_up);
                    }
                    else{
                        sign_pow_dn = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);
                        sign_FM = pow(-1.0, sign_pow_dn);
                    }




                    One_point_oprts[opr_no][site].value.push_back(((1.0*s_val)/6.0)*sign_FM*iota);
                    One_point_oprts[opr_no][site].rows.push_back(m_new);
                    One_point_oprts[opr_no][site].columns.push_back(m);

                    One_point_oprts[opr_no][site].value.push_back(((-1.0*s_val)/6.0)*sign_FM*iota);
                    One_point_oprts[opr_no][site].rows.push_back(m);
                    One_point_oprts[opr_no][site].columns.push_back(m_new);

                }



                gamma_p=2;
                for(gamma=0;gamma<=1;gamma++){


                    if(gamma==0){
                        factor_comp = iota*(2.0/6.0);
                    }
                    if(gamma==1){
                        factor_comp = one*(s_val*(2.0/6.0));

                    }

                    if(
                            (
                                (s_val==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                                )
                            ||
                            (
                                (s_val==-1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                                )

                            )
                    {


                        if(s_val==1){
                            D_dn = (int) (basis.D_dn_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) );
                        }
                        else{
                            D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) );
                        }

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);


                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        max_up = 2*basis.Length + basis.Length -1;
                        min_up = 0;
                        max_dn = 2*basis.Length + basis.Length -1;
                        min_dn = 0;

                        if(s_val==1){
                            sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                            if(l != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                            if(lp != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }
                        else{
                            sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                            if(lp != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                            if(l != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }





                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);


                        One_point_oprts[opr_no][site].value.push_back(factor_comp*sign_FM);
                        One_point_oprts[opr_no][site].rows.push_back(m_new);
                        One_point_oprts[opr_no][site].columns.push_back(m);

                        One_point_oprts[opr_no][site].value.push_back(conj(factor_comp*sign_FM));
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m_new);

                    }
                }

            }




        }




    }



    //n_{1by2_pm(sby2)}

    for(int opr_no=10;opr_no<=11;opr_no++){
        if(opr_no==10){
            s_val=-1;
        }
        else{
            s_val=1;
        }

        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;

            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                for(gamma=0;gamma<3;gamma++){

                    if(s_val==1){
                        if(gamma==2){
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);
                        }
                    }
                    else{
                        if(gamma==2){
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);
                        }
                    }

                    if(gamma==2){
                        factor=1.0/3.0;
                    }
                    else{
                        factor=1.0/3.0;
                    }

                    if(value!=0){
                        One_point_oprts[opr_no][site].value.push_back(factor*value*one);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }




                gamma=0;gamma_p=1;
                if(
                        (
                            (s_val==1)
                            &&
                            (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                            )
                        ||
                        (
                            (s_val==-1)
                            &&
                            (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                            )

                        )
                {


                    if(s_val==1){
                        D_dn = (int) (basis.D_dn_basis[i]);
                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                    }
                    else{
                        D_up = (int) (basis.D_up_basis[i]);
                        D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                    }


                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    l=gamma*basis.Length + site;
                    lp=gamma_p*basis.Length + site;

                    if(s_val==1){
                        sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);
                        sign_FM = pow(-1.0, sign_pow_up);
                    }
                    else{
                        sign_pow_dn = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);
                        sign_FM = pow(-1.0, sign_pow_dn);
                    }




                    One_point_oprts[opr_no][site].value.push_back(((1.0*s_val)/3.0)*sign_FM*iota);
                    One_point_oprts[opr_no][site].rows.push_back(m_new);
                    One_point_oprts[opr_no][site].columns.push_back(m);

                    One_point_oprts[opr_no][site].value.push_back(((-1.0*s_val)/3.0)*sign_FM*iota);
                    One_point_oprts[opr_no][site].rows.push_back(m);
                    One_point_oprts[opr_no][site].columns.push_back(m_new);

                }



                gamma_p=2;
                for(gamma=0;gamma<=1;gamma++){


                    if(gamma==0){
                        factor_comp = iota*(-1.0/3.0);
                    }
                    if(gamma==1){
                        factor_comp = one*(s_val*(-1.0/3.0));

                    }

                    if(
                            (
                                (s_val==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                                )
                            ||
                            (
                                (s_val==-1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                                )

                            )
                    {


                        if(s_val==1){
                            D_dn = (int) (basis.D_dn_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) );
                        }
                        else{
                            D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) );
                        }

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);


                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        max_up = 2*basis.Length + basis.Length -1;
                        min_up = 0;
                        max_dn = 2*basis.Length + basis.Length -1;
                        min_dn = 0;


                        if(s_val==1){
                            sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                            if(l != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                            if(lp != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }
                        else{
                            sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                            if(lp != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                            if(l != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }





                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);


                        One_point_oprts[opr_no][site].value.push_back(factor_comp*sign_FM);
                        One_point_oprts[opr_no][site].rows.push_back(m_new);
                        One_point_oprts[opr_no][site].columns.push_back(m);

                        One_point_oprts[opr_no][site].value.push_back(conj(factor_comp*sign_FM));
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m_new);

                    }
                }

            }




        }




    }



    //Del_{3by2_pm(sby2)_1by2_pm(sby2)}
    for(int opr_no=12;opr_no<=13;opr_no++){
        if(opr_no==12){
            s_val=1;
        }
        else{
            s_val=-1;
        }

        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;

            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                gamma=2;
                for(gamma_p=0;gamma_p<=1;gamma_p++){
                    if(gamma_p==1){
                        factor_comp = ((-1.0*s_val)/(3.0*sqrt(2.0)))*one;
                    }
                    if(gamma_p==0){
                        factor_comp = ((-1.0)/(3.0*sqrt(2.0)))*iota;
                    }


                    if(
                            (
                                (s_val==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                                )
                            ||
                            (
                                (s_val==-1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                                )

                            )
                    {

                        if(s_val==1){
                            D_dn = (int) (basis.D_dn_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) );
                        }
                        else{
                            D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) );
                        }

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);


                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        max_up = 2*basis.Length + basis.Length -1;
                        min_up = 0;
                        max_dn = 2*basis.Length + basis.Length -1;
                        min_dn = 0;


                        if(s_val==1){
                            sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                            if(l != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                            if(lp != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }
                        else{
                            sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                            if(lp != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                            if(l != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }

                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

                        One_point_oprts[opr_no][site].value.push_back(factor_comp*sign_FM);
                        One_point_oprts[opr_no][site].rows.push_back(m_new);
                        One_point_oprts[opr_no][site].columns.push_back(m);

                    }
                }


                gamma_p=2;
                for(gamma=0;gamma<=1;gamma++){
                    if(gamma==1){
                        factor_comp = ((2.0*s_val)/(3.0*sqrt(2.0)))*one;
                    }
                    if(gamma==0){
                        factor_comp = ((-2.0)/(3.0*sqrt(2.0)))*iota;
                    }
                    if(
                            (
                                (s_val==-1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                                )
                            ||
                            (
                                (s_val==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                                )

                            )
                    {

                        if(s_val==-1){
                            D_dn = (int) (basis.D_dn_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) );
                        }
                        else{
                            D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site) );
                            D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) );
                        }

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        max_up = 2*basis.Length + basis.Length -1;
                        min_up = 0;
                        max_dn = 2*basis.Length + basis.Length -1;
                        min_dn = 0;

                        if(s_val==-1){
                            sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                            if(l != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(lp, min_dn, basis.D_dn_basis[i]);
                            if(lp != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }
                        else{
                            sign_pow_up = one_bits_in_bw(max_up ,lp,basis.D_up_basis[i]) ;
                            if(lp != max_up){
                                sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                            }
                            sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                            if(l != min_dn){
                                sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                            }
                        }

                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

                        One_point_oprts[opr_no][site].value.push_back(factor_comp*sign_FM);
                        One_point_oprts[opr_no][site].rows.push_back(m_new);
                        One_point_oprts[opr_no][site].columns.push_back(m);

                    }
                }


                for(gamma=0;gamma<=1;gamma++){
                    if(gamma==0){
                        gamma_p=1;
                        factor_comp = iota*((1.0*s_val)/(sqrt(2.0)*3.0));

                    }
                    else{
                        gamma_p=0;
                        factor_comp = iota*((-1.0*s_val)/(sqrt(2.0)*3.0));
                    }

                    if(
                            (
                                (s_val==-1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                                )
                            ||
                            (
                                (s_val==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                                )

                            )
                    {

                        if(s_val==-1){
                            D_dn = (int) (basis.D_dn_basis[i]);
                            D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                        }
                        else{
                            D_up = (int) (basis.D_up_basis[i]);
                            D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site));
                        }

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                        l=gamma*basis.Length + site;
                        lp=gamma_p*basis.Length + site;

                        if(s_val==-1){
                            sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);
                            sign_FM = pow(-1.0, sign_pow_up);
                        }
                        else{
                            sign_pow_dn = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);
                            sign_FM = pow(-1.0, sign_pow_dn);
                        }


                        One_point_oprts[opr_no][site].value.push_back(factor_comp*sign_FM);
                        One_point_oprts[opr_no][site].rows.push_back(m_new);
                        One_point_oprts[opr_no][site].columns.push_back(m);

                    }
                }



                for(gamma=0;gamma<3;gamma++){

                    if(gamma==2){
                        factor_comp = one*((2.0)/(sqrt(2.0)*3.0));
                    }
                    else{
                        factor_comp = one*((-1.0)/(sqrt(2.0)*3.0));
                    }

                    if(s_val==-1){
                        if(gamma==2){
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);
                        }
                    }
                    else{
                        if(gamma==2){
                            value=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);}
                        else{
                            value=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);
                        }
                    }

                    if(value!=0){
                        One_point_oprts[opr_no][site].value.push_back(factor_comp*value*one);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }
            }
        }
    }


#endif


#ifdef USE_COMPLEX
    one_point_obs.resize(17);
    one_point_obs[14]="Lz";
    one_point_obs[15]="Lplus";
    one_point_obs[16]="Lminus";

    One_point_oprts.resize(17);

    for(int i=14;i<17;i++){
        One_point_oprts[i].resize(basis.Length);
    }

    cout<<"To construct operators n_{jm} eq-4 of PhysRevB.96.155111 is used"<<endl;

    /*
        int gamma, gamma_p;
        int D_up, D_dn;
        int l, lp;
        int sign_pow_up, sign_pow_dn;
        double sign_FM;
        double value;
        int m,j,m_new,i_new, j_new;
        int max_up, max_dn;
        int min_up, min_dn;
        complex<double> iota(0,1);
        */

    //Lz (14)
    //---------------------------Lz is constructed-----------------------------------//

    for(int site=0;site<basis.Length;site++){
        One_point_oprts[14][site].nrows = basis.D_up_basis.size();
        One_point_oprts[14][site].ncols = One_point_oprts[14][site].nrows;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            // 0=xz, 1=yz, 2=xy
            gamma=0;gamma_p=1;

            //for up spin
            if(
                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_dn = (int) (basis.D_dn_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[14][site].value.push_back(sign_FM*iota);
                One_point_oprts[14][site].rows.push_back(m_new);
                One_point_oprts[14][site].columns.push_back(m);

                One_point_oprts[14][site].value.push_back((-1.0)*sign_FM*iota);
                One_point_oprts[14][site].rows.push_back(m);
                One_point_oprts[14][site].columns.push_back(m_new);

            }



            //for down spin
            if(
                    (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_up = (int) (basis.D_up_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[14][site].value.push_back(sign_FM*iota);
                One_point_oprts[14][site].rows.push_back(m_new);
                One_point_oprts[14][site].columns.push_back(m);

                One_point_oprts[14][site].value.push_back((-1.0)*sign_FM*iota);
                One_point_oprts[14][site].rows.push_back(m);
                One_point_oprts[14][site].columns.push_back(m_new);

            }


        }

    }


    //------------------------------------------------------------------------------//


    //Lplus (15)
    //---------------------------Lplus is constructed-----------------------------------//

    for(int site=0;site<basis.Length;site++){
        One_point_oprts[15][site].nrows = basis.D_up_basis.size();
        One_point_oprts[15][site].ncols = One_point_oprts[15][site].nrows;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;


            // 0=xz, 1=yz, 2=xy

            gamma=2;gamma_p=1;

            //for up spin
            if(
                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_dn = (int) (basis.D_dn_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[15][site].value.push_back(sign_FM*one);
                One_point_oprts[15][site].rows.push_back(m_new);
                One_point_oprts[15][site].columns.push_back(m);

                One_point_oprts[15][site].value.push_back((-1.0)*sign_FM*one);
                One_point_oprts[15][site].rows.push_back(m);
                One_point_oprts[15][site].columns.push_back(m_new);

            }



            //for down spin
            if(
                    (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_up = (int) (basis.D_up_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[15][site].value.push_back(sign_FM*one);
                One_point_oprts[15][site].rows.push_back(m_new);
                One_point_oprts[15][site].columns.push_back(m);

                One_point_oprts[15][site].value.push_back((-1.0)*sign_FM*one);
                One_point_oprts[15][site].rows.push_back(m);
                One_point_oprts[15][site].columns.push_back(m_new);

            }



            gamma=2;gamma_p=0;

            //for up spin
            if(
                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_up = (int) (basis.D_up_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_dn = (int) (basis.D_dn_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_up_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[15][site].value.push_back(sign_FM*iota);
                One_point_oprts[15][site].rows.push_back(m_new);
                One_point_oprts[15][site].columns.push_back(m);

                One_point_oprts[15][site].value.push_back((-1.0)*sign_FM*iota);
                One_point_oprts[15][site].rows.push_back(m);
                One_point_oprts[15][site].columns.push_back(m_new);

            }



            //for down spin
            if(
                    (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==0)
                    )
            {

                D_dn = (int) (basis.D_dn_basis[i] - pow(2,gamma*basis.Length + site) + pow(2,gamma_p*basis.Length + site) );
                D_up = (int) (basis.D_up_basis[i]);

                if(basis.Restricted==true){
                    nup_temp = __builtin_popcount(D_up);
                    ndn_temp = __builtin_popcount(D_dn);
                    // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                    //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                    m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                }
                else{
                    i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                    j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                    m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                            basis.Nup_offsets[__builtin_popcount(D_up)].first;
                }

                //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                l=gamma*basis.Length + site;
                lp=gamma_p*basis.Length + site;

                sign_pow_up = one_bits_in_bw(lp ,l,basis.D_dn_basis[i]);

                sign_FM = pow(-1.0, sign_pow_up);


                One_point_oprts[15][site].value.push_back(sign_FM*iota);
                One_point_oprts[15][site].rows.push_back(m_new);
                One_point_oprts[15][site].columns.push_back(m);

                One_point_oprts[15][site].value.push_back((-1.0)*sign_FM*iota);
                One_point_oprts[15][site].rows.push_back(m);
                One_point_oprts[15][site].columns.push_back(m_new);

            }
        }

        One_point_oprts[16][site]=Dagger(One_point_oprts[15][site]);
    }


    //------------------------------------------------------------------------------//

    //Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2)
    one_point_obs.resize(33);
    one_point_obs[17]="Sz";
    one_point_obs[18]="Splus";
    one_point_obs[19]="Sminus";
    one_point_obs[20]="Jz_eff";
    one_point_obs[21]="Jplus_eff";
    one_point_obs[22]="Jminus_eff";

    one_point_obs[23]="Mz";
    one_point_obs[24]="Mplus";
    one_point_obs[25]="Mminus";

    one_point_obs[26]="Mz_eff";
    one_point_obs[27]="Mplus_eff";
    one_point_obs[28]="Mminus_eff";

    one_point_obs[29]="Jz";
    one_point_obs[30]="Jplus";
    one_point_obs[31]="Jminus";

    one_point_obs[32]="N_local";


    One_point_oprts.resize(33);

    for(int i=17;i<33;i++){
        One_point_oprts[i].resize(basis.Length);
    }



    //----------------------Sz---------------
    for(int site=0;site<basis.Length;site++){
        One_point_oprts[17][site].nrows = basis.D_up_basis.size();
        One_point_oprts[17][site].ncols = One_point_oprts[17][site].nrows;


        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            value=0;
            for(int gamma=0;gamma<3;gamma++){

                value+=0.5*( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                               bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                             );

            }

            if(value!=0){
                One_point_oprts[17][site].value.push_back(value*one);
                One_point_oprts[17][site].rows.push_back(m);
                One_point_oprts[17][site].columns.push_back(m);
            }

        }

    }


    //-------------------------------//


    //----------------------Splus---------------
    for(int site=0;site<basis.Length;site++){
        One_point_oprts[18][site].nrows = basis.D_up_basis.size();
        One_point_oprts[18][site].ncols = One_point_oprts[18][site].nrows;


        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;


            for(int gamma=0;gamma<3;gamma++){

                //Sp_site_gamma[site]
                //there have to be only down electron in gamma, site

                if(((bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0)
                    ))
                {

                    D_up = (int) (basis.D_up_basis[i]
                                  + pow(2,gamma*basis.Length + site) );
                    D_dn = (int) (basis.D_dn_basis[i]
                                  - pow(2,gamma*basis.Length + site) );

                    //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                    //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                    if(basis.Restricted==true){
                        nup_temp = __builtin_popcount(D_up);
                        ndn_temp = __builtin_popcount(D_dn);
                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                    }
                    else{
                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                    }

                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                    //m_new = basis.D_dn_basis.size()*i_new + j_new;


                    l=gamma*basis.Length + site;
                    //lp=gamma_p*basis.Length + site;
                    max_up = 2*basis.Length + basis.Length -1;
                    min_up = 0;
                    max_dn = 2*basis.Length + basis.Length -1;
                    min_dn = 0;


                    sign_pow_up = one_bits_in_bw(max_up ,l,basis.D_up_basis[i]) ;
                    if(l != max_up){
                        sign_pow_up += bit_value(basis.D_up_basis[i],max_up);
                    }
                    sign_pow_dn = one_bits_in_bw(l, min_dn, basis.D_dn_basis[i]);
                    if(l != min_dn){
                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                    }

                    //try this as well
                    /*
                                    sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                    sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                    */

                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);



                    //assert(m_new<m);

                    One_point_oprts[18][site].value.push_back(sign_FM*one);
                    One_point_oprts[18][site].rows.push_back(m_new);
                    One_point_oprts[18][site].columns.push_back(m);

                }
            }
        }


        One_point_oprts[19][site]=Dagger(One_point_oprts[18][site]);

        //14 15 16 17 18 19 20 21 22
        //Lz Lp Lm Sz Sp Sm Jzeff Jpeff Jmeff

        Sum(One_point_oprts[17][site] , One_point_oprts[14][site] , One_point_oprts[20][site], 1.0, -1.0);
        Sum(One_point_oprts[18][site] , One_point_oprts[15][site] , One_point_oprts[21][site], 1.0, -1.0);
        Sum(One_point_oprts[19][site] , One_point_oprts[16][site] , One_point_oprts[22][site], 1.0, -1.0);

        //23 24 25
        //Mz Mp Mm
        Sum(One_point_oprts[14][site] , One_point_oprts[17][site] , One_point_oprts[23][site], 1.0, 2.0);
        Sum(One_point_oprts[15][site] , One_point_oprts[18][site] , One_point_oprts[24][site], 1.0, 2.0);
        Sum(One_point_oprts[16][site] , One_point_oprts[19][site] , One_point_oprts[25][site], 1.0, 2.0);


        //26 27 28
        //Mzeff Mpeff Mmeff
        Sum(One_point_oprts[17][site] , One_point_oprts[14][site] , One_point_oprts[26][site], 1.0, -2.0);
        Sum(One_point_oprts[18][site] , One_point_oprts[15][site] , One_point_oprts[27][site], 1.0, -2.0);
        Sum(One_point_oprts[19][site] , One_point_oprts[16][site] , One_point_oprts[28][site], 1.0, -2.0);

        //29 30 31
        //Jz Jp Jm
        Sum(One_point_oprts[14][site] , One_point_oprts[17][site] , One_point_oprts[29][site], 1.0, 1.0);
        Sum(One_point_oprts[15][site] , One_point_oprts[18][site] , One_point_oprts[30][site], 1.0, 1.0);
        Sum(One_point_oprts[16][site] , One_point_oprts[19][site] , One_point_oprts[31][site], 1.0, 1.0);


        /*
        one_point_obs[0]="n_0_up";
        one_point_obs[1]="n_1_up";
        one_point_obs[2]="n_2_up";
        one_point_obs[3]="n_0_dn";
        one_point_obs[4]="n_1_dn";
        one_point_obs[5]="n_2_dn";
        */


        // Sum(one_point_obs[0][site], one_point_obs[1][site], one_point_obs[32][site], 1.0, 1.0);



        //N_local[site]---------------------------------------------------
        One_point_oprts[32][site].nrows = basis.D_up_basis.size();
        One_point_oprts[32][site].ncols = One_point_oprts[32][site].nrows;
        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            //n_orb_spin[site]:
            value=0.0;
            for(orb=0;orb<3;orb++){
                value += (1.0*bit_value(basis.D_up_basis[i],orb*basis.Length + site) + bit_value(basis.D_dn_basis[j],orb*basis.Length + site));
            }

            if(value!=0){
                One_point_oprts[32][site].value.push_back(value*one);
                One_point_oprts[32][site].rows.push_back(m);
                One_point_oprts[32][site].columns.push_back(m);
            }

        }
        //--------------------------------------------------------------------

    }


    //-------------------------------//



    one_point_obs.resize(38);
    one_point_obs[33]="Delta_3by2_1by2_1by2_m1by2";
    one_point_obs[34]="Delta_3by2_m1by2_1by2_1by2";
    one_point_obs[35]="Delta_3by2_1by2_1by2_1by2";
    one_point_obs[36]="Delta_3by2_m1by2_1by2_m1by2";
    one_point_obs[37]="Dummy";
    One_point_oprts.resize(38);

    for(int i=33;i<38;i++){
        One_point_oprts[i].resize(basis.Length);
    }


    Mat_2_doub AMat, AMat0; //[sigma,alpha][sigma_p,beta]
    AMat.resize(6);AMat0.resize(6);
    for(int ind=0;ind<6;ind++){
        AMat[ind].resize(6);
        AMat0[ind].resize(6);
    }
    Mat_1_doub Vec_0, Vec_1, Vec_2;
    Vec_0.resize(6); Vec_1.resize(6); Vec_2.resize(6);
    int SPIN_UP=0;
    int SPIN_DN=1;
    int XZ_,YZ_,XY_;
    XZ_=0;YZ_=1;XY_=2;

    int Jz_SPIN_, Jz_SPIN_P_;
    int SPIN_, SPIN_P_;
    int BAR_SPIN_, BAR_SPIN_P_;
    //Del_{3by2_pm(sby2)_1by2_pm(sby2)}
    for(int opr_no=33;opr_no<=37;opr_no++){

        if(opr_no==33){
            //Del_3by2_1by2_1by2_m1by2
            Jz_SPIN_=-1;
            Jz_SPIN_P_=1;
            SPIN_=SPIN_DN;
            SPIN_P_=SPIN_UP;
            BAR_SPIN_=SPIN_UP;
            BAR_SPIN_P_=SPIN_DN;

        }
        if(opr_no==34){
            //Del_3by2_m1by2_1by2_1by2
            Jz_SPIN_=1;
            Jz_SPIN_P_=-1;
            SPIN_=SPIN_UP;
            SPIN_P_=SPIN_DN;
            BAR_SPIN_=SPIN_DN;
            BAR_SPIN_P_=SPIN_UP;
        }
        if(opr_no==35){
            //Del_3by2_1by2_1by2_1by2
            Jz_SPIN_=1;
            Jz_SPIN_P_=1;
            SPIN_=SPIN_UP;
            SPIN_P_=SPIN_UP;
            BAR_SPIN_=SPIN_DN;
            BAR_SPIN_P_=SPIN_DN;
        }
        if(opr_no==36){
            //Del_3by2_m1by2_1by2_m1by2
            Jz_SPIN_=-1;
            Jz_SPIN_P_=-1;
            SPIN_=SPIN_DN;
            SPIN_P_=SPIN_DN;
            BAR_SPIN_=SPIN_UP;
            BAR_SPIN_P_=SPIN_UP;
        }
        if(opr_no==37){
            //Dummy, declared below
        }


        //index = spin*3 + orb
        Vec_1=Vec_0;
        Vec_1[(BAR_SPIN_P_*3) + YZ_] = complex<double> ((-1.0*Jz_SPIN_P_)/sqrt(6.0) , 0.0);
        Vec_1[(BAR_SPIN_P_*3) + XZ_] = complex<double> (0.0, (-1.0)/sqrt(6.0));
        Vec_1[(SPIN_P_*3) + XY_] = complex<double> ((2.0)/sqrt(6.0), 0.0);
        //        Vec_1[ (SPIN_UP*3) + YZ_] = complex<double>(0.0, -1.0/sqrt(2.0)); //*iota_comp;
        //        Vec_1[ (SPIN_UP*3) + XZ_] = complex<double>(1.0/sqrt(2.0), 0.0);

        //        Vec_1[ (SPIN_UP*3) + YZ_] = complex<double>(1.0, 0.0); //*iota_comp;
        //        Vec_1[ (SPIN_UP*3) + XZ_] = complex<double>(0.0, 0.0);

        Vec_2=Vec_0;
        //        Vec_2[ (SPIN_UP*3) + YZ_] = complex<double>(0.0, 1.0/sqrt(2.0));
        //        Vec_2[ (SPIN_UP*3) + XZ_] = complex<double>(1.0/sqrt(2.0), 0.0);
        //        Vec_2[ (SPIN_UP*3) + YZ_] = complex<double>(0.0, 0.0);
        //        Vec_2[ (SPIN_UP*3) + XZ_] = complex<double>(1.0, 0.0);

        Vec_2[(BAR_SPIN_*3) + YZ_] = complex<double> ((1.0*Jz_SPIN_)/sqrt(3.0) , 0.0);
        Vec_2[(BAR_SPIN_*3) + XZ_] = complex<double> (0.0, (-1.0)/sqrt(3.0));
        Vec_2[(SPIN_*3) + XY_] = complex<double> ((1.0)/sqrt(3.0), 0.0);

        for(int sigma=0;sigma<2;sigma++){
            for(int alpha=0;alpha<3;alpha++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int alpha_p=0;alpha_p<3;alpha_p++){

                        AMat[(sigma*3) + alpha][(sigma_p*3) + alpha_p] = Vec_1[(sigma*3) + alpha]*
                                Vec_2[(sigma_p*3) + alpha_p];
                    }
                }
            }
        }

        if(opr_no==37){
            //Dummy = Splus + Sz [Weird :)]
            AMat=AMat0;
            AMat[0][3]=(1.0*one);AMat[1][4]=(1.0*one);AMat[2][5]=(1.0*one);
            AMat[0][0]=(0.5*one);AMat[1][1]=(0.5*one);AMat[2][2]=(0.5*one);
            AMat[3][3]=(-0.5*one);AMat[4][4]=(-0.5*one);AMat[5][5]=(-0.5*one);
        }

        for(int site=0;site<basis.Length;site++){
            Get_CdaggerC_type_Opr(AMat, One_point_oprts[opr_no][site], site); //Create this
        }

    }


#endif


}

template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site){


    int i_new, j_new, m_new;
    int l,lp, max_up, min_up, min_dn, max_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    int nup_temp, ndn_temp;

    OPR.value.clear();
    OPR.rows.clear();
    OPR.columns.clear();
    OPR.nrows = basis.D_up_basis.size();
    OPR.ncols = OPR.nrows;

    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;

    int value_;
    double_type value_diagonal;


    //Diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;

        for(int sigma=0;sigma<2;sigma++){
            for(int gamma=0;gamma<3;gamma++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int gamma_p=0;gamma_p<3;gamma_p++){

                        //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma)*c(sigma_p,gamma_p)
                        if(AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p] != zero){

                            if( ((sigma*3) + gamma) == ((sigma_p*3) + gamma_p) ){
                                assert(sigma == sigma_p);
                                assert(gamma == gamma_p);
                                if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);}
                                if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);}
                                if(value_ != 0){
                                    value_diagonal += AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p]*one;
                                }
                            }

                        }

                    }
                }
            }
        }

        OPR.value.push_back(value_diagonal*one);
        OPR.rows.push_back(i);
        OPR.columns.push_back(i);
    }




    //OFF-diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){

        for(int sigma=0;sigma<2;sigma++){
            for(int gamma=0;gamma<3;gamma++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int gamma_p=0;gamma_p<3;gamma_p++){

                        //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma)*c(sigma_p,gamma_p)
                        if(AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p] != zero){

                            if( ((sigma*3) + gamma) != ((sigma_p*3) + gamma_p) ){
                                if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                    check=(bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==1);
                                    check = (check &&
                                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                    check=(bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==1);
                                    check = (check &&
                                             (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                    check=(bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site)==1);
                                    check = (check &&
                                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                    check=(bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==1);
                                    check = (check &&
                                             (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==0));
                                }


                                if(check)
                                {

                                    D_up = (int) (basis.D_up_basis[i]
                                                  + ((1-sigma)*pow(2,gamma*basis.Length + site))
                                                  - ((1-sigma_p)*pow(2,gamma_p*basis.Length + site)) );
                                    D_dn = (int) (basis.D_dn_basis[i]
                                                  + (sigma*pow(2,gamma*basis.Length + site))
                                                  - (sigma_p*pow(2,gamma_p*basis.Length + site)) );

                                    //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                    if(basis.Restricted==true){
                                        nup_temp = __builtin_popcount(D_up);
                                        ndn_temp = __builtin_popcount(D_dn);
                                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    }
                                    else{
                                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                    }

                                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                    //m_new = basis.D_dn_basis.size()*i_new + j_new;


                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site;

                                    max_up = 3*basis.Length -1;
                                    min_up = 0;
                                    max_dn = 3*basis.Length -1;
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

                                    //try this as well
                                    /*
                                    sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                    sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                    */

                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

                                    //assert(m_new<m);

                                    OPR.value.push_back(sign_FM*AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p]*one);
                                    OPR.rows.push_back(m_new);
                                    OPR.columns.push_back(i);

                                }
                            }

                        }

                    }
                }
            }
        }
    }


}


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site, int site_p){

    //sum_{sigma,sigma_p,alpha,alpha_p} AMat[sigma,alpha][sigma_p,alpha_p]
    //c_{sigma,alpha}^{dagger,site}*c_{sigma_p,alpha_p,site_p}
    int i_new, j_new, m_new;
    int l,lp, max_up, min_up, min_dn, max_dn;
    int sign_pow_up, sign_pow_dn;
    double sign_FM;
    int D_up, D_dn;
    int nup_temp, ndn_temp;

    OPR.value.clear();
    OPR.rows.clear();
    OPR.columns.clear();
    OPR.nrows = basis.D_up_basis.size();
    OPR.ncols = OPR.nrows;

    bool check;
    int SPIN_UP=0;
    int SPIN_DN=1;

    int value_;
    double_type value_diagonal;


    //Diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){
        value_diagonal=zero;

        for(int sigma=0;sigma<2;sigma++){
            for(int gamma=0;gamma<3;gamma++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int gamma_p=0;gamma_p<3;gamma_p++){

                        //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma)*c(sigma_p,gamma_p)
                        if( AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p] != zero){

                            if( ( ((sigma*3) + gamma) == ((sigma_p*3) + gamma_p) ) &&
                                    (site == site_p)
                                    ){
                                assert(sigma == sigma_p);
                                assert(gamma == gamma_p);
                                if(sigma==SPIN_UP){value_=bit_value(basis.D_up_basis[i],gamma*basis.Length + site);}
                                if(sigma==SPIN_DN){value_=bit_value(basis.D_dn_basis[i],gamma*basis.Length + site);}
                                if(value_ != 0){
                                    value_diagonal += AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p]*one;
                                }
                            }

                        }

                    }
                }
            }
        }

        OPR.value.push_back(value_diagonal*one);
        OPR.rows.push_back(i);
        OPR.columns.push_back(i);
    }




    //OFF-diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){

        for(int sigma=0;sigma<2;sigma++){
            for(int gamma=0;gamma<3;gamma++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int gamma_p=0;gamma_p<3;gamma_p++){

                        //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma,site)*c(sigma_p,gamma_p,site_p)
                        if(AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p] != zero){

                            if( ( ((sigma*3) + gamma) != ((sigma_p*3) + gamma_p) ) ||
                                    (site != site_p)
                                    ){
                                if(sigma==SPIN_UP && sigma_p==SPIN_UP){
                                    check=(bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site_p)==1);
                                    check = (check &&
                                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_DN && sigma_p==SPIN_DN){
                                    check=(bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site_p)==1);
                                    check = (check &&
                                             (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_UP && sigma_p==SPIN_DN){
                                    check=(bit_value(basis.D_dn_basis[i],gamma_p*basis.Length + site_p)==1);
                                    check = (check &&
                                             (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0));
                                }
                                if(sigma==SPIN_DN && sigma_p==SPIN_UP){
                                    check=(bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site_p)==1);
                                    check = (check &&
                                             (bit_value(basis.D_dn_basis[i],gamma*basis.Length + site)==0));
                                }


                                if(check)
                                {

                                    D_up = (int) (basis.D_up_basis[i]
                                                  + ((1-sigma)*pow(2,gamma*basis.Length + site))
                                                  - ((1-sigma_p)*pow(2,gamma_p*basis.Length + site_p)) );
                                    D_dn = (int) (basis.D_dn_basis[i]
                                                  + (sigma*pow(2,gamma*basis.Length + site))
                                                  - (sigma_p*pow(2,gamma_p*basis.Length + site_p)) );

                                    //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                    if(basis.Restricted==true){
                                        nup_temp = __builtin_popcount(D_up);
                                        ndn_temp = __builtin_popcount(D_dn);
                                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    }
                                    else{
                                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                    }

                                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                    //m_new = basis.D_dn_basis.size()*i_new + j_new;


                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    max_up = 3*basis.Length -1;
                                    min_up = 0;
                                    max_dn = 3*basis.Length -1;
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

                                    //try this as well
                                    /*
                                    sign_pow_up = one_bits_in_bw(l,min_up, basis.D_up_basis[i]) + bit_value(basis.D_up_basis[i],min_up);
                                    sign_pow_dn = one_bits_in_bw(max_dn,lp, basis.D_up_basis[i])+ bit_value(basis.D_dn_basis[i],max_dn);
                                    */

                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

                                    //assert(m_new<m);

                                    OPR.value.push_back(sign_FM*AMat[(sigma*3) + gamma][(sigma_p*3) + gamma_p]*one);
                                    OPR.rows.push_back(m_new);
                                    OPR.columns.push_back(i);

                                }
                            }

                        }

                    }
                }
            }
        }
    }


}


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_ExcitonCoherence_Length(Mat_1_doub &vector_used){

#ifdef USE_COMPLEX
    Matrix_COO OPR_;
    Mat_1_doub temp_vec;
    double_type Num_, Den_, value_;

    Mat_2_doub AMat, AMat0; //[sigma,alpha][sigma_p,beta]
    AMat.resize(6);AMat0.resize(6);
    for(int ind=0;ind<6;ind++){
        AMat[ind].resize(6);
        AMat0[ind].resize(6);
    }
    Mat_1_doub Vec_0, Vec_1, Vec_2;
    Vec_0.resize(6); Vec_1.resize(6); Vec_2.resize(6);
    int SPIN_UP=0;
    int SPIN_DN=1;
    int XZ_,YZ_,XY_;
    XZ_=0;YZ_=1;XY_=2;

    int Jz_SPIN_, Jz_SPIN_P_;
    int SPIN_, SPIN_P_;
    int BAR_SPIN_, BAR_SPIN_P_;

    for(int type=0;type<2;type++){
        if(type==0){
            //a^{dagger}_{1/2, 1/2, site} a_{3/2, 1/2, site_p}
            Jz_SPIN_=1;
            Jz_SPIN_P_=1;
            SPIN_=SPIN_UP;
            SPIN_P_=SPIN_UP;
            BAR_SPIN_=SPIN_DN;
            BAR_SPIN_P_=SPIN_DN;
        }
        if(type==1){
            //a^{dagger}_{1/2, -1/2, site} a_{3/2, -1/2, site_p}
            Jz_SPIN_=-1;
            Jz_SPIN_P_=-1;
            SPIN_=SPIN_DN;
            SPIN_P_=SPIN_DN;
            BAR_SPIN_=SPIN_UP;
            BAR_SPIN_P_=SPIN_UP;
        }

        //index = spin*3 + orb
        Vec_1=Vec_0;
        Vec_2=Vec_0;


        Vec_1[(BAR_SPIN_*3) + YZ_] = complex<double> ((1.0*Jz_SPIN_)/sqrt(3.0) , 0.0);
        Vec_1[(BAR_SPIN_*3) + XZ_] = complex<double> (0.0, (1.0)/sqrt(3.0));
        Vec_1[(SPIN_*3) + XY_] = complex<double> ((1.0)/sqrt(3.0), 0.0);

        Vec_2[(BAR_SPIN_P_*3) + YZ_] = complex<double> ((-1.0*Jz_SPIN_P_)/sqrt(6.0) , 0.0);
        Vec_2[(BAR_SPIN_P_*3) + XZ_] = complex<double> (0.0, (1.0)/sqrt(6.0));
        Vec_2[(SPIN_P_*3) + XY_] = complex<double> ((2.0)/sqrt(6.0), 0.0);


        /*
        Vec_1[(BAR_SPIN_P_*3) + YZ_] = complex<double> ((-1.0*Jz_SPIN_P_)/sqrt(6.0) , 0.0);
        Vec_1[(BAR_SPIN_P_*3) + XZ_] = complex<double> (0.0, (-1.0)/sqrt(6.0));
        Vec_1[(SPIN_P_*3) + XY_] = complex<double> ((2.0)/sqrt(6.0), 0.0);

        Vec_2[(BAR_SPIN_*3) + YZ_] = complex<double> ((1.0*Jz_SPIN_)/sqrt(3.0) , 0.0);
        Vec_2[(BAR_SPIN_*3) + XZ_] = complex<double> (0.0, (-1.0)/sqrt(3.0));
        Vec_2[(SPIN_*3) + XY_] = complex<double> ((1.0)/sqrt(3.0), 0.0);
        */

        for(int sigma=0;sigma<2;sigma++){
            for(int alpha=0;alpha<3;alpha++){

                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    for(int alpha_p=0;alpha_p<3;alpha_p++){

                        AMat[(sigma*3) + alpha][(sigma_p*3) + alpha_p] = Vec_1[(sigma*3) + alpha]*
                                Vec_2[(sigma_p*3) + alpha_p];
                    }
                }
            }
        }


        cout<<"a^{dagger}(1/2,m)a(3/2,m):----------------------------------"<<endl;
        Num_=zero;
        Den_=zero;
        for(int site=0;site<basis.Length;site++){
            for(int site_p=0;site_p<basis.Length;site_p++){
                Get_CdaggerC_type_Opr(AMat, OPR_, site, site_p);

                Matrix_COO_vector_multiplication("cx", OPR_, vector_used, temp_vec);
                value_ = dot_product(temp_vec,vector_used);

                Num_ += (pow( abs(site-site_p),2 )*1.0)*(value_*conj(value_));
                Den_ += value_*conj(value_);

                cout<<value_<<" ";

            }
            cout<<endl;
        }
        cout<<"-------------------------------------------------------------"<<endl;

        cout<<"Coherence Length, r_coh(m=";
        if(Jz_SPIN_==1){
            cout<<"1/2";
        }
        if(Jz_SPIN_==-1){
            cout<<"-1/2";
        }
        //cout<<") = "<< sqrt(Num_.real()/ Den_.real())<<endl;
        cout<<") = "<< Num_.real()<<"  "<< Den_.real()<<"   "<< sqrt(Num_.real()/ Den_.real())<<endl;
    }

#endif
}


template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_Delta_Matrix(LANCZOS & lanczos){

    double_type value =zero;
    Mat_1_doub vec_out, vec_in;
    Matrix_COO OPR_;
    int site_=0;

    cout<<scientific<<setprecision(4);

    for(int opr_no=33;opr_no<=36;opr_no++){
        cout<<"Matrix for "<<one_point_obs[opr_no]<<"[site="<<site_<<"] in eigenstates:"<<endl;

        cout<<"------------------------------------------------------------------"<<endl;
        for(int row=0;row<lanczos.Eig_vecs.size();row++){
            for(int col=0;col<lanczos.Eig_vecs.size();col++){

                OPR_ =  One_point_oprts[opr_no][site_];
                vec_in = lanczos.Eig_vecs[col];
                Matrix_COO_vector_multiplication("XX", OPR_, vec_in , vec_out);
                value=dot_product(vec_out, lanczos.Eig_vecs[row]);
                if(abs(value)<0.0001){
                    value = zero;
                }
                cout<<value<<" ";
            }
            cout<<endl;
        }

        cout<<"------------------------------------------------------------------"<<endl;

    }

}

template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Initialize_two_point_to_calculate(){
    two_point_obs.resize(4);
    two_point_obs[0]="SzSz";
    two_point_obs[1]="SpSm";
    two_point_obs[2]="SmSp";
    two_point_obs[3]="a^{dagger}{1by2,1by2}a_{3by2,1by2}";
    Two_point_oprts.resize(4);


    int T_no_oprs=3;
    int nup_temp, ndn_temp;



    for(int i=0;i<T_no_oprs;i++){
        Two_point_oprts[i].resize(basis.Length);
        for(int j=0;j<basis.Length;j++){
            Two_point_oprts[i][j].resize(basis.Length);
        }
    }




    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        if(two_point_obs[opr_no]=="SzSz"){




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

                        m=i;
                        j=i;

                        value=0;
                        for(int gamma=0;gamma<3;gamma++){
                            for(int gamma_p=0;gamma_p<3;gamma_p++){
                                value+=0.25*( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                              ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2) -
                                                bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2) )
                                              );

                            }
                        }



                        if(value!=0){
                            Two_point_oprts[opr_no][site][site2].value.push_back(value*one);
                            Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                            Two_point_oprts[opr_no][site][site2].columns.push_back(m);
                        }

                    }

                }

            }

        }


        if(two_point_obs[opr_no]=="SpSm"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){

                        m=i;
                        j=i;


                        for(int gamma=0;gamma<3;gamma++){
                            for(int gamma_p=0;gamma_p<3;gamma_p++){
                                //Sp_site_gamma[site]*Sm_site_gamma_p[site2]  Hunds coupling:
                                //there have to be ony up electron in gamma_p, site2
                                //there have to be only down electron in gamma, site

                                if(((bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)
                                    &&
                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==0)
                                    )
                                        &&
                                        ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==1)
                                         &&
                                         (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==0)
                                         ))
                                {

                                    D_up = (int) (basis.D_up_basis[i] - pow(2,gamma_p*basis.Length + site2)
                                                  + pow(2,gamma*basis.Length + site) );
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site2)
                                                  - pow(2,gamma*basis.Length + site) );

                                    //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                                    if(basis.Restricted==true){
                                        nup_temp = __builtin_popcount(D_up);
                                        ndn_temp = __builtin_popcount(D_dn);
                                        //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                        //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    }
                                    else{
                                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                    }

                                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                    //m_new = basis.D_dn_basis.size()*i_new + j_new;

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site2;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                    //assert(m_new<m);

                                    Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM*one);
                                    Two_point_oprts[opr_no][site][site2].rows.push_back(m_new);
                                    Two_point_oprts[opr_no][site][site2].columns.push_back(m);


                                }

                                if((site==site2) && (gamma==gamma_p)){


                                    if(
                                            ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==1)
                                             &&
                                             (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==0)
                                             )
                                            )
                                    {
                                        Two_point_oprts[opr_no][site][site2].value.push_back(one);
                                        Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                                        Two_point_oprts[opr_no][site][site2].columns.push_back(m);

                                    }



                                }

                            }
                        }

                    }

                }

            }

        }




        if(two_point_obs[opr_no]=="SmSp"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){

                        m=i;
                        j=i;


                        for(int gamma=0;gamma<3;gamma++){
                            for(int gamma_p=0;gamma_p<3;gamma_p++){
                                //Sm_site_gamma[site]*Sp_site_gamma_p[site2]  Hunds coupling:
                                //there have to be ony up electron in gamma_p, site2
                                //there have to be only down electron in gamma, site

                                if(((bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==0)
                                    &&
                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                                    )
                                        &&
                                        ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==0)
                                         &&
                                         (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==1)
                                         ))
                                {

                                    D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site2)
                                                  - pow(2,gamma*basis.Length + site) );
                                    D_dn = (int) (basis.D_dn_basis[j] - pow(2,gamma_p*basis.Length + site2)
                                                  + pow(2,gamma*basis.Length + site) );

                                    //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                    //m_new = basis.D_dn_basis.size()*i_new + j_new;
                                    if(basis.Restricted==true){
                                        nup_temp = __builtin_popcount(D_up);
                                        ndn_temp = __builtin_popcount(D_dn);
                                        // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                                        //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                                        m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                                    }
                                    else{
                                        i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                                        j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                                        m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                                basis.Nup_offsets[__builtin_popcount(D_up)].first;
                                    }

                                    //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                                    //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site2;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                    //assert(m_new<m);

                                    Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM*one);
                                    Two_point_oprts[opr_no][site][site2].rows.push_back(m_new);
                                    Two_point_oprts[opr_no][site][site2].columns.push_back(m);


                                }

                                if((site==site2) && (gamma==gamma_p)){


                                    if(
                                            ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==0)
                                             &&
                                             (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==1)
                                             )
                                            )
                                    {
                                        Two_point_oprts[opr_no][site][site2].value.push_back(one);
                                        Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                                        Two_point_oprts[opr_no][site][site2].columns.push_back(m);

                                    }



                                }

                            }
                        }

                    }

                }

            }

        }
        //HERE::
        //     if(two_point_obs[opr_no]=="a^{\dagger}{1by2,1by2}a_{3by2,1by2}"){

        //   }



    }


}



template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Initialize_Opr_for_Dynamics(LANCZOS &lanczos_GS){

    int nup_temp, ndn_temp;
    if(Dyn_Momentum_Resolved){


        Hamiltonian_1_COO Oprs_local;
        Oprs_local.resize(basis.Length);



#ifdef USE_COMPLEX
        complex<double> iota_temp(0.0,1.0);

        if(Dyn_opr_string == "Sz"){


            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].nrows = basis.D_up_basis.size();
                Oprs_local[site].ncols = Oprs_local[site].nrows;
            }

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;
            double value;
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

                for (int i=0;i<basis.D_up_basis.size();i++){

                    m=i;

                    //Sz_total[site]:
                    value=0;
                    for (int orb=0;orb<3;orb++){
                        value = value + 0.5*(bit_value(basis.D_up_basis[i],orb*basis.Length + site)
                                             - bit_value(basis.D_dn_basis[j],orb*basis.Length + site) );
                    }

                    if(value!=0){
                        Oprs_local[site].value.push_back(value*one);
                        Oprs_local[site].rows.push_back(m);
                        Oprs_local[site].columns.push_back(m);
                    }

                }

            }
            //local operators created ----------------------------------------------

            //In Momentum space--only OBC---------------------------------------------------

            Matrix_COO temp;
            temp=Oprs_local[0];
            complex<double> value1, value2;
            for(int site=0;site<basis.Length-1;site++){
                //value2=sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                value2 = exp((site+1)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                if(site==0){
                    value1=exp((site)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                    //value1=sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    Sum(temp, Oprs_local[site+1], temp, value1, value2);}
                else{
                    Sum(temp, Oprs_local[site+1], temp, 1.0, value2);
                }

            }

            Dyn_opr=temp;

            //----------------------------------------------------------------------


            temp.value.clear();
            temp.rows.clear();
            temp.columns.clear();
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

            }





        }



        if(Dyn_opr_string == "Splus"){


            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].nrows = basis.D_up_basis.size();
                Oprs_local[site].ncols = Oprs_local[site].nrows;
            }

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;
            double value;
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

                Oprs_local[site] = One_point_oprts[18][site];

            }
            //local operators created ----------------------------------------------

            //In Momentum space--only OBC---------------------------------------------------

            Matrix_COO temp;
            temp=Oprs_local[0];
            complex<double> value1, value2;
            for(int site=0;site<basis.Length-1;site++){
                //value2=sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                value2 = exp((site+1)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                if(site==0){
                    value1=exp((site)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                    //value1=sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    Sum(temp, Oprs_local[site+1], temp, value1, value2);}
                else{
                    Sum(temp, Oprs_local[site+1], temp, 1.0, value2);
                }

            }

            Dyn_opr=temp;

            //----------------------------------------------------------------------


            temp.value.clear();
            temp.rows.clear();
            temp.columns.clear();
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

            }

        }

        if(Dyn_opr_string == "Sminus"){


            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].nrows = basis.D_up_basis.size();
                Oprs_local[site].ncols = Oprs_local[site].nrows;
            }

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>
            int m,j;
            double value;
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

                Oprs_local[site] = One_point_oprts[19][site];

            }
            //local operators created ----------------------------------------------

            //In Momentum space--only OBC---------------------------------------------------

            Matrix_COO temp;
            temp=Oprs_local[0];
            complex<double> value1, value2;
            for(int site=0;site<basis.Length-1;site++){
                //value2=sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                value2 = exp((site+1)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                if(site==0){
                    value1=exp((site)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                    //value1=sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    Sum(temp, Oprs_local[site+1], temp, value1, value2);}
                else{
                    Sum(temp, Oprs_local[site+1], temp, 1.0, value2);
                }

            }

            Dyn_opr=temp;

            //----------------------------------------------------------------------


            temp.value.clear();
            temp.rows.clear();
            temp.columns.clear();
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

            }

        }



        cout<<"here"<<endl;

        Matrix_COO Iden_temp;

        if(Dyn_opr_string == "Delta_SOC"){


            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>

            for(int site=0;site<basis.Length;site++){
                Iden_temp=Identity_COO(One_point_oprts[12][site].nrows,One_point_oprts[12][site].ncols);
                Sum(One_point_oprts[12][site],Iden_temp , Oprs_local[site], one, (-1.0*lanczos_GS.One_point_observables_values[12][site]));
                //                Oprs_local[site] = One_point_oprts[12][site];


            }
            //local operators created ----------------------------------------------

            //In Momentum space--only PBC---------------------------------------------------




            Matrix_COO temp;
            temp=Oprs_local[0];
            complex<double> value1, value2;
            for(int site=0;site<basis.Length-1;site++){
                value2=exp((site+1)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                if(site==0){
                    value1=exp((site)*Dyn_Momentum*PI*iota_temp)*sqrt(1.0/(basis.Length));
                    Sum(temp, Oprs_local[site+1], temp, value1, value2);}
                else{
                    Sum(temp, Oprs_local[site+1], temp, one, value2);
                }

            }

            Dyn_opr=temp;


            //----------------------------------------------------------------------


            temp.value.clear();
            temp.rows.clear();
            temp.columns.clear();
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

            }


        }



        if(Dyn_opr_string == "Delta_SOC_Local"){

            int Local_site=int((basis.Length/2.0) + 0.5);

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>

            for(int site=0;site<basis.Length;site++){
                //Iden_temp=Identity_COO(One_point_oprts[12][site].nrows,One_point_oprts[12][site].ncols);
                //Sum(One_point_oprts[12][site],Iden_temp , Oprs_local[site], one, (-1.0*lanczos_GS.One_point_observables_values[12][site]));
                Oprs_local[site] = One_point_oprts[12][site];


            }
            //local operators created ----------------------------------------------

            Dyn_opr=Oprs_local[Local_site];


            //----------------------------------------------------------------------


        }


        if(Dyn_opr_string == "Sz_Local"){

            int Local_site=int((basis.Length/2.0) - 0.5);
            cout <<"LOCAL_SITE = "<<Local_site<<endl;

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>

            for(int site=0;site<basis.Length;site++){
                //Iden_temp=Identity_COO(One_point_oprts[12][site].nrows,One_point_oprts[12][site].ncols);
                //Sum(One_point_oprts[12][site],Iden_temp , Oprs_local[site], one, (-1.0*lanczos_GS.One_point_observables_values[12][site]));
                Oprs_local[site] = One_point_oprts[17][site];


            }
            //local operators created ----------------------------------------------

            Dyn_opr=Oprs_local[Local_site];


            //----------------------------------------------------------------------


        }

        if(Dyn_opr_string == "Splus_Local"){

            int Local_site=int((basis.Length/2.0) - 0.5);
            cout <<"LOCAL_SITE = "<<Local_site<<endl;

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>

            for(int site=0;site<basis.Length;site++){
                //Iden_temp=Identity_COO(One_point_oprts[12][site].nrows,One_point_oprts[12][site].ncols);
                //Sum(One_point_oprts[12][site],Iden_temp , Oprs_local[site], one, (-1.0*lanczos_GS.One_point_observables_values[12][site]));
                Oprs_local[site] = One_point_oprts[18][site];


            }
            //local operators created ----------------------------------------------

            Dyn_opr=Oprs_local[Local_site];


            //----------------------------------------------------------------------


        }


#endif




    }

    else{
        int gamma;
        if(Dyn_opr_string == "J_0"){
            gamma=0;
        }
        if(Dyn_opr_string == "J_1"){
            gamma=1;
        }
        if(Dyn_opr_string == "J_2"){
            gamma=2;
        }

        if(Dyn_opr_string == "J_0" || Dyn_opr_string == "J_1" || Dyn_opr_string == "J_2"){
            Matrix_COO Opr;
            Opr.nrows=basis.D_up_basis.size();
            Opr.ncols=Opr.nrows;
            Opr.value.clear();
            Opr.columns.clear();
            Opr.rows.clear();


            double value;
            int m,j;
            int D_up,D_dn;
            int i_new,j_new;
            int m_new;
            double sign_FM;
            int sign_pow_up, sign_pow_dn;
            int l,lp;
            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                value=0;


                for(int site=0;site<basis.Length-1 ;site++){

                    int site_p= site +1;


                    //---------------Hopping for up electrons-------------------//
                    //there have to be one up electron in gamma_p, site
                    //there have to be no up electron in gamma, site
                    if(
                            (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_up_basis[i],gamma*basis.Length + site_p)==0)
                            )
                    {

                        D_up = (int) (basis.D_up_basis[i] + pow(2,gamma*basis.Length + site_p)
                                      - pow(2,gamma*basis.Length + site) );
                        D_dn = basis.D_dn_basis[j];


                        //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                        //j_new = j;

                        //m_new = basis.D_dn_basis.size()*i_new + j_new;
                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            //m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //      basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                        l=gamma*basis.Length + site;
                        lp=gamma*basis.Length + site_p;

                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                        sign_FM = pow(-1.0, sign_pow_up);



                        if((Hopping_mat_NN[gamma][gamma])!=0){

                            Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[gamma][gamma])*one);
                            Opr.rows.push_back((m_new));
                            Opr.columns.push_back((m));
                            Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma][gamma])*one);
                            Opr.rows.push_back((m));
                            Opr.columns.push_back((m_new));

                        }

                    } // if up hopping possible


                    //---------------Hopping for dn electrons-------------------//
                    //there have to be one dn electron in gamma, site
                    //there have to be no dn electron in gamma, site
                    if(
                            (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site_p)==0)
                            )
                    {

                        D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma*basis.Length + site_p)
                                      - pow(2,gamma*basis.Length + site) );
                        D_up = basis.D_up_basis[i];

                        //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                        //i_new = i;

                        //m_new = basis.D_dn_basis.size()*i_new + j_new;

                        if(basis.Restricted==true){
                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            //  m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //        basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                        }
                        else{
                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;
                        }

                        //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                        //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                        l=gamma*basis.Length + site;
                        lp=gamma*basis.Length + site_p;

                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                        sign_FM = pow(-1.0, sign_pow_dn);



                        if((Hopping_mat_NN[gamma][gamma])!=0){

                            Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[gamma][gamma])*one);
                            Opr.rows.push_back((m_new));
                            Opr.columns.push_back((m));
                            Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma][gamma])*one);
                            Opr.rows.push_back((m));
                            Opr.columns.push_back((m_new));


                        }

                    } // if up hopping possible






                } // site

            } // "i" i.e up_decimals



            Dyn_opr=Opr;
            Opr.value.clear();
            Opr.columns.clear();
            Opr.rows.clear();


        }












    }

}



template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_c_on_GS(LANCZOS & lanczos, BASIS_3_orb_Hubb_chain_GC & basis_Nm1,
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
        orb_val = TRIO_VEC[n].orb_;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        for (int i=0;i<basis.D_up_basis.size();i++){


            if(spin_val==0){
                //For c_up|GS>
                if(bit_value(basis.D_up_basis[i],orb_val*basis.Length + site_val)==1){
                    l = orb_val*basis.Length + site_val;

                    D_up_new = (int) (basis.D_up_basis[i] - pow(2,orb_val*basis.Length + site_val) );
                    D_dn_new = basis.D_dn_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Nm1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Nm1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_up = 3*basis.Length -1;
                    min_up =0 ;
                    sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                    if(l != min_up){
                        sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                    }

                    sign_FM = pow(-1.0, sign_pow_up);

                    value = sign_FM*lanczos.Eig_vec[i]*value_in;

                    State_c_on_GS[i_new] += value;


                }
            }

            if(spin_val==1){
                //For c_dn|GS>
                if(bit_value(basis.D_dn_basis[i],orb_val*basis.Length + site_val)==1){
                    l = orb_val*basis.Length + site_val;

                    D_dn_new = (int) (basis.D_dn_basis[i] - pow(2,orb_val*basis.Length + site_val) );
                    D_up_new = basis.D_up_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Nm1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Nm1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Nm1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_dn = 3*basis.Length -1;
                    min_dn=0;
                    sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[i]) ;
                    if(l != min_dn){
                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                    }
                    sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                    sign_FM = pow(-1.0, sign_pow_dn);

                    value = sign_FM*lanczos.Eig_vec[i]*value_in;

                    State_c_on_GS[i_new] += value;


                }
            }



        }


    }



}



template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Get_cdagger_on_GS(LANCZOS & lanczos, BASIS_3_orb_Hubb_chain_GC & basis_Np1,
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
        orb_val = TRIO_VEC[n].orb_;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        for (int i=0;i<basis.D_up_basis.size();i++){


            if(spin_val==0){
                //For c_dagger_up|GS>
                if(bit_value(basis.D_up_basis[i],orb_val*basis.Length + site_val)==0){
                    l = orb_val*basis.Length + site_val;

                    D_up_new = (int) (basis.D_up_basis[i] + pow(2,orb_val*basis.Length + site_val) );
                    D_dn_new = basis.D_dn_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Np1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Np1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_up = 3*basis.Length -1;
                    min_up=0;
                    sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                    if(l != min_up){
                        sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                    }

                    sign_FM = pow(-1.0, sign_pow_up);

#ifdef USE_COMPLEX
                    value = sign_FM*lanczos.Eig_vec[i]*conj(value_in);
#endif
#ifndef USE_COMPLEX
                    value = sign_FM*lanczos.Eig_vec[i]*(value_in);
#endif


                    State_cdagger_on_GS[i_new] += value;


                }
            }

            if(spin_val==1){
                //For c_dn|GS>
                if(bit_value(basis.D_dn_basis[i],orb_val*basis.Length + site_val)==0){
                    l = orb_val*basis.Length + site_val;

                    D_dn_new = (int) (basis.D_dn_basis[i] + pow(2,orb_val*basis.Length + site_val) );
                    D_up_new = basis.D_up_basis[i];

                    i_up = Find_int_in_intarray(D_up_new,basis_Np1.Canonical_partition_up[__builtin_popcount(D_up_new)]);
                    i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)]);
                    i_new = (basis_Np1.Canonical_partition_dn[__builtin_popcount(D_up_new)].size()*i_up + i_dn) +
                            basis_Np1.Nup_offsets[__builtin_popcount(D_up_new)].first;

                    max_dn = 3*basis.Length -1;
                    min_dn=0;
                    sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[i]) ;
                    if(l != min_dn){
                        sign_pow_dn += bit_value(basis.D_dn_basis[i],min_dn);
                    }
                    sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                    sign_FM = pow(-1.0, sign_pow_dn);


#ifdef USE_COMPLEX
                    value = sign_FM*lanczos.Eig_vec[i]*conj(value_in);
#endif
#ifndef USE_COMPLEX
                    value = sign_FM*lanczos.Eig_vec[i]*(value_in);
#endif

                    State_cdagger_on_GS[i_new] += value;


                }
            }



        }


    }



}

template <typename Basis_type>
void MODEL_3_orb_Hubb_chain_GC<Basis_type>::Calculate_Local_Obs_for_States_to_Look(LANCZOS & lanczos){

    bool grouping_by_orb2 =true;




    if(lanczos.calculate_local_obs_for_states_to_look == true){
        div_t divresult;

        int nup,ndn,temp_d;





        for(int Ts=0;Ts<lanczos.states_to_look.size();Ts++){

            Mat_1_int nup_2_old,ndn_2_old;
            nup_2_old.resize(basis.Length);ndn_2_old.resize(basis.Length);
            for(int site=0;site<basis.Length;site++){
                nup_2_old[site]=0;ndn_2_old[site]=0;
            }
            Mat_1_int nup_2,ndn_2;
            nup_2.resize(basis.Length);ndn_2.resize(basis.Length);

            Mat_2_int nup_2_group, ndn_2_group;
            nup_2_group.clear();ndn_2_group.clear();

            Mat_1_doub group_weight;
            group_weight.clear();



            string out0 = lanczos.file_Loc_obs_in_basis_of_states + NumberToString (lanczos.states_to_look[Ts])+ ".txt";
            ofstream file_out0(out0.c_str());

            for(int bi=0;bi<lanczos.no_basis_to_check;bi++){
                file_out0<<"#For Basis no = "<<lanczos.Overlaps[Ts][bi].second<<"["<<bi<<"]"<<endl;
                // m=basis.D_dn_basis.size()*i + j;
                divresult = div (lanczos.Overlaps[Ts][bi].second,basis.D_dn_basis.size());

                int bi_up = divresult.quot;
                int bi_dn = divresult.rem;


                for(int gamma=2;gamma>=0;gamma--){
                    file_out0<<"#";
                    for(int site=0;site<basis.Length;site++){
                        nup = bit_value(basis.D_up_basis[bi_up],gamma*basis.Length + site);
                        ndn = bit_value(basis.D_dn_basis[bi_dn],gamma*basis.Length + site);
                        if(nup != ndn){
                            temp_d=nup - ndn;
                        }
                        else{
                            temp_d = nup + ndn;
                        }

                        file_out0<<temp_d<<" ";

                        // file_out0<<"nup"<<gamma<<"["<<site<<"] = "<<bit_value(basis.D_up_basis[bi_up],gamma*basis.Length + site)<<endl;
                        // file_out0<<"ndn"<<gamma<<"["<<site<<"] = "<<bit_value(basis.D_dn_basis[bi_dn],gamma*basis.Length + site)<<endl;
                    }
                    file_out0<<endl;
                }

                file_out0<<lanczos.Overlaps[Ts][bi].first<<endl;
            }




            if(grouping_by_orb2==true){



                string out2 = lanczos.file_Loc_obs_in_basis_of_states + NumberToString (lanczos.states_to_look[Ts])+ "_groupwise.txt";
                ofstream file_out2(out2.c_str());

                for(int bi=0;bi<lanczos.no_basis_to_check;bi++){

                    divresult = div (lanczos.Overlaps[Ts][bi].second,basis.D_dn_basis.size());

                    int bi_up = divresult.quot;
                    int bi_dn = divresult.rem;

                    for(int site=0;site<basis.Length;site++){
                        nup_2[site] = bit_value(basis.D_up_basis[bi_up],2*basis.Length + site);
                        ndn_2[site] = bit_value(basis.D_dn_basis[bi_dn],2*basis.Length + site);
                    }

                    int pos;
                    if( present_before(nup_2, ndn_2, nup_2_group, ndn_2_group, pos) == false ) //if new group
                    {
                        group_weight.push_back(lanczos.Overlaps[Ts][bi].first);
                        nup_2_group.push_back(nup_2);
                        ndn_2_group.push_back(ndn_2);
                    }
                    else{
                        bool test = present_before(nup_2, ndn_2, nup_2_group, ndn_2_group, pos);
                        group_weight[pos] = group_weight[pos]  + lanczos.Overlaps[Ts][bi].first;
                    }




                }


                file_out2<<"#group_index  group_weight group_configuration"<<endl;
                for(int gi=0;gi<group_weight.size();gi++){
                    file_out2<<gi<<"  "<<group_weight[gi]<<" ";
                    for(int site=0;site<basis.Length;site++){

                        if(nup_2_group[gi][site] != ndn_2_group[gi][site]){
                            temp_d=nup_2_group[gi][site] - ndn_2_group[gi][site];
                        }
                        else{
                            temp_d = nup_2_group[gi][site] + ndn_2_group[gi][site];
                        }
                        file_out2<<temp_d<<" ";
                    }
                    file_out2<<endl;

                }



            }




        }

    }

}





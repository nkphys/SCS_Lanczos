/*
This class includes the Model for which Lanczos is being done
*/
//#ifndef USE_COMPLEX
#include "Model_2_orb_Hubbard_chain.h"
#include <stdlib.h>
#include <string>
using namespace std;
#define PI 3.14159265
/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ] [_______________________  _  ]
    site----->[012....................(L-1)] [012....................(L-1)]
    orbital-->[.........orb - "0"..........] [.........orb - "1"..........]

2)  similarly for "down spin" basis

3)  For total
    m=basis.D_dn_basis.size()*i + j;
*/



void MODEL_2_orb_Hubb_chain::Act_Hamil(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    cout<<"NOT WORKING AT PRESENT"<<endl;

}

void MODEL_2_orb_Hubb_chain::Add_diagonal_terms(BASIS_2_orb_Hubb_chain &basis){

    cout<<"Started Hamiltonian construction: Diagonal"<<endl;

    Hamil.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    Hamil.ncols = Hamil.nrows;



    //Remember H[l][m]=<l|H|m>
    int m;
    double value;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;
            //intra-orbital coulomb repulsion:
            value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);

            //inter-orbital coulomb repulsion:
            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
                    for(int site=0;site<basis.Length;site++){
                        value+=(U_p - (J_H*0.5))*
                                ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                                    bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                  ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site) +
                                    bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site) )
                                  );
                    }
                }
            }

            //SzSz Hunds coupling:
            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
                    for(int site=0;site<basis.Length;site++){
                        value+=0.25*(-J_H*2.0)*
                                ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                    bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                  ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site) -
                                    bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site) )
                                  );
                    }
                }
            }


            //Anisotropy Dz_Anisotropy * (Sz(i,0) + Sz(i,1))^2   == Sz(i,0)^2 + Sz(i,1)^2 + 2*Sz(i,0)*Sz(i,1)
            //-------------------------------------------------------------------------------------------------

            //Sz(i,0)^2 + Sz(i,1)^2
            for(int gamma=0;gamma<2;gamma++){
                for(int site=0;site<basis.Length;site++){
                    value += 0.25*Dz_Anisotropy*
                            (( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                               bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                             ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                               bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                             );

                }
            }

            //2*Sz(i,0)*Sz(i,1)
            for(int site=0;site<basis.Length;site++){
                value += 0.25*2.0*Dz_Anisotropy*
                        (( bit_value(basis.D_up_basis[i],0*basis.Length + site) -
                           bit_value(basis.D_dn_basis[j],0*basis.Length + site) )*
                         ( bit_value(basis.D_up_basis[i],1*basis.Length + site) -
                           bit_value(basis.D_dn_basis[j],1*basis.Length + site) )
                         );
            }

            //-----------------------------------------------------------------------------------------------------


            //Crystal Field Splitting (CFE):
            for(int gamma=0;gamma<2;gamma++){
                for(int site=0;site<basis.Length;site++){
                    value+=(CFS[gamma])*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                              );
                }
            }



            //magnetic Field
            for(int gamma=0;gamma<2;gamma++){
                for(int site=0;site<basis.Length;site++){
                    value+=0.5*(H_field[site])*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                              );
                }
            }





            if(value!=0){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m);
                Hamil.columns.push_back(m);
            }


            if(m%1000 ==1){
                //cout<<"done "<<m<<" basis"<<endl;
            }

        }
    }


    cout<<"Done Hamiltonian construction: Diagonal"<<endl;

}

void MODEL_2_orb_Hubb_chain::Add_non_diagonal_terms(BASIS_2_orb_Hubb_chain &basis){


    bool PAIRHOPPINGINCLUDED = true;
    if(!PAIRHOPPINGINCLUDED){
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
        cout<<"  ---------- PAIR HOPPING is switched OFF --------- "<<endl;
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    }

    cout<<"Started Hamiltonian construction: Non Diagonal"<<endl;

    Hamil.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    Hamil.ncols = Hamil.nrows;

    double value;
    int m;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM ;
    int sign_pow_up, sign_pow_dn ;
    int l,lp;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;


            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
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

                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                            assert(m_new<m);

                            Hamil.value.push_back(sign_FM*(0.5*(-J_H*2.0))*one);
                            Hamil.rows.push_back(m_new);
                            Hamil.columns.push_back(m);

                        } // if SpSm possible


                        //Pair hopping: P_i\gamma = c_i\gamma\dn c_i\gamma\dn
                        //there have to be pair present at i,gamma_p
                        //there have to nothing at i,gamma

                        if(PAIRHOPPINGINCLUDED){
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

                                i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l=gamma*basis.Length + site;
                                lp=gamma_p*basis.Length + site;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);



                                assert(m_new<m);
                                Hamil.value.push_back(sign_FM*J_H*one);
                                Hamil.rows.push_back((m_new));
                                Hamil.columns.push_back((m));


                            } //Pair-Hopping

                        }

                    } // site
                } //gamma_p
            } //gamma

            if(m%1000 ==1){
                //cout<<"done "<<m<<" basis"<<endl;
            }

        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals


    cout<<"Done Hamiltonian construction: Non Diagonal"<<endl;

}

void MODEL_2_orb_Hubb_chain::Add_connections(BASIS_2_orb_Hubb_chain &basis){


    cout<<"Started Hamiltonian construction: Connections"<<endl;

    Hamil.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    Hamil.ncols = Hamil.nrows;

    double value;
    double_type hopp_val;
    int m;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int l,lp;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;

            for(int site=0;site<basis.Length ;site++){

                for(int gamma=0;gamma<2;gamma++){

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
                                        PBC
                                        )

                                    )
                            { // nearest neighbour

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


                                    i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    j_new = j;

                                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up);



                                    if(abs(Hopping_mat_NN[gamma_p][gamma])>0.000001){
                                        assert(m_new<m);
                                        if(neigh==-1 || neigh==basis.Length-1){
                                            hopp_val=Hopping_mat_NN[gamma_p][gamma];
                                        }
                                        else{
                                            hopp_val=conjugate(Hopping_mat_NN[gamma_p][gamma]);
                                        }
                                        Hamil.value.push_back(-1.0*sign_FM*(hopp_val)*one);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));
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

                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site_p)
                                                  - pow(2,gamma*basis.Length + site) );


                                    j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                                    i_new = i;

                                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn);



                                    if(abs(Hopping_mat_NN[gamma_p][gamma])>0.00001){
                                        assert(m_new<m);

                                        if(neigh==-1 || neigh==basis.Length-1){
                                            hopp_val=Hopping_mat_NN[gamma_p][gamma];
                                        }
                                        else{
                                            hopp_val=conjugate(Hopping_mat_NN[gamma_p][gamma]);
                                        }

                                        //conjugate(hopp_val) is used to preserve time-reversal symmetry
                                        Hamil.value.push_back(-1.0*sign_FM*(conjugate(hopp_val))*one);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));}

                                } // if dn hopping possible


                            }//nearest neighbour
                        } //gamma_p

                    }//site_p

                } //gamma


            } // site

            if(m%1000 ==1){
                //cout<<"done "<<m<<" basis"<<endl;
            }

        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals


    cout<<"Done Hamiltonian construction: Connections"<<endl;

}


void MODEL_2_orb_Hubb_chain::Create_states_with_hole_hopping(BASIS_2_orb_Hubb_chain &basis,
                                                             Mat_1_doub Alpha_, Mat_1_doub &vec){



    //    //This routine should be used for only PBC
    //    assert(PBC);

    //    Mat_1_doub vec_new;
    //    vec_new.resize(vec.size());

    //    int gamma_p;
    //    double distance;
    //    int site_pp;
    //    double value;
    //    int m;
    //    int D_up,D_dn;
    //    int i_new,j_new;
    //    int m_new;
    //    double sign_FM;
    //    int sign_pow_up, sign_pow_dn;
    //    int l,lp;
    //    int neigh;



    //    for (int i=0;i<basis.D_up_basis.size();i++){
    //        for (int j=0;j<basis.D_dn_basis.size();j++){
    //            m=basis.D_dn_basis.size()*i + j;

    //            value=0;


    //            for(int site=0;site<basis.Length ;site++){

    //                for(int site_p=0;site_p<basis.Length ;site_p++){
    //                    neigh = site_p - site;
    //                    if(neigh==(basis.Length - 1)){
    //                        neigh=-1;
    //                    }
    //                    if(neigh==-(basis.Length - 1)){
    //                        neigh=1;
    //                    }

    //                    if(abs(neigh)==1){

    //                        for(int gamma=0;gamma<2;gamma++){

    //                            /*
    //                            //hopping type(a) in Elbio's Picture
    //                            if(gamma==1 && neigh ==1){
    //                                value=Alpha_[0];
    //                            }
    //                            //hopping type(b) in Elbio's Picture
    //                            if(gamma==1 && neigh ==-1){
    //                                value=Alpha_[1];
    //                            }
    //                            //hopping type(c) in Elbio's Picture
    //                            if(gamma==0 && neigh ==-1){
    //                                value=Alpha_[2];
    //                            }
    //                            //hopping type(d) in Elbio's Picture
    //                            if(gamma==0 && neigh ==1){
    //                                value=Alpha_[3];
    //                            }
    //                            */

    //                            //---------------Hopping for up electrons-------------------//
    //                            //there have to be one up electron in gamma, site
    //                            //gamma_p, site_p site have to be fully empty
    //                            if(
    //                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
    //                                    &&
    //                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site_p)==0)
    //                                    &&
    //                                    (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site_p)==0)
    //                                    )
    //                            {

    //                                gamma_p = 1 - gamma;
    //                                site_pp = site_p + neigh;
    //                                site_pp = (site_pp + basis.Length)%basis.Length;

    //                                //another hole is sitting on different orbital, and leads to holes on same rung [distance =sqrt(2), distance_new =1].
    //                                if(
    //                                        (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
    //                                        &&
    //                                        (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site)==0)
    //                                        ){
    //                                    value=Alpha_[0];
    //                                    // cout<<"site = "<<site<<"site_p = "<<site_p<<endl;


    //                                }
    //                                //another hole is sitting on different orbital, and leads to holes on next nearest diagonal [distance_new = sqrt(5)].
    //                                else if(
    //                                        (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site_pp)==0)
    //                                        &&
    //                                        (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site_pp)==0)
    //                                        ){
    //                                    value=Alpha_[1];
    //                                    //cout<<"site_pp = "<<site_pp<<"site_p = "<<site_p<<endl;


    //                                }
    //                                else{
    //                                    if(!(vec[m]==0)){
    //                                        cout<<"INITIAL STATE DID NOT HAVE HOLES ALONG DIAGONAL"<<endl;
    //                                        assert((vec[m]==0));
    //                                    }

    //                                }


    //                                D_up = (int) (basis.D_up_basis[i] + pow(2,gamma*basis.Length + site_p)
    //                                              - pow(2,gamma*basis.Length + site) );


    //                                i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
    //                                j_new = j;

    //                                m_new = basis.D_dn_basis.size()*i_new + j_new;

    //                                l=gamma*basis.Length + site;
    //                                lp=gamma*basis.Length + site_p;

    //                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

    //                                sign_FM = pow(-1.0, sign_pow_up);


    //                                vec_new[m_new] +=sign_FM*one*vec[m]*value;


    //                            } // if up hopping possible


    //                            //---------------Hopping for dn electrons-------------------//
    //                            //there have to be one dn electron in gamma, site
    //                            //there have to be no dn electron in gamma_p, site_p
    //                            if(
    //                                    (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)
    //                                    &&
    //                                    (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site_p)==0)
    //                                    &&
    //                                    (bit_value(basis.D_up_basis[i],gamma*basis.Length + site_p)==0)
    //                                    )
    //                            {


    //                                gamma_p = 1 - gamma;
    //                                site_pp = site_p + neigh;
    //                                site_pp = (site_pp + basis.Length)%basis.Length;

    //                                //another hole is sitting on different orbital, and leads to holes on same rung [distance_new =1].
    //                                if(
    //                                        (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site)==0)
    //                                        &&
    //                                        (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site)==0)
    //                                        ){
    //                                    value=Alpha_[0];




    //                                }
    //                                //another hole is sitting on different orbital, and leads to holes on next nearest diagonal [distance_new = sqrt(5)].
    //                                else if(
    //                                        (bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site_pp)==0)
    //                                        &&
    //                                        (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site_pp)==0)
    //                                        ){
    //                                    value=Alpha_[1];


    //                                }
    //                                else{
    //                                    if(!(vec[m]==0)){
    //                                        cout<<"INITIAL STATE DID NOT HAVE HOLES ALONG DIAGONAL"<<endl;
    //                                        assert((vec[m]==0));
    //                                    }

    //                                }



    //                                D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma*basis.Length + site_p)
    //                                              - pow(2,gamma*basis.Length + site) );


    //                                j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
    //                                i_new = i;

    //                                m_new = basis.D_dn_basis.size()*i_new + j_new;

    //                                l=gamma*basis.Length + site;
    //                                lp=gamma*basis.Length + site_p;

    //                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

    //                                sign_FM = pow(-1.0, sign_pow_dn);

    //                                vec_new[m_new] +=sign_FM*one*vec[m]*value;

    //                            } // if dn hopping possible


    //                        } //gamma


    //                    }

    //                }//site_p


    //            } // site

    //        }// "j" i.e dn_decimals
    //    } // "i" i.e up_decimals





    //    vec = vec_new;




}




void MODEL_2_orb_Hubb_chain::Discard_double_occupancies(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub &vec){


    int m;

    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            for(int site=0;site<basis.Length ;site++){

                for(int gamma=0;gamma<2;gamma++){

                    if(
                            (bit_value(basis.D_up_basis[i],gamma*basis.Length + site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[j],gamma*basis.Length + site)==1)

                            ){
                        vec[m] = zero;

                    }




                }
            }}}

}


void MODEL_2_orb_Hubb_chain::Check_orbital_symmetry(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_){

    Mat_1_doub vec1_;
    double norm_;
    double norm1_;
    double_type overlap_temp,overlap_temp2 ;
    double eps_ = 0.00001;


    norm_ = abs(dot_product(vec_,vec_));
    for(int i=0;i<vec_.size();i++){
        vec_[i] = vec_[i]/sqrt(norm_);
    }

    vec1_ = Act_Orbital_Exchange(basis, vec_);

    norm1_ = abs(dot_product(vec1_,vec1_));
    for(int i=0;i<vec_.size();i++){
        vec1_[i] = vec1_[i]/sqrt(norm1_);
    }

    overlap_temp = dot_product(vec_,vec1_);
    overlap_temp2 = overlap_temp*conjugate(overlap_temp);

    // if(  abs((overlap_temp2) - 1.0) > eps_ ){
    //   cout <<"This state is not an eigenstate of Orbital Symmetry Operator"<<endl;
    //}

    //else{
    cout <<"This state is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;
    // }



}


Mat_1_doub MODEL_2_orb_Hubb_chain::Act_Translational_operator(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_){



    //bit_value(basis.D_dn_basis[j],o0*basis.Length + site0)
    Mat_1_doub vec1_;
    vec1_.resize(vec_.size());


    int index_no, index_no_new;
    int i_new, j_new, l,l_new;
    int gamma, site_new;
    int sign_pow_dn;
    double sign_FM_dn;
    int D_dn_new;
    int sign_pow_up;
    double sign_FM_up;
    int D_up_new;
    int site;


    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            index_no=basis.D_dn_basis.size()*i + j;

            D_dn_new=0;
            D_up_new=0;


            sign_FM_dn=1.0;
            sign_FM_up=1.0;

            for( gamma=0;gamma<2;gamma++){
                for( site=0;site<basis.Length;site++){
                    site_new = site + 1;
                    site_new = (site_new)%basis.Length;

                    l= gamma*basis.Length + site;
                    l_new = gamma*basis.Length + site_new;


                    //--moving dn spin--
                    if(site == basis.Length -1){
                        sign_pow_dn = one_bits_in_bw( gamma*basis.Length, gamma*basis.Length + basis.Length -1,
                                                      basis.D_dn_basis[j]) +
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length);
                        if(bit_value(basis.D_dn_basis[j],l)==1){
                            sign_FM_dn = sign_FM_dn*pow(-1.0, sign_pow_dn);
                        }
                    }


                    D_dn_new += pow(2,l_new)*bit_value(basis.D_dn_basis[j],l);
                    //----------------------

                    //--moving up spin--
                    if(site == basis.Length -1){
                        sign_pow_up = one_bits_in_bw( gamma*basis.Length, gamma*basis.Length + basis.Length -1,
                                                      basis.D_up_basis[i]) +
                                bit_value(basis.D_up_basis[i],gamma*basis.Length);
                        if(bit_value(basis.D_up_basis[i],l)==1){
                            sign_FM_up = sign_FM_up*pow(-1.0, sign_pow_up);
                        }
                    }


                    D_up_new += pow(2,l_new)*bit_value(basis.D_up_basis[i],l);
                    //----------------------


                }
            }

            j_new = Find_int_in_intarray(D_dn_new,basis.D_dn_basis);
            i_new = Find_int_in_intarray(D_up_new,basis.D_up_basis);
            index_no_new=basis.D_dn_basis.size()*i_new + j_new;

            vec1_[index_no_new] = sign_FM_up*sign_FM_dn*vec_[index_no];

            //if(vec_[index_no] !=0){
            //  cout<<index_no<<"  "<<vec_[index_no]<<"  "<<vec1_[index_no_new]<<"  "<<sign_FM_up<<"  "<<sign_FM_dn<<endl;
            //}





        }}


    return vec1_;

}

Mat_1_doub MODEL_2_orb_Hubb_chain::Act_Orbital_Exchange(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_){



    //bit_value(basis.D_dn_basis[j],o0*basis.Length + site0)
    Mat_1_doub vec1_;
    vec1_.resize(vec_.size());


    int index_no, index_no_new;
    int i_new, j_new, l,l_new;
    int gamma, gamma_new;
    int sign_pow_dn;
    double sign_FM_dn;
    int D_dn_new;
    int sign_pow_up;
    double sign_FM_up;
    int D_up_new;
    int site;


    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            index_no=basis.D_dn_basis.size()*i + j;

            sign_FM_dn=1.0;
            sign_FM_up=1.0;

            D_dn_new=0;
            D_up_new=0;

            for( gamma=0;gamma<2;gamma++){
                for( site=0;site<basis.Length;site++){
                    gamma_new = abs(gamma -1);
                    l= gamma*basis.Length + site;
                    l_new = gamma_new*basis.Length + site;


                    //--orbital exchange dn spin--
                    sign_pow_dn = one_bits_in_bw(l,l_new,basis.D_dn_basis[j]);
                    if(bit_value(basis.D_dn_basis[j],l) != bit_value(basis.D_dn_basis[j],l_new)){
                        sign_FM_dn = sign_FM_dn*pow(-1.0, sign_pow_dn);
                    }
                    else if(bit_value(basis.D_dn_basis[j],l)==0){
                        assert(bit_value(basis.D_dn_basis[j],l)==0);
                        assert(bit_value(basis.D_dn_basis[j],l_new)==0);
                        sign_FM_dn = sign_FM_dn*1.0;
                    }
                    else{
                        assert(bit_value(basis.D_dn_basis[j],l)==1);
                        assert(bit_value(basis.D_dn_basis[j],l_new)==1);
                        sign_FM_dn = sign_FM_dn*(-1.0);
                    }

                    D_dn_new += pow(2,l_new)*bit_value(basis.D_dn_basis[j],l);
                    //----------------------

                    //--orbital exchange up spin--
                    sign_pow_up = one_bits_in_bw(l,l_new,basis.D_up_basis[i]);
                    if(bit_value(basis.D_up_basis[i],l) != bit_value(basis.D_up_basis[i],l_new)){
                        sign_FM_up = sign_FM_up*pow(-1.0, sign_pow_up);
                    }
                    else if(bit_value(basis.D_up_basis[i],l)==0){
                        assert(bit_value(basis.D_up_basis[i],l)==0);
                        assert(bit_value(basis.D_up_basis[i],l_new)==0);
                        sign_FM_up = sign_FM_up*1.0;
                    }
                    else{
                        assert(bit_value(basis.D_up_basis[i],l)==1);
                        assert(bit_value(basis.D_up_basis[i],l_new)==1);
                        sign_FM_up = sign_FM_up*(-1.0);
                    }

                    D_up_new += pow(2,l_new)*bit_value(basis.D_up_basis[i],l);
                    //----------------------


                }
            }

            j_new = Find_int_in_intarray(D_dn_new,basis.D_dn_basis);
            i_new = Find_int_in_intarray(D_up_new,basis.D_up_basis);
            index_no_new=basis.D_dn_basis.size()*i_new + j_new;

            vec1_[index_no_new] = sign_FM_up*sign_FM_dn*vec_[index_no];





        }}


    return vec1_;
}


void MODEL_2_orb_Hubb_chain::Get_Pair_Operator_Matrix(BASIS_2_orb_Hubb_chain &basis, BASIS_2_orb_Hubb_chain &basis_nm2,
                                                      trio_int trio_0,trio_int trio_1, double value_out_ ){

    //Convention of state  (---upupupup---,--dndndndn-----)|0>

    int o1,o0,site0,site1,s0,s1;

    double value=0.0;

    //c_{trio_1}c_{trio_0}
    o0=trio_0.orb_;site0=trio_0.site_;s0=trio_0.spin_;
    o1=trio_1.orb_;site1=trio_1.site_;s1=trio_1.spin_;
    double value_out = value_out_;


    int column_no, row_no;

    int D_up_new, D_dn_new ,i_new, j_new, l0, l1;
    double sign_pow_dn, sign_pow_up;


    Pair_Annihilation.ncols=basis.D_up_basis.size()*basis.D_dn_basis.size();
    Pair_Annihilation.nrows=basis_nm2.D_up_basis.size()*basis_nm2.D_dn_basis.size();
    Pair_Annihilation.value.clear();
    Pair_Annihilation.columns.clear();Pair_Annihilation.rows.clear();




    //c_[o1][site1][s1]c_[o0][site0][s0]
    assert(s1 != s0);
    if(s0==0){
        swap(s0,s1);
        swap(o0,o1);
        swap(site0,site1);
        value_out = -1.0*value_out_;
    }
    else{
        value_out=1.0*value_out_;
    }

    //assuming s0=1=dn, s1=0=up
    assert(s0==1);


    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            column_no=basis.D_dn_basis.size()*i + j;


            value=1.0;

            if(
                    (bit_value(basis.D_dn_basis[j],o0*basis.Length + site0)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],o1*basis.Length + site1)==1)
                    )
            {


                D_up_new = (int) (basis.D_up_basis[i] - pow(2,o1*basis.Length + site1) );
                D_dn_new = (int) (basis.D_dn_basis[j] - pow(2,o0*basis.Length + site0) );

                i_new = Find_int_in_intarray(D_up_new,basis_nm2.D_up_basis);
                j_new = Find_int_in_intarray(D_dn_new,basis_nm2.D_dn_basis);

                row_no = basis_nm2.D_dn_basis.size()*i_new + j_new;


                l0=o0*basis.Length + site0;
                l1=o1*basis.Length + site1;

                /*
                value = value*(pow(-1.0,basis.Nup));
                sign_pow_dn = one_bits_in_bw(2*basis.Length-1,l0,basis.D_dn_basis[j]);
                sign_pow_up = one_bits_in_bw(2*basis.Length-1,l1,basis.D_up_basis[i]);
                if(l0 != 2*basis.Length-1){
                    sign_pow_dn += bit_value(basis.D_dn_basis[j],2*basis.Length-1);
                }
                if(l1 != 2*basis.Length-1){
                    sign_pow_up +=bit_value(basis.D_up_basis[i],2*basis.Length-1);
                }
                */


                value = value*(pow(-1.0,basis.Nup));
                sign_pow_dn = one_bits_in_bw(l0,0,basis.D_dn_basis[j]);
                sign_pow_up = one_bits_in_bw(l1,0,basis.D_up_basis[i]);
                if(l0 != 0){
                    sign_pow_dn += bit_value(basis.D_dn_basis[j],0);
                }
                if(l1 != 0){
                    sign_pow_up +=bit_value(basis.D_up_basis[i],0);
                }


                value = value*(pow(-1.0,(sign_pow_dn+sign_pow_up)));

                Pair_Annihilation.value.push_back(value*one*value_out);
                Pair_Annihilation.rows.push_back(row_no);
                Pair_Annihilation.columns.push_back(column_no);

            }



        }
    }


}


void MODEL_2_orb_Hubb_chain::Variational_state_optimization(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){

    //    double dp=0.005;
    //    double eps_= 0.000001;
    //    double VALUE_CHECK, VALUE_CHECK_OLD ;
    //    double overlap_with_GS;
    //    double Norm_Ansatz;





    //    if(basis.Length==6){
    //        vector<double> alpha_final;
    //        vector<double> alpha;
    //        vector<double> alpha_min, alpha_max;
    //        alpha.resize(8);
    //        alpha_final.resize(8);
    //        alpha_min.resize(8); alpha_max.resize(8);
    //        VALUE_CHECK_OLD=0.0;

    //        /*
    //alpha[0] = -1.1833269165e-02
    //alpha[1] = 4.3848768798e-02
    //alpha[2] = -4.3848768798e-02
    //alpha[3] = 4.3848768798e-02
    //alpha[4] = 4.3848768798e-02
    //alpha[5] = -4.3848768798e-02
    //alpha[6] = 4.3848768798e-02
    //alpha[7] = -1.1833269165e-02
    //overlap_with_GS = 6.1889488366e-01
    //         */

    //        /*
    //alpha[0] = -1.2526611734e-02
    //alpha[1] = 3.5832057778e-02
    //alpha[2] = -4.1125021016e-02
    //alpha[3] = 5.7003910729e-02
    //alpha[4] = 3.5832057778e-02
    //alpha[5] = -4.1125021016e-02
    //alpha[6] = 4.6417984253e-02
    //alpha[7] = -2.3112538210e-02
    //overlap with GS = 6.1985191582e-01
    //         */

    //        alpha_min[0] = dot_product(GS_,BASIS_STATES_ANSATZ[0]);
    //        alpha_min[1] = dot_product(GS_,BASIS_STATES_ANSATZ[2]);
    //        alpha_min[2] = dot_product(GS_,BASIS_STATES_ANSATZ[14]);
    //        alpha_min[3] = dot_product(GS_,BASIS_STATES_ANSATZ[26]);
    //        alpha_min[4] = dot_product(GS_,BASIS_STATES_ANSATZ[28]);
    //        alpha_min[5] = dot_product(GS_,BASIS_STATES_ANSATZ[34]);
    //        alpha_min[6] = dot_product(GS_,BASIS_STATES_ANSATZ[40]);
    //        alpha_min[7] = dot_product(GS_,BASIS_STATES_ANSATZ[46]);

    //        alpha_max[0] = dot_product(GS_,BASIS_STATES_ANSATZ[0]);
    //        alpha_max[1] = dot_product(GS_,BASIS_STATES_ANSATZ[2]);
    //        alpha_max[2] = dot_product(GS_,BASIS_STATES_ANSATZ[14]);
    //        alpha_max[3] = dot_product(GS_,BASIS_STATES_ANSATZ[26]);
    //        alpha_max[4] = dot_product(GS_,BASIS_STATES_ANSATZ[28]);
    //        alpha_max[5] = dot_product(GS_,BASIS_STATES_ANSATZ[34]);
    //        alpha_max[6] = dot_product(GS_,BASIS_STATES_ANSATZ[40]);
    //        alpha_max[7] = dot_product(GS_,BASIS_STATES_ANSATZ[46]);





    //        //cout<<"here 1"<<endl;

    //        alpha[0]=alpha_min[0];
    //        while(alpha[0]<=alpha_max[0]){

    //            alpha[1]=alpha_min[1];
    //            //alpha[1]=-1.0*alpha[0];
    //            while(alpha[1]<=alpha_max[1]){

    //                alpha[2]=alpha_min[2];
    //                //alpha[2]=alpha[0];
    //                while(alpha[2]<=alpha_max[2]){

    //                    alpha[3]=alpha_min[3];
    //                    //alpha[3]=-1.0*alpha[0];
    //                    while(alpha[3]<=alpha_max[3]){

    //                        alpha[4]=alpha_min[4];
    //                        //alpha[4]=-1.0*alpha[0];
    //                        while(alpha[4]<=alpha_max[4]){

    //                            alpha[5]=alpha_min[5];
    //                            //alpha[5]=alpha[0];
    //                            while(alpha[5]<=alpha_max[5]){

    //                                alpha[6]=alpha_min[6];
    //                                //alpha[6]=-1.0*alpha[0];
    //                                while(alpha[6]<=alpha_max[6]){

    //                                    alpha[7]=alpha_min[7];
    //                                    //alpha[7]=alpha[0];
    //                                    while(alpha[7]<=alpha_max[7]){

    //                                        //cout<<"here 2"<<endl;

    //                                        State_.clear();
    //                                        State_.resize(BASIS_STATES_ANSATZ[0].size());


    //                                        for(int m=0;m<=1;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[0]*(BASIS_STATES_ANSATZ[m][i]);

    //                                            }
    //                                        }
    //                                        //cout<<"1"<<endl;
    //                                        for(int m=2;m<=13;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[1]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"2"<<endl;
    //                                        for(int m=14;m<=25;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[2]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"3"<<endl;
    //                                        for(int m=26;m<=27;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[3]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"4"<<endl;
    //                                        for(int m=28;m<=33;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[4]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"5"<<endl;
    //                                        for(int m=34;m<=39;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[5]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"6"<<endl;
    //                                        for(int m=40;m<=45;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[6]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"7"<<endl;
    //                                        for(int m=52;m<=57;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[6]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"8"<<endl;
    //                                        for(int m=46;m<=51;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[7]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"9"<<endl;
    //                                        for(int m=58;m<=63;m++){
    //                                            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                                                State_[i] += alpha[7]*(BASIS_STATES_ANSATZ[m][i]);
    //                                            }
    //                                        }
    //                                        //cout<<"10"<<endl;



    //                                        Norm_Ansatz=dot_product(State_,State_);
    //                                        if(Norm_Ansatz <= eps_){
    //                                            VALUE_CHECK=0.0;
    //                                        }
    //                                        else{
    //                                            overlap_with_GS=dot_product(GS_,State_);
    //                                            overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //                                            VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //                                        }


    //                                        if( VALUE_CHECK > VALUE_CHECK_OLD){
    //                                            VALUE_CHECK_OLD=VALUE_CHECK;
    //                                            alpha_final=alpha;
    //                                        }

    //                                        for(int n=0;n<8;n++){
    //                                            cout<<"alpha["<<n<<"] = "<<alpha[n]<<"    ";
    //                                        }
    //                                        cout<< "is done"<<endl;

    //                                        alpha[7]+=dp*1;

    //                                    }
    //                                    alpha[6]+=dp*1;
    //                                }
    //                                alpha[5]+=dp*1;
    //                            }
    //                            alpha[4]+=dp*1;
    //                        }
    //                        alpha[3]+=dp*1;
    //                    }
    //                    alpha[2]+=dp*1;
    //                }

    //                alpha[1]+=dp*1;
    //            }

    //            alpha[0]+=dp;
    //        }


    //        State_.clear();
    //        State_.resize(BASIS_STATES_ANSATZ[0].size());
    //        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){

    //            for(int m=0;m<=1;m++){
    //                State_[i] += alpha_final[0]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=2;m<=13;m++){
    //                State_[i] += alpha_final[1]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=14;m<=25;m++){
    //                State_[i] += alpha_final[2]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=26;m<=27;m++){
    //                State_[i] += alpha_final[3]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=28;m<=33;m++){
    //                State_[i] += alpha_final[4]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=34;m<=39;m++){
    //                State_[i] += alpha_final[5]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=40;m<=45;m++){
    //                State_[i] += alpha_final[6]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=52;m<=57;m++){
    //                State_[i] += alpha_final[6]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=46;m<=51;m++){
    //                State_[i] += alpha_final[7]*(BASIS_STATES_ANSATZ[m][i]);
    //            }
    //            for(int m=58;m<=63;m++){
    //                State_[i] += alpha_final[7]*(BASIS_STATES_ANSATZ[m][i]);
    //            }

    //        }
    //        Norm_Ansatz=dot_product(State_,State_);
    //        overlap_with_GS=dot_product(GS_,State_);
    //        overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //        VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //        cout<<"XXXXXXXX-----VARIATIONAL STATE OPTIMIZATION-------XXXXXXXXXX"<<endl;
    //        for(int n=0;n<8;n++){
    //            cout<<"alpha["<<n<<"] = "<<alpha_final[n]/sqrt(Norm_Ansatz)<<endl;
    //        }
    //        cout<<"(overlap_with_GS)^2 = "<<(overlap_with_GS/Norm_Ansatz)<<endl;
    //        cout<<"Norm_Ansatz = "<<Norm_Ansatz<<endl;
    //        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;



    //    }


    //    if(basis.Length==4){
    //        double alpha_final, beta_final, gamma_final, delta_final;
    //        double alpha, beta, gamma, delta;

    //        VALUE_CHECK_OLD=0.0;
    //        double alpha_min, alpha_max, beta_min, beta_max;
    //        double gamma_min, gamma_max, delta_min, delta_max;

    //        /*
    //    alpha_min=-1; alpha_max=1;
    //    beta_min=-1; beta_max=1;
    //    gamma_min=-1; gamma_max=1;
    //    delta_min=-1; delta_max=1;
    //    */

    //        /*//9.0501610389e-01
    //        alpha_min=-9.6964354728e-02; alpha_max=-9.6964354728e-02;
    //        beta_min=-1.2120544341e-01; beta_max=-1.2120544341e-01;
    //        gamma_min=1.0908489907e-01; gamma_max=1.0908489907e-01;
    //        delta_min=-1.0908489907e-01; delta_max=-1.0908489907e-01;
    //        */


    //        /*
    //    alpha = -1.0022660134e-01
    //    beta = -1.1796553740e-01
    //    gamma = 1.0909606937e-01
    //    delta = -1.0909606937e-01
    //    overlap_with_GS = 9.0503998632e-01
    //    */

    //        /*
    //    alpha = -9.9917982671e-02
    //    beta = -1.1827233816e-01
    //    gamma = 1.0909516042e-01
    //    delta = -1.0909516042e-01
    //    overlap_with_GS = 9.0504018577e-01
    //*/

    //        /*
    //        alpha = -9.9961035979e-02
    //        beta = -1.1831778857e-01
    //        gamma = 1.0907940444e-01
    //        delta = -1.0908240483e-01
    //        overlap_with_GS = 9.0504021161e-01
    //        */

    //        /*
    //    alpha_min=-1; alpha_max=1;
    //    beta_min=-1; beta_max=1;
    //    gamma_min=-1; gamma_max=1;
    //    delta_min=-1; delta_max=1;
    //*/


    //        alpha_min=dot_product(GS_,BASIS_STATES_ANSATZ[0]); alpha_max=dot_product(GS_,BASIS_STATES_ANSATZ[0]);
    //        beta_min=dot_product(GS_,BASIS_STATES_ANSATZ[2]); beta_max=dot_product(GS_,BASIS_STATES_ANSATZ[2]);
    //        gamma_min=dot_product(GS_,BASIS_STATES_ANSATZ[4]); gamma_max=dot_product(GS_,BASIS_STATES_ANSATZ[4]);
    //        delta_min=dot_product(GS_,BASIS_STATES_ANSATZ[12]); delta_max=dot_product(GS_,BASIS_STATES_ANSATZ[12]);


    //        /*
    //        alpha_min=-9.9961035979e-02; alpha_max=-9.9961035979e-02;
    //        beta_min=-1.1831778857e-01; beta_max=-1.1831778857e-01;
    //        gamma_min=1.0907940444e-01; gamma_max=1.0907940444e-01;
    //        delta_min=-1.0908240483e-01; delta_max=-1.0908240483e-01;
    //*/


    //        /*
    //        alpha_min=-1.0; alpha_max=1.0;
    //        beta_min=-1.0; beta_max=1.0;
    //        gamma_min=-1.0; gamma_max=1.0;
    //        delta_min=-1.0; delta_max=1.0;
    //        */
    //        /*
    //        alpha_min=-6.6295939532e-02; alpha_max=-6.6295939532e-02;
    //        beta_min= -1.5745285639e-01; beta_max= -1.5745285639e-01;
    //        gamma_min=1.0773090174e-01; gamma_max=1.0773090174e-01;
    //        delta_min=-1.0773090174e-01; delta_max=-1.0773090174e-01;

    //*/

    //        for(alpha=alpha_min-0.0;alpha<=alpha_max+0.0;alpha=alpha+dp){
    //            for(beta=beta_min-0.0;beta<=beta_max+0.0;beta=beta+dp){
    //                for(gamma=gamma_min-0.0;gamma<=gamma_max+0.0;gamma=gamma+dp){
    //                    for(delta=delta_min-0.0;delta<=delta_max+0.0;delta=delta+dp){


    //                        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                            State_[i] = alpha*(BASIS_STATES_ANSATZ[0][i] + BASIS_STATES_ANSATZ[1][i]) +
    //                                    beta*(BASIS_STATES_ANSATZ[2][i] + BASIS_STATES_ANSATZ[3][i]) +
    //                                    gamma*(BASIS_STATES_ANSATZ[4][i] + BASIS_STATES_ANSATZ[5][i]
    //                                    + BASIS_STATES_ANSATZ[6][i] + BASIS_STATES_ANSATZ[7][i]
    //                                    + BASIS_STATES_ANSATZ[9][i] + BASIS_STATES_ANSATZ[8][i]
    //                                    + BASIS_STATES_ANSATZ[10][i] + BASIS_STATES_ANSATZ[11][i]) +
    //                                    delta*(BASIS_STATES_ANSATZ[12][i] + BASIS_STATES_ANSATZ[13][i]
    //                                    + BASIS_STATES_ANSATZ[14][i] + BASIS_STATES_ANSATZ[15][i]);
    //                        }
    //                        Norm_Ansatz=dot_product(State_,State_);
    //                        if(Norm_Ansatz <= eps_){
    //                            VALUE_CHECK=0.0;
    //                        }
    //                        else{
    //                            overlap_with_GS=dot_product(GS_,State_);
    //                            overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //                            VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //                        }


    //                        if( VALUE_CHECK > VALUE_CHECK_OLD){
    //                            VALUE_CHECK_OLD=VALUE_CHECK;
    //                            alpha_final=alpha;beta_final=beta;
    //                            gamma_final=gamma;delta_final=delta;
    //                        }

    //                    }
    //                }

    //            }
    //            cout<<"alpha = "<<alpha<<" check done"<<endl;

    //        }


    //        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //            State_[i] = alpha_final*(BASIS_STATES_ANSATZ[0][i] + BASIS_STATES_ANSATZ[1][i]) +
    //                    beta_final*(BASIS_STATES_ANSATZ[2][i] + BASIS_STATES_ANSATZ[3][i]) +
    //                    gamma_final*(BASIS_STATES_ANSATZ[4][i] + BASIS_STATES_ANSATZ[5][i]
    //                    + BASIS_STATES_ANSATZ[6][i] + BASIS_STATES_ANSATZ[7][i]
    //                    + BASIS_STATES_ANSATZ[9][i] + BASIS_STATES_ANSATZ[8][i]
    //                    + BASIS_STATES_ANSATZ[10][i] + BASIS_STATES_ANSATZ[11][i]) +
    //                    delta_final*(BASIS_STATES_ANSATZ[12][i] + BASIS_STATES_ANSATZ[13][i]
    //                    + BASIS_STATES_ANSATZ[14][i] + BASIS_STATES_ANSATZ[15][i]);
    //        }
    //        Norm_Ansatz=dot_product(State_,State_);
    //        overlap_with_GS=dot_product(GS_,State_);
    //        overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //        VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //        cout<<"XXXXXXXX-----VARIATIONAL STATE OPTIMIZATION-------XXXXXXXXXX"<<endl;
    //        cout<<"alpha = "<<alpha_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"beta = "<<beta_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"gamma = "<<gamma_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"delta = "<<delta_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"(overlap_with_GS)^2 = "<<overlap_with_GS/Norm_Ansatz<<endl;
    //        cout<<"Norm_Ansatz = "<<Norm_Ansatz<<endl;
    //        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;



    //    }



    //    if(basis.Length==2){
    //        double alpha_final, beta_final;
    //        double alpha, beta;
    //        double alpha_min, alpha_max, beta_min, beta_max;


    //        alpha_min=0.5773502691896; alpha_max=0.5773502691896;
    //        beta_min=-0.5773502691896; beta_max=-0.5773502691896;

    //        /*
    //        alpha_min=1.0; alpha_max=1.0;
    //        beta_min=-1.0; beta_max=-1.0;
    //        */

    //        VALUE_CHECK_OLD=0.0;
    //        for(alpha=alpha_min;alpha<=alpha_max;alpha=alpha+dp){
    //            for(beta=beta_min;beta<=beta_max;beta=beta+dp){

    //                for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //                    State_[i] = alpha*(BASIS_STATES_ANSATZ[0][i]) +
    //                            beta*(BASIS_STATES_ANSATZ[1][i]) ;
    //                }
    //                Norm_Ansatz=dot_product(State_,State_);
    //                if(Norm_Ansatz <= eps_){
    //                    VALUE_CHECK=0.0;
    //                }
    //                else{
    //                    overlap_with_GS=dot_product(GS_,State_);
    //                    overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //                    VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //                }


    //                if( VALUE_CHECK > VALUE_CHECK_OLD){
    //                    VALUE_CHECK_OLD=VALUE_CHECK;
    //                    alpha_final=alpha;beta_final=beta;
    //                }



    //            }

    //        }


    //        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
    //            State_[i] = alpha_final*(BASIS_STATES_ANSATZ[0][i]) +
    //                    beta_final*(BASIS_STATES_ANSATZ[1][i]) ;
    //        }
    //        Norm_Ansatz=dot_product(State_,State_);
    //        overlap_with_GS=dot_product(GS_,State_);
    //        overlap_with_GS = overlap_with_GS*overlap_with_GS;
    //        VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

    //        cout<<"XXXXXXXX-----VARIATIONAL STATE OPTIMIZATION-------XXXXXXXXXX"<<endl;
    //        cout<<"alpha = "<<alpha_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"beta = "<<beta_final/sqrt(Norm_Ansatz)<<endl;
    //        cout<<"overlap_with_GS = "<<VALUE_CHECK<<endl;
    //        //cout<<"Norm_Ansatz = "<<Norm_Ansatz<<endl;
    //        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;



    //    }

}



//void MODEL_2_orb_Hubb_chain::Create_Anzatz_basis_2Holes_new(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){


//    Variational_state_pair_Geometry.resize(basis.Length-1);
//    Variational_state_pair_spin_symmetry.resize(basis.Length-1);
//    Variational_state_pair_orbital_symmetry.resize(basis.Length-1);
//    Variational_state_pair_sites.resize(basis.Length-1);

//    Mat_1_doub State_temp1, State_temp2;
//    int state_count;
//    double temp1, temp2;



//    int No_of_classes,No_of_states;
//    Mat_1_int States_in_class;
//    Ansatz_2holes_Basis_Overlap_with_GS.clear();


//    if(basis.Length==6){
//    No_of_classes = 36;
//    No_of_states = No_of_classes*basis.Length*2;

//    States_in_class.resize(No_of_classes);
//    variational_2holes_classes_contruction_strings.resize(No_of_classes);

//    States_in_class[0]=12;
//    variational_2holes_classes_contruction_strings[0]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
//    States_in_class[1]=12;
//    variational_2holes_classes_contruction_strings[1]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
//    States_in_class[2]=12;
//    variational_2holes_classes_contruction_strings[2]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
//    States_in_class[3]=12;
//    variational_2holes_classes_contruction_strings[3]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
//    States_in_class[4]=12;
//    variational_2holes_classes_contruction_strings[4]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
//    States_in_class[5]=12;
//    variational_2holes_classes_contruction_strings[5]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
//    States_in_class[6]=12;
//    variational_2holes_classes_contruction_strings[6]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
//    States_in_class[7]=12;
//    variational_2holes_classes_contruction_strings[7]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
//    States_in_class[8]=12;
//    variational_2holes_classes_contruction_strings[8]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
//    States_in_class[9]=12;
//    variational_2holes_classes_contruction_strings[9]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
//    States_in_class[10]=12;
//    variational_2holes_classes_contruction_strings[10]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
//    States_in_class[11]=12;
//    variational_2holes_classes_contruction_strings[11]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
//    States_in_class[12]=12;
//    variational_2holes_classes_contruction_strings[12]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
//    States_in_class[13]=12;
//    variational_2holes_classes_contruction_strings[13]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
//    States_in_class[14]=12;
//    variational_2holes_classes_contruction_strings[14]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
//    States_in_class[15]=12;
//    variational_2holes_classes_contruction_strings[15]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
//    States_in_class[16]=12;
//    variational_2holes_classes_contruction_strings[16]="AC S NOS2 0 1 AC S NOS1 1 3 AC S NOS2 2 4 D S NOS1 4 5 AC S NOS1 5 0";
//    States_in_class[17]=12;
//    variational_2holes_classes_contruction_strings[17]="AC S NOS2 0 1 AC S NOS1 1 3 AC S NOS2 2 4 AC S NOS1 4 5 D S NOS2 5 0";
//    States_in_class[18]=12;
//    variational_2holes_classes_contruction_strings[18]="D S NOS1 0 1 AC S NOS1 1 3 AC S NOS2 2 4 D S NOS1 4 5 D S NOS1 5 0";
//    States_in_class[19]=12;
//    variational_2holes_classes_contruction_strings[19]="D S NOS1 0 1 AC S NOS1 1 3 AC S NOS2 2 4 AC S NOS1 4 5 AC S NOS2 5 0";
//    States_in_class[20]=12;
//    variational_2holes_classes_contruction_strings[20]="AC S NOS2 0 1 D S NOS1 1 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
//    States_in_class[21]=12;
//    variational_2holes_classes_contruction_strings[21]="AC S NOS2 0 1 D S NOS1 1 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
//    States_in_class[22]=12;
//    variational_2holes_classes_contruction_strings[22]="AC S NOS2 0 1 D S NOS1 1 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
//    States_in_class[23]=12;
//    variational_2holes_classes_contruction_strings[23]="AC S NOS2 0 1 D S NOS1 1 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
//    States_in_class[24]=12;
//    variational_2holes_classes_contruction_strings[24]="D S NOS1 0 1 D S NOS1 1 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
//    States_in_class[25]=12;
//    variational_2holes_classes_contruction_strings[25]="D S NOS1 0 1 D S NOS1 1 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
//    States_in_class[26]=12;
//    variational_2holes_classes_contruction_strings[26]="D S NOS1 0 1 D S NOS1 1 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
//    States_in_class[27]=12;
//    variational_2holes_classes_contruction_strings[27]="D S NOS1 0 1 D S NOS1 1 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
//    States_in_class[28]=12;
//    variational_2holes_classes_contruction_strings[28]="D S NOS2 0 1 AC S NOS2 1 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
//    States_in_class[29]=12;
//    variational_2holes_classes_contruction_strings[29]="D S NOS2 0 1 AC S NOS2 1 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
//    States_in_class[30]=12;
//    variational_2holes_classes_contruction_strings[30]="D S NOS2 0 1 AC S NOS2 1 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
//    States_in_class[31]=12;
//    variational_2holes_classes_contruction_strings[31]="D S NOS2 0 1 AC S NOS2 1 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
//    States_in_class[32]=12;
//    variational_2holes_classes_contruction_strings[32]="AC S NOS1 0 1 AC S NOS2 1 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
//    States_in_class[33]=12;
//    variational_2holes_classes_contruction_strings[33]="AC S NOS1 0 1 AC S NOS2 1 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
//    States_in_class[34]=12;
//    variational_2holes_classes_contruction_strings[34]="AC S NOS1 0 1 AC S NOS2 1 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
//    States_in_class[35]=12;
//    variational_2holes_classes_contruction_strings[35]="AC S NOS1 0 1 AC S NOS2 1 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";

//        }



//    if(basis.Length==4){
//    No_of_classes = 7;
//    No_of_states = No_of_classes*basis.Length*2;

//    States_in_class.resize(No_of_classes);
//    variational_2holes_classes_contruction_strings.resize(No_of_classes);

//    States_in_class[0]=16;
//    variational_2holes_classes_contruction_strings[0]="D S NOS1 0 2 D S NOS1 2 3 D S NOS1 3 0";
//    States_in_class[1]=16;
//    variational_2holes_classes_contruction_strings[1]="D S NOS1 0 2 AC S NOS1 2 3 AC S NOS2 3 0";
//    States_in_class[2]=16;
//    variational_2holes_classes_contruction_strings[2]="AC S NOS2 0 2 D S NOS1 2 3 AC S NOS1 3 0";
//    States_in_class[3]=16;
//    variational_2holes_classes_contruction_strings[3]="D S NOS1 0 1 D S NOS1 2 3 D S NOS2 3 0";
//    States_in_class[4]=16;
//    variational_2holes_classes_contruction_strings[4]="AC S NOS2 0 1 D S NOS1 2 3 AC S NOS1 3 0";
//    States_in_class[5]=16;
//    variational_2holes_classes_contruction_strings[5]="AC S NOS2 0 1 AC S NOS1 2 3 D S NOS2 3 0";
//    States_in_class[6]=16;
//    variational_2holes_classes_contruction_strings[6]="AC S NOS1 0 2 AC S NOS2 1 3 D S NOS1 3 0";
//        }

//    int temp_int;

//    BASIS_STATES_ANSATZ_2HOLES.clear();


//    string variational_state_contruction_;


//    string file_out_overlap_GS = "RVB_2holes_GSNm2_overlap.txt";

//    ofstream outfile_out_overlap_GS(file_out_overlap_GS.c_str());
//    outfile_out_overlap_GS<<"# state_no   overlap_with_GS(N-2)"<<endl;



//    for(int i=0;i<No_of_classes;i++){
//        ostringstream ss_int;
//        ss_int << i;

//        variational_state_contruction_=variational_2holes_classes_contruction_strings[i];
//        //------------UPDATING Variational_state_pair******--------
//        stringstream variational_state_contruction_stream(variational_state_contruction_);
//        for(int n=0;n<basis.Length-1;n++){
//            variational_state_contruction_stream >> Variational_state_pair_Geometry[n];
//            assert(Variational_state_pair_Geometry[n]=="D" ||
//                   Variational_state_pair_Geometry[n]=="AC");

//            variational_state_contruction_stream >> Variational_state_pair_spin_symmetry[n];
//            assert(Variational_state_pair_spin_symmetry[n]=="S" ||
//                   Variational_state_pair_spin_symmetry[n]=="T");

//            variational_state_contruction_stream >> Variational_state_pair_orbital_symmetry[n];
//            assert(Variational_state_pair_orbital_symmetry[n] =="OA" ||
//                   Variational_state_pair_orbital_symmetry[n] =="OS" ||
//                   Variational_state_pair_orbital_symmetry[n] =="NOS1" ||
//                   Variational_state_pair_orbital_symmetry[n] =="NOS2");

//            variational_state_contruction_stream >> temp_int;
//            Variational_state_pair_sites[n].first =temp_int;
//            assert(temp_int>=0 && temp_int<basis.Length);

//            variational_state_contruction_stream >> temp_int;
//            Variational_state_pair_sites[n].second =temp_int;
//            assert(temp_int>=0 && temp_int<basis.Length);
//        }


//        Get_Variational_State(basis,basis.Length-1);
//        BASIS_STATES_ANSATZ_2HOLES.push_back(State_);
//        temp1=dot_product(GS_,State_);
//        Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
//        state_count=0;

//        for(int Trns_n=1;Trns_n<=basis.Length;Trns_n++){
//            for(int Orb_Ex_n=1;Orb_Ex_n=1<=2;Orb_Ex_n++){
//                for(int Ref_n=1;Ref_n=1<=2;Ref_n++){
//HERE!!
//                }
//            }
//        }



//        cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
//              ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
//        cout<<temp1<<endl;
//        outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
//        ostringstream ssn_int;
//        ssn_int << state_count;
//        string file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
//                "_state_" + ssn_int.str() +".txt";

//        ofstream outfile_Variational_states(file_out_Variational_states.c_str());

//        outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
//                                 <<variational_state_contruction_<<endl;

//        for(int j=0;j<State_.size();j++){
//            if(State_[j]!=0.0){
//                outfile_Variational_states<<j<<"\t"<<State_[j]<<endl;
//            }
//        }
//        outfile_Variational_states.close();


//        //Orbital_Exchange
//        State_temp1 = Act_Orbital_Exchange(basis, State_);
//        OEn++;
//        BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp1);
//        temp1=dot_product(GS_,State_temp1);
//        Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
//        state_count++;
//        cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
//              ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
//        cout<<temp1<<endl;
//        outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
//        ssn_int.str("");
//        ssn_int.clear();
//        ssn_int << state_count;
//        file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
//                "_state_" + ssn_int.str() +".txt";
//        outfile_Variational_states.open(file_out_Variational_states.c_str());
//        outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
//                                 <<variational_state_contruction_<<endl;

//        for(int j=0;j<State_temp1.size();j++){
//            if(State_temp1[j]!=0.0){
//                outfile_Variational_states<<j<<"\t"<<State_temp1[j]<<endl;
//            }
//        }
//        outfile_Variational_states.close();

//        //TRANSLATIONS
//        for(int T_no=1;T_no<=basis.Length-1;T_no++){
//            State_temp2 = Act_Translational_operator(basis, State_temp1);
//            Tn++;
//            BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp2);
//            temp1=dot_product(GS_,State_temp2);
//            Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
//            state_count++;
//            cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
//                  ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
//            cout<<temp1<<endl;
//            outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
//            ssn_int.str("");
//            ssn_int.clear();
//            ssn_int << state_count;
//            file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
//                    "_state_" + ssn_int.str() +".txt";
//            outfile_Variational_states.open(file_out_Variational_states.c_str());
//            outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
//                                     <<variational_state_contruction_<<endl;

//            for(int j=0;j<State_temp2.size();j++){
//                if(State_temp2[j]!=0.0){
//                    outfile_Variational_states<<j<<"\t"<<State_temp2[j]<<endl;
//                }
//            }
//            outfile_Variational_states.close();


//            //Orbital_Exchange
//            State_temp1 = Act_Orbital_Exchange(basis, State_temp2);
//            OEn++;
//            BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp1);
//            temp1=dot_product(GS_,State_temp1);
//            Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
//            state_count++;
//            cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
//                  ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
//            cout<<temp1<<endl;
//            outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
//            ssn_int.str("");
//            ssn_int.clear();
//            ssn_int << state_count;
//            file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
//                    "_state_" + ssn_int.str() +".txt";
//            outfile_Variational_states.open(file_out_Variational_states.c_str());
//            outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
//                                     <<variational_state_contruction_<<endl;

//            for(int j=0;j<State_temp1.size();j++){
//                if(State_temp1[j]!=0.0){
//                    outfile_Variational_states<<j<<"\t"<<State_temp1[j]<<endl;
//                }
//            }
//            outfile_Variational_states.close();

//        }

//        assert(state_count+1 == States_in_class[i]);





//        //-----------------------------




//    }



//    string file_out_overlaps = "overlap_bw_RVB_2holes_states.txt";

//    ofstream outfile_out_overlaps(file_out_overlaps.c_str());
//    outfile_out_overlaps<<"# state_no(i)   state_no(j)  <i|j> "<<endl;

//    OVERLAP_MATRIX_2holes_Basis.resize(BASIS_STATES_ANSATZ_2HOLES.size());
//    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
//        OVERLAP_MATRIX_2holes_Basis[i].resize(BASIS_STATES_ANSATZ_2HOLES.size());
//    }

//    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
//        for(int j=0;j<BASIS_STATES_ANSATZ_2HOLES.size();j++){

//            OVERLAP_MATRIX_2holes_Basis[i][j]=dot_product(BASIS_STATES_ANSATZ_2HOLES[i],BASIS_STATES_ANSATZ_2HOLES[j]);
//            outfile_out_overlaps<<i<<"  "<<j<< "   "<<OVERLAP_MATRIX_2holes_Basis[i][j]<<endl;
//            //cout<<"i="<<i<<", j="<<j<<", overlap="<<OVERLAP_MATRIX_2holes_Basis[i][j]<<endl;
//        }
//    }

//    /*
//    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
//        for(int j=0;j<BASIS_STATES_ANSATZ_2HOLES.size();j++){
//            cout<<dot_product(BASIS_STATES_ANSATZ_2HOLES[i],BASIS_STATES_ANSATZ_2HOLES[j])<<"  ";
//        }
//        cout<<endl;
//    }
//    */


//}


void MODEL_2_orb_Hubb_chain::Create_Anzatz_basis_2Holes(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){


    //    Variational_state_pair_Geometry.resize(basis.Length-1);
    //    Variational_state_pair_spin_symmetry.resize(basis.Length-1);
    //    Variational_state_pair_orbital_symmetry.resize(basis.Length-1);
    //    Variational_state_pair_sites.resize(basis.Length-1);

    //    Mat_1_doub State_temp1, State_temp2;
    //    int state_count;
    //    double temp1, temp2;



    //    int No_of_classes,No_of_states;
    //    Mat_1_int States_in_class;
    //    Ansatz_2holes_Basis_Overlap_with_GS.clear();


    //    if(basis.Length==6){
    //        No_of_classes = 36;
    //        No_of_states = No_of_classes*basis.Length*2;

    //        States_in_class.resize(No_of_classes);
    //        variational_2holes_classes_contruction_strings.resize(No_of_classes);

    //        States_in_class[0]=12;
    //        variational_2holes_classes_contruction_strings[0]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //        States_in_class[1]=12;
    //        variational_2holes_classes_contruction_strings[1]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //        States_in_class[2]=12;
    //        variational_2holes_classes_contruction_strings[2]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //        States_in_class[3]=12;
    //        variational_2holes_classes_contruction_strings[3]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //        States_in_class[4]=12;
    //        variational_2holes_classes_contruction_strings[4]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //        States_in_class[5]=12;
    //        variational_2holes_classes_contruction_strings[5]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //        States_in_class[6]=12;
    //        variational_2holes_classes_contruction_strings[6]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //        States_in_class[7]=12;
    //        variational_2holes_classes_contruction_strings[7]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //        States_in_class[8]=12;
    //        variational_2holes_classes_contruction_strings[8]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //        States_in_class[9]=12;
    //        variational_2holes_classes_contruction_strings[9]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //        States_in_class[10]=12;
    //        variational_2holes_classes_contruction_strings[10]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //        States_in_class[11]=12;
    //        variational_2holes_classes_contruction_strings[11]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //        States_in_class[12]=12;
    //        variational_2holes_classes_contruction_strings[12]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //        States_in_class[13]=12;
    //        variational_2holes_classes_contruction_strings[13]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //        States_in_class[14]=12;
    //        variational_2holes_classes_contruction_strings[14]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //        States_in_class[15]=12;
    //        variational_2holes_classes_contruction_strings[15]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //        States_in_class[16]=12;
    //        variational_2holes_classes_contruction_strings[16]="AC S NOS2 0 1 AC S NOS1 1 3 AC S NOS2 2 4 D S NOS1 4 5 AC S NOS1 5 0";
    //        States_in_class[17]=12;
    //        variational_2holes_classes_contruction_strings[17]="AC S NOS2 0 1 AC S NOS1 1 3 AC S NOS2 2 4 AC S NOS1 4 5 D S NOS2 5 0";
    //        States_in_class[18]=12;
    //        variational_2holes_classes_contruction_strings[18]="D S NOS1 0 1 AC S NOS1 1 3 AC S NOS2 2 4 D S NOS1 4 5 D S NOS1 5 0";
    //        States_in_class[19]=12;
    //        variational_2holes_classes_contruction_strings[19]="D S NOS1 0 1 AC S NOS1 1 3 AC S NOS2 2 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //        States_in_class[20]=12;
    //        variational_2holes_classes_contruction_strings[20]="AC S NOS2 0 1 D S NOS1 1 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //        States_in_class[21]=12;
    //        variational_2holes_classes_contruction_strings[21]="AC S NOS2 0 1 D S NOS1 1 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //        States_in_class[22]=12;
    //        variational_2holes_classes_contruction_strings[22]="AC S NOS2 0 1 D S NOS1 1 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //        States_in_class[23]=12;
    //        variational_2holes_classes_contruction_strings[23]="AC S NOS2 0 1 D S NOS1 1 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //        States_in_class[24]=12;
    //        variational_2holes_classes_contruction_strings[24]="D S NOS1 0 1 D S NOS1 1 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //        States_in_class[25]=12;
    //        variational_2holes_classes_contruction_strings[25]="D S NOS1 0 1 D S NOS1 1 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //        States_in_class[26]=12;
    //        variational_2holes_classes_contruction_strings[26]="D S NOS1 0 1 D S NOS1 1 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //        States_in_class[27]=12;
    //        variational_2holes_classes_contruction_strings[27]="D S NOS1 0 1 D S NOS1 1 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //        States_in_class[28]=12;
    //        variational_2holes_classes_contruction_strings[28]="D S NOS2 0 1 AC S NOS2 1 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //        States_in_class[29]=12;
    //        variational_2holes_classes_contruction_strings[29]="D S NOS2 0 1 AC S NOS2 1 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //        States_in_class[30]=12;
    //        variational_2holes_classes_contruction_strings[30]="D S NOS2 0 1 AC S NOS2 1 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //        States_in_class[31]=12;
    //        variational_2holes_classes_contruction_strings[31]="D S NOS2 0 1 AC S NOS2 1 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //        States_in_class[32]=12;
    //        variational_2holes_classes_contruction_strings[32]="AC S NOS1 0 1 AC S NOS2 1 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //        States_in_class[33]=12;
    //        variational_2holes_classes_contruction_strings[33]="AC S NOS1 0 1 AC S NOS2 1 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //        States_in_class[34]=12;
    //        variational_2holes_classes_contruction_strings[34]="AC S NOS1 0 1 AC S NOS2 1 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //        States_in_class[35]=12;
    //        variational_2holes_classes_contruction_strings[35]="AC S NOS1 0 1 AC S NOS2 1 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";

    //    }



    //    if(basis.Length==4){
    //        No_of_classes = 7;
    //        No_of_states = No_of_classes*basis.Length*2;

    //        States_in_class.resize(No_of_classes);
    //        variational_2holes_classes_contruction_strings.resize(No_of_classes);

    //        States_in_class[0]=16;
    //        variational_2holes_classes_contruction_strings[0]="D S NOS1 0 2 D S NOS1 2 3 D S NOS1 3 0";
    //        States_in_class[1]=16;
    //        variational_2holes_classes_contruction_strings[1]="D S NOS1 0 2 AC S NOS1 2 3 AC S NOS2 3 0";
    //        States_in_class[2]=16;
    //        variational_2holes_classes_contruction_strings[2]="AC S NOS2 0 2 D S NOS1 2 3 AC S NOS1 3 0";
    //        States_in_class[3]=16;
    //        variational_2holes_classes_contruction_strings[3]="D S NOS1 0 1 D S NOS1 2 3 D S NOS2 3 0";
    //        States_in_class[4]=16;
    //        variational_2holes_classes_contruction_strings[4]="AC S NOS2 0 1 D S NOS1 2 3 AC S NOS1 3 0";
    //        States_in_class[5]=16;
    //        variational_2holes_classes_contruction_strings[5]="AC S NOS2 0 1 AC S NOS1 2 3 D S NOS2 3 0";
    //        States_in_class[6]=16;
    //        variational_2holes_classes_contruction_strings[6]="AC S NOS1 0 2 AC S NOS2 1 3 D S NOS1 3 0";
    //    }

    //    int temp_int;

    //    BASIS_STATES_ANSATZ_2HOLES.clear();


    //    string variational_state_contruction_;


    //    string file_out_overlap_GS = "RVB_2holes_GSNm2_overlap.txt";

    //    ofstream outfile_out_overlap_GS(file_out_overlap_GS.c_str());
    //    outfile_out_overlap_GS<<"# state_no   overlap_with_GS(N-2)"<<endl;



    //    for(int i=0;i<No_of_classes;i++){
    //        ostringstream ss_int;
    //        ss_int << i;



    //        variational_state_contruction_=variational_2holes_classes_contruction_strings[i];


    //        //------------UPDATING Variational_state_pair******--------
    //        stringstream variational_state_contruction_stream(variational_state_contruction_);
    //        for(int n=0;n<basis.Length-1;n++){
    //            variational_state_contruction_stream >> Variational_state_pair_Geometry[n];
    //            assert(Variational_state_pair_Geometry[n]=="D" ||
    //                   Variational_state_pair_Geometry[n]=="AC");

    //            variational_state_contruction_stream >> Variational_state_pair_spin_symmetry[n];
    //            assert(Variational_state_pair_spin_symmetry[n]=="S" ||
    //                   Variational_state_pair_spin_symmetry[n]=="T");

    //            variational_state_contruction_stream >> Variational_state_pair_orbital_symmetry[n];
    //            assert(Variational_state_pair_orbital_symmetry[n] =="OA" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="OS" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="NOS1" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="NOS2");

    //            variational_state_contruction_stream >> temp_int;
    //            Variational_state_pair_sites[n].first =temp_int;
    //            assert(temp_int>=0 && temp_int<basis.Length);


    //            variational_state_contruction_stream >> temp_int;
    //            Variational_state_pair_sites[n].second =temp_int;
    //            assert(temp_int>=0 && temp_int<basis.Length);
    //        }


    //        int Tn, OEn, Refn;
    //        Tn=0;OEn=0;Refn=0;

    //        Get_Variational_State(basis,basis.Length-1);
    //        BASIS_STATES_ANSATZ_2HOLES.push_back(State_);
    //        temp1=dot_product(GS_,State_);
    //        Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
    //        state_count=0;
    //        cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
    //              ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
    //        cout<<temp1<<endl;
    //        outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
    //        ostringstream ssn_int;
    //        ssn_int << state_count;
    //        string file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
    //                "_state_" + ssn_int.str() +".txt";

    //        ofstream outfile_Variational_states(file_out_Variational_states.c_str());

    //        outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
    //                                 <<variational_state_contruction_<<endl;

    //        for(int j=0;j<State_.size();j++){
    //            if(State_[j]!=0.0){
    //                outfile_Variational_states<<j<<"\t"<<State_[j]<<endl;
    //            }
    //        }
    //        outfile_Variational_states.close();


    //        //Orbital_Exchange
    //        State_temp1 = Act_Orbital_Exchange(basis, State_);
    //        OEn++;
    //        BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp1);
    //        temp1=dot_product(GS_,State_temp1);
    //        Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
    //        state_count++;
    //        cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
    //              ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
    //        cout<<temp1<<endl;
    //        outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
    //        ssn_int.str("");
    //        ssn_int.clear();
    //        ssn_int << state_count;
    //        file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
    //                "_state_" + ssn_int.str() +".txt";
    //        outfile_Variational_states.open(file_out_Variational_states.c_str());
    //        outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
    //                                 <<variational_state_contruction_<<endl;

    //        for(int j=0;j<State_temp1.size();j++){
    //            if(State_temp1[j]!=0.0){
    //                outfile_Variational_states<<j<<"\t"<<State_temp1[j]<<endl;
    //            }
    //        }
    //        outfile_Variational_states.close();

    //        //TRANSLATIONS
    //        for(int T_no=1;T_no<=basis.Length-1;T_no++){
    //            State_temp2 = Act_Translational_operator(basis, State_temp1);
    //            Tn++;
    //            BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp2);
    //            temp1=dot_product(GS_,State_temp2);
    //            Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
    //            state_count++;
    //            cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
    //                  ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
    //            cout<<temp1<<endl;
    //            outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
    //            ssn_int.str("");
    //            ssn_int.clear();
    //            ssn_int << state_count;
    //            file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
    //                    "_state_" + ssn_int.str() +".txt";
    //            outfile_Variational_states.open(file_out_Variational_states.c_str());
    //            outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
    //                                     <<variational_state_contruction_<<endl;

    //            for(int j=0;j<State_temp2.size();j++){
    //                if(State_temp2[j]!=0.0){
    //                    outfile_Variational_states<<j<<"\t"<<State_temp2[j]<<endl;
    //                }
    //            }
    //            outfile_Variational_states.close();


    //            //Orbital_Exchange
    //            State_temp1 = Act_Orbital_Exchange(basis, State_temp2);
    //            OEn++;
    //            BASIS_STATES_ANSATZ_2HOLES.push_back(State_temp1);
    //            temp1=dot_product(GS_,State_temp1);
    //            Ansatz_2holes_Basis_Overlap_with_GS.push_back(temp1);
    //            state_count++;
    //            cout<<"CLASS NO. "<<i<<", STATE NO. "<<state_count<<
    //                  ", State type = (Trans)"<<Tn<<"_(Orb_Exc)"<<OEn<<"("<<variational_state_contruction_<<")"<<endl;
    //            cout<<temp1<<endl;
    //            outfile_out_overlap_GS<<(i*12)+state_count<<"   "<<temp1<<endl;
    //            ssn_int.str("");
    //            ssn_int.clear();
    //            ssn_int << state_count;
    //            file_out_Variational_states = "Variational_state_2holes_Class_" + ss_int.str() +
    //                    "_state_" + ssn_int.str() +".txt";
    //            outfile_Variational_states.open(file_out_Variational_states.c_str());
    //            outfile_Variational_states<<"#VARIATIONAL CLASS NO. "<<i<<", STATE NO. "<<state_count<<", "
    //                                     <<variational_state_contruction_<<endl;

    //            for(int j=0;j<State_temp1.size();j++){
    //                if(State_temp1[j]!=0.0){
    //                    outfile_Variational_states<<j<<"\t"<<State_temp1[j]<<endl;
    //                }
    //            }
    //            outfile_Variational_states.close();

    //        }

    //        assert(state_count+1 == States_in_class[i]);





    //        //-----------------------------




    //    }



    //    string file_out_overlaps = "overlap_bw_RVB_2holes_states.txt";

    //    ofstream outfile_out_overlaps(file_out_overlaps.c_str());
    //    outfile_out_overlaps<<"# state_no(i)   state_no(j)  <i|j> "<<endl;

    //    OVERLAP_MATRIX_2holes_Basis.resize(BASIS_STATES_ANSATZ_2HOLES.size());
    //    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
    //        OVERLAP_MATRIX_2holes_Basis[i].resize(BASIS_STATES_ANSATZ_2HOLES.size());
    //    }

    //    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
    //        for(int j=0;j<BASIS_STATES_ANSATZ_2HOLES.size();j++){

    //            OVERLAP_MATRIX_2holes_Basis[i][j]=dot_product(BASIS_STATES_ANSATZ_2HOLES[i],BASIS_STATES_ANSATZ_2HOLES[j]);
    //            outfile_out_overlaps<<i<<"  "<<j<< "   "<<OVERLAP_MATRIX_2holes_Basis[i][j]<<endl;
    //            //cout<<"i="<<i<<", j="<<j<<", overlap="<<OVERLAP_MATRIX_2holes_Basis[i][j]<<endl;
    //        }
    //    }

    //    /*
    //    for(int i=0;i<BASIS_STATES_ANSATZ_2HOLES.size();i++){
    //        for(int j=0;j<BASIS_STATES_ANSATZ_2HOLES.size();j++){
    //            cout<<dot_product(BASIS_STATES_ANSATZ_2HOLES[i],BASIS_STATES_ANSATZ_2HOLES[j])<<"  ";
    //        }
    //        cout<<endl;
    //    }
    //    */





}

void MODEL_2_orb_Hubb_chain::Get_overlap_matrix_for_Anzatz_basis(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){

    //    int No_of_basis;


    //    if(basis.Length==6){
    //        No_of_basis=64;
    //        variational_state_contruction_strings.resize(No_of_basis);
    //        for(int i=0;i<64;i++){
    //            if(i==0){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==1){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==2){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==3){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==4){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==5){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==6){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==7){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==8){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==9){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==10){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==11){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==12){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==13){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==14){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==15){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==16){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==17){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==18){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==19){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==20){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==21){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==22){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==23){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==24){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==25){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==26){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==27){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==28){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==29){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==30){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==31){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==32){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==33){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==34){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==35){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==36){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==37){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==38){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==39){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==40){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==41){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==42){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==43){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==44){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==45){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==46){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==47){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==48){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==49){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==50){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==51){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
    //            }
    //            if(i==52){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==53){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==54){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==55){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==56){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==57){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==58){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //            if(i==59){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
    //            }
    //            if(i==60){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
    //            }
    //            if(i==61){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
    //            }
    //            if(i==62){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
    //            }
    //            if(i==63){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
    //            }
    //        }

    //    }


    //    if(basis.Length==4){
    //        No_of_basis=16;
    //        variational_state_contruction_strings.resize(No_of_basis);
    //        for(int i=0;i<16;i++){
    //            if(i==0){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 0";
    //            }
    //            if(i==1){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 0";
    //            }
    //            if(i==2){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 0";
    //            }
    //            if(i==3){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 0";
    //            }
    //            if(i==4){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 0";
    //            }
    //            if(i==5){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 0";
    //            }
    //            if(i==6){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 0";
    //            }
    //            if(i==7){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 0";
    //            }
    //            if(i==8){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 0";
    //            }
    //            if(i==9){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 0";
    //            }
    //            if(i==10){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 0";
    //            }
    //            if(i==11){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 0";
    //            }
    //            if(i==12){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 0";
    //            }
    //            if(i==13){
    //                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 0";
    //            }
    //            if(i==14){
    //                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 0";
    //            }
    //            if(i==15){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 0";
    //            }
    //        }

    //    }

    //    if(basis.Length==2){
    //        No_of_basis=2;
    //        variational_state_contruction_strings.resize(No_of_basis);

    //        for(int i=0;i<2;i++){
    //            if(i==0){
    //                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 0";
    //            }
    //            if(i==1){
    //                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 0";
    //            }

    //        }

    //    }


    //    int temp_int;

    //    BASIS_STATES_ANSATZ.resize(No_of_basis);

    //    string variational_state_contruction_;

    //    cout<<"-------OVERLAPS OF ANSATZ BASIS WITH GS--------------"<<endl;






    //    Ansatz_Basis_Overlap_with_GS.resize(No_of_basis);
    //    for(int i=0;i<No_of_basis;i++){
    //        ostringstream ss_int;
    //        ss_int << i;

    //        string file_out_Variational_states = "Variational_state" + ss_int.str() + ".txt";
    //        ofstream outfile_Variational_states(file_out_Variational_states.c_str());

    //        variational_state_contruction_=variational_state_contruction_strings[i];


    //        //------------UPDATING Variational_state_pair******--------

    //        stringstream variational_state_contruction_stream(variational_state_contruction_);
    //        for(int n=0;n<basis.Length;n++){
    //            variational_state_contruction_stream >> Variational_state_pair_Geometry[n];
    //            assert(Variational_state_pair_Geometry[n]=="D" ||
    //                   Variational_state_pair_Geometry[n]=="AC");

    //            variational_state_contruction_stream >> Variational_state_pair_spin_symmetry[n];
    //            assert(Variational_state_pair_spin_symmetry[n]=="S" ||
    //                   Variational_state_pair_spin_symmetry[n]=="T");

    //            variational_state_contruction_stream >> Variational_state_pair_orbital_symmetry[n];
    //            assert(Variational_state_pair_orbital_symmetry[n] =="OA" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="OS" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="NOS1" ||
    //                   Variational_state_pair_orbital_symmetry[n] =="NOS2");

    //            variational_state_contruction_stream >> temp_int;
    //            Variational_state_pair_sites[n].first =temp_int;
    //            assert(temp_int>=0 && temp_int<basis.Length);


    //            variational_state_contruction_stream >> temp_int;
    //            Variational_state_pair_sites[n].second =temp_int;
    //            assert(temp_int>=0 && temp_int<basis.Length);


    //        }

    //        //-----------------------------




    //        Get_Variational_State(basis,basis.Length);
    //        BASIS_STATES_ANSATZ[i]=State_;

    //        Ansatz_Basis_Overlap_with_GS[i]=dot_product(GS_,BASIS_STATES_ANSATZ[i]);
    //        cout<<i<<" ,State type = "<<variational_state_contruction_<<endl;
    //        cout<<Ansatz_Basis_Overlap_with_GS[i]<<endl;

    //        outfile_Variational_states<<"#VARIATIONAL STATE NO. "<<i<<",  "<<variational_state_contruction_<<endl;
    //        for(int j=0;j<BASIS_STATES_ANSATZ[i].size();j++){
    //            if(BASIS_STATES_ANSATZ[i][j]!=0.0){
    //                outfile_Variational_states<<j<<"\t"<<BASIS_STATES_ANSATZ[i][j]<<endl;
    //            }
    //        }


    //    }
    //    cout<<"---------------------------------------------------"<<endl;


    //    overlap_matrix_for_Anzatz_basis.resize(No_of_basis);
    //    for(int i=0;i<No_of_basis;i++){
    //        overlap_matrix_for_Anzatz_basis[i].resize(No_of_basis);
    //        for(int j=0;j<No_of_basis;j++){
    //            overlap_matrix_for_Anzatz_basis[i][j]=dot_product(BASIS_STATES_ANSATZ[i],BASIS_STATES_ANSATZ[j]);
    //        }
    //    }


}



void MODEL_2_orb_Hubb_chain::Create_OS_TS_states_by_reading(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){

    //    int No_of_basis, No_of_classes;
    //    if(basis.Length==2){
    //        No_of_basis=2;
    //        No_of_classes=2;
    //    }
    //    if(basis.Length==4){
    //        No_of_basis=16;
    //        No_of_classes=4;
    //    }
    //    if(basis.Length==6){
    //        No_of_basis=64;
    //        No_of_classes=8;
    //    }


    //    STATES_OS_TS.clear();
    //    STATES_OS_TS.resize(No_of_classes);
    //    assert(BASIS_STATES_ANSATZ.size() == No_of_basis);

    //    Mat_2_int Classes_;
    //    Classes_.resize(No_of_classes);

    //    if(basis.Length==6){
    //        for(int i=0;i<No_of_classes;i++){
    //            Classes_[i].clear();

    //            if(i==0){
    //                for(int j=0;j<2;j++){ //change to j=0;j<2
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==1){
    //                for(int j=2;j<14;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==2){
    //                for(int j=14;j<26;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==3){
    //                for(int j=26;j<28;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==4){
    //                for(int j=28;j<34;j++){  //change to j=28;j<34
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==5){
    //                for(int j=34;j<40;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==6){
    //                for(int j=40;j<46;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //                for(int j=52;j<58;j++){
    //                    Classes_[i].push_back(j);
    //                }

    //            }
    //            if(i==7){
    //                for(int j=46;j<52;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //                for(int j=58;j<64;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }

    //        }
    //    }


    //    if(basis.Length==4){
    //        for(int i=0;i<No_of_classes;i++){
    //            Classes_[i].clear();

    //            if(i==0){
    //                for(int j=0;j<2;j++){ //change to j=0;j<2
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==1){
    //                for(int j=2;j<4;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==2){
    //                for(int j=4;j<12;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //            if(i==3){
    //                for(int j=12;j<16;j++){
    //                    Classes_[i].push_back(j);
    //                }
    //            }
    //        }
    //    }




    //    for(int i=0;i<No_of_classes;i++){

    //        STATES_OS_TS[i].clear();
    //        STATES_OS_TS[i].resize(BASIS_STATES_ANSATZ[0].size());

    //        for(int j_ind=0;j_ind<Classes_[i].size();j_ind++){
    //            for(int m=0;m<STATES_OS_TS[i].size();m++){
    //                STATES_OS_TS[i][m] += (BASIS_STATES_ANSATZ[Classes_[i][j_ind]][m]);
    //            }
    //        }


    //        cout<<"overlap of class "<<i<<" with GS = "<<dot_product(GS_,STATES_OS_TS[i])/(Classes_[i].size())<<", No of states = "<<Classes_[i].size()<<endl;


    //    }

    //    //cout<<"here 0"<<endl;


}


void MODEL_2_orb_Hubb_chain::Create_STATES_OS_TS_DIAGONAL_HOLES(BASIS_2_orb_Hubb_chain &basis,
                                                                BASIS_2_orb_Hubb_chain &basis_nm2,
                                                                Mat_2_doub &STATES_OS_TS_){



    //    assert( (basis_nm2.Length -1) == basis_nm2.Ndn);
    //    assert( (basis_nm2.Length -1) == basis_nm2.Nup);



    //    STATES_OS_TS_DIAGONAL_HOLES.resize(STATES_OS_TS_.size());


    //    //Creating [Delta_{D}|V_{i}(N)>] states

    //    double alpha_S, alpha_O;
    //    Mat_1_doub vec1, vec2, vec_total;
    //    trio_int trio_0, trio_1;


    //    for(int ci=0;ci<STATES_OS_TS_DIAGONAL_HOLES.size();ci++){

    //        vec_total.clear();
    //        vec_total.resize(basis_nm2.D_dn_basis.size()*basis_nm2.D_up_basis.size());

    //        for(int site=0;site<basis_nm2.Length;site++){
    //            // for(int site=0;site<1;site++){
    //            int site_p1= site +1;
    //            if(site_p1==basis_nm2.Length){
    //                site_p1=0;
    //            }
    //            //----------Spin-Singlet, Orbital Antisymmetric pair Anhilation From Variational states------------//
    //            alpha_S =-1.0;alpha_O=-1.0;

    //            //ALONG THE DIAGONAL
    //            trio_0.orb_ = 1; trio_0.spin_=1; trio_0.site_=site_p1;
    //            trio_1.orb_ = 0; trio_1.spin_=0; trio_1.site_=site;
    //            Get_Pair_Operator_Matrix(basis, basis_nm2, trio_0, trio_1, 1.0*alpha_S);
    //            Matrix_COO_vector_multiplication("full",Pair_Annihilation , STATES_OS_TS_[ci], vec1);

    //            trio_0.orb_ = 1; trio_0.spin_=0; trio_0.site_=site_p1;
    //            trio_1.orb_ = 0; trio_1.spin_=1; trio_1.site_=site;
    //            Get_Pair_Operator_Matrix(basis, basis_nm2, trio_0, trio_1, 1.0*1.0);
    //            Matrix_COO_vector_multiplication("full",Pair_Annihilation , STATES_OS_TS_[ci], vec2);
    //            Subtract(vec1, -1.0, vec2);


    //            trio_0.orb_ = 0; trio_0.spin_=1; trio_0.site_=site_p1;
    //            trio_1.orb_ = 1; trio_1.spin_=0; trio_1.site_=site;
    //            Get_Pair_Operator_Matrix(basis, basis_nm2, trio_0, trio_1, 1.0*alpha_S*alpha_O);  //change to 1.0*alp--
    //            Matrix_COO_vector_multiplication("full",Pair_Annihilation , STATES_OS_TS_[ci], vec2);
    //            Subtract(vec1, -1.0, vec2);

    //            trio_0.orb_ = 0; trio_0.spin_=0; trio_0.site_=site_p1;
    //            trio_1.orb_ = 1; trio_1.spin_=1; trio_1.site_=site;
    //            Get_Pair_Operator_Matrix(basis, basis_nm2, trio_0, trio_1, 1.0*alpha_O);
    //            Matrix_COO_vector_multiplication("full",Pair_Annihilation , STATES_OS_TS_[ci], vec2);  //change to 1.0*alpha--
    //            Subtract(vec1, -1.0, vec2);

    //            //----------------------------------------//

    //            Subtract(vec_total, -1.0, vec1);

    //        }

    //        STATES_OS_TS_DIAGONAL_HOLES[ci]=vec_total;


    //    }
    //    //--------------Created [Delta_{D}|V_{i}(N)>] states--------



}


void MODEL_2_orb_Hubb_chain::Optimize_Anzatz_basis_2Holes(BASIS_2_orb_Hubb_chain &basis_nm2, Mat_1_doub &OPT_VEC,
                                                          string read_overlaps){


    //    int index,index_i,index_j;
    //    double val_;
    //    double overlap_;
    //    double norm_;
    //    double VALUE_CHECK, VALUE_CHECK_OLD;
    //    Mat_1_doub Alpha_;
    //    Mat_1_doub Alpha_Class;
    //    string sLine;

    //    if(read_overlaps == "READ_OVERLAPS"){


    //        BASIS_STATES_ANSATZ_2HOLES.clear();
    //        BASIS_STATES_ANSATZ_2HOLES.resize(12*36);
    //        for(int cl_no=0;cl_no<36;cl_no++){
    //            ostringstream cl_int;
    //            cl_int << cl_no;

    //            for(int s_no=0;s_no<12;s_no++){
    //                BASIS_STATES_ANSATZ_2HOLES[s_no + (cl_no*12)].clear();
    //                BASIS_STATES_ANSATZ_2HOLES[s_no + (cl_no*12)].resize(basis_nm2.D_up_basis.size()*basis_nm2.D_dn_basis.size());
    //                ostringstream s_int;
    //                s_int << s_no;
    //                ifstream infile_BASIS;
    //                string file_in_BASIS = "Variational_state_2holes_Class_"+ cl_int.str() +"_state_"+s_int.str() + ".txt";
    //                infile_BASIS.open(file_in_BASIS.c_str());
    //                cout<<"using "<<file_in_BASIS<<endl;

    //                getline(infile_BASIS, sLine);
    //                while(!infile_BASIS.eof()){
    //                    getline(infile_BASIS, sLine);
    //                    stringstream ssLine;
    //                    ssLine << sLine;
    //                    ssLine>>index;
    //                    ssLine>>val_;
    //                    BASIS_STATES_ANSATZ_2HOLES[s_no + (cl_no*12)][index]=val_;
    //                    // cout<<index<<"  "<<val_<<endl;

    //                }

    //                cout<<"state: class="<<cl_no<<", s_no="<<s_no<<" read"<<endl;


    //            }
    //        }





    //        Ansatz_2holes_Basis_Overlap_with_GS.clear();


    //        cout<<"OVERLAPS for 2-holes ansatz are read"<<endl;

    //        ifstream outfile_out_overlap_GS;
    //        string file_out_overlap_GS = "RVB_2holes_GSNm2_overlap.txt";
    //        outfile_out_overlap_GS.open(file_out_overlap_GS.c_str());

    //        Ansatz_2holes_Basis_Overlap_with_GS.resize(432);
    //        getline(outfile_out_overlap_GS, sLine);
    //        //cout<<sLine<<endl;
    //        while(!outfile_out_overlap_GS.eof()){
    //            getline(outfile_out_overlap_GS, sLine);
    //            //cout<<sLine<<endl;
    //            stringstream ssLine;
    //            ssLine << sLine;
    //            ssLine>>index;
    //            ssLine>>val_;
    //            Ansatz_2holes_Basis_Overlap_with_GS[index]=val_;
    //        }

    //        //cout<<"Ansatz_2holes_Basis_Overlap_with_GS.size = "<<Ansatz_2holes_Basis_Overlap_with_GS.size()<<endl;
    //        //assert(Ansatz_2holes_Basis_Overlap_with_GS.size()==36*12);
    //        // cout<<"here3"<<endl;


    //        OVERLAP_MATRIX_2holes_Basis.clear();
    //        OVERLAP_MATRIX_2holes_Basis.resize(432);
    //        for(int i=0;i<432;i++){
    //            OVERLAP_MATRIX_2holes_Basis[i].resize(432);
    //        }
    //        ifstream outfile_out_overlaps;
    //        string file_out_overlaps = "overlap_bw_RVB_2holes_states.txt";
    //        outfile_out_overlaps.open(file_out_overlaps.c_str());

    //        getline(outfile_out_overlaps, sLine);
    //        while(!outfile_out_overlaps.eof()){
    //            getline(outfile_out_overlaps, sLine);
    //            stringstream ssLine;
    //            ssLine << sLine;
    //            ssLine>>index_i;
    //            ssLine>>index_j;
    //            ssLine>>val_;
    //            OVERLAP_MATRIX_2holes_Basis[index_i][index_j] = val_;
    //        }

    //    }
    //    else{
    //        assert(Ansatz_2holes_Basis_Overlap_with_GS.size()==36*12);
    //    }


    //    Alpha_Class.resize(36);

    //    Alpha_.resize(Ansatz_2holes_Basis_Overlap_with_GS.size());
    //    assert(Alpha_.size() == (36*12));

    //    Mat_1_doub Alpha_min;
    //    Alpha_min.resize(Ansatz_2holes_Basis_Overlap_with_GS.size());

    //    Mat_1_doub Alpha_Final;
    //    Alpha_Final.resize(Ansatz_2holes_Basis_Overlap_with_GS.size());

    //    Mat_1_doub Alpha_max;
    //    Alpha_max.resize(Ansatz_2holes_Basis_Overlap_with_GS.size());

    //    Mat_1_doub Dalpha_;
    //    Dalpha_.resize(Ansatz_2holes_Basis_Overlap_with_GS.size());



    //    Mat_1_doub Alpha_min_Class;
    //    Alpha_min_Class.resize(36);

    //    Mat_1_doub Alpha_Final_Class;
    //    Alpha_Final_Class.resize(36);

    //    Mat_1_doub Alpha_max_Class;
    //    Alpha_max_Class.resize(36);

    //    Mat_1_doub Dalpha_Class;
    //    Dalpha_Class.resize(36);




    //    for(int i=0;i<16;i++){
    //        Dalpha_Class[i]=0.1;
    //        Alpha_min_Class[i]=-1.0;
    //        Alpha_max_Class[i]=1.0;
    //    }

    //    for(int i=16;i<20;i++){
    //        Dalpha_Class[i]=0.05;
    //        Alpha_min_Class[i]=0.0;
    //        Alpha_max_Class[i]=0.0;
    //    }



    //    for(int i=20;i<28;i++){
    //        if(i==23){
    //            Dalpha_Class[i]=0.1;
    //            Alpha_min_Class[i]=0.0;
    //            Alpha_max_Class[i]=0.0;
    //        }
    //        else{
    //            Dalpha_Class[i]=0.1;
    //            Alpha_min_Class[i]=0.0;
    //            Alpha_max_Class[i]=0.0;
    //        }
    //    }
    //    for(int i=28;i<36;i++){
    //        Dalpha_Class[i]=0.1;
    //        Alpha_min_Class[i]=0.0;
    //        Alpha_max_Class[i]=0.0;
    //    }



    //    /*
    //    Alpha_min_Class[0]=0.0863;Alpha_max_Class[0]=0.0863;
    //    Alpha_min_Class[1]=0.076;Alpha_max_Class[1]=0.076;
    //    Alpha_min_Class[2]=0.0866;Alpha_max_Class[2]=0.0866;
    //    Alpha_min_Class[3]=0.0743;Alpha_max_Class[3]=0.0743;
    //    Alpha_min_Class[4]=0.0827;Alpha_max_Class[4]=0.0827;
    //    Alpha_min_Class[5]=0.0728;Alpha_max_Class[5]=0.0728;
    //    Alpha_min_Class[6]=0.0736;Alpha_max_Class[6]=0.0736;
    //    Alpha_min_Class[8]=0.0647;Alpha_max_Class[8]=0.0647;
    //    Alpha_min_Class[9]=0.0655;Alpha_max_Class[9]=0.0655;
    //    Alpha_min_Class[12]=0.0642;Alpha_max_Class[12]=0.0642;
    //    */
    //    Alpha_min_Class[0]=4.65597022195399709577e-02;Alpha_max_Class[0]=4.65597022195399709577e-02;
    //    Alpha_min_Class[1]=3.32569301568142669523e-02;Alpha_max_Class[1]=3.32569301568142669523e-02;
    //    Alpha_min_Class[5]=2.66055441254514184191e-02;Alpha_max_Class[5]=2.66055441254514184191e-02;
    //    Alpha_min_Class[8]=1.99541580940885664164e-02;Alpha_max_Class[8]=1.99541580940885664164e-02;


    //    /*
    //    Alpha_min_Class[16]=0.0751;Alpha_max_Class[16]=0.0751;
    //    Alpha_min_Class[17]=0.0763;Alpha_max_Class[17]=0.0763;
    //    Alpha_min_Class[18]=0.0721;Alpha_max_Class[18]=0.0721;
    //*/
    //    //Optimized Class-set 2 [-ive signs are added internally]
    //    Alpha_min_Class[16]=1.00335009313597803282e-01;Alpha_max_Class[16]=1.00335009313597803282e-01;
    //    Alpha_min_Class[17]=1.00335009313597803282e-01;Alpha_max_Class[17]=1.00335009313597803282e-01;
    //    Alpha_min_Class[18]=8.36125077613315004221e-02;Alpha_max_Class[18]=8.36125077613315004221e-02;


    //    /*
    //    Alpha_min_Class[20]=0.0580;Alpha_max_Class[20]=0.0580;
    //    Alpha_min_Class[21]=0.0576;Alpha_max_Class[21]=0.0576;
    //    Alpha_min_Class[23]=0.0603;Alpha_max_Class[23]=0.0603;
    //    Alpha_min_Class[24]=0.0557;Alpha_max_Class[24]=0.0557;
    //    Alpha_min_Class[26]=0.0573;Alpha_max_Class[26]=0.0573;
    //*/



    //    //Optimized Class-set 3 [-ive signs are added internally]
    //    Alpha_min_Class[20]=7.13589189776528964870e-02;Alpha_max_Class[20]=7.13589189776528964870e-02;
    //    Alpha_min_Class[21]=2.67595946166198465910e-02;Alpha_max_Class[21]=2.67595946166198465910e-02;
    //    Alpha_min_Class[23]=4.45993243610330256099e-02;Alpha_max_Class[23]=4.45993243610330256099e-02;
    //    Alpha_min_Class[24]=1.33797973083099371733e-02;Alpha_max_Class[24]=1.33797973083099371733e-02;
    //    Alpha_min_Class[26]=7.58188514137562025175e-02;Alpha_max_Class[26]=7.58188514137562025175e-02;


    //    /*
    //    Alpha_min_Class[28]=0.0580;Alpha_max_Class[28]=0.0580;
    //    Alpha_min_Class[29]=0.0572;Alpha_max_Class[29]=0.0572;
    //    Alpha_min_Class[30]=0.0600;Alpha_max_Class[30]=0.0600;
    //    Alpha_min_Class[33]=0.0607;Alpha_max_Class[33]=0.0607;
    //    */

    //    //Optimized Class-set 4 [-ive signs are added internally]
    //    Alpha_min_Class[28]=4.06353010052327642820e-02;Alpha_max_Class[28]=4.06353010052327642820e-02;
    //    Alpha_min_Class[29]=4.87623612062793101996e-02;Alpha_max_Class[29]=4.87623612062793101996e-02;
    //    Alpha_min_Class[30]=6.50164816083724089735e-02;Alpha_max_Class[30]=6.50164816083724089735e-02;
    //    Alpha_min_Class[33]=6.50164816083724089735e-02;Alpha_max_Class[33]=6.50164816083724089735e-02;





    //    /*
    //    Alpha_min_Class[0]=0.0863;Alpha_max_Class[0]=0.0863;
    //    Alpha_min_Class[1]=-0.076;Alpha_max_Class[0]=-0.076;
    //    Alpha_min_Class[2]=-0.0866;Alpha_max_Class[0]=-0.0866;
    //    Alpha_min_Class[3]=0.0743;Alpha_max_Class[0]=0.0743;
    //    Alpha_min_Class[4]=-0.0827;Alpha_max_Class[0]=-0.0827;
    //    Alpha_min_Class[5]=0.0728;Alpha_max_Class[0]=0.0728;
    //    Alpha_min_Class[6]=-0.0736;Alpha_max_Class[0]=-0.0736;
    //    Alpha_min_Class[8]=-0.0647;Alpha_max_Class[0]=-0.0647;
    //    Alpha_min_Class[9]=0.0655;Alpha_max_Class[0]=0.0655;
    //    Alpha_min_Class[12]=-0.0642;Alpha_max_Class[0]=-0.0642;
    //    */






    //    double a_CS1, a_CS2, a_CS3, a_CS4;
    //    double a_CS1_final, a_CS2_final, a_CS3_final, a_CS4_final;

    //    double a_CS1_min, a_CS2_min, a_CS3_min, a_CS4_min;
    //    double a_CS1_max, a_CS2_max, a_CS3_max, a_CS4_max;

    //    double d_a_CS;


    //    //a_CS1=-9.00000000000000022204e-01  a_CS2=-5.00000000000000111022e-01
    //    //a_CS3=-3.00000000000000155431e-01  a_CS4=-4.00000000000000133227e-01
    //    //-----------------
    //    a_CS1_min=-9.00000000000000022204e-01;a_CS1_max=-9.00000000000000022204e-01;
    //    a_CS2_min=-5.00000000000000111022e-01;a_CS2_max=-5.00000000000000111022e-01;
    //    a_CS3_min=-3.00000000000000155431e-01;a_CS3_max=-3.00000000000000155431e-01;
    //    a_CS4_min=-4.00000000000000133227e-01;a_CS4_max=-4.00000000000000133227e-01;
    //    d_a_CS = 0.1;

    //    //------------------

    //    /*a_CS1=-6.00000000000000088818e-01
    //     a_CS2=-9.00000000000000022204e-01
    //            a_CS3=-3.00000000000000155431e-01
    //            a_CS4=-4.00000000000000133227e-01
    //            */

    //    VALUE_CHECK_OLD=0;
    //    for(a_CS1=a_CS1_min;a_CS1<=a_CS1_max;a_CS1=a_CS1+d_a_CS){

    //        for(a_CS2=a_CS2_min;a_CS2<=a_CS2_max;a_CS2=a_CS2+d_a_CS){

    //            for(a_CS3=a_CS3_min;a_CS3<=a_CS3_max;a_CS3=a_CS3+d_a_CS){

    //                for(a_CS4=a_CS4_min;a_CS4<=a_CS4_max;a_CS4=a_CS4+d_a_CS){

    //                    Alpha_Class[0]=Alpha_min_Class[0];
    //                    while(Alpha_Class[0]<=Alpha_max_Class[0]){
    //                        Alpha_Class[7]=Alpha_Class[0];
    //                        Alpha_Class[2]=Alpha_Class[0];
    //                        Alpha_Class[4]=Alpha_Class[0];


    //                        Alpha_Class[1]=Alpha_min_Class[1];
    //                        while(Alpha_Class[1]<=Alpha_max_Class[1]){
    //                            Alpha_Class[11]=Alpha_Class[1];
    //                            Alpha_Class[3]=Alpha_Class[1];
    //                            Alpha_Class[13]=Alpha_Class[1];

    //                            //   Alpha_Class[2]=Alpha_min_Class[2];
    //                            //  while(Alpha_Class[2]<=Alpha_max_Class[2]){

    //                            //   Alpha_Class[3]=Alpha_min_Class[3];
    //                            //  while(Alpha_Class[3]<=Alpha_max_Class[3]){
    //                            //     Alpha_Class[13]=Alpha_Class[3];

    //                            //    Alpha_Class[4]=Alpha_min_Class[4];
    //                            //   while(Alpha_Class[4]<=Alpha_max_Class[4]){

    //                            Alpha_Class[5]=Alpha_min_Class[5];
    //                            while(Alpha_Class[5]<=Alpha_max_Class[5]){
    //                                Alpha_Class[10]=Alpha_Class[5];
    //                                Alpha_Class[6]=Alpha_Class[5];
    //                                Alpha_Class[15]=Alpha_Class[5];

    //                                //     Alpha_Class[6]=Alpha_min_Class[6];
    //                                //    while(Alpha_Class[6]<=Alpha_max_Class[6]){
    //                                //       Alpha_Class[15]=Alpha_Class[6];


    //                                Alpha_Class[8]=Alpha_min_Class[8];
    //                                while(Alpha_Class[8]<=Alpha_max_Class[8]){
    //                                    Alpha_Class[12]=Alpha_Class[8];
    //                                    Alpha_Class[9]=Alpha_Class[8];
    //                                    Alpha_Class[14]=Alpha_Class[8];


    //                                    //      Alpha_Class[9]=Alpha_min_Class[9];
    //                                    //     while(Alpha_Class[9]<=Alpha_max_Class[9]){
    //                                    //        Alpha_Class[14]=Alpha_Class[9];

    //                                    //      Alpha_Class[12]=Alpha_min_Class[12];
    //                                    //     while(Alpha_Class[12]<=Alpha_max_Class[12]){



    //                                    Alpha_Class[16]=Alpha_min_Class[16];
    //                                    while(Alpha_Class[16]<=Alpha_max_Class[16]){
    //                                        Alpha_Class[19]=Alpha_Class[16];


    //                                        Alpha_Class[17]=Alpha_min_Class[17];
    //                                        while(Alpha_Class[17]<=Alpha_max_Class[17]){
    //                                            Alpha_Class[18]=Alpha_min_Class[18];
    //                                            while(Alpha_Class[18]<=Alpha_max_Class[18]){


    //                                                Alpha_Class[20]=Alpha_min_Class[20];
    //                                                while(Alpha_Class[20]<=Alpha_max_Class[20]){
    //                                                    Alpha_Class[22]=Alpha_Class[20];
    //                                                    Alpha_Class[25]=Alpha_Class[20];

    //                                                    Alpha_Class[21]=Alpha_min_Class[21];
    //                                                    while(Alpha_Class[21]<=Alpha_max_Class[21]){
    //                                                        Alpha_Class[27]=Alpha_Class[21];


    //                                                        Alpha_Class[23]=Alpha_min_Class[23];
    //                                                        while(Alpha_Class[23]<=Alpha_max_Class[23]){

    //                                                            Alpha_Class[24]=Alpha_min_Class[24];
    //                                                            while(Alpha_Class[24]<=Alpha_max_Class[24]){


    //                                                                Alpha_Class[26]=Alpha_min_Class[26];
    //                                                                while(Alpha_Class[26]<=Alpha_max_Class[26]){




    //                                                                    Alpha_Class[28]=Alpha_min_Class[28];
    //                                                                    while(Alpha_Class[28]<=Alpha_max_Class[28]){
    //                                                                        Alpha_Class[35]=1.0*Alpha_Class[28];


    //                                                                        Alpha_Class[29]=Alpha_min_Class[29];
    //                                                                        while(Alpha_Class[29]<=Alpha_max_Class[29]){
    //                                                                            Alpha_Class[31]=1.0*Alpha_Class[29];


    //                                                                            Alpha_Class[30]=Alpha_min_Class[30];
    //                                                                            while(Alpha_Class[30]<=Alpha_max_Class[30]){
    //                                                                                Alpha_Class[32]=1.0*Alpha_Class[30];


    //                                                                                Alpha_Class[33]=Alpha_min_Class[33];
    //                                                                                while(Alpha_Class[33]<=Alpha_max_Class[33]){
    //                                                                                    Alpha_Class[34]=1.0*Alpha_Class[33];



    //                                                                                    for(int i=0;i<35;i++){
    //                                                                                        for(int j=i*12;j<((i*12) + 12);j++){
    //                                                                                            if(i>=0 && i<16){
    //                                                                                                Alpha_[j]= a_CS1*Alpha_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[j]);
    //                                                                                            }
    //                                                                                            if(i>=16 && i<20){
    //                                                                                                Alpha_[j]= a_CS2*Alpha_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[j]);
    //                                                                                            }
    //                                                                                            if(i>=20 && i<28){
    //                                                                                                Alpha_[j]= a_CS3*Alpha_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[j]);
    //                                                                                            }
    //                                                                                            if(i>=28 && i<36){
    //                                                                                                Alpha_[j]= a_CS4*Alpha_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[j]);
    //                                                                                            }
    //                                                                                        }

    //                                                                                    }



    //                                                                                    /*
    //                    for(int i=192;i<240;i++){
    //                        Alpha_[i]=a_CS2;
    //                        Alpha_[i]= Alpha_[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[i]);
    //                    }
    //                    for(int i=240;i<336;i++){
    //                        Alpha_[i]=a_CS3;
    //                        Alpha_[i]= Alpha_[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[i]);
    //                    }
    //                    for(int i=336;i<432;i++){
    //                        Alpha_[i]=a_CS4;
    //                        Alpha_[i]= Alpha_[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[i]);
    //                    }
    //                    */






    //                                                                                    overlap_=0.0;
    //                                                                                    for(int i=0;i<432;i++){
    //                                                                                        overlap_ += Ansatz_2holes_Basis_Overlap_with_GS[i]*Alpha_[i];
    //                                                                                    }



    //                                                                                    norm_=0.0;
    //                                                                                    for(int i=0;i<432;i++){
    //                                                                                        for(int j=0;j<432;j++){
    //                                                                                            norm_ += OVERLAP_MATRIX_2holes_Basis[i][j]*Alpha_[i]*Alpha_[j];
    //                                                                                        }
    //                                                                                    }



    //                                                                                    VALUE_CHECK = (overlap_*overlap_)/(norm_);

    //                                                                                    if( VALUE_CHECK > VALUE_CHECK_OLD){
    //                                                                                        VALUE_CHECK_OLD=VALUE_CHECK;

    //                                                                                        Alpha_Final=Alpha_;
    //                                                                                        Alpha_Final_Class =Alpha_Class;
    //                                                                                        a_CS1_final = a_CS1; a_CS2_final = a_CS2;
    //                                                                                        a_CS3_final = a_CS3; a_CS4_final = a_CS4;

    //                                                                                    }

    //                                                                                    cout<<scientific<<setprecision(4);
    //                                                                                    cout<<"a_CS1="<<a_CS1<<"  "<<"a_CS2="<<a_CS2<<"  "<<
    //                                                                                          "a_CS3="<<a_CS3<<"  "<<"a_CS4="<<a_CS4<<endl;

    //                                                                                    for(int i=0;i<36;i++){
    //                                                                                        cout<<Alpha_Class[i]<<"   ";
    //                                                                                    }

    //                                                                                    cout<<scientific<<setprecision(20);
    //                                                                                    cout <<"  value="<<VALUE_CHECK<<endl;




    //                                                                                    Alpha_Class[33]=Alpha_Class[33]+ Dalpha_Class[33];
    //                                                                                }


    //                                                                                Alpha_Class[30]=Alpha_Class[30]+ Dalpha_Class[30];
    //                                                                            }
    //                                                                            Alpha_Class[29]=Alpha_Class[29]+ Dalpha_Class[29];
    //                                                                        }
    //                                                                        Alpha_Class[28]=Alpha_Class[28]+ Dalpha_Class[28];
    //                                                                    }

    //                                                                    Alpha_Class[26]=Alpha_Class[26]+ Dalpha_Class[26];
    //                                                                }

    //                                                                Alpha_Class[24]=Alpha_Class[24]+ Dalpha_Class[24];
    //                                                            }
    //                                                            Alpha_Class[23]=Alpha_Class[23]+ Dalpha_Class[23];
    //                                                        }

    //                                                        Alpha_Class[21]=Alpha_Class[21]+ Dalpha_Class[21];
    //                                                    }
    //                                                    Alpha_Class[20]=Alpha_Class[20]+ Dalpha_Class[20];
    //                                                }

    //                                                Alpha_Class[18]=Alpha_Class[18]+ Dalpha_Class[18];
    //                                            }
    //                                            Alpha_Class[17]=Alpha_Class[17]+ Dalpha_Class[17];
    //                                        }
    //                                        Alpha_Class[16]=Alpha_Class[16]+ Dalpha_Class[16];
    //                                    }



    //                                    // Alpha_Class[12]=Alpha_Class[12]+ Dalpha_Class[12];
    //                                    // }


    //                                    //Alpha_Class[9]=Alpha_Class[9]+ Dalpha_Class[9];
    //                                    // }
    //                                    Alpha_Class[8]=Alpha_Class[8]+ Dalpha_Class[8];
    //                                }

    //                                //  Alpha_Class[6]=Alpha_Class[6]+ Dalpha_Class[6];
    //                                // }
    //                                Alpha_Class[5]=Alpha_Class[5]+ Dalpha_Class[5];
    //                            }
    //                            // Alpha_Class[4]=Alpha_Class[4]+ Dalpha_Class[4];
    //                            //}
    //                            // Alpha_Class[3]=Alpha_Class[3]+ Dalpha_Class[3];
    //                            // }
    //                            //Alpha_Class[2]=Alpha_Class[2]+ Dalpha_Class[2];
    //                            //}
    //                            Alpha_Class[1]=Alpha_Class[1]+ Dalpha_Class[1];
    //                        }
    //                        Alpha_Class[0]=Alpha_Class[0]+ Dalpha_Class[0];
    //                    }


    //                }
    //            }

    //        }


    //    }

    //    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;







    //    overlap_=0.0;
    //    for(int i=0;i<432;i++){
    //        overlap_ += Ansatz_2holes_Basis_Overlap_with_GS[i]*Alpha_Final[i];
    //    }
    //    norm_=0.0;
    //    for(int i=0;i<432;i++){
    //        for(int j=0;j<432;j++){
    //            norm_ += OVERLAP_MATRIX_2holes_Basis[i][j]*Alpha_Final[i]*Alpha_Final[j];
    //        }
    //    }


    //    cout<<"overlap/sqrt(norm) = "<< overlap_/sqrt(norm_)<<endl;
    //    cout <<"norm = "<<norm_<<endl;


    //    cout<<"OPTIMIZED VALUES FOR NORM=1:"<<endl;
    //    cout<<"a_CS1="<<a_CS1_final<<"  "<<"a_CS2="<<a_CS2_final<<"  "<<
    //          "a_CS3="<<a_CS3_final<<"  "<<"a_CS4="<<a_CS4_final<<endl;

    //    double temp_alpha;
    //    for(int i=0;i<36;i++){

    //        cout<<"["<<i<<"]";
    //        if(i>=0 && i<16){
    //            temp_alpha= a_CS1_final*Alpha_Final_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[(i*12)]);
    //        }
    //        if(i>=16 && i<20){
    //            temp_alpha= a_CS2_final*Alpha_Final_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[(i*12)]);
    //        }
    //        if(i>=20 && i<28){
    //            temp_alpha= a_CS3_final*Alpha_Final_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[(i*12)]);
    //        }
    //        if(i>=28 && i<36){
    //            temp_alpha= a_CS4_final*Alpha_Final_Class[i]*sign_of_double(Ansatz_2holes_Basis_Overlap_with_GS[(i*12)]);
    //        }
    //        cout<<temp_alpha/sqrt(norm_)<<"   ";


    //    }



    //    OPT_VEC.clear();
    //    OPT_VEC.resize(basis_nm2.D_dn_basis.size()*basis_nm2.D_up_basis.size());
    //    for(int i=0;i<432;i++){
    //        Sum(OPT_VEC, 1.0,BASIS_STATES_ANSATZ_2HOLES[i],Alpha_Final[i] );
    //    }

    //    double norm_temp = dot_product(OPT_VEC,OPT_VEC);
    //    for(int i =0;i<OPT_VEC.size();i++){
    //        OPT_VEC[i]=OPT_VEC[i]/sqrt(norm_temp);
    //    }




}

void MODEL_2_orb_Hubb_chain::Optimize_overlap_Nm2_Variational_State(BASIS_2_orb_Hubb_chain &basis,
                                                                    BASIS_2_orb_Hubb_chain &basis_nm2,
                                                                    Mat_2_doub &STATES_OS_TS_,
                                                                    Mat_1_doub &Eig_vec,
                                                                    Mat_1_doub &vec_final){




    //IMPROVING sum_{i}{i=0..7} a^{i}_{0}*Holes(dis=1)[Delta_{D}|V_{i}(N)>] + a^{i}_{1}*Holes(dis=sqrt(5))[Delta_{D}|V_{i}(N)>] + beta{i}[Delta_{D}|V_{i}(N)>]----//

    //    Mat_1_doub vec_total_updated;
    //    Mat_1_doub vec_temp;

    //    Mat_1_doub alpha_used;
    //    alpha_used.resize(2);


    //    vec_temp.clear();
    //    vec_total_updated.clear();

    //    double  VALUE_CHECK_OLD, VALUE_CHECK;

    //    double norm_updated_, overlap_;

    //    double dalpha_=0.5;
    //    double dbeta_=100;
    //    double alpha_min=(1.0/3.0) - 1.0;
    //    double alpha_max=(1.0/3.0) + 1.0;

    //    Mat_1_doub beta_min;
    //    Mat_1_doub beta_max;

    //    beta_min.resize(STATES_OS_TS_DIAGONAL_HOLES.size());
    //    beta_max.resize(STATES_OS_TS_DIAGONAL_HOLES.size());

    //    /*
    //overlap of class 0 with GS = -3.13043859422097336953e-01, No of states = 2
    //overlap of class 1 with GS = 3.25409978887851336626e-01, No of states = 12
    //overlap of class 2 with GS = -3.39092014449482725968e-01, No of states = 12
    //overlap of class 3 with GS = 3.54256365525317862275e-01, No of states = 2
    //overlap of class 4 with GS = 3.24857861824018656627e-01, No of states = 6
    //overlap of class 5 with GS = -3.38179544421070676918e-01, No of states = 6
    //overlap of class 6 with GS = 3.38299670520067541357e-01, No of states = 12
    //overlap of class 7 with GS = -3.24811522677815101012e-01, No of states = 12
    //     */

    //    if(basis.Length==6){


    //        beta_min[0]=-3.13043859422097336953e-01  ; beta_max[0]=-3.13043859422097336953e-01  ;
    //        beta_min[1]=3.25409978887851336626e-01  ; beta_max[1]=3.25409978887851336626e-01  ;
    //        beta_min[2]=-3.39092014449482725968e-01  ; beta_max[2]= -3.39092014449482725968e-01 ;
    //        beta_min[3]=3.54256365525317862275e-01 ; beta_max[3]= 3.54256365525317862275e-01 ;
    //        beta_min[4]=3.24857861824018656627e-01  ; beta_max[4]= 3.24857861824018656627e-01 ;
    //        beta_min[5]= -3.38179544421070676918e-01 ; beta_max[5]= -3.38179544421070676918e-01 ;
    //        beta_min[6]=3.38299670520067541357e-01  ; beta_max[6]= 3.38299670520067541357e-01 ;
    //        beta_min[7]=-3.24811522677815101012e-01  ; beta_max[7]=-3.24811522677815101012e-01  ;


    //        /* for(int i=0;i<8;i++){
    //        beta_min[i]=-1.0;beta_max[i]=1.0;
    //        }*/
    //        /*
    //        beta_min[1]=0;beta_max[1]=0;
    //        beta_min[2]=0;beta_max[2]=0;
    //        beta_min[3]=0;beta_max[3]=0;
    //        beta_min[4]=0;beta_max[4]=0;
    //        beta_min[5]=0;beta_max[5]=0;
    //        beta_min[6]=0;beta_max[6]=0;
    //        beta_min[7]=0;beta_max[7]=0;
    //        */


    //    }



    //    Mat_2_doub Alpha_;
    //    Alpha_.resize(8);
    //    for(int i=0;i<8;i++){
    //        Alpha_[i].resize(2);
    //    }

    //    Mat_2_doub _Alpha_final_;
    //    _Alpha_final_.resize(8);
    //    for(int i=0;i<8;i++){
    //        _Alpha_final_[i].resize(2);
    //    }

    //    Mat_1_doub Beta_;
    //    Beta_.resize(8);


    //    Mat_1_doub Beta_final;
    //    Beta_final.resize(8);


    //    VALUE_CHECK_OLD=0.0;


    //    Alpha_[0][0]=alpha_min;
    //    while(Alpha_[0][0]<=alpha_max){

    //        Alpha_[1][0]=alpha_min;
    //        while(Alpha_[1][0]<=alpha_max){

    //            Alpha_[2][0]=alpha_min;
    //            while(Alpha_[2][0]<=alpha_max){

    //                Alpha_[3][0]=alpha_min;
    //                while(Alpha_[3][0]<=alpha_max){

    //                    Alpha_[4][0]=alpha_min;
    //                    while(Alpha_[4][0]<=alpha_max){

    //                        Alpha_[5][0]=alpha_min;
    //                        while(Alpha_[5][0]<=alpha_max){

    //                            Alpha_[6][0]=alpha_min;
    //                            while(Alpha_[6][0]<=alpha_max){

    //                                Alpha_[7][0]=alpha_min;
    //                                while(Alpha_[7][0]<=alpha_max){

    //                                    Alpha_[0][1]=alpha_min;
    //                                    while(Alpha_[0][1]<=alpha_max){

    //                                        Alpha_[1][1]=alpha_min;
    //                                        while(Alpha_[1][1]<=alpha_max){

    //                                            Alpha_[2][1]=alpha_min;
    //                                            while(Alpha_[2][1]<=alpha_max){

    //                                                Alpha_[3][1]=alpha_min;
    //                                                while(Alpha_[3][1]<=alpha_max){

    //                                                    Alpha_[4][1]=alpha_min;
    //                                                    while(Alpha_[4][1]<=alpha_max){

    //                                                        Alpha_[5][1]=alpha_min;
    //                                                        while(Alpha_[5][1]<=alpha_max){

    //                                                            Alpha_[6][1]=alpha_min;
    //                                                            while(Alpha_[6][1]<=alpha_max){

    //                                                                Alpha_[7][1]=alpha_min;
    //                                                                while(Alpha_[7][1]<=alpha_max){

    //                                                                    Beta_[0]=beta_min[0];
    //                                                                    while(Beta_[0]<=beta_max[0]){

    //                                                                        Beta_[1]=beta_min[1];
    //                                                                        while(Beta_[1]<=beta_max[1]){

    //                                                                            Beta_[2]=beta_min[2];
    //                                                                            while(Beta_[2]<=beta_max[2]){

    //                                                                                Beta_[3]=beta_min[3];
    //                                                                                while(Beta_[3]<=beta_max[3]){

    //                                                                                    Beta_[4]=beta_min[4];
    //                                                                                    while(Beta_[4]<=beta_max[4]){

    //                                                                                        Beta_[5]=beta_min[5];
    //                                                                                        while(Beta_[5]<=beta_max[5]){

    //                                                                                            Beta_[6]=beta_min[6];
    //                                                                                            while(Beta_[6]<=beta_max[6]){

    //                                                                                                Beta_[7]=beta_min[7];
    //                                                                                                while(Beta_[7]<=beta_max[7]){

    //                                                                                                    if( true

    //                                                                                                            //(Alpha_[3]==Alpha_[1] )
    //                                                                                                            //&&
    //                                                                                                            //(Alpha_[2]==Alpha_[0] )

    //                                                                                                            ){

    //                                                                                                        vec_total_updated.clear();
    //                                                                                                        vec_total_updated.resize(basis_nm2.D_dn_basis.size()*basis_nm2.D_up_basis.size());

    //                                                                                                        for(int ci=0;ci<STATES_OS_TS_DIAGONAL_HOLES.size();ci++){

    //                                                                                                            // Alpha_[ci][0]=Beta_[ci];
    //                                                                                                            // Alpha_[ci][1]=Beta_[ci];
    //                                                                                                            alpha_used[0]=Alpha_[ci][0]*Beta_[ci];
    //                                                                                                            alpha_used[1]=Alpha_[ci][1]*Beta_[ci];

    //                                                                                                            vec_temp=STATES_OS_TS_DIAGONAL_HOLES[ci];
    //                                                                                                            Create_states_with_hole_hopping(basis_nm2, alpha_used, vec_temp);
    //                                                                                                            Sum(vec_total_updated, 1.0, vec_temp, 1.0);

    //                                                                                                            vec_temp=STATES_OS_TS_DIAGONAL_HOLES[ci];
    //                                                                                                            Sum(vec_total_updated, 1.0, vec_temp, Beta_[ci]);


    //                                                                                                        }





    //                                                                                                        norm_updated_=dot_product(vec_total_updated, vec_total_updated);
    //                                                                                                        // cout <<"Norm(H[Delta|Ansatz(N)>])  = "<<norm_updated_<<endl;
    //                                                                                                        overlap_ = dot_product(Eig_vec, vec_total_updated);
    //                                                                                                        // cout <<"<GS(N-2)|H [Delta|Ansatz(N)>]  = "<<(overlap_)<<endl;
    //                                                                                                        VALUE_CHECK = (overlap_*overlap_)/(norm_updated_);

    //                                                                                                        if( VALUE_CHECK > VALUE_CHECK_OLD){
    //                                                                                                            VALUE_CHECK_OLD=VALUE_CHECK;
    //                                                                                                            _Alpha_final_=Alpha_;
    //                                                                                                            Beta_final = Beta_;
    //                                                                                                        }

    //                                                                                                        for(int i=0;i<Alpha_.size();i++){
    //                                                                                                            for(int j=0;j<2;j++){
    //                                                                                                                cout<<"Alpha_["<<i<<"]["<<j<<"] = "<<Alpha_[i][j]<<", ";
    //                                                                                                            }
    //                                                                                                        }
    //                                                                                                        for(int i=0;i<Beta_.size();i++){
    //                                                                                                            cout<<"Beta_["<<i<<"] = "<<Beta_[i]<<", ";
    //                                                                                                        }
    //                                                                                                        cout<<endl;
    //                                                                                                        cout<<", value = "<<VALUE_CHECK<<endl;



    //                                                                                                    }

    //                                                                                                    Beta_[7] +=dbeta_;
    //                                                                                                }
    //                                                                                                Beta_[6] +=dbeta_;
    //                                                                                            }
    //                                                                                            Beta_[5] +=dbeta_;
    //                                                                                        }
    //                                                                                        Beta_[4] +=dbeta_;
    //                                                                                    }
    //                                                                                    Beta_[3] +=dbeta_;
    //                                                                                }
    //                                                                                Beta_[2] +=dbeta_;
    //                                                                            }
    //                                                                            Beta_[1] +=dbeta_;
    //                                                                        }
    //                                                                        Beta_[0] +=dbeta_;
    //                                                                    }
    //                                                                    Alpha_[7][1] +=dalpha_;
    //                                                                }
    //                                                                Alpha_[6][1] +=dalpha_;
    //                                                            }
    //                                                            Alpha_[5][1] +=dalpha_;
    //                                                        }
    //                                                        Alpha_[4][1] +=dalpha_;
    //                                                    }
    //                                                    Alpha_[3][1] +=dalpha_;
    //                                                }
    //                                                Alpha_[2][1] +=dalpha_;
    //                                            }
    //                                            Alpha_[1][1] +=dalpha_;
    //                                        }
    //                                        Alpha_[0][1] +=dalpha_;
    //                                    }
    //                                    Alpha_[7][0] +=dalpha_;
    //                                }
    //                                Alpha_[6][0] +=dalpha_;
    //                            }
    //                            Alpha_[5][0] +=dalpha_;
    //                        }
    //                        Alpha_[4][0] +=dalpha_;
    //                    }
    //                    Alpha_[3][0] +=dalpha_;
    //                }
    //                Alpha_[2][0] +=dalpha_;
    //            }
    //            Alpha_[1][0] +=dalpha_;
    //        }
    //        Alpha_[0][0] +=dalpha_;
    //    }


    //    for(int i=0;i<_Alpha_final_.size();i++){
    //        for(int j=0;j<2;j++){
    //            cout<<"Alpha_final["<<i<<"]["<<j<<"] = "<<_Alpha_final_[i][j]<<", ";
    //        }
    //    }
    //    for(int i=0;i<Beta_final.size();i++){
    //        cout<<"Beta_final["<<i<<"] = "<<Beta_final[i]<<", ";
    //    }
    //    cout<<endl;

    //    //WRITING FINAL SOLUTION-----------------//
    //    vec_total_updated.clear();
    //    vec_total_updated.resize(basis_nm2.D_dn_basis.size()*basis_nm2.D_up_basis.size());

    //    for(int ci=0;ci<STATES_OS_TS_DIAGONAL_HOLES.size();ci++){

    //        // _Alpha_final_[ci][0]=Beta_final[ci];
    //        // _Alpha_final_[ci][1]=Beta_final[ci];

    //        alpha_used[0]=_Alpha_final_[ci][0]*Beta_final[ci];
    //        alpha_used[1]=_Alpha_final_[ci][1]*Beta_final[ci];

    //        vec_temp=STATES_OS_TS_DIAGONAL_HOLES[ci];
    //        Create_states_with_hole_hopping(basis_nm2, alpha_used, vec_temp);
    //        Sum(vec_total_updated, 1.0, vec_temp, 1.0);

    //        vec_temp=STATES_OS_TS_DIAGONAL_HOLES[ci];
    //        Sum(vec_total_updated, 1.0, vec_temp, Beta_final[ci]);


    //    }

    //    norm_updated_=dot_product(vec_total_updated, vec_total_updated);
    //    cout <<"Norm(|2-hole VAR>)  = "<<norm_updated_<<endl;
    //    overlap_ = dot_product(Eig_vec, vec_total_updated);
    //    cout <<"<GS(N-2)|2-hole VAR>]  = "<<(overlap_)<<endl;


    //    vec_final=vec_total_updated;

}

void MODEL_2_orb_Hubb_chain::Read_Anzatz_basis(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_){

    //    int No_of_basis;
    //    if(basis.Length==2){
    //        No_of_basis=2;
    //    }
    //    if(basis.Length==4){
    //        No_of_basis=16;
    //    }
    //    if(basis.Length==6){
    //        No_of_basis=64;
    //    }

    //    string state_type;
    //    BASIS_STATES_ANSATZ.clear();
    //    BASIS_STATES_ANSATZ.resize(No_of_basis);

    //    Ansatz_Basis_Overlap_with_GS.clear();
    //    Ansatz_Basis_Overlap_with_GS.resize(No_of_basis);

    //    int index;
    //    double value;

    //    ostringstream ss_length;
    //    ss_length << basis.Length;



    //    for(int i=0;i<No_of_basis;i++){
    //        ostringstream ss_int;
    //        ss_int << i;
    //        string sLine;

    //        ifstream infile_Variational_states;

    //        string file_in_Variational_states = "Variational_states_L"+ ss_length.str() +"/Variational_state" + ss_int.str() + ".txt";
    //        infile_Variational_states.open(file_in_Variational_states.c_str());
    //        getline(infile_Variational_states, state_type);


    //        BASIS_STATES_ANSATZ[i].clear();
    //        BASIS_STATES_ANSATZ[i].resize(basis.D_dn_basis.size()*basis.D_up_basis.size());

    //        while (!infile_Variational_states.eof())
    //        {
    //            getline(infile_Variational_states, sLine);
    //            stringstream ssLine;
    //            ssLine << sLine;
    //            ssLine>>index;
    //            ssLine>>value;

    //            BASIS_STATES_ANSATZ[i][index]=value;

    //        }


    //        infile_Variational_states.close();


    //        Ansatz_Basis_Overlap_with_GS[i]=dot_product(GS_,BASIS_STATES_ANSATZ[i]);
    //        cout<<i<<" ,State type = "<<state_type<<endl;
    //        cout<<Ansatz_Basis_Overlap_with_GS[i]<<endl;


    //    }


}


double MODEL_2_orb_Hubb_chain::Get_Holes_Projected_state_probability(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec){

    //    double prob_num, prob_den, prob_conditional;

    //    Mat_1_doub vec1, vec2;
    //    int bit_VALUE_up, bit_VALUE_dn;
    //    int index_no;



    //    int orb, spin, site;


    //    //FOR DENOMINATOR
    //    vec1=vec;
    //    for (int i=0;i<basis.D_up_basis.size();i++){
    //        for (int j=0;j<basis.D_dn_basis.size();j++){
    //            index_no=basis.D_dn_basis.size()*i + j;


    //            orb= Hole_positions[0].orb_;
    //            site= Hole_positions[0].site_;

    //            bit_VALUE_dn = bit_value(basis.D_dn_basis[j],orb*basis.Length + site);
    //            bit_VALUE_up = bit_value(basis.D_up_basis[i],orb*basis.Length + site);


    //            if(
    //                    (bit_VALUE_up==1 || bit_VALUE_dn==1 )
    //                    ){

    //                vec1[index_no] = 0.0;

    //            }

    //        }}



    //    //FOR NUMERATOR
    //    vec2=vec;
    //    for(int n=0;n<Hole_positions.size();n++){
    //        for (int i=0;i<basis.D_up_basis.size();i++){
    //            for (int j=0;j<basis.D_dn_basis.size();j++){
    //                index_no=basis.D_dn_basis.size()*i + j;


    //                orb= Hole_positions[n].orb_;

    //                site= Hole_positions[n].site_;

    //                bit_VALUE_dn = bit_value(basis.D_dn_basis[j],orb*basis.Length + site);

    //                bit_VALUE_up = bit_value(basis.D_up_basis[i],orb*basis.Length + site);


    //                if(
    //                        (bit_VALUE_up==1 || bit_VALUE_dn==1)
    //                        ){

    //                    vec2[index_no] = 0.0;

    //                }

    //            }}
    //    }


    //    prob_num = dot_product(vec, vec2);
    //    prob_den = dot_product(vec, vec1);

    //    prob_conditional = prob_num/prob_den;



    //    return prob_conditional;

}



Mat_1_doub MODEL_2_orb_Hubb_chain::Get_2Holes_Projected_state(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec){


    //    Mat_1_doub vec1;
    //    int bit_VALUE_up_hole_0, bit_VALUE_dn_hole_0;
    //    int bit_VALUE_up_hole_1, bit_VALUE_dn_hole_1;
    //    int index_no;
    //    double norm_;



    //    int orb0, site0;
    //    int orb1, site1;



    //    vec1=vec;
    //    for (int i=0;i<basis.D_up_basis.size();i++){
    //        for (int j=0;j<basis.D_dn_basis.size();j++){
    //            index_no=basis.D_dn_basis.size()*i + j;


    //            orb0= Hole_positions[0].orb_;
    //            site0= Hole_positions[0].site_;
    //            orb1= Hole_positions[1].orb_;
    //            site1= Hole_positions[1].site_;

    //            bit_VALUE_dn_hole_0 = bit_value(basis.D_dn_basis[j],orb0*basis.Length + site0);
    //            bit_VALUE_up_hole_0 = bit_value(basis.D_up_basis[i],orb0*basis.Length + site0);

    //            bit_VALUE_dn_hole_1 = bit_value(basis.D_dn_basis[j],orb1*basis.Length + site1);
    //            bit_VALUE_up_hole_1 = bit_value(basis.D_up_basis[i],orb1*basis.Length + site1);


    //            if(
    //                    ( (bit_VALUE_up_hole_0==1 || bit_VALUE_dn_hole_0==1 )
    //                      ||
    //                      (bit_VALUE_up_hole_1==1 || bit_VALUE_dn_hole_1==1 )
    //                      )
    //                    ){

    //                vec1[index_no] = 0.0;

    //            }

    //        }}


    //    norm_=dot_product(vec1,vec1);

    //    for(int i=0;i<vec1.size();i++){
    //        vec1[i] = vec1[i]/sqrt(norm_);
    //    }

    //    return vec1;
}


void MODEL_2_orb_Hubb_chain::Get_Holes_Projected_state(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec){

    //    Mat_1_doub vec1;
    //    vec1=vec;

    //    int bit_VALUE_up, bit_VALUE_dn;
    //    int index_no;

    //    int orb, site;
    //    double norm_;


    //    for (int i=0;i<basis.D_up_basis.size();i++){
    //        for (int j=0;j<basis.D_dn_basis.size();j++){
    //            index_no=basis.D_dn_basis.size()*i + j;


    //            orb= Hole_positions[0].orb_;
    //            site= Hole_positions[0].site_;

    //            bit_VALUE_dn = bit_value(basis.D_dn_basis[j],orb*basis.Length + site);
    //            bit_VALUE_up = bit_value(basis.D_up_basis[i],orb*basis.Length + site);


    //            if(
    //                    (bit_VALUE_up==1 || bit_VALUE_dn==1 )
    //                    ){

    //                vec1[index_no] = 0.0;

    //            }

    //        }}


    //    //Normalizing vec1;
    //    norm_ = dot_product(vec1,vec1);

    //    for(int i=0;i<vec1.size();i++){
    //        vec1[i] = vec1[i]/sqrt(norm_);
    //    }

    //    vec = vec1;

}

void MODEL_2_orb_Hubb_chain::Read_Mat_2_trio(Mat_2_trio_int &MAT_TEMP, Mat_1_doub &VALUES_TEMP,
                                             int pair_no){
    //Already having superpositions; see notes
    MAT_TEMP.resize(4);
    for(int i=0;i<4;i++){
        MAT_TEMP[i].resize(2);
    }
    VALUES_TEMP.resize(4);
    double VALUE_TYPE_1, VALUE_TYPE_2;


    if(Variational_state_pair_orbital_symmetry[pair_no]=="NOS1"){
        VALUE_TYPE_1 = 1.0;
        VALUE_TYPE_2 = 0.0;
    }
    else if(Variational_state_pair_orbital_symmetry[pair_no]=="NOS2"){
        VALUE_TYPE_1 = 0.0;
        VALUE_TYPE_2 = 1.0;
    }
    else if(Variational_state_pair_orbital_symmetry[pair_no]=="OS"){
        VALUE_TYPE_1 = 1.0;
        VALUE_TYPE_2 = 1.0;
    }
    else if(Variational_state_pair_orbital_symmetry[pair_no]=="OA"){
        VALUE_TYPE_1 = 1.0;
        VALUE_TYPE_2 = -1.0;
    }


    if(Variational_state_pair_Geometry[pair_no]=="D"
            ){

        //1.0*c^{\dagger}_{a,i,up}  c^{\dagger}_{b,i+1,dn}*VALUE_TYPE_1
        MAT_TEMP[0][0].orb_=0;MAT_TEMP[0][0].spin_=0;
        MAT_TEMP[0][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[0][1].orb_=1;MAT_TEMP[0][1].spin_=1;
        MAT_TEMP[0][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[0]=1.0*VALUE_TYPE_1;

        // +\-1.0*c^{\dagger}_{a,i,dn}  c^{\dagger}_{b,i+1,up}*VALUE_TYPE_1
        MAT_TEMP[1][0].orb_=0;MAT_TEMP[1][0].spin_=1;
        MAT_TEMP[1][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[1][1].orb_=1;MAT_TEMP[1][1].spin_=0;
        MAT_TEMP[1][1].site_=Variational_state_pair_sites[pair_no].second;
        if(Variational_state_pair_spin_symmetry[pair_no]=="S"){
            VALUES_TEMP[1]=-1.0*VALUE_TYPE_1;
        }
        else{
            assert(Variational_state_pair_spin_symmetry[pair_no]=="T");
            VALUES_TEMP[1]=1.0*VALUE_TYPE_1;
        }

        //For next two terms

        //value_oa*c^{\dagger}_{b,i,up}  c^{\dagger}_{a,i+1,dn}
        MAT_TEMP[2][0].orb_=1;MAT_TEMP[2][0].spin_=0;
        MAT_TEMP[2][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[2][1].orb_=0;MAT_TEMP[2][1].spin_=1;
        MAT_TEMP[2][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[2]=1.0*VALUE_TYPE_2;

        // +\-1.0*value_oa*c^{\dagger}_{b,i,dn}  c^{\dagger}_{a,i+1,up}
        MAT_TEMP[3][0].orb_=1;MAT_TEMP[3][0].spin_=1;
        MAT_TEMP[3][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[3][1].orb_=0;MAT_TEMP[3][1].spin_=0;
        MAT_TEMP[3][1].site_=Variational_state_pair_sites[pair_no].second;
        if(Variational_state_pair_spin_symmetry[pair_no]=="S"){
            VALUES_TEMP[3]=-1.0*VALUE_TYPE_2;
        }
        else{
            assert(Variational_state_pair_spin_symmetry[pair_no]=="T");
            VALUES_TEMP[3]=1.0*VALUE_TYPE_2;
        }

    }
    else{

        assert(Variational_state_pair_Geometry[pair_no]=="AC");

        //1.0*c^{\dagger}_{a,i,up}  c^{\dagger}_{a,i+1,dn}
        MAT_TEMP[0][0].orb_=0;MAT_TEMP[0][0].spin_=0;
        MAT_TEMP[0][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[0][1].orb_=0;MAT_TEMP[0][1].spin_=1;
        MAT_TEMP[0][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[0]=1.0*VALUE_TYPE_1;

        // +\-1.0*c^{\dagger}_{a,i,dn}  c^{\dagger}_{a,i+1,up}
        MAT_TEMP[1][0].orb_=0;MAT_TEMP[1][0].spin_=1;
        MAT_TEMP[1][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[1][1].orb_=0;MAT_TEMP[1][1].spin_=0;
        MAT_TEMP[1][1].site_=Variational_state_pair_sites[pair_no].second;
        if(Variational_state_pair_spin_symmetry[pair_no]=="S"){
            VALUES_TEMP[1]=-1.0*VALUE_TYPE_1;
        }
        else{
            assert(Variational_state_pair_spin_symmetry[pair_no]=="T");
            VALUES_TEMP[1]=1.0*VALUE_TYPE_1;
        }

        //For next two terms

        //value_oa*c^{\dagger}_{b,i,up}  c^{\dagger}_{b,i+1,dn}
        MAT_TEMP[2][0].orb_=1;MAT_TEMP[2][0].spin_=0;
        MAT_TEMP[2][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[2][1].orb_=1;MAT_TEMP[2][1].spin_=1;
        MAT_TEMP[2][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[2]=1.0*VALUE_TYPE_2;

        // +\-1.0*value_oa*c^{\dagger}_{b,i,dn}  c^{\dagger}_{b,i+1,up}
        MAT_TEMP[3][0].orb_=1;MAT_TEMP[3][0].spin_=1;
        MAT_TEMP[3][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[3][1].orb_=1;MAT_TEMP[3][1].spin_=0;
        MAT_TEMP[3][1].site_=Variational_state_pair_sites[pair_no].second;
        if(Variational_state_pair_spin_symmetry[pair_no]=="S"){
            VALUES_TEMP[3]=-1.0*VALUE_TYPE_2;
        }
        else{
            assert(Variational_state_pair_spin_symmetry[pair_no]=="T");
            VALUES_TEMP[3]=1.0*VALUE_TYPE_2;
        }

    }

    //---------------------------------------//
    //--------------------------------------//



}

void MODEL_2_orb_Hubb_chain::Get_Variational_State(BASIS_2_orb_Hubb_chain &basis, int no_of_pairs){


    //    Mat_2_trio_int MAT_FINAL;
    //    Mat_1_doub VALUES_FINAL;

    //    Mat_2_trio_int MAT_TEMP;
    //    Mat_1_doub VALUES_TEMP;

    //    Read_Mat_2_trio(MAT_FINAL, VALUES_FINAL,0);

    //    for(int i=1;i<no_of_pairs;i++){
    //        Read_Mat_2_trio(MAT_TEMP, VALUES_TEMP,i);

    //        Direct_product_of_Mat_2_trio_int(MAT_FINAL, VALUES_FINAL,
    //                                         MAT_TEMP, VALUES_TEMP,
    //                                         MAT_FINAL, VALUES_FINAL);
    //    }


    //    Mat_1_doub State_temp;
    //    State_.resize(basis.D_dn_basis.size()*basis.D_up_basis.size());
    //    State_temp.resize(basis.D_dn_basis.size()*basis.D_up_basis.size());
    //    for(int i=0;i<State_.size();i++){
    //        State_[i]=0.0;
    //    }


    //    for(int i=0;i<MAT_FINAL.size();i++){

    //        if(VALUES_FINAL[i]!=0.0){
    //            Get_State_by_acting_creation_operators(basis, MAT_FINAL[i],State_temp);

    //            //Add VALUE_TEMP[i]*State_temp in State_
    //            Subtract( State_, -1.0*VALUES_FINAL[i], State_temp);
    //        }

    //    }

    //    //Normalize State_;
    //    double tmpnrm_type_double,tmpnrm;
    //    tmpnrm_type_double=dot_product(State_,State_);
    //    tmpnrm=sqrt(tmpnrm_type_double);
    //    if(tmpnrm !=0){
    //        for(int i=0;i<State_.size();i++){
    //            State_[i] = (State_[i]/(tmpnrm));
    //        }
    //    }


}


void MODEL_2_orb_Hubb_chain::Get_State_by_acting_creation_operators(BASIS_2_orb_Hubb_chain &basis,
                                                                    Mat_1_trio_int MAT_,
                                                                    Mat_1_doub &State_temp){


    //    double FM_sign_up_down, FM_sign_up, FM_sign_down;
    //    int up_down_int, up_int, down_int;
    //    int offset;

    //    int Nup_check, Ndown_check;
    //    int Dup_temp, Ddn_temp;

    //    bool up_repeat, down_repeat;

    //    State_temp.resize(basis.D_dn_basis.size()*basis.D_up_basis.size());
    //    for(int i=0;i<State_temp.size();i++){
    //        State_temp[i]=0.0;
    //    }
    //    Mat_1_int positions_up;
    //    Mat_1_int positions_down;


    //    positions_up.clear();
    //    positions_down.clear();

    //    up_down_int=0;
    //    offset=0;
    //    Nup_check=0;
    //    Ndown_check=0;
    //    for(int i=0;i<MAT_.size();i++){

    //        //up
    //        if(MAT_[i].spin_==0){
    //            positions_up.push_back( MAT_[i].site_ + (basis.Length*MAT_[i].orb_) );
    //            offset++;
    //            Nup_check++;
    //        }
    //        else{
    //            //down
    //            assert(MAT_[i].spin_==1);
    //            positions_down.push_back( MAT_[i].site_ + (basis.Length*MAT_[i].orb_) );

    //            up_down_int += offset;
    //            Ndown_check++;
    //        }


    //    }


    //    up_repeat=false;
    //    for (int i = 0; i < positions_up.size() - 1; i++){
    //        for (int j = i + 1;j < positions_up.size(); j++){
    //            if (positions_up[i] == positions_up[j]){
    //                up_repeat=true;
    //            } // then this is a duplicate
    //        }
    //    }

    //    down_repeat=false;
    //    for (int i = 0; i < positions_down.size() - 1; i++){
    //        for (int j = i + 1;j < positions_down.size(); j++){
    //            if (positions_down[i] == positions_down[j]){
    //                down_repeat=true;
    //            } // then this is a duplicate
    //        }
    //    }



    //    if(up_repeat==false &&
    //            down_repeat==false){
    //        assert(Nup_check == basis.Nup);
    //        assert(Ndown_check == basis.Ndn);

    //        FM_sign_up_down = pow(-1.0, 1.0*up_down_int);


    //        up_int=minSwaps(positions_up, positions_up.size());
    //        FM_sign_up=pow(-1.0, 1.0*up_int);

    //        down_int=minSwaps(positions_down, positions_down.size());
    //        FM_sign_down=pow(-1.0, 1.0*down_int);

    //        Dup_temp=0;
    //        for(int i=0;i<positions_up.size();i++){
    //            Dup_temp +=pow(2,positions_up[i]);
    //        }

    //        Ddn_temp=0;
    //        for(int i=0;i<positions_down.size();i++){
    //            Ddn_temp +=pow(2,positions_down[i]);
    //        }


    //        int i_new, j_new, m_new;

    //        i_new = Find_int_in_intarray(Dup_temp,basis.D_up_basis);
    //        j_new = Find_int_in_intarray(Ddn_temp,basis.D_dn_basis);

    //        m_new = basis.D_dn_basis.size()*i_new + j_new;


    //        State_temp[m_new]=FM_sign_up_down*FM_sign_up*FM_sign_down;

    //    }



}


void MODEL_2_orb_Hubb_chain::Get_c_on_GS(Mat_1_doub & EigVec_, BASIS_2_orb_Hubb_chain & basis_Nm1, BASIS_2_orb_Hubb_chain & basis,
                                         Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


    int site_val ;
    int orb_val ;
    int spin_val ;
    int D_dn_new, D_up_new;
    int m_new, m, i_up,i_dn;
    int max_up, max_dn, min_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    int l;
    double sign_FM;
    double_type value_in, value;

    State_c_on_GS.clear();
    State_c_on_GS.resize(basis_Nm1.D_up_basis.size()*basis_Nm1.D_dn_basis.size());

    assert(TRIO_VEC.size() == values.size());

    for(int n=0;n<TRIO_VEC.size();n++){
        site_val = TRIO_VEC[n].site_;
        orb_val = TRIO_VEC[n].orb_;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        if(spin_val==0){
            //For c_up|GS>
            assert( (basis_Nm1.Ndn==basis.Ndn) &
                    (basis_Nm1.Nup==basis.Nup-1)
                    );
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++)
                {
                    m = basis.D_dn_basis.size()*i + j;


                    if(bit_value(basis.D_up_basis[i], orb_val*basis.Length + site_val)==1){
                        l = orb_val*basis.Length + site_val;

                        D_up_new = (int) (basis.D_up_basis[i] - pow(2,orb_val*basis.Length + site_val) );
                        D_dn_new = basis.D_dn_basis[j];


                        i_up = Find_int_in_intarray(D_up_new,basis_Nm1.D_up_basis);
                        i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.D_dn_basis);
                        //i_up = Find_int_in_intarray_smartly(D_up_new,basis_Nm1.D_up_basis,basis_Nm1.partitions_up,basis_Nm1.Dup_val_at_partitions);
                        //i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Nm1.D_dn_basis,basis_Nm1.partitions_dn,basis_Nm1.Ddn_val_at_partitions);

                        m_new = basis_Nm1.D_dn_basis.size()*i_up + i_dn;

                        max_up = 2*basis.Length -1;
                        min_up =0 ;
                        sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                        if(l != min_up){
                            sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                        }

                        sign_FM = pow(-1.0, sign_pow_up);

                        value = sign_FM*EigVec_[m]*value_in;

                        State_c_on_GS[m_new] += value;
                    }
                }
            }
        }


        //For c_dn|GS>
        if(spin_val==1){
            assert( (basis_Nm1.Ndn==basis.Ndn-1) &
                    (basis_Nm1.Nup==basis.Nup)
                    );
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++)
                {
                    m = basis.D_dn_basis.size()*i + j;

                    if(bit_value(basis.D_dn_basis[j], orb_val*basis.Length + site_val)==1){
                        l = orb_val*basis.Length + site_val;

                        D_dn_new = (int) (basis.D_dn_basis[j] - pow(2, orb_val*basis.Length + site_val));
                        D_up_new = basis.D_up_basis[i];

                        i_up = Find_int_in_intarray(D_up_new,basis_Nm1.D_up_basis);
                        i_dn = Find_int_in_intarray(D_dn_new,basis_Nm1.D_dn_basis);
                        //  i_up = Find_int_in_intarray_smartly(D_up_new,basis_Nm1.D_up_basis,basis_Nm1.partitions_up,basis_Nm1.Dup_val_at_partitions);
                        //  i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Nm1.D_dn_basis,basis_Nm1.partitions_dn,basis_Nm1.Ddn_val_at_partitions);

                        m_new = basis_Nm1.D_dn_basis.size()*i_up + i_dn;

                        max_dn = 2*basis.Length -1;
                        min_dn=0;
                        sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[j]) ;
                        if(l != min_dn){
                            sign_pow_dn += bit_value(basis.D_dn_basis[j],min_dn);
                        }
                        sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                        sign_FM = pow(-1.0, sign_pow_dn);

                        value = sign_FM*EigVec_[m]*value_in;

                        State_c_on_GS[m_new] += value;

                    }
                }
            }
        }


    }

}



void MODEL_2_orb_Hubb_chain::Get_cdagger_on_GS(Mat_1_doub & EigVec_, BASIS_2_orb_Hubb_chain & basis_Np1, BASIS_2_orb_Hubb_chain & basis,
                                               Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


    int site_val ;
    int orb_val ;
    int spin_val ;
    int D_dn_new, D_up_new;
    int m_new,m, i_up,i_dn;
    int max_up,max_dn, min_up, min_dn;
    int sign_pow_up, sign_pow_dn;
    int l;
    double sign_FM;
    double_type value_in, value;

    State_cdagger_on_GS.clear();
    State_cdagger_on_GS.resize(basis_Np1.D_up_basis.size()*basis_Np1.D_dn_basis.size());

    assert(TRIO_VEC.size() == values.size());

    for(int n=0;n<TRIO_VEC.size();n++){
        site_val = TRIO_VEC[n].site_;
        orb_val = TRIO_VEC[n].orb_;
        spin_val = TRIO_VEC[n].spin_;
        value_in = values[n];

        //For c_dagger_up|GS>
        if(spin_val==0){
            assert( (basis_Np1.Ndn==basis.Ndn) &
                    (basis_Np1.Nup==basis.Nup+1)
                    );
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++)
                {
                    m = basis.D_dn_basis.size()*i + j;


                    if(bit_value(basis.D_up_basis[i],orb_val*basis.Length +site_val)==0){
                        l = orb_val*basis.Length +site_val;

                        D_up_new = (int) (basis.D_up_basis[i] + pow(2, orb_val*basis.Length +site_val) );
                        D_dn_new = basis.D_dn_basis[j];

                        i_up = Find_int_in_intarray(D_up_new,basis_Np1.D_up_basis);
                        i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.D_dn_basis);
                        // i_up = Find_int_in_intarray_smartly(D_up_new,basis_Np1.D_up_basis,basis_Np1.partitions_up,basis_Np1.Dup_val_at_partitions);
                        // i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Np1.D_dn_basis,basis_Np1.partitions_dn,basis_Np1.Ddn_val_at_partitions);

                        m_new = basis_Np1.D_dn_basis.size()*i_up + i_dn;

                        max_up = 2*basis.Length -1;
                        min_up=0;
                        sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                        if(l != min_up){
                            sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                        }

                        sign_FM = pow(-1.0, sign_pow_up);

                        value = sign_FM*EigVec_[m]*(conjugate(value_in));

                        State_cdagger_on_GS[m_new] += value;

                    }
                }
            }
        }

        //For c_dagger_dn|GS>
        if(spin_val==1){
            assert( (basis_Np1.Ndn==basis.Ndn+1) &
                    (basis_Np1.Nup==basis.Nup)
                    );
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++)
                {
                    m = basis.D_dn_basis.size()*i + j;

                    if(bit_value(basis.D_dn_basis[j],orb_val*basis.Length +site_val)==0){
                        l = orb_val*basis.Length +site_val;

                        D_dn_new = (int) (basis.D_dn_basis[j] + pow(2, orb_val*basis.Length +site_val) );
                        D_up_new = basis.D_up_basis[i];

                        i_up = Find_int_in_intarray(D_up_new,basis_Np1.D_up_basis);
                        i_dn = Find_int_in_intarray(D_dn_new,basis_Np1.D_dn_basis);
                        //i_up = Find_int_in_intarray_smartly(D_up_new,basis_Np1.D_up_basis,basis_Np1.partitions_up,basis_Np1.Dup_val_at_partitions);
                        //i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Np1.D_dn_basis,basis_Np1.partitions_dn,basis_Np1.Ddn_val_at_partitions);

                        m_new = basis_Np1.D_dn_basis.size()*i_up + i_dn;

                        max_dn = 2*basis.Length -1;
                        min_dn=0;
                        sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[j]) ;
                        if(l != min_dn){
                            sign_pow_dn += bit_value(basis.D_dn_basis[j],min_dn);
                        }
                        sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                        sign_FM = pow(-1.0, sign_pow_dn);

                        value = sign_FM*EigVec_[m]*(conjugate(value_in));

                        State_cdagger_on_GS[m_new] += value;
                    }
                }
            }
        }


    }



}


void MODEL_2_orb_Hubb_chain::Read_parameters_for_variational_state
(BASIS_2_orb_Hubb_chain &basis, string filename){
    //    string filepath = filename;
    //    string variational_state_contruction_;
    //    string Variational_State_Contruction_ = "Variational_state_contruction = ";


    //    int offset;
    //    string line;
    //    ifstream inputfile(filepath.c_str());


    //    if(inputfile.is_open())
    //    {
    //        while(!inputfile.eof())
    //        {
    //            getline(inputfile,line);

    //            if ((offset = line.find(Variational_State_Contruction_, 0)) != string::npos) {
    //                variational_state_contruction_ = line.substr (offset+Variational_State_Contruction_.length());				}

    //        }
    //        inputfile.close();
    //    }
    //    else
    //    {cout<<"Unable to open input file while in the Model class."<<endl;}





    //    Variational_state_pair_Geometry.resize(basis.Length);
    //    Variational_state_pair_spin_symmetry.resize(basis.Length);
    //    Variational_state_pair_orbital_symmetry.resize(basis.Length);
    //    Variational_state_pair_sites.resize(basis.Length);

    //    int temp_int;


    //    stringstream variational_state_contruction_stream(variational_state_contruction_);

    //    for(int n=0;n<basis.Length;n++){
    //        variational_state_contruction_stream >> Variational_state_pair_Geometry[n];
    //        assert(Variational_state_pair_Geometry[n]=="D" ||
    //               Variational_state_pair_Geometry[n]=="AC");

    //        variational_state_contruction_stream >> Variational_state_pair_spin_symmetry[n];
    //        assert(Variational_state_pair_spin_symmetry[n]=="S" ||
    //               Variational_state_pair_spin_symmetry[n]=="T");

    //        variational_state_contruction_stream >> Variational_state_pair_orbital_symmetry[n];
    //        assert(Variational_state_pair_orbital_symmetry[n] =="OA" ||
    //               Variational_state_pair_orbital_symmetry[n] =="OS" ||
    //               Variational_state_pair_orbital_symmetry[n] =="NOS1" ||
    //               Variational_state_pair_orbital_symmetry[n] =="NOS2");

    //        variational_state_contruction_stream >> temp_int;
    //        Variational_state_pair_sites[n].first =temp_int;
    //        assert(temp_int>=0 && temp_int<basis.Length);


    //        variational_state_contruction_stream >> temp_int;
    //        Variational_state_pair_sites[n].second =temp_int;
    //        assert(temp_int>=0 && temp_int<basis.Length);


    //    }




}

void MODEL_2_orb_Hubb_chain::Read_parameters(BASIS_2_orb_Hubb_chain &basis, string filename){


    string filepath = filename;


    double temp_val;
    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string ucoul, Ucoul = "U = ";
    string jhund, Jhund = "JHund = ";
    string upcoul, Upcoul = "Uprime = ";
    string dz_anisotropy_, Dz_Anisotropy_ = "Dz_Anisotropy = ";
    string hmag, Hmag = "H_mag = ";
    string cfs_, CFS_ = "CFS = ";
    string hopp0_, Hopp0_ = "Hopping_mat[0][orb] = ";
    string hopp1_, Hopp1_ = "Hopping_mat[1][orb] = ";




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

            if ((offset = line.find(Length, 0)) != string::npos) {
                length = line.substr (offset + Length.length());		}

            if ((offset = line.find(Ndn, 0)) != string::npos) {
                ndn = line.substr (offset + Ndn.length());		}

            if ((offset = line.find(Nup, 0)) != string::npos) {
                nup= line.substr (offset + Nup.length());		}

            if ((offset = line.find(Ucoul, 0)) != string::npos) {
                ucoul= line.substr (offset + Ucoul.length());		}

            if ((offset = line.find(Jhund, 0)) != string::npos) {
                jhund = line.substr (offset + Jhund.length());		}

            if ((offset = line.find(Dz_Anisotropy_, 0)) != string::npos) {
                dz_anisotropy_ = line.substr (offset + Dz_Anisotropy_.length());		}

            if ((offset = line.find(Upcoul, 0)) != string::npos) {
                upcoul = line.substr (offset + Upcoul.length());		}

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(CFS_, 0)) != string::npos) {
                cfs_ = line.substr (offset+CFS_.length());				}

            if ((offset = line.find(Hopp0_, 0)) != string::npos) {
                hopp0_ = line.substr (offset+Hopp0_.length());				}

            if ((offset = line.find(Hopp1_, 0)) != string::npos) {
                hopp1_ = line.substr (offset+Hopp1_.length());				}




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
    basis.Ndn=atoi(ndn.c_str());
    basis.Nup=atoi(nup.c_str());


    U=atof(ucoul.c_str());
    J_H=atof(jhund.c_str());
    U_p=atof(upcoul.c_str());
    Dz_Anisotropy=atof(dz_anisotropy_.c_str());

    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }


    stringstream cfs_stream;
    cfs_stream<<cfs_;

    CFS.clear();
    CFS.resize(2);


    for(int n=0;n<2;n++){
        cfs_stream >> temp_val;
        CFS[n]=temp_val;

    }



    Hopping_mat_NN.clear();
    Hopping_mat_NN.resize(2);
    for (int i=0;i<2;i++){
        Hopping_mat_NN[i].resize(2);
    }

    //Hopping_mat_NN[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}


    stringstream hopp0_stream(hopp0_);

    for(int n=0;n<2;n++){
        hopp0_stream >> Hopping_mat_NN[0][n];
    }
    stringstream hopp1_stream(hopp1_);

    for(int n=0;n<2;n++){
        hopp1_stream >> Hopping_mat_NN[1][n];
    }



    Momentum_values.resize(basis.Length);
    if(PBC==true){
        for(int n=0;n<basis.Length;n++){
            Momentum_values[n]=(2.0*n)/(basis.Length*1.0);
        }
    }
    else{
        for(int n=0;n<basis.Length;n++){
            Momentum_values[n]=((n+1)*1.0)/((basis.Length + 1)*1.0);
        }
    }


}


void MODEL_2_orb_Hubb_chain::Read_parameters_for_dynamics(string filename){

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



void MODEL_2_orb_Hubb_chain::Initialize_macro_oprs_to_calculate(BASIS_2_orb_Hubb_chain &basis){
    macro_obs.resize(6);

    macro_obs[0]="U_intra";
    macro_obs[1]="U_inter";
    macro_obs[2]="Hunds_z";
    macro_obs[3]="Hunds_pm";
    macro_obs[4]="U_Pair";
    macro_obs[5]="H_KE";

    Macro_oprts.resize(6);



    int T_no_oprs=6;
    double_type hopp_val;




    for(int i=0;i<T_no_oprs;i++){
        Macro_oprts[i].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
        Macro_oprts[i].ncols = Macro_oprts[i].nrows;
    }


    //Remember H[l][m]=<l|H|m>
    int m;

    double value1, value2, value3, value4;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value1=0;
            //intra-orbital coulomb repulsion:
            value1+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);

            if(value1!=0){
                Macro_oprts[0].value.push_back(value1);
                Macro_oprts[0].rows.push_back(m);
                Macro_oprts[0].columns.push_back(m);
            }


            value2=0;
            //inter-orbital coulomb repulsion:
            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
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


            if(value2!=0){
                Macro_oprts[1].value.push_back(value2);
                Macro_oprts[1].rows.push_back(m);
                Macro_oprts[1].columns.push_back(m);
            }


            value3=0;
            //SzSz Hunds coupling:
            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
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

            if(value3!=0){
                Macro_oprts[2].value.push_back(value3);
                Macro_oprts[2].rows.push_back(m);
                Macro_oprts[2].columns.push_back(m);
            }

            //CFS
            value4=0;
            for(int gamma=0;gamma<2;gamma++){
                for(int site=0;site<basis.Length;site++){
                    value4+=(CFS[gamma])*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                              );
                }
            }

            if(value4!=0){
                Macro_oprts[5].value.push_back(value4);
                Macro_oprts[5].rows.push_back(m);
                Macro_oprts[5].columns.push_back(m);
            }


        }
    }



    double value;

    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM ;
    int sign_pow_up, sign_pow_dn ;
    int l,lp;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;


            for(int gamma=0;gamma<2;gamma++){
                for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
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

                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                            assert(m_new<m);

                            Macro_oprts[3].value.push_back(sign_FM*(0.5*(-J_H*2.0)));
                            Macro_oprts[3].rows.push_back(m_new);
                            Macro_oprts[3].columns.push_back(m);

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

                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);



                            assert(m_new<m);
                            Macro_oprts[4].value.push_back(sign_FM*J_H);
                            Macro_oprts[4].rows.push_back((m_new));
                            Macro_oprts[4].columns.push_back((m));


                        } //Pair-Hopping


                    } // site
                } //gamma_p
            } //gamma

        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals




    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;


            for(int site=0;site<basis.Length ;site++){

                for(int gamma=0;gamma<2;gamma++){

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
                                        PBC
                                        )

                                    )
                            { // nearest neighbour

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


                                    i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                    j_new = j;

                                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up);



                                    if(abs(Hopping_mat_NN[gamma_p][gamma])>0.000001){
                                        assert(m_new<m);
                                        if(neigh==-1 || neigh==basis.Length-1){
                                            hopp_val=Hopping_mat_NN[gamma_p][gamma];
                                        }
                                        else{
                                            hopp_val=conjugate(Hopping_mat_NN[gamma_p][gamma]);
                                        }

                                        Macro_oprts[5].value.push_back(-1.0*sign_FM*(hopp_val));
                                        Macro_oprts[5].rows.push_back((m_new));
                                        Macro_oprts[5].columns.push_back((m));
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

                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site_p)
                                                  - pow(2,gamma*basis.Length + site) );


                                    j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                                    i_new = i;

                                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn);



                                    if(abs(Hopping_mat_NN[gamma_p][gamma])>0.000001){
                                        assert(m_new<m);
                                        if(neigh==-1 || neigh==basis.Length-1){
                                            hopp_val=Hopping_mat_NN[gamma_p][gamma];
                                        }
                                        else{
                                            hopp_val=conjugate(Hopping_mat_NN[gamma_p][gamma]);
                                        }
                                        Macro_oprts[5].value.push_back(-1.0*sign_FM*(conjugate(hopp_val)));
                                        Macro_oprts[5].rows.push_back((m_new));
                                        Macro_oprts[5].columns.push_back((m));}

                                } // if up hopping possible


                            }//nearest neighbour
                        } //gamma_p

                    }//site_p

                } //gamma


            } // site
        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals






}




void MODEL_2_orb_Hubb_chain::Initialize_one_point_to_calculate(BASIS_2_orb_Hubb_chain &basis){
    one_point_obs.resize(7);
    one_point_obs[0]="n_0_up";
    one_point_obs[1]="n_1_up";

    one_point_obs[2]="n_0_dn";
    one_point_obs[3]="n_1_dn";

    one_point_obs[4]="n_i"; //write now n_k is n_i

    one_point_obs[5]="n_0_upn_0_dn";
    one_point_obs[6]="n_1_upn_1_dn";


    int T_no_oprs=one_point_obs.size();

    One_point_oprts.resize(T_no_oprs);



    int orb;
    int spin;


    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }




    for(int opr_no=0;opr_no<4;opr_no++){


        if(one_point_obs[opr_no]=="n_0_up" || one_point_obs[opr_no]=="n_0_dn"){
            orb=0;
        }
        if(one_point_obs[opr_no]=="n_1_up" || one_point_obs[opr_no]=="n_1_dn"){
            orb=1;
        }

        if(one_point_obs[opr_no]=="n_0_up" || one_point_obs[opr_no]=="n_1_up"){
            spin=0;
        }
        else{
            spin=1;
        }


        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }


        //Remember OPR[l][m]=<l|OPR|m>
        int m;
        double value;


        for(int site=0;site<basis.Length;site++){

            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++){
                    m=basis.D_dn_basis.size()*i + j;

                    //n_orb_spin[site]:
                    if(spin==0){
                        value=bit_value(basis.D_up_basis[i],orb*basis.Length + site);
                    }
                    else{
                        value=bit_value(basis.D_dn_basis[j],orb*basis.Length + site);
                    }





                    if(value!=0){
                        One_point_oprts[opr_no][site].value.push_back(value);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }
            }

        }



    }




    //for n_i------------------------------------------------------------------
    for(int site=0;site<basis.Length;site++){
        One_point_oprts[4][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
        One_point_oprts[4][site].ncols = One_point_oprts[4][site].nrows;
    }



    Hamiltonian_1_COO n_i;
    n_i.resize(basis.Length);
    Matrix_COO temp;
    for(int site=0;site<basis.Length;site++){
        temp = One_point_oprts[0][site];
        for(int dof=1;dof<4;dof++){
            Sum(temp, One_point_oprts[dof][site], temp, 1.0, 1.0);
        }

        n_i[site]=temp;
        One_point_oprts[4][site]=temp;
    }


    for(int site=0;site<basis.Length;site++){
        for(int opr_no=5;opr_no<7;opr_no++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }
    }


    int m;
    double value;
    for(int opr_no=5;opr_no<7;opr_no++){
        if(opr_no==5){
            orb=0;
        }
        else if(opr_no==6){
            orb=1;
        }

        for(int site=0;site<basis.Length;site++){
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++){
                    m=basis.D_dn_basis.size()*i + j;

                    value=bit_value(basis.D_up_basis[i],orb*basis.Length + site)*
                            bit_value(basis.D_dn_basis[j],orb*basis.Length + site);

                    if(value!=0){
                        One_point_oprts[opr_no][site].value.push_back(value);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }
            }

        }
    }


    /*

    for(int k_n=0;k_n<Momentum_values.size();k_n++){


        double value1, value2;

        temp = n_i[0];

        for(int site=0;site<basis.Length-1;site++){
            if(PBC==true){
                value2=cos((site+2)*Momentum_values[k_n]*PI)*sqrt(1.0/(basis.Length));
            }
            else{
                value2=sin((site+2)*Momentum_values[k_n]*PI)*sqrt(2.0/(basis.Length +1));
            }
            if(site==0){

                if(PBC==true){
                    value1=cos((site+2)*Momentum_values[k_n]*PI)*sqrt(1.0/(basis.Length));
                }
                else{
                    value1=sin((site+1)*Momentum_values[k_n]*PI)*sqrt(2.0/(basis.Length +1));
                }
                Sum(temp, n_i[site+1], temp, value1, value2);}
            else{
                Sum(temp, n_i[site+1], temp, 1.0, value2);
            }

        }

        One_point_oprts[6][k_n]=temp;


    }
*/
    temp.value.clear();
    temp.rows.clear();
    temp.columns.clear();

    for(int site=0;site<basis.Length;site++){
        n_i[site].value.clear();
        n_i[site].rows.clear();
        n_i[site].columns.clear();
    }



    // n_k is done-----------------------------------------------------------------


}


void MODEL_2_orb_Hubb_chain::Initialize_two_point_operator_sites_specific(string type , Matrix_COO &OPR_, int site, int site2, BASIS_2_orb_Hubb_chain &basis){

    OPR_.columns.clear();
    OPR_.rows.clear();
    OPR_.value.clear();
    OPR_.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    OPR_.ncols = basis.D_up_basis.size()*basis.D_dn_basis.size();

    int m;
    double value;
    int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
    double sign_FM;

    if(type == "SzSz"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;

                value=0;
                for(int gamma=0;gamma<2;gamma++){
                    for(int gamma_p=0;gamma_p<2;gamma_p++){
                        value+=0.25*( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                        bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                      ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2) -
                                        bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2) )
                                      );

                    }
                }

                if(value!=0){
                    OPR_.value.push_back(value);
                    OPR_.rows.push_back(m);
                    OPR_.columns.push_back(m);
                }
            }
        }
    }

    if(type == "SpSm"){

        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;


                for(int gamma=0;gamma<2;gamma++){
                    for(int gamma_p=0;gamma_p<2;gamma_p++){
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

                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                            //assert(m_new<m);

                            OPR_.value.push_back(sign_FM);
                            OPR_.rows.push_back(m_new);
                            OPR_.columns.push_back(m);


                        }

                        if((site==site2) && (gamma==gamma_p)){
                            if(
                                    ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==1)
                                     &&
                                     (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==0)
                                     )
                                    )
                            {
                                OPR_.value.push_back(1.0);
                                OPR_.rows.push_back(m);
                                OPR_.columns.push_back(m);
                            }
                        }
                    }
                }
            }
        }
    }

    if(type == "SmSp"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;

                for(int gamma=0;gamma<2;gamma++){
                    for(int gamma_p=0;gamma_p<2;gamma_p++){
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

                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                            //assert(m_new<m);
                            OPR_.value.push_back(sign_FM);
                            OPR_.rows.push_back(m_new);
                            OPR_.columns.push_back(m);
                        }

                        if((site==site2) && (gamma==gamma_p)){
                            if(
                                    ((bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2)==0)
                                     &&
                                     (bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2)==1)
                                     )
                                    )
                            {
                                OPR_.value.push_back(1.0);
                                OPR_.rows.push_back(m);
                                OPR_.columns.push_back(m);
                            }
                        }
                    }
                }
            }
        }
    }


}


void MODEL_2_orb_Hubb_chain::Initialize_two_point_to_calculate(BASIS_2_orb_Hubb_chain &basis){
    two_point_obs.resize(15);
    two_point_obs[0]="SzSz";
    two_point_obs[1]="SpSm";
    two_point_obs[2]="SmSp";
    two_point_obs[3]="SzSz_00";
    two_point_obs[4]="SzSz_11";
    two_point_obs[5]="SzSz_01";
    two_point_obs[6]="SzSz_10";
    two_point_obs[7]="SpSm_00";
    two_point_obs[8]="SpSm_11";
    two_point_obs[9]="SpSm_01";
    two_point_obs[10]="SpSm_10";
    two_point_obs[11]="SmSp_00";
    two_point_obs[12]="SmSp_11";
    two_point_obs[13]="SmSp_01";
    two_point_obs[14]="SmSp_10";
    Two_point_oprts.resize(15);


    int T_no_oprs=15;



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
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;
            double value;


            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;

                            value=0;
                            for(int gamma=0;gamma<2;gamma++){
                                for(int gamma_p=0;gamma_p<2;gamma_p++){
                                    value+=0.25*( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                                    bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                                  ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2) -
                                                    bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2) )
                                                  );

                                }
                            }



                            if(value!=0){
                                Two_point_oprts[opr_no][site][site2].value.push_back(value);
                                Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                                Two_point_oprts[opr_no][site][site2].columns.push_back(m);
                            }
                        }
                    }

                }

            }

        }


        if(two_point_obs[opr_no]=="SzSz_00" || two_point_obs[opr_no]=="SzSz_01" || two_point_obs[opr_no]=="SzSz_10" || two_point_obs[opr_no]=="SzSz_11"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;
            double value;


            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;

                            value=0;
                            int gamma,gamma_p;
                            if(two_point_obs[opr_no]=="SzSz_00"){
                                gamma=0;gamma_p=0;
                            }
                            else if(two_point_obs[opr_no]=="SzSz_01"){
                                gamma=0;gamma_p=1;
                            }
                            else if(two_point_obs[opr_no]=="SzSz_10"){
                                gamma=1;gamma_p=0;
                            }
                            else{
                                assert( (two_point_obs[opr_no]=="SzSz_11"));
                                gamma=1;gamma_p=1;

                            }
                            value+=0.25*( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) -
                                            bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )*
                                          ( bit_value(basis.D_up_basis[i],gamma_p*basis.Length + site2) -
                                            bit_value(basis.D_dn_basis[j],gamma_p*basis.Length + site2) )
                                          );






                            if(value!=0){
                                Two_point_oprts[opr_no][site][site2].value.push_back(value);
                                Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                                Two_point_oprts[opr_no][site][site2].columns.push_back(m);
                            }
                        }
                    }

                }

            }

        }


        if(two_point_obs[opr_no]=="SpSm"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;


                            for(int gamma=0;gamma<2;gamma++){
                                for(int gamma_p=0;gamma_p<2;gamma_p++){
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

                                        i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                        j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                        m_new = basis.D_dn_basis.size()*i_new + j_new;

                                        l=gamma*basis.Length + site;
                                        lp=gamma_p*basis.Length + site2;

                                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                        //assert(m_new<m);

                                        Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM);
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
                                            Two_point_oprts[opr_no][site][site2].value.push_back(1.0);
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

        }



        if(two_point_obs[opr_no]=="SpSm_00" || two_point_obs[opr_no]=="SpSm_11" || two_point_obs[opr_no]=="SpSm_01"|| two_point_obs[opr_no]=="SpSm_10"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;

                            int gamma,gamma_p;
                            if(two_point_obs[opr_no]=="SpSm_00"){
                                gamma=0;gamma_p=0;
                            }
                            else if(two_point_obs[opr_no]=="SpSm_01"){
                                gamma=0;gamma_p=1;
                            }
                            else if(two_point_obs[opr_no]=="SpSm_10"){
                                gamma=1;gamma_p=0;
                            }
                            else{
                                assert( (two_point_obs[opr_no]=="SpSm_11"));
                                gamma=1;gamma_p=1;

                            }

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

                                i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l=gamma*basis.Length + site;
                                lp=gamma_p*basis.Length + site2;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                //assert(m_new<m);

                                Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM);
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
                                    Two_point_oprts[opr_no][site][site2].value.push_back(1.0);
                                    Two_point_oprts[opr_no][site][site2].rows.push_back(m);
                                    Two_point_oprts[opr_no][site][site2].columns.push_back(m);

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
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;


                            for(int gamma=0;gamma<2;gamma++){
                                for(int gamma_p=0;gamma_p<2;gamma_p++){
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

                                        i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                        j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                        m_new = basis.D_dn_basis.size()*i_new + j_new;

                                        l=gamma*basis.Length + site;
                                        lp=gamma_p*basis.Length + site2;

                                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                        //assert(m_new<m);

                                        Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM);
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
                                            Two_point_oprts[opr_no][site][site2].value.push_back(1.0);
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

        }


        if(two_point_obs[opr_no]=="SmSp_00" || two_point_obs[opr_no]=="SmSp_11" || two_point_obs[opr_no]=="SmSp_01" || two_point_obs[opr_no]=="SmSp_10"){




            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;


            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;


                            int gamma,gamma_p;
                            if(two_point_obs[opr_no]=="SmSp_00"){
                                gamma=0;gamma_p=0;
                            }
                            else if(two_point_obs[opr_no]=="SmSp_01"){
                                gamma=0;gamma_p=1;
                            }
                            else if(two_point_obs[opr_no]=="SmSp_10"){
                                gamma=1;gamma_p=0;
                            }
                            else{
                                assert( (two_point_obs[opr_no]=="SmSp_11"));
                                gamma=1;gamma_p=1;

                            }
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

                                i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l=gamma*basis.Length + site;
                                lp=gamma_p*basis.Length + site2;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                //assert(m_new<m);

                                Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM);
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
                                    Two_point_oprts[opr_no][site][site2].value.push_back(1.0);
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


}


void MODEL_2_orb_Hubb_chain::Initialize_Opr_for_Dynamics(BASIS_2_orb_Hubb_chain &basis){

    if(Dyn_Momentum_Resolved){


        Hamiltonian_1_COO Oprs_local;
        Oprs_local.resize(basis.Length);

        if(Dyn_opr_string == "Sz"){


            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                Oprs_local[site].ncols = Oprs_local[site].nrows;
            }

            //creating Local opeartors-------------------------------------------
            //Remember OPR[l][m]=<l|OPR|m>
            int m;
            double value;
            for(int site=0;site<basis.Length;site++){
                Oprs_local[site].value.clear();
                Oprs_local[site].rows.clear();
                Oprs_local[site].columns.clear();

                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //Sz_total[site]:
                        value=0;
                        for (int orb=0;orb<2;orb++){
                            value = value + 0.5*(bit_value(basis.D_up_basis[i],orb*basis.Length + site)
                                                 - bit_value(basis.D_dn_basis[j],orb*basis.Length + site) );
                        }

                        if(value!=0){
                            Oprs_local[site].value.push_back(value);
                            Oprs_local[site].rows.push_back(m);
                            Oprs_local[site].columns.push_back(m);
                        }
                    }
                }

            }
            //local operators created ----------------------------------------------

            //In Momentum space--------For PBC use cosine, for OBC use sin----------

            Matrix_COO temp;
            temp=Oprs_local[0];
            double value1, value2;
            for(int site=0;site<basis.Length-1;site++){
                if(PBC==false){
                    value2=sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                }
                else{
                    value2=cos((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                }
                if(site==0){
                    if(PBC==false){
                        value1=sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    }
                    else{
                        value1=cos((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    }
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







    }

    else{
        int gamma;
        if(Dyn_opr_string == "J_0"){
            gamma=0;
        }
        if(Dyn_opr_string == "J_1"){
            gamma=1;
        }


        if(Dyn_opr_string == "J_0" || Dyn_opr_string == "J_1"){
            Matrix_COO Opr;
            Opr.nrows=basis.D_up_basis.size()*basis.D_dn_basis.size();
            Opr.ncols=Opr.nrows;
            Opr.value.clear();
            Opr.columns.clear();
            Opr.rows.clear();


            double value;
            int m;
            int D_up,D_dn;
            int i_new,j_new;
            int m_new;
            double sign_FM;
            int sign_pow_up, sign_pow_dn;
            int l,lp;
            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++){
                    m=basis.D_dn_basis.size()*i + j;

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


                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma*basis.Length + site_p;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);



                            if(abs(Hopping_mat_NN[gamma][gamma])>0.000001){

                                Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[gamma][gamma]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma][gamma]));
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


                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                            i_new = i;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l=gamma*basis.Length + site;
                            lp=gamma*basis.Length + site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);



                            if(abs(Hopping_mat_NN[gamma][gamma])>0.000001){

                                Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[gamma][gamma]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma][gamma]));
                                Opr.rows.push_back((m));
                                Opr.columns.push_back((m_new));


                            }

                        } // if up hopping possible






                    } // site
                }// "j" i.e dn_decimals
            } // "i" i.e up_decimals



            Dyn_opr=Opr;
            Opr.value.clear();
            Opr.columns.clear();
            Opr.rows.clear();


        }












    }

}



void MODEL_2_orb_Hubb_chain::Calculate_Local_Obs_for_States_to_Look(bool calculate_local_obs_for_states_to_look,
                                                                    Mat_1_int & states_to_look,
                                                                    string file_Loc_obs_in_basis_of_states,
                                                                    int no_basis_to_check,
                                                                    Mat_2_pair_realint &Overlaps,
                                                                    BASIS_2_orb_Hubb_chain & basis){

    bool grouping_by_orb2 =true;




    if(calculate_local_obs_for_states_to_look == true){
        div_t divresult;

        int nup,ndn,temp_d;


        for(int Ts=0;Ts<states_to_look.size();Ts++){

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



            string out0 = file_Loc_obs_in_basis_of_states + NumberToString (states_to_look[Ts])+ ".txt";
            ofstream file_out0(out0.c_str());

            for(int bi=0;bi<no_basis_to_check;bi++){
                file_out0<<"#For Basis no = "<<Overlaps[Ts][bi].second<<"["<<bi<<"]"<<endl;
                // m=basis.D_dn_basis.size()*i + j;
                divresult = div (Overlaps[Ts][bi].second,basis.D_dn_basis.size());

                int bi_up = divresult.quot;
                int bi_dn = divresult.rem;


                for(int gamma=1;gamma>=0;gamma--){
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

                file_out0<<Overlaps[Ts][bi].first<<endl;
            }




            if(grouping_by_orb2==true){

                string out2 = file_Loc_obs_in_basis_of_states + NumberToString (states_to_look[Ts])+ "_groupwise.txt";
                ofstream file_out2(out2.c_str());

                for(int bi=0;bi<no_basis_to_check;bi++){

                    divresult = div (Overlaps[Ts][bi].second,basis.D_dn_basis.size());

                    int bi_up = divresult.quot;
                    int bi_dn = divresult.rem;

                    for(int site=0;site<basis.Length;site++){
                        nup_2[site] = bit_value(basis.D_up_basis[bi_up],2*basis.Length + site);
                        ndn_2[site] = bit_value(basis.D_dn_basis[bi_dn],2*basis.Length + site);
                    }

                    int pos;
                    if( present_before(nup_2, ndn_2, nup_2_group, ndn_2_group, pos) == false ) //if new group
                    {
                        group_weight.push_back(Overlaps[Ts][bi].first);
                        nup_2_group.push_back(nup_2);
                        ndn_2_group.push_back(ndn_2);
                    }
                    else{
                        bool test = present_before(nup_2, ndn_2, nup_2_group, ndn_2_group, pos);
                        group_weight[pos] = group_weight[pos]  + Overlaps[Ts][bi].first;
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
//#endif

/*
This class includes the Model for which Lanczos is being done
*/
//#ifndef USE_COMPLEX
#include "Model_2_orb_Hubbard_chain_KSector.h"
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

void MODEL_2_orb_Hubb_chain_KSector::Add_diagonal_terms(BASIS_2_orb_Hubb_chain_KSector &basis){


    cout<<"Started Hamiltonian construction: Diagonal"<<endl;

    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());

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


    cout<<"Done Hamiltonian construction: Diagonal"<<endl;

}

void MODEL_2_orb_Hubb_chain_KSector::Add_non_diagonal_terms(BASIS_2_orb_Hubb_chain_KSector &basis){


    cout<<"Started Hamiltonian construction: Non Diagonal"<<endl;
    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;


    double value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
    int l,lp;
    int range_min, range_max;
    bool row_found_;
    double_type phase_;
    int Inv_Trnsltns_;
    complex<double> iota_(0.0,1.0);
    int D_up_temp ,D_dn_temp;

    bool repeating_rows;
    int row_counter;
    int check_min, check_max;
    //cout<<"here 2"<<endl;



    for (int i=0;i<basis.D_up_basis.size();i++){

        /******TO REMOVE*********************
        cout<<"up("<<i<<") : ";
        print_binary_of_decimal(basis.D_up_basis[i]);
        cout<<"dn("<<i<<") : ";
        print_binary_of_decimal(basis.D_dn_basis[i]);
        cout<<endl;
        ************************************/

        m=i;
        j=i;
        value=0;
        row_counter=0;
        repeating_rows=false;
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

                        sign_pow_up=0;
                        sign_pow_dn=0;

                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );
                        D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site)
                                      - pow(2,gamma*basis.Length + site) );

                        /******TO REMOVE**********************
                        cout<<"After acting operator:"<<endl;
                        cout<<"up : ";
                        print_binary_of_decimal(D_up);
                        cout<<"dn : ";
                        print_binary_of_decimal(D_dn);
                        cout<<endl;
                        ***********************************/


                        D_up_temp=D_up;
                        D_dn_temp=D_dn;
                        row_found_=false;
                        for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

                            if(inv_trnsltns>0){
                                //Inv Translation on orb-0,spin_dn
                                sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                                        1*bit_value(D_dn_temp,0);
                                if(bit_value(D_dn_temp,basis.Length -1)==1){
                                    sign_pow_dn += 1*sign_pow_dn_orb0;
                                }
                                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                                //Inv Translation on orb-1,spin_dn
                                sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                                        1*bit_value(D_dn_temp,basis.Length);
                                if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                                    sign_pow_dn += 1*sign_pow_dn_orb1;
                                }

                                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                                //Inv Translation on orb-0,spin_up
                                sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                                        1*bit_value(D_up_temp,0);
                                if(bit_value(D_up_temp,basis.Length -1)==1){
                                    sign_pow_up += 1*sign_pow_up_orb0;
                                }
                                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                                //Inv Translation on orb-1,spin_up
                                sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                                        1*bit_value(D_up_temp,basis.Length);
                                if(bit_value(D_up_temp,2*basis.Length -1)==1){
                                    sign_pow_up += 1*sign_pow_up_orb1;
                                }
                                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

                            }
                            else{
                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                            }

                            /******TO REMOVE********************
                            cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
                            cout<<"up : ";
                            print_binary_of_decimal(D_up_temp);
                            cout<<"dn : ";
                            print_binary_of_decimal(D_dn_temp);
                            cout<<endl;
                            ***********************************/

                            if(D_up_temp>=basis.Dup_Range.size())
                            {
                                row_found_=false;
                            }
                            else
                            {
                                assert(D_up_temp<basis.Dup_Range.size());
                                range_min=basis.Dup_Range[D_up_temp].first;
                                range_max=basis.Dup_Range[D_up_temp].second;
                                if(range_min==-1)
                                {
                                    row_found_=false;
                                    assert(range_max==-1);
                                }
                                else
                                {
                                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                    if(i_new==-1){
                                        row_found_=false;
                                    }
                                    else{
                                        row_found_=true;
                                        Inv_Trnsltns_=inv_trnsltns;
                                        break;
                                    }
                                }

                            }
                        }

                        if(row_found_==true){
                            m_new = i_new;

#ifdef USE_COMPLEX
                            phase_=exp(-1.0*iota_*( (2.0*PI*(1.0*basis.Momentum_n))/(basis.Length) )*(1.0*Inv_Trnsltns_) )*
                                    sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif
#ifndef USE_COMPLEX
                            if(basis.Momentum_n!=0){
                                cout<<"ONLY K=0 is allowed in real space calculations"<<endl;
                            }
                            assert(basis.Momentum_n==0);
                            phase_=one*sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif


                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site;

                            sign_pow_up += one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn += one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                            assert(m_new<m);


                            check_min=Hamil.rows.size()-1;
                            check_max=(Hamil.rows.size()-1)-row_counter;
                            for(int check_=check_min;check_>check_max;check_--){
                                if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                    Hamil.value[check_] +=sign_FM*(0.5*(-J_H*2.0))*one*phase_;
                                    repeating_rows=true;
                                    break;
                                }
                            }

                            if(!repeating_rows){
                                Hamil.value.push_back(sign_FM*(0.5*(-J_H*2.0))*one*phase_);
                                Hamil.rows.push_back(m_new);
                                Hamil.columns.push_back(m);
                                row_counter++;
                            }

                            //cout<<m<<"=====>"<<m_new<<endl;

                        }

                    } // if SpSm possible

                } // site
            } //gamma_p
        } //gamma


        //Pair hopping: P_i\gamma = c_i\gamma\dn c_i\gamma\dn
        //there have to be pair present at i,gamma_p
        //there have to nothing at i,gamma
        repeating_rows=false;
        for(int gamma=0;gamma<2;gamma++){
            for(int gamma_p=gamma+1;gamma_p<2;gamma_p++){
                for(int site=0;site<basis.Length;site++){
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

                        sign_pow_up=0;
                        sign_pow_dn=0;

                        D_up = (int) (basis.D_up_basis[i] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );
                        D_dn = (int) (basis.D_dn_basis[j] - pow(2,gamma_p*basis.Length + site)
                                      + pow(2,gamma*basis.Length + site) );

                        D_up_temp=D_up;
                        D_dn_temp=D_dn;
                        row_found_=false;
                        for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

                            if(inv_trnsltns>0){
                                //Inv Translation on orb-0,spin_dn
                                sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                                        1*bit_value(D_dn_temp,0);
                                if(bit_value(D_dn_temp,basis.Length -1)==1){
                                    sign_pow_dn += 1*sign_pow_dn_orb0;
                                }
                                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                                //Inv Translation on orb-1,spin_dn
                                sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                                        1*bit_value(D_dn_temp,basis.Length);
                                if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                                    sign_pow_dn += 1*sign_pow_dn_orb1;
                                }

                                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                                //Inv Translation on orb-0,spin_up
                                sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                                        1*bit_value(D_up_temp,0);
                                if(bit_value(D_up_temp,basis.Length -1)==1){
                                    sign_pow_up += 1*sign_pow_up_orb0;
                                }
                                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                                //Inv Translation on orb-1,spin_up
                                sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                                        1*bit_value(D_up_temp,basis.Length);
                                if(bit_value(D_up_temp,2*basis.Length -1)==1){
                                    sign_pow_up += 1*sign_pow_up_orb1;
                                }
                                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

                            }
                            else{
                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                            }

                            /******TO REMOVE********************
                            cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
                            cout<<"up : ";
                            print_binary_of_decimal(D_up_temp);
                            cout<<"dn : ";
                            print_binary_of_decimal(D_dn_temp);
                            cout<<endl;
                            ***********************************/

                            if(D_up_temp>=basis.Dup_Range.size())
                            {
                                row_found_=false;
                            }
                            else
                            {
                                assert(D_up_temp<basis.Dup_Range.size());
                                range_min=basis.Dup_Range[D_up_temp].first;
                                range_max=basis.Dup_Range[D_up_temp].second;
                                if(range_min==-1)
                                {
                                    row_found_=false;
                                    assert(range_max==-1);
                                }
                                else
                                {
                                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                    if(i_new==-1){
                                        row_found_=false;
                                    }
                                    else{
                                        row_found_=true;
                                        Inv_Trnsltns_=inv_trnsltns;
                                        break;
                                    }
                                }

                            }
                        }


                        if(row_found_==true){
                            m_new = i_new;

#ifdef USE_COMPLEX
                            phase_=exp(-1.0*iota_*( (2.0*PI*basis.Momentum_n)/(basis.Length) )*(1.0*Inv_Trnsltns_) )*
                                    sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif
#ifndef USE_COMPLEX
                            if(basis.Momentum_n!=0){
                                cout<<"ONLY K=0 is allowed in real space calculations"<<endl;
                            }
                            assert(basis.Momentum_n==0);
                            phase_=one*sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif


                            l=gamma*basis.Length + site;
                            lp=gamma_p*basis.Length + site;

                            sign_pow_up += one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn += one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);

                            assert(m_new<m);


                            check_min=Hamil.rows.size()-1;
                            check_max=(Hamil.rows.size()-1)-row_counter;
                            for(int check_=check_min;check_>check_max;check_--){
                                if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                    Hamil.value[check_] +=sign_FM*J_H*one*phase_;
                                    repeating_rows=true;
                                    break;
                                }
                            }

                            if(!repeating_rows){
                                Hamil.value.push_back(sign_FM*J_H*one*phase_);
                                Hamil.rows.push_back(m_new);
                                Hamil.columns.push_back(m);
                                row_counter++;
                            }
                        }


                    } //Pair-Hopping
                } // site
            } //gamma_p
        } //gamma

        if(m%1000 ==1){
            //cout<<"done "<<m<<" basis"<<endl;
        }


    } // "i" i.e up_decimals, dn_decimals


    cout<<"Done Hamiltonian construction: Non Diagonal"<<endl;

}

void MODEL_2_orb_Hubb_chain_KSector::Add_connections(BASIS_2_orb_Hubb_chain_KSector &basis){


    cout<<"Started Hamiltonian construction: Connections"<<endl;
    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;

    double value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM, sign_FM_by_Translations;
    int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
    int sign_pow_up, sign_pow_dn;
    int l,lp;
    int range_min, range_max;
    bool row_found_;
    double_type phase_;
    int Inv_Trnsltns_;
    complex<double> iota_(0.0,1.0);
    int D_up_temp ,D_dn_temp;
    bool repeating_rows;
    int row_counter;
    int check_min, check_max;

    for (int i=0;i<basis.D_up_basis.size();i++){

        /******TO REMOVE*********************
        cout<<"up("<<i<<") : ";
        print_binary_of_decimal(basis.D_up_basis[i]);
        cout<<"dn("<<i<<") : ";
        print_binary_of_decimal(basis.D_dn_basis[i]);
        cout<<endl;
        ************************************/

        m=i;
        j=i;

        value=0;

        row_counter=0;
        for(int site=0;site<basis.Length ;site++){
            for(int gamma=0;gamma<2;gamma++){
                for(int site_p=0;site_p<basis.Length ;site_p++){
                    int neigh =site-site_p;

                    for(int gamma_p=0;gamma_p<2;gamma_p++){

                        if( (abs(neigh)==1)
                                ||
                                (
                                    ((abs(neigh)==(basis.Length - 1) ))
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
                                    &&
                                    ((Hopping_mat_NN[gamma_p][gamma])!=0)

                                    )
                            {

                                sign_pow_up=0;
                                sign_pow_dn=0;

                                D_up = (int) (basis.D_up_basis[i] + pow(2,gamma_p*basis.Length + site_p)
                                              - pow(2,gamma*basis.Length + site) );
                                D_dn = basis.D_dn_basis[m] ;

                                /******TO REMOVE***********************
                                cout<<"After up hopping:"<<endl;
                                cout<<"up : ";
                                print_binary_of_decimal(D_up);
                                cout<<"dn : ";
                                print_binary_of_decimal(D_dn);
                                cout<<endl;
                                ***********************************/

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;


                                for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

                                    if(inv_trnsltns>0){
                                        //Inv Translation on orb-0,spin_dn
                                        sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                                                1*bit_value(D_dn_temp,0);
                                        if(bit_value(D_dn_temp,basis.Length -1)==1){
                                            sign_pow_dn += 1*sign_pow_dn_orb0;
                                        }
                                        D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                                        //Inv Translation on orb-1,spin_dn
                                        sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                                                1*bit_value(D_dn_temp,basis.Length);
                                        if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                                            sign_pow_dn += 1*sign_pow_dn_orb1;
                                        }

                                        D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                                        //Inv Translation on orb-0,spin_up
                                        sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                                                1*bit_value(D_up_temp,0);
                                        if(bit_value(D_up_temp,basis.Length -1)==1){
                                            sign_pow_up += 1*sign_pow_up_orb0;
                                        }
                                        D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                                        //Inv Translation on orb-1,spin_up
                                        sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                                                1*bit_value(D_up_temp,basis.Length);
                                        if(bit_value(D_up_temp,2*basis.Length -1)==1){
                                            sign_pow_up += 1*sign_pow_up_orb1;
                                        }
                                        D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

                                    }
                                    else{
                                        D_up_temp=D_up;
                                        D_dn_temp=D_dn;
                                    }

                                    /******TO REMOVE********************
                                    cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
                                    cout<<"up : ";
                                    print_binary_of_decimal(D_up_temp);
                                    cout<<"dn : ";
                                    print_binary_of_decimal(D_dn_temp);
                                    cout<<endl;
                                    ***********************************/

                                    if(D_up_temp>=basis.Dup_Range.size())
                                    {
                                        row_found_=false;
                                    }
                                    else
                                    {
                                        assert(D_up_temp<basis.Dup_Range.size());
                                        range_min=basis.Dup_Range[D_up_temp].first;
                                        range_max=basis.Dup_Range[D_up_temp].second;
                                        if(range_min==-1)
                                        {
                                            row_found_=false;
                                            assert(range_max==-1);
                                        }
                                        else
                                        {
                                            i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                            if(i_new==-1){
                                                row_found_=false;
                                            }
                                            else{
                                                row_found_=true;
                                                Inv_Trnsltns_=inv_trnsltns;
                                                break;
                                            }
                                        }

                                    }
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( (2.0*PI*basis.Momentum_n)/(basis.Length) )*(1.0*Inv_Trnsltns_) )
                                            *sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif
#ifndef USE_COMPLEX
                                    if(basis.Momentum_n!=0){
                                        cout<<"ONLY K=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_n==0);
                                    phase_=one*sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif


                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_up += one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);

                                    if(m_new<=m)
                                    {
                                        if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                            repeating_rows=false;
                                            check_min=Hamil.rows.size()-1;
                                            check_max=(Hamil.rows.size()-1)-row_counter;
                                            // cout<<check_min<<endl;
                                            // cout<<check_max<<endl;
                                            for(int check_=check_min;check_>check_max;check_--){
                                                if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                                    Hamil.value[check_] +=-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*phase_;
                                                    repeating_rows=true;
                                                    break;
                                                }
                                            }

                                            if(!repeating_rows){
                                                Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*phase_);
                                                Hamil.rows.push_back(m_new);
                                                Hamil.columns.push_back(m);
                                                row_counter++;
                                            }
                                        }
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
                                    &&
                                    ((Hopping_mat_NN[gamma_p][gamma])!=0)
                                    )
                            {

                                sign_pow_dn=0;
                                sign_pow_up=0;

                                D_up = basis.D_up_basis[m] ;
                                D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma_p*basis.Length + site_p)
                                              - pow(2,gamma*basis.Length + site) );

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;

                                /******TO REMOVE***********************
                                cout<<"After down hopping:"<<endl;
                                cout<<"up : ";
                                print_binary_of_decimal(D_up);
                                cout<<"dn : ";
                                print_binary_of_decimal(D_dn);
                                cout<<endl;
                                ***********************************/

                                row_found_=false;
                                for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

                                    if(inv_trnsltns>0){
                                        //Inv Translation on orb-0,spin_dn
                                        sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                                                1*bit_value(D_dn_temp,0);
                                        if(bit_value(D_dn_temp,basis.Length -1)==1){
                                            sign_pow_dn += 1*sign_pow_dn_orb0;
                                        }
                                        D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                                        //Inv Translation on orb-1,spin_dn
                                        sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                                                1*bit_value(D_dn_temp,basis.Length);
                                        if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                                            sign_pow_dn += 1*sign_pow_dn_orb1;
                                        }

                                        D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                                        //Inv Translation on orb-0,spin_up
                                        sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                                                1*bit_value(D_up_temp,0);
                                        if(bit_value(D_up_temp,basis.Length -1)==1){
                                            sign_pow_up += 1*sign_pow_up_orb0;
                                        }
                                        D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                                        //Inv Translation on orb-1,spin_up
                                        sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                                                1*bit_value(D_up_temp,basis.Length);
                                        if(bit_value(D_up_temp,2*basis.Length -1)==1){
                                            sign_pow_up += 1*sign_pow_up_orb1;
                                        }
                                        D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

                                    }
                                    else{
                                        D_up_temp=D_up;
                                        D_dn_temp=D_dn;
                                    }

                                    /******TO REMOVE********************
                                    cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
                                    cout<<"up : ";
                                    print_binary_of_decimal(D_up_temp);
                                    cout<<"dn : ";
                                    print_binary_of_decimal(D_dn_temp);
                                    cout<<endl;
                                    ***********************************/

                                    if(D_up_temp>=basis.Dup_Range.size())
                                    {
                                        row_found_=false;
                                    }
                                    else
                                    {
                                        assert(D_up_temp<basis.Dup_Range.size());
                                        range_min=basis.Dup_Range[D_up_temp].first;
                                        range_max=basis.Dup_Range[D_up_temp].second;
                                        if(range_min==-1)
                                        {
                                            row_found_=false;
                                            assert(range_max==-1);
                                        }
                                        else
                                        {
                                            i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                            if(i_new==-1){
                                                row_found_=false;
                                            }
                                            else{
                                                row_found_=true;
                                                Inv_Trnsltns_=inv_trnsltns;
                                                break;
                                            }
                                        }

                                    }

                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( (2.0*PI*basis.Momentum_n)/(basis.Length) )*(1.0*Inv_Trnsltns_) )
                                            *sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif
#ifndef USE_COMPLEX
                                    if(basis.Momentum_n!=0){
                                        cout<<"ONLY K=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_n==0);
                                    phase_=one*sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif

                                    l=gamma*basis.Length + site;
                                    lp=gamma_p*basis.Length + site_p;

                                    sign_pow_dn += one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn+sign_pow_up);


                                    if(m_new<=m)
                                    {
                                        if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                            repeating_rows=false;
                                            check_min=Hamil.rows.size()-1;
                                            check_max=(Hamil.rows.size()-1)-row_counter;
                                            for(int check_=check_min;check_>check_max;check_--){
                                                if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                                    Hamil.value[check_] +=-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*phase_;
                                                    repeating_rows=true;
                                                    break;
                                                }
                                            }

                                            if(!repeating_rows){
                                                Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma])*one*phase_);
                                                Hamil.rows.push_back(m_new);
                                                Hamil.columns.push_back(m);
                                                row_counter++;
                                            }
                                        }
                                    }

                                }
                            } // if dn hopping possible


                        }//nearest neighbour
                    } //gamma_p

                }//site_p

            } //gamma


        } // site

        if(m%1000 ==1){
            //cout<<"done "<<m<<" basis"<<endl;
        }

    } // "i" i.e up_decimals


    cout<<"Done Hamiltonian construction: Connections"<<endl;

}


Mat_1_doub MODEL_2_orb_Hubb_chain_KSector::Act_Orbital_Exchange(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub vec_){



    //bit_value(basis.D_dn_basis[j],o0*basis.Length + site0)
    Mat_1_doub vec1_;
    vec1_.resize(vec_.size());


    int index_no, index_no_new;
    int i_new, j_new, l,l_new;
    int gamma, gamma_new;
    int sign_pow_dn;
    double sign_FM_dn;
    int D_dn_new;
    int D_up_temp, D_dn_temp;
    int sign_pow_up;
    double sign_FM_up;
    int D_up_new;
    int site;
    int j;
    bool row_found_;
    int range_min, range_max;
    int sign_pow_dn_orb0 ,sign_pow_dn_orb1,sign_pow_up_orb0, sign_pow_up_orb1;
    int Inv_Trnsltns_;


    for (int i=0;i<basis.D_up_basis.size();i++){
        j=i;
        index_no=i;

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



        D_up_temp=D_up_new;
        D_dn_temp=D_dn_new;
        row_found_=false;
        sign_pow_dn=0;
        sign_pow_up=0;
        for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

            if(inv_trnsltns>0){
                //Inv Translation on orb-0,spin_dn
                sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,0);
                if(bit_value(D_dn_temp,basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb0;
                }
                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_dn
                sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,basis.Length);
                if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb1;
                }

                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                //Inv Translation on orb-0,spin_up
                sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,0);
                if(bit_value(D_up_temp,basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb0;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_up
                sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,basis.Length);
                if(bit_value(D_up_temp,2*basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb1;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

            }
            else{
                D_up_temp=D_up_new;
                D_dn_temp=D_dn_new;
            }

            /******TO REMOVE********************
            cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
            cout<<"up : ";
            print_binary_of_decimal(D_up_temp);
            cout<<"dn : ";
            print_binary_of_decimal(D_dn_temp);
            cout<<endl;
            ***********************************/

            if(D_up_temp>=basis.Dup_Range.size())
            {
                row_found_=false;
            }
            else
            {
                assert(D_up_temp<basis.Dup_Range.size());
                range_min=basis.Dup_Range[D_up_temp].first;
                range_max=basis.Dup_Range[D_up_temp].second;
                if(range_min==-1)
                {
                    row_found_=false;
                    assert(range_max==-1);
                }
                else
                {
                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                    if(i_new==-1){
                        row_found_=false;
                    }
                    else{
                        row_found_=true;
                        Inv_Trnsltns_=inv_trnsltns;
                        break;
                    }
                }

            }
        }





        //        j_new = Find_int_in_intarray(D_dn_new,basis.D_dn_basis);
        //      i_new = Find_int_in_intarray(D_up_new,basis.D_up_basis);
        //     index_no_new=basis.D_dn_basis.size()*i_new + j_new;

        if(row_found_==true){
            index_no_new=i_new;
            sign_FM_up = sign_FM_up*pow(-1.0,sign_pow_up);
            sign_FM_dn = sign_FM_dn*pow(-1.0,sign_pow_dn);
            vec1_[index_no_new] = sign_FM_up*sign_FM_dn*vec_[index_no];
        }
        else{
            cout<<"issue : not found "<<D_up_new<<", "<<D_dn_new<<endl;
            assert(row_found_);
        }





    }


    return vec1_;
}


Mat_1_doub MODEL_2_orb_Hubb_chain_KSector::Act_Reflection_about_Central_site(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub vec_){



    //bit_value(basis.D_dn_basis[j],o0*basis.Length + site0)
    Mat_1_doub vec1_;
    vec1_.resize(vec_.size());


    int index_no, index_no_new;
    int i_new, j_new, l,l_new;
    int gamma, gamma_new;
    int sign_pow_dn;
    double sign_FM_dn;
    int D_dn_new;
    int D_up_temp, D_dn_temp;
    int sign_pow_up;
    double sign_FM_up;
    int D_up_new;
    int site, site_new;
    int j;
    bool row_found_;
    int range_min, range_max;
    int sign_pow_dn_orb0 ,sign_pow_dn_orb1,sign_pow_up_orb0, sign_pow_up_orb1;
    int Inv_Trnsltns_;


    for (int i=0;i<basis.D_up_basis.size();i++){
        j=i;
        index_no=i;

        sign_FM_dn=1.0;
        sign_FM_up=1.0;

        D_dn_new=0;
        D_up_new=0;

        for( gamma=0;gamma<2;gamma++){
            for( site=0;site<basis.Length;site++){
                site_new = basis.Length -1 -site;
                l= gamma*basis.Length + site;
                l_new = gamma*basis.Length + site_new;


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



        D_up_temp=D_up_new;
        D_dn_temp=D_dn_new;
        row_found_=false;
        sign_pow_dn=0;
        sign_pow_up=0;
        for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

            if(inv_trnsltns>0){
                //Inv Translation on orb-0,spin_dn
                sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,0);
                if(bit_value(D_dn_temp,basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb0;
                }
                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_dn
                sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,basis.Length);
                if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb1;
                }

                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                //Inv Translation on orb-0,spin_up
                sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,0);
                if(bit_value(D_up_temp,basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb0;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_up
                sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,basis.Length);
                if(bit_value(D_up_temp,2*basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb1;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

            }
            else{
                D_up_temp=D_up_new;
                D_dn_temp=D_dn_new;
            }

            /******TO REMOVE********************
            cout<<"Translation("<<inv_trnsltns<<"): "<<endl;
            cout<<"up : ";
            print_binary_of_decimal(D_up_temp);
            cout<<"dn : ";
            print_binary_of_decimal(D_dn_temp);
            cout<<endl;
            ***********************************/

            if(D_up_temp>=basis.Dup_Range.size())
            {
                row_found_=false;
            }
            else
            {
                assert(D_up_temp<basis.Dup_Range.size());
                range_min=basis.Dup_Range[D_up_temp].first;
                range_max=basis.Dup_Range[D_up_temp].second;
                if(range_min==-1)
                {
                    row_found_=false;
                    assert(range_max==-1);
                }
                else
                {
                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                    if(i_new==-1){
                        row_found_=false;
                    }
                    else{
                        row_found_=true;
                        Inv_Trnsltns_=inv_trnsltns;
                        break;
                    }
                }

            }
        }





        //        j_new = Find_int_in_intarray(D_dn_new,basis.D_dn_basis);
        //      i_new = Find_int_in_intarray(D_up_new,basis.D_up_basis);
        //     index_no_new=basis.D_dn_basis.size()*i_new + j_new;

        if(row_found_==true){
            index_no_new=i_new;
            sign_FM_up = sign_FM_up*pow(-1.0,sign_pow_up);
            sign_FM_dn = sign_FM_dn*pow(-1.0,sign_pow_dn);
            vec1_[index_no_new] = sign_FM_up*sign_FM_dn*vec_[index_no];
        }
        else{
            cout<<"issue : not found "<<D_up_new<<", "<<D_dn_new<<endl;
            assert(row_found_);
        }





    }


    return vec1_;
}


void MODEL_2_orb_Hubb_chain_KSector::Initialize_two_point_operator_sites_orbital_specific(string type , Matrix_COO &OPR_, int site1, int gamma1,  int site2, int gamma2, BASIS_2_orb_Hubb_chain_KSector &basis){

    OPR_.columns.clear();
    OPR_.rows.clear();
    OPR_.value.clear();
    OPR_.nrows = basis.D_up_basis.size();
    OPR_.ncols = basis.D_up_basis.size();

    int m,j;
    double value;
    int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
    int D_up_temp, D_dn_temp, sign_pow_dn_orb0, sign_pow_up_orb0, sign_pow_dn_orb1, sign_pow_up_orb1;
    bool row_found_;
    complex<double> iota_(0.0,1.0);
    int range_min, range_max, Inv_Trnsltns_;
    int check_min, check_max;
    double_type phase_;
    double sign_FM;
    int index;
    bool first_time;
    int site1_offseted, site2_offseted;
    int row_counter;
    bool repeating_rows;

    if(type == "SzSz"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            m=i;
            j=i;
            first_time=true;
            for(int site_offset=0;site_offset<basis.Length;site_offset++){
                site1_offseted = (site1 + site_offset)%basis.Length;
                site2_offseted = (site2 + site_offset)%basis.Length;
                value=0.25*( ( bit_value(basis.D_up_basis[i],gamma1*basis.Length + site1_offseted) -
                               bit_value(basis.D_dn_basis[j],gamma1*basis.Length + site1_offseted) )*
                             ( bit_value(basis.D_up_basis[i],gamma2*basis.Length + site2_offseted) -
                               bit_value(basis.D_dn_basis[j],gamma2*basis.Length + site2_offseted) )
                             );

                if(value!=0 && first_time){
                    OPR_.value.push_back((value*one)*(1.0/basis.Length));
                    OPR_.rows.push_back(m);
                    OPR_.columns.push_back(m);
                    first_time=false;
                }
                else if(value!=0 && !first_time){
                    index = OPR_.value.size();
                    OPR_.value[index - 1] += (value*one)*(1.0/basis.Length);
                }
            }
        }
    }

    if(type == "SpSm"){

        for (int i=0;i<basis.D_up_basis.size();i++){
            m=i;
            j=i;
            value=0;
            row_counter=0;
            repeating_rows=false;

            for(int site_offset=0;site_offset<basis.Length;site_offset++){
                //Sp_site_gamma1[site1]*Sm_site_gamma2[site2]  Hunds coupling:
                //there have to be ony up electron in gamma_p, site2
                //there have to be only down electron in gamma, site
                site1_offseted = (site1 + site_offset)%basis.Length;
                site2_offseted = (site2 + site_offset)%basis.Length;

                if(((bit_value(basis.D_dn_basis[j],gamma1*basis.Length + site1_offseted)==1)
                    &&
                    (bit_value(basis.D_up_basis[i],gamma1*basis.Length + site1_offseted)==0)
                    )
                        &&
                        ((bit_value(basis.D_up_basis[i],gamma2*basis.Length + site2_offseted)==1)
                         &&
                         (bit_value(basis.D_dn_basis[j],gamma2*basis.Length + site2_offseted)==0)
                         ))
                {


                    sign_pow_up=0;
                    sign_pow_dn=0;

                    D_up = (int) (basis.D_up_basis[i] - pow(2,gamma2*basis.Length + site2_offseted)
                                  + pow(2,gamma1*basis.Length + site1_offseted) );
                    D_dn = (int) (basis.D_dn_basis[j] + pow(2,gamma2*basis.Length + site2_offseted)
                                  - pow(2,gamma1*basis.Length + site1_offseted) );


                    D_up_temp=D_up;
                    D_dn_temp=D_dn;
                    row_found_=false;
                    for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

                        if(inv_trnsltns>0){
                            //Inv Translation on orb-0,spin_dn
                            sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                                    1*bit_value(D_dn_temp,0);
                            if(bit_value(D_dn_temp,basis.Length -1)==1){
                                sign_pow_dn += 1*sign_pow_dn_orb0;
                            }
                            D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                            //Inv Translation on orb-1,spin_dn
                            sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                                    1*bit_value(D_dn_temp,basis.Length);
                            if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                                sign_pow_dn += 1*sign_pow_dn_orb1;
                            }

                            D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                            //Inv Translation on orb-0,spin_up
                            sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                                    1*bit_value(D_up_temp,0);
                            if(bit_value(D_up_temp,basis.Length -1)==1){
                                sign_pow_up += 1*sign_pow_up_orb0;
                            }
                            D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                            //Inv Translation on orb-1,spin_up
                            sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                                    1*bit_value(D_up_temp,basis.Length);
                            if(bit_value(D_up_temp,2*basis.Length -1)==1){
                                sign_pow_up += 1*sign_pow_up_orb1;
                            }
                            D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

                        }
                        else{
                            D_up_temp=D_up;
                            D_dn_temp=D_dn;
                        }


                        if(D_up_temp>=basis.Dup_Range.size())
                        {
                            row_found_=false;
                        }
                        else
                        {
                            assert(D_up_temp<basis.Dup_Range.size());
                            range_min=basis.Dup_Range[D_up_temp].first;
                            range_max=basis.Dup_Range[D_up_temp].second;
                            if(range_min==-1)
                            {
                                row_found_=false;
                                assert(range_max==-1);
                            }
                            else
                            {
                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                if(i_new==-1){
                                    row_found_=false;
                                }
                                else{
                                    row_found_=true;
                                    Inv_Trnsltns_=inv_trnsltns;
                                    break;
                                }
                            }

                        }
                    }

                    if(row_found_==true){
                        m_new = i_new;

#ifdef USE_COMPLEX
                        phase_=exp(-1.0*iota_*( (2.0*PI*(1.0*basis.Momentum_n))/(basis.Length) )*(1.0*Inv_Trnsltns_) )*
                                sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));

#endif
#ifndef USE_COMPLEX
                        if(basis.Momentum_n!=0){
                            cout<<"ONLY K=0 is allowed in real space calculations"<<endl;
                        }
                        assert(basis.Momentum_n==0);
                        phase_=one*sqrt((1.0*basis.D_Period[m])/(1.0*basis.D_Period[m_new]));
#endif


                        l=gamma1*basis.Length + site1_offseted;
                        lp=gamma2*basis.Length + site2_offseted;

                        sign_pow_up += one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                        sign_pow_dn += one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                        //assert(m_new<m);
                        repeating_rows=false;
                        check_min=OPR_.rows.size()-1;
                        check_max=(OPR_.rows.size()-1)-row_counter;
                        for(int check_=check_min;check_>check_max;check_--){
                            if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                OPR_.value[check_] +=(1.0*sign_FM)*one*phase_*(1.0/basis.Length);
                                repeating_rows=true;
                                break;
                            }
                        }
                        if(!repeating_rows){
                            OPR_.value.push_back((1.0*sign_FM)*one*phase_*(1.0/basis.Length));
                            OPR_.rows.push_back(m_new);
                            OPR_.columns.push_back(m);
                            row_counter++;
                        }



                    }

                }

                if((site1_offseted==site2_offseted) && (gamma1==gamma2)){
                    if(
                            ((bit_value(basis.D_up_basis[i],gamma2*basis.Length + site2_offseted)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j],gamma2*basis.Length + site2_offseted)==0)
                             )
                            )
                    {


                        repeating_rows=false;
                        check_min=OPR_.rows.size()-1;
                        check_max=(OPR_.rows.size()-1)-row_counter;
                        for(int check_=check_min;check_>check_max;check_--){
                            if(OPR_.rows[check_]==m && OPR_.columns[check_]==m){
                                OPR_.value[check_] +=one*(1.0/basis.Length);
                                repeating_rows=true;
                                break;
                            }
                        }
                        if(!repeating_rows){
                            OPR_.value.push_back(one*(1.0/basis.Length));
                            OPR_.rows.push_back(m);
                            OPR_.columns.push_back(m);
                            row_counter++;
                        }

                    }
                }

            }
        }
    }

    if(type == "SmSp"){
        //nothing
    }


}


string MODEL_2_orb_Hubb_chain_KSector::Bond_type_int_tag_to_string(int int_temp){

    string bond_type;
    if(int_temp==0){
        bond_type="D S NOS1";
    }
    if(int_temp==1){
        bond_type="D S NOS2";
    }
    if(int_temp==2){
        bond_type="AC S NOS1";
    }
    if(int_temp==3){
        bond_type="AC S NOS2";
    }

    return bond_type;
}



void MODEL_2_orb_Hubb_chain_KSector::Create_variational_state_strings_2_holes(BASIS_2_orb_Hubb_chain_KSector &basis){


    string Variational_String_;
    int bond1_type,bond2_type,bond3_type,bond4_type,bond5_type,bond6_type,bond7_type;
    string Variational_String0_, Variational_String1_, Variational_String2_, Variational_String3_;
    string Variational_String4_, Variational_String5_, Variational_String6_, Variational_String7_;
    Mat_2_int Allowed_bonds_left_constrained;
    Allowed_bonds_left_constrained.resize(4);

    //D S NOS1 (0)  --->D S NOS1 (0),  AC S NOS1 (2)
    Allowed_bonds_left_constrained[0].push_back(0);Allowed_bonds_left_constrained[0].push_back(2);

    //D S NOS2 (1)  --->D S NOS2 (1),  AC S NOS2 (3)
    Allowed_bonds_left_constrained[1].push_back(1);Allowed_bonds_left_constrained[1].push_back(3);

    //AC S NOS1 (2)  --->D S NOS2 (1),  AC S NOS2 (3)
    Allowed_bonds_left_constrained[2].push_back(1);Allowed_bonds_left_constrained[2].push_back(3);

    //AC S NOS2 (3)  --->D S NOS1 (0),  AC S NOS1 (2)
    Allowed_bonds_left_constrained[3].push_back(0);Allowed_bonds_left_constrained[3].push_back(2);


    if(basis.Length==8){

        //For holes sitting on diagonal on sites-2,3(/) [Only nearest neighbour singlets]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                if(  (bond1_type==0) || (bond1_type==3) ){
                    Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                    Variational_String1_ += " 1 2 ";

                    for(bond3_type=0;bond3_type<4;bond3_type++){
                        if(  (bond3_type==0) || (bond3_type==2) ){
                            Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                            Variational_String3_ += " 3 4 ";

                            for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                                bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                                Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                                Variational_String4_ += " 4 5 ";

                                for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                    bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                    Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                    Variational_String5_ += " 5 6 ";

                                    for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                        bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                        Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                        Variational_String6_ += " 6 7 ";

                                        for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                            bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                            Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                            Variational_String7_ += " 7 0";

                                            if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                    bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                                Variational_String_= Variational_String0_ + Variational_String1_ +
                                                        Variational_String3_ + Variational_String4_ + Variational_String5_
                                                        + Variational_String6_  + Variational_String7_;
                                                variational_state_contruction_strings.push_back(Variational_String_);

                                            }
                                        }
                                    }
                                }
                            }
                        }}
                }}
        }


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;

        //For holes sitting on diagonal on sites-2,3(\) [Only nearest neighbour singlets]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                if(  (bond1_type==1) || (bond1_type==2) ){
                    Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                    Variational_String1_ += " 1 2 ";

                    for(bond3_type=0;bond3_type<4;bond3_type++){
                        if(  (bond3_type==3) || (bond3_type==1) ){
                            Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                            Variational_String3_ += " 3 4 ";

                            for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                                bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                                Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                                Variational_String4_ += " 4 5 ";

                                for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                    bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                    Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                    Variational_String5_ += " 5 6 ";

                                    for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                        bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                        Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                        Variational_String6_ += " 6 7 ";

                                        for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                            bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                            Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                            Variational_String7_ += " 7 0";

                                            if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                    bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                                Variational_String_= Variational_String0_ + Variational_String1_ +
                                                        Variational_String3_ + Variational_String4_ + Variational_String5_
                                                        + Variational_String6_  + Variational_String7_;
                                                variational_state_contruction_strings.push_back(Variational_String_);

                                            }
                                        }
                                    }

                                }
                            }
                        }}
                }}
        }


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on diagonal on sites-2,3 (/) [Next nearest neighbour AC singlets across holes]
        for(int bond0_type=0;bond0_type<4;bond0_type++){

            if((bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=2;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                bond3_type=3;
                Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                Variational_String3_ += " 2 4 ";

                for(bond4_type=0;bond4_type<4;bond4_type++){
                    if((bond4_type==0) || (bond4_type==2) ){
                        Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                        Variational_String4_ += " 4 5 ";

                        for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                            bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                            Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                            Variational_String5_ += " 5 6 ";

                            for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                Variational_String6_ += " 6 7 ";

                                for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                    bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                    Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                    Variational_String7_ += " 7 0";

                                    if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                            bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                        Variational_String_= Variational_String0_ + Variational_String1_ +
                                                Variational_String3_ + Variational_String4_ + Variational_String5_
                                                + Variational_String6_  + Variational_String7_;
                                        variational_state_contruction_strings.push_back(Variational_String_);

                                    }
                                }
                            }

                        }
                    }}
            }}


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on diagonal on sites-2,3 (\) [Next nearest neighbour AC singlets across holes]
        for(int bond0_type=0;bond0_type<4;bond0_type++){

            if((bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=3;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                bond3_type=2;
                Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                Variational_String3_ += " 2 4 ";

                for(bond4_type=0;bond4_type<4;bond4_type++){
                    if((bond4_type==1) || (bond4_type==3) ){
                        Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                        Variational_String4_ += " 4 5 ";

                        for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                            bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                            Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                            Variational_String5_ += " 5 6 ";

                            for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                Variational_String6_ += " 6 7 ";

                                for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                    bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                    Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                    Variational_String7_ += " 7 0";

                                    if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                            bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                        Variational_String_= Variational_String0_ + Variational_String1_ +
                                                Variational_String3_ + Variational_String4_ + Variational_String5_
                                                + Variational_String6_  + Variational_String7_;
                                        variational_state_contruction_strings.push_back(Variational_String_);

                                    }
                                }
                            }

                        }
                    }}
            }}


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on rung on site-2 [Diagonal next nearest neighbour singlet across the holes (/)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=0;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==0) || (bond3_type==2) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 6 ";

                                for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                    bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                    Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                    Variational_String6_ += " 6 7 ";

                                    for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                        bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                        Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                        Variational_String7_ += " 7 0";

                                        if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                            Variational_String_= Variational_String0_ + Variational_String1_ +
                                                    Variational_String3_ + Variational_String4_ + Variational_String5_
                                                    + Variational_String6_  + Variational_String7_;
                                            variational_state_contruction_strings.push_back(Variational_String_);

                                        }
                                    }
                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;



        //For holes sitting on rung on site-2 [Diagonal next nearest neighbour singlet across the holes (\)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=1;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==1) || (bond3_type==3) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 6 ";

                                for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                    bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                    Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                    Variational_String6_ += " 6 7 ";

                                    for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                        bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                        Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                        Variational_String7_ += " 7 0";

                                        if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                            Variational_String_= Variational_String0_ + Variational_String1_ +
                                                    Variational_String3_ + Variational_String4_ + Variational_String5_
                                                    + Variational_String6_  + Variational_String7_;
                                            variational_state_contruction_strings.push_back(Variational_String_);

                                        }
                                    }
                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;

        //For holes sitting on rung on site-2 [Along chain next nearest neighbour singlet across the holes (upper orb.)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=3;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==0) || (bond3_type==2) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 6 ";

                                for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                    bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                    Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                    Variational_String6_ += " 6 7 ";

                                    for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                        bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                        Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                        Variational_String7_ += " 7 0";

                                        if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                            Variational_String_= Variational_String0_ + Variational_String1_ +
                                                    Variational_String3_ + Variational_String4_ + Variational_String5_
                                                    + Variational_String6_  + Variational_String7_;
                                            variational_state_contruction_strings.push_back(Variational_String_);

                                        }
                                    }
                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on rung on site-2 [Along chain next nearest neighbour singlet across the holes (lower orb.)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=2;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==1) || (bond3_type==3) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 6 ";

                                for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                    bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                    Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                    Variational_String6_ += " 6 7 ";

                                    for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                        bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                        Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                        Variational_String7_ += " 7 0";

                                        if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                            Variational_String_= Variational_String0_ + Variational_String1_ +
                                                    Variational_String3_ + Variational_String4_ + Variational_String5_
                                                    + Variational_String6_  + Variational_String7_;
                                            variational_state_contruction_strings.push_back(Variational_String_);

                                        }
                                    }
                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


    }


    if(basis.Length==6){


        //For holes sitting on diagonal on sites-2,3(/) [Only nearest neighbour singlets]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                if(  (bond1_type==0) || (bond1_type==3) ){
                    Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                    Variational_String1_ += " 1 2 ";

                    for(bond3_type=0;bond3_type<4;bond3_type++){
                        if(  (bond3_type==0) || (bond3_type==2) ){
                            Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                            Variational_String3_ += " 3 4 ";

                            for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                                bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                                Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                                Variational_String4_ += " 4 5 ";

                                for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                    bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                    Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                    Variational_String5_ += " 5 0";

                                    if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                            bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                        Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                                Variational_String3_ + Variational_String4_ + Variational_String5_;
                                        variational_state_contruction_strings.push_back(Variational_String_);

                                    }

                                }
                            }
                        }}
                }}
        }

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;

        //For holes sitting on diagonal on sites-2,3(\) [Only nearest neighbour singlets]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                if(  (bond1_type==1) || (bond1_type==2) ){
                    Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                    Variational_String1_ += " 1 2 ";

                    for(bond3_type=0;bond3_type<4;bond3_type++){
                        if(  (bond3_type==3) || (bond3_type==1) ){
                            Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                            Variational_String3_ += " 3 4 ";

                            for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                                bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                                Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                                Variational_String4_ += " 4 5 ";

                                for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                    bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                    Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                    Variational_String5_ += " 5 0";

                                    if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                            bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                        Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                                Variational_String3_ + Variational_String4_ + Variational_String5_;
                                        variational_state_contruction_strings.push_back(Variational_String_);

                                    }

                                }
                            }
                        }}
                }}
        }


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on diagonal on sites-2,3 (/) [Next nearest neighbour AC singlets across holes]
        for(int bond0_type=0;bond0_type<4;bond0_type++){

            if((bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=2;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                bond3_type=3;
                Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                Variational_String3_ += " 2 4 ";

                for(bond4_type=0;bond4_type<4;bond4_type++){
                    if((bond4_type==0) || (bond4_type==2) ){
                        Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                        Variational_String4_ += " 4 5 ";

                        for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                            bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                            Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                            Variational_String5_ += " 5 0";

                            if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                    bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                        Variational_String3_ + Variational_String4_ + Variational_String5_;
                                variational_state_contruction_strings.push_back(Variational_String_);

                            }

                        }
                    }}
            }}


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on diagonal on sites-2,3 (\) [Next nearest neighbour AC singlets across holes]
        for(int bond0_type=0;bond0_type<4;bond0_type++){

            if((bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=3;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                bond3_type=2;
                Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                Variational_String3_ += " 2 4 ";

                for(bond4_type=0;bond4_type<4;bond4_type++){
                    if((bond4_type==1) || (bond4_type==3) ){
                        Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                        Variational_String4_ += " 4 5 ";

                        for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                            bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                            Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                            Variational_String5_ += " 5 0";

                            if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                    bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                        Variational_String3_ + Variational_String4_ + Variational_String5_;
                                variational_state_contruction_strings.push_back(Variational_String_);

                            }

                        }
                    }}
            }}


        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on rung on site-2 [Diagonal next nearest neighbour singlet across the holes (/)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=0;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==0) || (bond3_type==2) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 0";

                                if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                        bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                    Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                            Variational_String3_ + Variational_String4_ + Variational_String5_;
                                    variational_state_contruction_strings.push_back(Variational_String_);

                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;



        //For holes sitting on rung on site-2 [Diagonal next nearest neighbour singlet across the holes (\)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=1;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==1) || (bond3_type==3) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 0";

                                if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                        bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                    Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                            Variational_String3_ + Variational_String4_ + Variational_String5_;
                                    variational_state_contruction_strings.push_back(Variational_String_);

                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;

        //For holes sitting on rung on site-2 [Along chain next nearest neighbour singlet across the holes (upper orb.)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==1) || (bond0_type==2) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=3;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==0) || (bond3_type==2) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 0";

                                if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                        bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                    Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                            Variational_String3_ + Variational_String4_ + Variational_String5_;
                                    variational_state_contruction_strings.push_back(Variational_String_);

                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


        //For holes sitting on rung on site-2 [Along chain next nearest neighbour singlet across the holes (lower orb.)]
        for(int bond0_type=0;bond0_type<4;bond0_type++){
            if(  (bond0_type==0) || (bond0_type==3) ){
                Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
                Variational_String0_ += " 0 1 ";

                bond1_type=2;
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 3 ";

                for(bond3_type=0;bond3_type<4;bond3_type++){
                    if(  (bond3_type==1) || (bond3_type==3) ){
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 0";

                                if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                        bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                    Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                            Variational_String3_ + Variational_String4_ + Variational_String5_;
                                    variational_state_contruction_strings.push_back(Variational_String_);

                                }

                            }
                        }
                    }}
            }}

        cout<<"---------STATES UPTO NOW = "<<variational_state_contruction_strings.size()<<endl;


    }

    if(basis.Length==4 && false){

        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 2 ";

                for(int index2=0;index2<Allowed_bonds_left_constrained[bond1_type].size();index2++){
                    bond2_type=Allowed_bonds_left_constrained[bond1_type][index2];
                    Variational_String2_ =Bond_type_int_tag_to_string(bond2_type);
                    Variational_String2_ += " 2 3 ";

                    for(int index3=0;index3<Allowed_bonds_left_constrained[bond2_type].size();index3++){
                        bond3_type=Allowed_bonds_left_constrained[bond2_type][index3];
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 0";

                        if(bond0_type==Allowed_bonds_left_constrained[bond3_type][0] ||
                                bond0_type==Allowed_bonds_left_constrained[bond3_type][1]    ){
                            Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                    Variational_String3_;
                            variational_state_contruction_strings.push_back(Variational_String_);

                        }
                    }
                }
            }
        }
    }

}

void MODEL_2_orb_Hubb_chain_KSector::Create_variational_state_strings_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis){


    string Variational_String_;
    int bond1_type,bond2_type,bond3_type,bond4_type,bond5_type,bond6_type,bond7_type;
    string Variational_String0_, Variational_String1_, Variational_String2_, Variational_String3_;
    string Variational_String4_, Variational_String5_, Variational_String6_, Variational_String7_;
    Mat_2_int Allowed_bonds_left_constrained;
    Allowed_bonds_left_constrained.resize(4);

    //D S NOS1 (0)  --->D S NOS1 (0),  AC S NOS1 (2)
    Allowed_bonds_left_constrained[0].push_back(0);Allowed_bonds_left_constrained[0].push_back(2);

    //D S NOS2 (1)  --->D S NOS2 (1),  AC S NOS2 (3)
    Allowed_bonds_left_constrained[1].push_back(1);Allowed_bonds_left_constrained[1].push_back(3);

    //AC S NOS1 (2)  --->D S NOS2 (1),  AC S NOS2 (3)
    Allowed_bonds_left_constrained[2].push_back(1);Allowed_bonds_left_constrained[2].push_back(3);

    //AC S NOS2 (3)  --->D S NOS1 (0),  AC S NOS1 (2)
    Allowed_bonds_left_constrained[3].push_back(0);Allowed_bonds_left_constrained[3].push_back(2);


    if(basis.Length==8){

        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 2 ";

                for(int index2=0;index2<Allowed_bonds_left_constrained[bond1_type].size();index2++){
                    bond2_type=Allowed_bonds_left_constrained[bond1_type][index2];
                    Variational_String2_ =Bond_type_int_tag_to_string(bond2_type);
                    Variational_String2_ += " 2 3 ";

                    for(int index3=0;index3<Allowed_bonds_left_constrained[bond2_type].size();index3++){
                        bond3_type=Allowed_bonds_left_constrained[bond2_type][index3];
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 6 ";

                                for(int index6=0;index6<Allowed_bonds_left_constrained[bond5_type].size();index6++){
                                    bond6_type=Allowed_bonds_left_constrained[bond5_type][index6];
                                    Variational_String6_ =Bond_type_int_tag_to_string(bond6_type);
                                    Variational_String6_ += " 6 7 ";

                                    for(int index7=0;index7<Allowed_bonds_left_constrained[bond6_type].size();index7++){
                                        bond7_type=Allowed_bonds_left_constrained[bond6_type][index7];
                                        Variational_String7_ =Bond_type_int_tag_to_string(bond7_type);
                                        Variational_String7_ += " 7 0";

                                        if(bond0_type==Allowed_bonds_left_constrained[bond7_type][0] ||
                                                bond0_type==Allowed_bonds_left_constrained[bond7_type][1]    ){
                                            Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                                    Variational_String3_ + Variational_String4_ + Variational_String5_ +
                                                    Variational_String6_ + Variational_String7_;
                                            variational_state_contruction_strings.push_back(Variational_String_);

                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    if(basis.Length==6){

        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 2 ";

                for(int index2=0;index2<Allowed_bonds_left_constrained[bond1_type].size();index2++){
                    bond2_type=Allowed_bonds_left_constrained[bond1_type][index2];
                    Variational_String2_ =Bond_type_int_tag_to_string(bond2_type);
                    Variational_String2_ += " 2 3 ";

                    for(int index3=0;index3<Allowed_bonds_left_constrained[bond2_type].size();index3++){
                        bond3_type=Allowed_bonds_left_constrained[bond2_type][index3];
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 4 ";

                        for(int index4=0;index4<Allowed_bonds_left_constrained[bond3_type].size();index4++){
                            bond4_type=Allowed_bonds_left_constrained[bond3_type][index4];
                            Variational_String4_ =Bond_type_int_tag_to_string(bond4_type);
                            Variational_String4_ += " 4 5 ";

                            for(int index5=0;index5<Allowed_bonds_left_constrained[bond4_type].size();index5++){
                                bond5_type=Allowed_bonds_left_constrained[bond4_type][index5];
                                Variational_String5_ =Bond_type_int_tag_to_string(bond5_type);
                                Variational_String5_ += " 5 0";

                                if(bond0_type==Allowed_bonds_left_constrained[bond5_type][0] ||
                                        bond0_type==Allowed_bonds_left_constrained[bond5_type][1]    ){
                                    Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                            Variational_String3_ + Variational_String4_ + Variational_String5_;
                                    variational_state_contruction_strings.push_back(Variational_String_);

                                }

                            }
                        }
                    }
                }
            }
        }
    }

    if(basis.Length==4){

        for(int bond0_type=0;bond0_type<4;bond0_type++){
            Variational_String0_=Bond_type_int_tag_to_string(bond0_type);
            Variational_String0_ += " 0 1 ";

            for(int index1=0;index1<Allowed_bonds_left_constrained[bond0_type].size();index1++){
                bond1_type=Allowed_bonds_left_constrained[bond0_type][index1];
                Variational_String1_ =Bond_type_int_tag_to_string(bond1_type);
                Variational_String1_ += " 1 2 ";

                for(int index2=0;index2<Allowed_bonds_left_constrained[bond1_type].size();index2++){
                    bond2_type=Allowed_bonds_left_constrained[bond1_type][index2];
                    Variational_String2_ =Bond_type_int_tag_to_string(bond2_type);
                    Variational_String2_ += " 2 3 ";

                    for(int index3=0;index3<Allowed_bonds_left_constrained[bond2_type].size();index3++){
                        bond3_type=Allowed_bonds_left_constrained[bond2_type][index3];
                        Variational_String3_ =Bond_type_int_tag_to_string(bond3_type);
                        Variational_String3_ += " 3 0";

                        if(bond0_type==Allowed_bonds_left_constrained[bond3_type][0] ||
                                bond0_type==Allowed_bonds_left_constrained[bond3_type][1]    ){
                            Variational_String_= Variational_String0_ + Variational_String1_ + Variational_String2_ +
                                    Variational_String3_;
                            variational_state_contruction_strings.push_back(Variational_String_);

                        }
                    }
                }
            }
        }
    }

}

void MODEL_2_orb_Hubb_chain_KSector::Read_Anzatz_basis(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_){

    int No_of_basis;
    if(basis.Length==2){
        No_of_basis=2;
    }
    if(basis.Length==4){
        No_of_basis=16;
    }
    if(basis.Length==6){
        No_of_basis=64;
    }
    if(basis.Length==8){
        No_of_basis=256;
    }

    string state_type;
    BASIS_STATES_ANSATZ.clear();
    BASIS_STATES_ANSATZ.resize(No_of_basis);

    Ansatz_Basis_Overlap_with_GS.clear();
    Ansatz_Basis_Overlap_with_GS.resize(No_of_basis);

    int index;
    double value;

    ostringstream ss_length;
    ss_length << basis.Length;



    for(int i=0;i<No_of_basis;i++){
        ostringstream ss_int;
        ss_int << i;
        string sLine;

        ifstream infile_Variational_states;

        string file_in_Variational_states = "Variational_states_L"+ ss_length.str() +"/Variational_state" + ss_int.str() + ".txt";
        infile_Variational_states.open(file_in_Variational_states.c_str());
        getline(infile_Variational_states, state_type);


        BASIS_STATES_ANSATZ[i].clear();
        BASIS_STATES_ANSATZ[i].resize(basis.D_up_basis.size());

        while (!infile_Variational_states.eof())
        {
            getline(infile_Variational_states, sLine);
            stringstream ssLine;
            ssLine << sLine;
            ssLine>>index;
            ssLine>>value;

            BASIS_STATES_ANSATZ[i][index]=value;

        }


        infile_Variational_states.close();


        Ansatz_Basis_Overlap_with_GS[i]=dot_product(GS_,BASIS_STATES_ANSATZ[i]);
        cout<<i<<" ,State type = "<<state_type<<endl;
        cout<<Ansatz_Basis_Overlap_with_GS[i]<<endl;

    }

    //Get degenerate overlaps XXXXXXXXXXXXXXXXXXXXXXXXXXXX//
    bool present_before=false;
    double eps_=0.0001;
    Distinct_overlaps.clear();
    Distinct_overlaps.push_back(Ansatz_Basis_Overlap_with_GS[0]);
    for(int i=1;i<Ansatz_Basis_Overlap_with_GS.size();i++){

        present_before=false;
        for(int j=0;j<Distinct_overlaps.size();j++){
            if(abs(Ansatz_Basis_Overlap_with_GS[i] - Distinct_overlaps[j])<eps_){
                present_before=true;
            }
        }

        if(!present_before){
            Distinct_overlaps.push_back(Ansatz_Basis_Overlap_with_GS[i]);
        }

    }

    cout<<"Distinct overlapsXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout<<"-----------------------------------------------------"<<endl;
    Degenerate_states.clear();
    Degenerate_states.resize(Distinct_overlaps.size());
    int check_total_states=0;
    for(int i=0;i<Distinct_overlaps.size();i++){
        cout<<"Overlap = "<< Distinct_overlaps[i]<<", for VB states :";
        for(int j=0;j<Ansatz_Basis_Overlap_with_GS.size();j++){
            if(abs(Ansatz_Basis_Overlap_with_GS[j] - Distinct_overlaps[i])<eps_){
                Degenerate_states[i].push_back(j);
                cout<<"  "<<j;
            }
        }
        cout<<endl;
        check_total_states +=Degenerate_states[i].size();
    }

    cout<<"check_total_states = "<<check_total_states <<endl;
    assert(check_total_states==Ansatz_Basis_Overlap_with_GS.size());
    cout<<"-----------------------------------------------------"<<endl;
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

    //Get degenerate overlaps DONE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//



    cout<<"---------------------------------------------------"<<endl;


    string file_in_Overlap_bw_ansatzs = "Variational_states_L"+ ss_length.str() +"/Overlap_bw_ansatzs.txt";
    ifstream infile_Overlap_bw_ansatzs(file_in_Overlap_bw_ansatzs.c_str());

    cout<<"Overlaps b/w anzatzs states are being read from : "<<file_in_Overlap_bw_ansatzs<<endl;

    overlap_matrix_for_Anzatz_basis.resize(No_of_basis);
    for(int i=0;i<No_of_basis;i++){
        overlap_matrix_for_Anzatz_basis[i].resize(No_of_basis);
    }

    int temp_i, temp_j;
    double temp_double_;
    for(int i=0;i<No_of_basis;i++){
        for(int j=0;j<No_of_basis;j++){
            infile_Overlap_bw_ansatzs>>temp_i>>temp_j>>temp_double_;
            assert(i==temp_i);
            assert(j==temp_j);
            overlap_matrix_for_Anzatz_basis[temp_i][temp_j]=temp_double_;
        }
    }



}

void MODEL_2_orb_Hubb_chain_KSector::Get_overlap_matrix_for_Anzatz_basis(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &GS_){

    int No_of_basis;
    int No_of_pairs;

    if(true){
        variational_state_contruction_strings.clear();

        if( (basis.Ndn == basis.Nup)  && (basis.Ndn == basis.Length) ){
            Create_variational_state_strings_half_filling(basis);
        }
        if( (basis.Ndn == basis.Nup)  && (basis.Ndn == basis.Length-1) ){
            Create_variational_state_strings_2_holes(basis);
        }
        cout<<"Total No of RVB states used = "<<variational_state_contruction_strings.size()<<endl;
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
        for(int i=0;i<variational_state_contruction_strings.size();i++){
            cout<<variational_state_contruction_strings[i]<<endl;
        }
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
        No_of_basis=variational_state_contruction_strings.size();
    }


    if(false){
        No_of_basis=64;
        variational_state_contruction_strings.resize(No_of_basis);
        for(int i=0;i<64;i++){
            if(i==0){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==1){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==2){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==3){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==4){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==5){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==6){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==7){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==8){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==9){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==10){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==11){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==12){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==13){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==14){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==15){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==16){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==17){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==18){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==19){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==20){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==21){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==22){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==23){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==24){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==25){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==26){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==27){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==28){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==29){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==30){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==31){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==32){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==33){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==34){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==35){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==36){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==37){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==38){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==39){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==40){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==41){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==42){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==43){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==44){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==45){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==46){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==47){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==48){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==49){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==50){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==51){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 D S NOS1 4 5 D S NOS1 5 0";
            }
            if(i==52){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 AC S NOS2 4 5 AC S NOS1 5 0";
            }
            if(i==53){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 AC S NOS2 5 0";
            }
            if(i==54){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==55){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==56){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==57){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==58){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
            if(i==59){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 4 AC S NOS1 4 5 D S NOS2 5 0";
            }
            if(i==60){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 4 D S NOS1 4 5 AC S NOS1 5 0";
            }
            if(i==61){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 4 AC S NOS2 4 5 D S NOS1 5 0";
            }
            if(i==62){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 4 D S NOS2 4 5 AC S NOS2 5 0";
            }
            if(i==63){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 4 D S NOS2 4 5 D S NOS2 5 0";
            }
        }

    }


    if(false){
        No_of_basis=16;
        variational_state_contruction_strings.resize(No_of_basis);
        for(int i=0;i<16;i++){
            if(i==0){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 D S NOS1 2 3 D S NOS1 3 0";
            }
            if(i==1){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 D S NOS2 2 3 D S NOS2 3 0";
            }
            if(i==2){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 AC S NOS1 2 3 AC S NOS2 3 0";
            }
            if(i==3){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 AC S NOS2 2 3 AC S NOS1 3 0";
            }
            if(i==4){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 2 AC S NOS1 2 3 AC S NOS2 3 0";
            }
            if(i==5){
                variational_state_contruction_strings[i]="D S NOS2 0 1 D S NOS2 1 2 AC S NOS2 2 3 AC S NOS1 3 0";
            }
            if(i==6){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 D S NOS1 2 3 AC S NOS1 3 0";
            }
            if(i==7){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 D S NOS2 2 3 AC S NOS2 3 0";
            }
            if(i==8){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 2 D S NOS1 2 3 D S NOS1 3 0";
            }
            if(i==9){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 AC S NOS1 1 2 D S NOS2 2 3 D S NOS2 3 0";
            }
            if(i==10){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 AC S NOS2 2 3 D S NOS1 3 0";
            }
            if(i==11){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 AC S NOS1 2 3 D S NOS2 3 0";
            }
            if(i==12){
                variational_state_contruction_strings[i]="D S NOS1 0 1 AC S NOS1 1 2 D S NOS2 2 3 AC S NOS2 3 0";
            }
            if(i==13){
                variational_state_contruction_strings[i]="D S NOS2 0 1 AC S NOS2 1 2 D S NOS1 2 3 AC S NOS1 3 0";
            }
            if(i==14){
                variational_state_contruction_strings[i]="AC S NOS2 0 1 D S NOS1 1 2 AC S NOS1 2 3 D S NOS2 3 0";
            }
            if(i==15){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 D S NOS2 1 2 AC S NOS2 2 3 D S NOS1 3 0";
            }
        }

    }

    if(basis.Length==2){
        No_of_basis=2;
        variational_state_contruction_strings.resize(No_of_basis);

        for(int i=0;i<2;i++){
            if(i==0){
                variational_state_contruction_strings[i]="AC S NOS1 0 1 AC S NOS2 1 0";
            }
            if(i==1){
                variational_state_contruction_strings[i]="D S NOS1 0 1 D S NOS1 1 0";
            }

        }

    }


    int temp_int;

    No_of_pairs=basis.Ndn;assert(basis.Ndn==basis.Nup);

    BASIS_STATES_ANSATZ.resize(No_of_basis);

    string variational_state_contruction_;

    cout<<"-------OVERLAPS OF ANSATZ BASIS WITH GS--------------"<<endl;
    Variational_state_pair_Geometry.resize(basis.Length);
    Variational_state_pair_spin_symmetry.resize(basis.Length);
    Variational_state_pair_orbital_symmetry.resize(basis.Length);
    Variational_state_pair_sites.resize(basis.Length);


    Ansatz_Basis_Overlap_with_GS.resize(No_of_basis);
    for(int i=0;i<No_of_basis;i++){
        ostringstream ss_int;
        ss_int << i;

        string file_out_Variational_states = "Variational_state" + ss_int.str() + ".txt";
        ofstream outfile_Variational_states(file_out_Variational_states.c_str());

        variational_state_contruction_=variational_state_contruction_strings[i];


        //------------UPDATING Variational_state_pair******--------

        stringstream variational_state_contruction_stream(variational_state_contruction_);
        for(int n=0;n<No_of_pairs;n++){
            variational_state_contruction_stream >> Variational_state_pair_Geometry[n];
            assert(Variational_state_pair_Geometry[n]=="D" ||
                   Variational_state_pair_Geometry[n]=="AC");

            variational_state_contruction_stream >> Variational_state_pair_spin_symmetry[n];
            assert(Variational_state_pair_spin_symmetry[n]=="S" ||
                   Variational_state_pair_spin_symmetry[n]=="T");

            variational_state_contruction_stream >> Variational_state_pair_orbital_symmetry[n];
            assert(Variational_state_pair_orbital_symmetry[n] =="OA" ||
                   Variational_state_pair_orbital_symmetry[n] =="OS" ||
                   Variational_state_pair_orbital_symmetry[n] =="NOS1" ||
                   Variational_state_pair_orbital_symmetry[n] =="NOS2");

            variational_state_contruction_stream >> temp_int;
            Variational_state_pair_sites[n].first =temp_int;
            assert(temp_int>=0 && temp_int<basis.Length);


            variational_state_contruction_stream >> temp_int;
            Variational_state_pair_sites[n].second =temp_int;
            assert(temp_int>=0 && temp_int<basis.Length);


        }

        //-----------------------------




        Get_Variational_State(basis,No_of_pairs);
        BASIS_STATES_ANSATZ[i]=State_;

        Ansatz_Basis_Overlap_with_GS[i]=dot_product(GS_,BASIS_STATES_ANSATZ[i]);
        cout<<i<<" ,State type = "<<variational_state_contruction_<<endl;
        cout<<Ansatz_Basis_Overlap_with_GS[i]<<endl;

        outfile_Variational_states<<"#VARIATIONAL STATE NO. "<<i<<",  "<<variational_state_contruction_<<endl;
        for(int j=0;j<BASIS_STATES_ANSATZ[i].size();j++){
            if(BASIS_STATES_ANSATZ[i][j]!=0.0){
                outfile_Variational_states<<j<<"\t"<<BASIS_STATES_ANSATZ[i][j]<<endl;
            }
        }


    }


    //Get degenerate overlaps XXXXXXXXXXXXXXXXXXXXXXXXXXXX//
    bool present_before=false;
    double eps_=0.000001;
    Distinct_overlaps.clear();
    Distinct_overlaps.push_back(Ansatz_Basis_Overlap_with_GS[0]);
    for(int i=1;i<Ansatz_Basis_Overlap_with_GS.size();i++){

        present_before=false;
        for(int j=0;j<Distinct_overlaps.size();j++){
            if(abs(Ansatz_Basis_Overlap_with_GS[i] - Distinct_overlaps[j])<eps_){
                present_before=true;
            }
        }

        if(!present_before){
            Distinct_overlaps.push_back(Ansatz_Basis_Overlap_with_GS[i]);
        }

    }

    cout<<"Distinct overlapsXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout<<"-----------------------------------------------------"<<endl;
    Degenerate_states.clear();
    Degenerate_states.resize(Distinct_overlaps.size());
    int check_total_states=0;
    for(int i=0;i<Distinct_overlaps.size();i++){
        cout<<"Overlap = "<< Distinct_overlaps[i]<<", for VB states :";
        for(int j=0;j<Ansatz_Basis_Overlap_with_GS.size();j++){
            if(abs(Ansatz_Basis_Overlap_with_GS[j] - Distinct_overlaps[i])<eps_){
                Degenerate_states[i].push_back(j);
                cout<<"  "<<j;
            }
        }
        cout<<endl;
        check_total_states +=Degenerate_states[i].size();
    }

    cout<<"check_total_states = "<<check_total_states <<endl;
    assert(check_total_states==Ansatz_Basis_Overlap_with_GS.size());
    cout<<"-----------------------------------------------------"<<endl;
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

    //Get degenerate overlaps DONE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//



    cout<<"---------------------------------------------------"<<endl;



    string file_out_Overlap_bw_ansatzs = "Overlap_bw_ansatzs.txt";
    ofstream outfile_Overlap_bw_ansatzs(file_out_Overlap_bw_ansatzs.c_str());

    overlap_matrix_for_Anzatz_basis.resize(No_of_basis);
    for(int i=0;i<No_of_basis;i++){
        overlap_matrix_for_Anzatz_basis[i].resize(No_of_basis);
        for(int j=0;j<No_of_basis;j++){
            overlap_matrix_for_Anzatz_basis[i][j]=dot_product(BASIS_STATES_ANSATZ[i],BASIS_STATES_ANSATZ[j]);
            outfile_Overlap_bw_ansatzs<<i<<"   "<<j<<"   "<<overlap_matrix_for_Anzatz_basis[i][j]<<endl;
        }
    }



    //JUST A CHECK-------------[REMOVE LATER]
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    /* Mat_1_doub vec_temp_28, vec_temp_29, vec_temp;
    double doub_temp;
    vec_temp_28 = BASIS_STATES_ANSATZ[28];
    vec_temp_29 = BASIS_STATES_ANSATZ[29];
    vec_temp=Act_Reflection_about_Central_site(basis,  vec_temp_28);
    doub_temp = dot_product(vec_temp_29, vec_temp);
    cout<<"<29|Ref.|28> = "<<doub_temp<<endl;*/

    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    //----------------------------


}

void MODEL_2_orb_Hubb_chain_KSector::Variational_state_optimization(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_){

    double dp=0.05;
    double eps_= 0.000001;
    double_type VALUE_CHECK, VALUE_CHECK_OLD;
    double_type overlap_with_GS;
    double Norm_Ansatz;
    double_type Norm_Ansatz_temp;
    int VB_state;
    int VB_state1, VB_state2;

    Mat_2_doub Combined_Basis;
    Combined_Basis.resize(Distinct_overlaps.size());


    for(int i=0;i<Distinct_overlaps.size();i++){
        Combined_Basis[i].resize(BASIS_STATES_ANSATZ[0].size());
        for(int j=0;j<Degenerate_states[i].size();j++){
            VB_state=Degenerate_states[i][j];

            for(int n=0;n<BASIS_STATES_ANSATZ[0].size();n++){
                Combined_Basis[i][n] += BASIS_STATES_ANSATZ[VB_state][n];
            }
        }
    }

    Mat_2_doub overlap_matrix_for_Combined_basis;
    overlap_matrix_for_Combined_basis.resize(Distinct_overlaps.size());
    for(int i=0;i<Distinct_overlaps.size();i++){
        overlap_matrix_for_Combined_basis[i].resize(Distinct_overlaps.size());
    }

    for(int i=0;i<Distinct_overlaps.size();i++){
        for(int j=0;j<Distinct_overlaps.size();j++){
            overlap_matrix_for_Combined_basis[i][j]=0.0;
            for(int l=0;l<Degenerate_states[i].size();l++){
                for(int m=0;m<Degenerate_states[j].size();m++){
                    VB_state1=Degenerate_states[i][l];
                    VB_state2=Degenerate_states[j][m];
                    overlap_matrix_for_Combined_basis[i][j] += overlap_matrix_for_Anzatz_basis[VB_state1][VB_state2];
                }
            }
        }
    }


    if(basis.Length==6 || basis.Length==8 || basis.Length==4){

        vector<double_type> alpha_final;
        vector<double_type> alpha;
        vector<double_type> alpha_min, alpha_max;
        alpha.resize(Distinct_overlaps.size());
        alpha_final.resize(Distinct_overlaps.size());
        alpha_min.resize(Distinct_overlaps.size()); alpha_max.resize(Distinct_overlaps.size());
        VALUE_CHECK_OLD=zero;

        /*
alpha[0] = -1.1833269165e-02
alpha[1] = 4.3848768798e-02
alpha[2] = -4.3848768798e-02
alpha[3] = 4.3848768798e-02
alpha[4] = 4.3848768798e-02
alpha[5] = -4.3848768798e-02
alpha[6] = 4.3848768798e-02
alpha[7] = -1.1833269165e-02
overlap_with_GS = 6.1889488366e-01
         */

        /*
alpha[0] = -1.2526611734e-02
alpha[1] = 3.5832057778e-02
alpha[2] = -4.1125021016e-02
alpha[3] = 5.7003910729e-02
alpha[4] = 3.5832057778e-02
alpha[5] = -4.1125021016e-02
alpha[6] = 4.6417984253e-02
alpha[7] = -2.3112538210e-02
overlap with GS = 6.1985191582e-01
         */

        for(int i=0;i<Distinct_overlaps.size();i++){
            alpha_min[i] = Distinct_overlaps[i]- (one*0.0);
            alpha_max[i] = Distinct_overlaps[i]+ (one*0.0);
            alpha_final[i] = Distinct_overlaps[i];
        }
        //  alpha_min[0] = Distinct_overlaps[0]-0.1;alpha_max[0] = Distinct_overlaps[0]+0.1;
        //  alpha_min[1] = Distinct_overlaps[1]-0.1;alpha_max[1] = Distinct_overlaps[1]+0.1;
        //  alpha_min[2] = Distinct_overlaps[2]-0.1;alpha_max[2] = Distinct_overlaps[2]+0.1;
        //  alpha_min[3] = Distinct_overlaps[3]-0.1;alpha_max[3] = Distinct_overlaps[3]+0.1;
        //  alpha_min[4] = Distinct_overlaps[4]-0.1;alpha_max[4] = Distinct_overlaps[4]+0.1;


        //cout<<"here 1"<<endl;

        //        if(basis.Length==8){
        //            alpha[0]=alpha_min[0];
        //            while(alpha[0]<=alpha_max[0]){

        //                alpha[1]=alpha_min[1];
        //                //alpha[1]=-1.0*alpha[0];
        //                while(alpha[1]<=alpha_max[1]){

        //                    alpha[2]=alpha_min[2];
        //                    //alpha[2]=alpha[0];
        //                    while(alpha[2]<=alpha_max[2]){

        //                        alpha[3]=alpha_min[3];
        //                        //alpha[3]=-1.0*alpha[0];
        //                        while(alpha[3]<=alpha_max[3]){

        //                            alpha[4]=alpha_min[4];
        //                            //alpha[4]=-1.0*alpha[0];
        //                            while(alpha[4]<=alpha_max[4]){

        //                                alpha[5]=alpha_min[5];
        //                                //alpha[5]=alpha[0];
        //                                while(alpha[5]<=alpha_max[5]){

        //                                    alpha[6]=alpha_min[6];
        //                                    //alpha[6]=-1.0*alpha[0];
        //                                    while(alpha[6]<=alpha_max[6]){

        //                                        alpha[7]=alpha_min[7];
        //                                        //alpha[7]=alpha[0];
        //                                        while(alpha[7]<=alpha_max[7]){

        //                                            alpha[8]=alpha_min[8];
        //                                            while(alpha[8]<=alpha_max[8]){

        //                                                alpha[9]=alpha_min[9];
        //                                                //alpha[9]=-1.0*alpha[0];
        //                                                while(alpha[9]<=alpha_max[9]){

        //                                                    alpha[10]=alpha_min[10];
        //                                                    //alpha[10]=alpha[0];
        //                                                    while(alpha[10]<=alpha_max[10]){

        //                                                        alpha[11]=alpha_min[11];
        //                                                        //alpha[11]=-1.0*alpha[0];
        //                                                        while(alpha[11]<=alpha_max[11]){

        //                                                            alpha[12]=alpha_min[12];
        //                                                            //alpha[12]=-1.0*alpha[0];
        //                                                            while(alpha[12]<=alpha_max[12]){

        //                                                                alpha[13]=alpha_min[13];
        //                                                                //alpha[13]=alpha[0];
        //                                                                while(alpha[13]<=alpha_max[13]){

        //                                                                    alpha[14]=alpha_min[14];
        //                                                                    //alpha[14]=-1.0*alpha[0];
        //                                                                    while(alpha[14]<=alpha_max[14]){

        //                                                                        alpha[15]=alpha_min[15];
        //                                                                        //alpha[15]=alpha[0];
        //                                                                        while(alpha[15]<=alpha_max[15]){

        //                                                                            alpha[16]=alpha_min[16];
        //                                                                            //alpha[16]=-1.0*alpha[0];
        //                                                                            while(alpha[16]<=alpha_max[16]){

        //                                                                                alpha[17]=alpha_min[17];
        //                                                                                //alpha[17]=alpha[0];
        //                                                                                while(alpha[17]<=alpha_max[17]){

        //                                                                                    //cout<<"here 2"<<endl;

        //                                                                                    //-------OLD COSTLY WAY------------//
        //                                                                                    /* State_.clear();
        //                                                                                    State_.resize(BASIS_STATES_ANSATZ[0].size());

        //                                                                                    for(int m=0;m<Distinct_overlaps.size();m++){
        //                                                                                        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
        //                                                                                            State_[i] += alpha[m]*(Combined_Basis[m][i]);
        //                                                                                        }
        //                                                                                    }

        //                                                                                    Norm_Ansatz=dot_product(State_,State_);
        //                                                                                    if(Norm_Ansatz <= eps_){
        //                                                                                        VALUE_CHECK=0.0;
        //                                                                                    }
        //                                                                                    else{
        //                                                                                        overlap_with_GS=dot_product(GS_,State_);
        //                                                                                        overlap_with_GS = overlap_with_GS*overlap_with_GS;
        //                                                                                        VALUE_CHECK=overlap_with_GS/Norm_Ansatz;

        //                                                                                    }*/
        //                                                                                    //---------------------------------//

        //                                                                                    //-------NEW FASTER WAY-----------//
        //                                                                                    Norm_Ansatz=0.0;
        //                                                                                    for(int i_=0;i_<Distinct_overlaps.size();i_++){
        //                                                                                        for(int j_=0;j_<Distinct_overlaps.size();j_++){
        //                                                                                            Norm_Ansatz +=alpha[i_]*alpha[j_]*overlap_matrix_for_Combined_basis[i_][j_];
        //                                                                                        }
        //                                                                                    }
        //                                                                                    if(Norm_Ansatz <= eps_){
        //                                                                                        VALUE_CHECK=0.0;
        //                                                                                    }
        //                                                                                    else{
        //                                                                                        overlap_with_GS=0.0;
        //                                                                                        for(int i_=0;i_<Distinct_overlaps.size();i_++){
        //                                                                                            overlap_with_GS +=  alpha[i_]*Distinct_overlaps[i_]*Degenerate_states[i_].size();
        //                                                                                        }
        //                                                                                        VALUE_CHECK=(overlap_with_GS*overlap_with_GS)/Norm_Ansatz;
        //                                                                                    }

        //                                                                                    //--------------------------------//



        //                                                                                    if( VALUE_CHECK > VALUE_CHECK_OLD){
        //                                                                                        VALUE_CHECK_OLD=VALUE_CHECK;
        //                                                                                        alpha_final=alpha;
        //                                                                                    }

        //                                                                                    for(int n=0;n<Distinct_overlaps.size();n++){
        //                                                                                        cout<<"alpha["<<n<<"] = "<<alpha[n]<<"    ";
        //                                                                                    }
        //                                                                                    cout<< "is done"<<endl;


        //                                                                                    //
        //                                                                                    alpha[17]+=dp*1;
        //                                                                                }
        //                                                                                alpha[16]+=dp*1;
        //                                                                            }
        //                                                                            alpha[15]+=dp*1;
        //                                                                        }
        //                                                                        alpha[14]+=dp*1;
        //                                                                    }
        //                                                                    alpha[13]+=dp*1;
        //                                                                }
        //                                                                alpha[12]+=dp*1;
        //                                                            }
        //                                                            alpha[11]+=dp*1;
        //                                                        }
        //                                                        alpha[10]+=dp*1;
        //                                                    }
        //                                                    alpha[9]+=dp*1;
        //                                                }
        //                                                alpha[8]+=dp;
        //                                            }
        //                                            alpha[7]+=dp*1;
        //                                        }
        //                                        alpha[6]+=dp*1;
        //                                    }
        //                                    alpha[5]+=dp*1;
        //                                }
        //                                alpha[4]+=dp*1;
        //                            }
        //                            alpha[3]+=dp*1;
        //                        }
        //                        alpha[2]+=dp*1;
        //                    }
        //                    alpha[1]+=dp*1;
        //                }
        //                alpha[0]+=dp;
        //            }
        //        }


        if(true){ //if State_ is needed
            State_.clear();
            State_.resize(BASIS_STATES_ANSATZ[0].size());

            for(int m=0;m<Distinct_overlaps.size();m++){
                for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
                    State_[i] += alpha_final[m]*(Combined_Basis[m][i]);
                }
            }
            Norm_Ansatz=abs(dot_product(State_,State_));
            for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
                State_[i] = State_[i]*(1.0/sqrt(Norm_Ansatz));
            }
            overlap_with_GS=dot_product(GS_,State_);
        }


        Norm_Ansatz_temp=zero;
        for(int i_=0;i_<Distinct_overlaps.size();i_++){
            for(int j_=0;j_<Distinct_overlaps.size();j_++){
                Norm_Ansatz_temp += alpha_final[i_]*alpha_final[j_]*overlap_matrix_for_Combined_basis[i_][j_];
            }
        }
        Norm_Ansatz=abs(Norm_Ansatz_temp);

        overlap_with_GS=zero;
        for(int i_=0;i_<Distinct_overlaps.size();i_++){
            overlap_with_GS +=  alpha_final[i_]*Distinct_overlaps[i_]*(1.0*Degenerate_states[i_].size());
        }


        cout<<"XXXXXXXX-----VARIATIONAL STATE OPTIMIZATION-------XXXXXXXXXX"<<endl;
        for(int n=0;n<Distinct_overlaps.size();n++){
            cout<<"alpha["<<n<<"] = "<<alpha_final[n]/sqrt(Norm_Ansatz)<<endl;
        }
        cout<<"(overlap_with_GS)/sqrt(Norm_Ansatz) = "<<overlap_with_GS*sqrt(1.0/Norm_Ansatz)<<endl;
        cout<<"Norm_Ansatz = "<<Norm_Ansatz<<endl;
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

    }



    if(basis.Length==5){
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
        //        beta_min=dot_product(GS_,BASIS_STATES_ANSATZ[11]); beta_max=dot_product(GS_,BASIS_STATES_ANSATZ[11]);
        //        gamma_min=dot_product(GS_,BASIS_STATES_ANSATZ[1]); gamma_max=dot_product(GS_,BASIS_STATES_ANSATZ[1]);
        //        delta_min=dot_product(GS_,BASIS_STATES_ANSATZ[2]); delta_max=dot_product(GS_,BASIS_STATES_ANSATZ[2]);


        //        /*
        //        alpha_min=9.9621e-02; alpha_max=9.9621e-02;
        //        beta_min=9.3288e-02; beta_max=9.3288e-02;
        //        gamma_min=-7.2306e-02; gamma_max=-7.2306e-02;
        //        delta_min=7.2054e-02; delta_max=7.2054e-02;
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

        //        for(alpha=alpha_min-0.1;alpha<=alpha_max+0.1;alpha=alpha+dp){
        //            for(beta=beta_min-0.1;beta<=beta_max+0.1;beta=beta+dp){
        //                for(gamma=gamma_min-0.1;gamma<=gamma_max+0.1;gamma=gamma+dp){
        //                    for(delta=delta_min-0.1;delta<=delta_max+0.1;delta=delta+dp){


        //                        for(int i=0;i<BASIS_STATES_ANSATZ[0].size();i++){
        //                            State_[i] = alpha*(BASIS_STATES_ANSATZ[0][i] + BASIS_STATES_ANSATZ[4][i]) +
        //                                    beta*(BASIS_STATES_ANSATZ[11][i] + BASIS_STATES_ANSATZ[15][i]) +
        //                                    gamma*(BASIS_STATES_ANSATZ[1][i] + BASIS_STATES_ANSATZ[3][i]
        //                                    + BASIS_STATES_ANSATZ[5][i] + BASIS_STATES_ANSATZ[7][i]
        //                                    + BASIS_STATES_ANSATZ[8][i] + BASIS_STATES_ANSATZ[10][i]
        //                                    + BASIS_STATES_ANSATZ[12][i] + BASIS_STATES_ANSATZ[14][i]) +
        //                                    delta*(BASIS_STATES_ANSATZ[2][i] + BASIS_STATES_ANSATZ[6][i]
        //                                    + BASIS_STATES_ANSATZ[9][i] + BASIS_STATES_ANSATZ[13][i]);
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
        //            State_[i] = alpha_final*(BASIS_STATES_ANSATZ[0][i] + BASIS_STATES_ANSATZ[4][i]) +
        //                    beta_final*(BASIS_STATES_ANSATZ[11][i] + BASIS_STATES_ANSATZ[15][i]) +
        //                    gamma_final*(BASIS_STATES_ANSATZ[1][i] + BASIS_STATES_ANSATZ[3][i]
        //                    + BASIS_STATES_ANSATZ[5][i] + BASIS_STATES_ANSATZ[7][i]
        //                    + BASIS_STATES_ANSATZ[8][i] + BASIS_STATES_ANSATZ[10][i]
        //                    + BASIS_STATES_ANSATZ[12][i] + BASIS_STATES_ANSATZ[14][i]) +
        //                    delta_final*(BASIS_STATES_ANSATZ[2][i] + BASIS_STATES_ANSATZ[6][i]
        //                    + BASIS_STATES_ANSATZ[9][i] + BASIS_STATES_ANSATZ[13][i]);
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



    }



    if(basis.Length==2){
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



    }

}


void MODEL_2_orb_Hubb_chain_KSector::Variational_state_optimization_old(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_){

    //    double dp=0.003;
    //    double eps_= 0.000001;
    //    double VALUE_CHECK, VALUE_CHECK_OLD;
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
    //        alpha_min=9.9621e-02; alpha_max=9.9621e-02;
    //        beta_min=9.3288e-02; beta_max=9.3288e-02;
    //        gamma_min=-7.2306e-02; gamma_max=-7.2306e-02;
    //        delta_min=7.2054e-02; delta_max=7.2054e-02;
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

void MODEL_2_orb_Hubb_chain_KSector::Read_Mat_2_trio(Mat_2_trio_int &MAT_TEMP, Mat_1_doub &VALUES_TEMP,
                                                     int pair_no){
    //Already having superpositions; see notes


    //Operator = VALUE_TEMP[0]*MAT_TEMP[0][0]*MAT_TEMP[0][1] +
    //           VALUE_TEMP[1]*MAT_TEMP[1][0]*MAT_TEMP[1][1] +
    //           VALUE_TEMP[2]*MAT_TEMP[2][0]*MAT_TEMP[2][1] +
    //           VALUE_TEMP[3]*MAT_TEMP[3][0]*MAT_TEMP[3][1]


    /*where,
    MAT_TEMP[0][0]*MAT_TEMP[0][1] = c^{\dagger}_{a,i,up}  c^{\dagger}_{b,i+1,dn}
    etc..
     */


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
    else if(Variational_state_pair_Geometry[pair_no]=="AC"){

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

    else {
        assert(Variational_state_pair_Geometry[pair_no]=="AR");
        assert(Variational_state_pair_spin_symmetry[pair_no]=="T");
        assert(Variational_state_pair_orbital_symmetry[pair_no]=="NOS");
        assert(Variational_state_pair_sites[pair_no].first==pair_no);
        assert(Variational_state_pair_sites[pair_no].second==pair_no);

        //1.0*c^{\dagger}_{a,i,up}  c^{\dagger}_{b,i,dn}*1.0
        MAT_TEMP[0][0].orb_=0;MAT_TEMP[0][0].spin_=0;
        MAT_TEMP[0][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[0][1].orb_=1;MAT_TEMP[0][1].spin_=1;
        MAT_TEMP[0][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[0]=1.0;

        // +1.0*c^{\dagger}_{a,i,dn}  c^{\dagger}_{b,i,up}
        MAT_TEMP[1][0].orb_=0;MAT_TEMP[1][0].spin_=1;
        MAT_TEMP[1][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[1][1].orb_=1;MAT_TEMP[1][1].spin_=0;
        MAT_TEMP[1][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[1]=1.0;

        //For next two terms [These are actually zeros]
        //c^{\dagger}_{b,i,up}  c^{\dagger}_{a,i+1,dn}
        MAT_TEMP[2][0].orb_=1;MAT_TEMP[2][0].spin_=0;
        MAT_TEMP[2][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[2][1].orb_=0;MAT_TEMP[2][1].spin_=1;
        MAT_TEMP[2][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[2]=0.0;

        // c^{\dagger}_{b,i,dn}  c^{\dagger}_{a,i+1,up}
        MAT_TEMP[3][0].orb_=1;MAT_TEMP[3][0].spin_=1;
        MAT_TEMP[3][0].site_=Variational_state_pair_sites[pair_no].first;
        MAT_TEMP[3][1].orb_=0;MAT_TEMP[3][1].spin_=0;
        MAT_TEMP[3][1].site_=Variational_state_pair_sites[pair_no].second;
        VALUES_TEMP[3]=0.0;

    }

    //---------------------------------------//
    //--------------------------------------//



}

void MODEL_2_orb_Hubb_chain_KSector::Get_Variational_State(BASIS_2_orb_Hubb_chain_KSector &basis, int no_of_pairs){


    Mat_2_trio_int MAT_FINAL;
    Mat_1_doub VALUES_FINAL;

    Mat_2_trio_int MAT_TEMP;
    Mat_1_doub VALUES_TEMP;

    Read_Mat_2_trio(MAT_FINAL, VALUES_FINAL,0);

    for(int i=1;i<no_of_pairs;i++){
        Read_Mat_2_trio(MAT_TEMP, VALUES_TEMP,i);

        Direct_product_of_Mat_2_trio_int(MAT_FINAL, VALUES_FINAL,
                                         MAT_TEMP, VALUES_TEMP,
                                         MAT_FINAL, VALUES_FINAL);
    }


    Mat_1_doub State_temp;
    State_.resize(basis.D_dn_basis.size());
    State_temp.resize(basis.D_dn_basis.size());
    for(int i=0;i<State_.size();i++){
        State_[i]=0.0;
    }


    for(int i=0;i<MAT_FINAL.size();i++){
        if(VALUES_FINAL[i]!=0.0){

            Get_State_by_acting_creation_operators(basis, MAT_FINAL[i],State_temp);

            //Add VALUE_TEMP[i]*State_temp in State_
            Subtract( State_, -1.0*VALUES_FINAL[i], State_temp);
        }

    }

    //Normalize State_;
    double_type tmpnrm_type_double;
    double tmpnrm;
    tmpnrm_type_double=dot_product(State_,State_);
    tmpnrm=sqrt(abs(tmpnrm_type_double));
    if(tmpnrm !=0){
        for(int i=0;i<State_.size();i++){
            State_[i] = (State_[i]*(1.0/tmpnrm));
        }
    }


}


void MODEL_2_orb_Hubb_chain_KSector::Get_State_by_acting_creation_operators(BASIS_2_orb_Hubb_chain_KSector &basis,
                                                                            Mat_1_trio_int MAT_,
                                                                            Mat_1_doub &State_temp){


    double FM_sign_up_down, FM_sign_up, FM_sign_down;
    int up_down_int, up_int, down_int;
    int offset;

    int Nup_check, Ndown_check;
    int Dup_temp, Ddn_temp;

    bool up_repeat, down_repeat;

    State_temp.resize(basis.D_up_basis.size());
    for(int i=0;i<State_temp.size();i++){
        State_temp[i]=0.0;
    }
    Mat_1_int positions_up;
    Mat_1_int positions_down;


    positions_up.clear();
    positions_down.clear();

    up_down_int=0;
    offset=0;
    Nup_check=0;
    Ndown_check=0;
    for(int i=0;i<MAT_.size();i++){

        //up
        if(MAT_[i].spin_==0){
            positions_up.push_back( MAT_[i].site_ + (basis.Length*MAT_[i].orb_) );
            offset++;
            Nup_check++;
        }
        else{
            //down
            assert(MAT_[i].spin_==1);
            positions_down.push_back( MAT_[i].site_ + (basis.Length*MAT_[i].orb_) );

            up_down_int += offset;
            Ndown_check++;
        }


    }


    up_repeat=false;
    for (int i = 0; i < positions_up.size() - 1; i++){
        for (int j = i + 1;j < positions_up.size(); j++){
            if (positions_up[i] == positions_up[j]){
                up_repeat=true;
            } // then this is a duplicate
        }
    }

    down_repeat=false;
    for (int i = 0; i < positions_down.size() - 1; i++){
        for (int j = i + 1;j < positions_down.size(); j++){
            if (positions_down[i] == positions_down[j]){
                down_repeat=true;
            } // then this is a duplicate
        }
    }



    if(up_repeat==false &&
            down_repeat==false){
        assert(Nup_check == basis.Nup);
        assert(Ndown_check == basis.Ndn);

        FM_sign_up_down = pow(-1.0, 1.0*up_down_int);

        up_int=minSwaps(positions_up, positions_up.size());
        FM_sign_up=pow(-1.0, 1.0*up_int);

        down_int=minSwaps(positions_down, positions_down.size());
        FM_sign_down=pow(-1.0, 1.0*down_int);

        Dup_temp=0;
        for(int i=0;i<positions_up.size();i++){
            Dup_temp +=pow(2,positions_up[i]);
        }

        Ddn_temp=0;
        for(int i=0;i<positions_down.size();i++){
            Ddn_temp +=pow(2,positions_down[i]);
        }


        int i_new, m_new;
        int D_up,D_dn;
        D_up=Dup_temp;
        D_dn=Ddn_temp;

        int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
        int sign_pow_up, sign_pow_dn;
        int range_min, range_max;
        int Inv_Trnsltns_;

        int D_up_temp=D_up;
        int D_dn_temp=D_dn;
        bool row_found_=false;

        sign_pow_up=0;sign_pow_dn=0;
        double sign_up_trans, sign_dn_trans;

        for(int inv_trnsltns=0;inv_trnsltns<basis.Length;inv_trnsltns++){

            if(inv_trnsltns>0){
                //Inv Translation on orb-0,spin_dn
                sign_pow_dn_orb0 = one_bits_in_bw(0, basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,0);
                if(bit_value(D_dn_temp,basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb0;
                }
                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_dn
                sign_pow_dn_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_dn_temp) +
                        1*bit_value(D_dn_temp,basis.Length);
                if(bit_value(D_dn_temp,2*basis.Length -1)==1){
                    sign_pow_dn += 1*sign_pow_dn_orb1;
                }

                D_dn_temp = Act_Translation_assuming_PBC(D_dn_temp,basis.Length,basis.n_orb*basis.Length-1);


                //Inv Translation on orb-0,spin_up
                sign_pow_up_orb0 = one_bits_in_bw(0, basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,0);
                if(bit_value(D_up_temp,basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb0;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,0,basis.Length-1);


                //Inv Translation on orb-1,spin_up
                sign_pow_up_orb1 = one_bits_in_bw(basis.Length, 2*basis.Length -1, D_up_temp) +
                        1*bit_value(D_up_temp,basis.Length);
                if(bit_value(D_up_temp,2*basis.Length -1)==1){
                    sign_pow_up += 1*sign_pow_up_orb1;
                }
                D_up_temp = Act_Translation_assuming_PBC(D_up_temp,basis.Length,basis.n_orb*basis.Length-1);

            }
            else{
                D_up_temp=D_up;
                D_dn_temp=D_dn;
            }

            if(D_up_temp>=basis.Dup_Range.size())
            {
                row_found_=false;
            }
            else
            {
                assert(D_up_temp<basis.Dup_Range.size());
                range_min=basis.Dup_Range[D_up_temp].first;
                range_max=basis.Dup_Range[D_up_temp].second;
                if(range_min==-1)
                {
                    row_found_=false;
                    assert(range_max==-1);
                }
                else
                {
                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                    if(i_new==-1){
                        row_found_=false;
                    }
                    else{
                        row_found_=true;
                        Inv_Trnsltns_=inv_trnsltns;
                        break;
                    }
                }

            }
        }

        //______________________________________________

        if(row_found_){
            m_new = i_new;
            sign_up_trans=pow(-1.0,sign_pow_up*1.0);
            sign_dn_trans=pow(-1.0,sign_pow_dn*1.0);
            State_temp[m_new]=FM_sign_up_down*FM_sign_up*FM_sign_down*sign_up_trans*
                    sign_dn_trans*
                    sqrt(basis.D_Period[m_new]/(1.0*basis.Length*basis.Length))*
                    ((1.0*basis.Length)/(1.0*basis.D_Period[m_new]));
        }


    }

}



void MODEL_2_orb_Hubb_chain_KSector::Perform_RVB_Analysis_at_2_hole_doped(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec){

    if( (basis.Ndn==basis.Ndn)
            &&
            (basis.Nup== basis.Length -1)
            )
    {

        if(true){
            Get_overlap_matrix_for_Anzatz_basis(basis,Eig_vec);}
        else{
            Read_Anzatz_basis(basis,Eig_vec);
        }
        cout<<scientific<<setprecision(4)<<endl;
        //cout<<"+++++++++OVERLAP MATRIX+++++++++++++"<<endl;
        //Print_Matrix(overlap_matrix_for_Anzatz_basis);
        //cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
        //    _MODEL.Create_OS_TS_states_by_reading(_BASIS,_LANCZOS.Eig_vec);

        //assert(false);
        cout<<scientific<<setprecision(10)<<endl;
        Variational_state_optimization(basis,Eig_vec);


        Mat_1_doub vec1_ ;
        double_type norm1_temp, overlap_temp;
        double norm1_, overlap_temp2;

        vec1_ = Act_Orbital_Exchange(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|RVB> is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;



        vec1_ = Act_Reflection_about_Central_site(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*(1.0/sqrt(norm1_));
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|RVB> is an eigenstate of Parity Operator [Ref. about plane cutting through center] with eigenvalue = "<< overlap_temp<<endl;


        //     cout<<"XXXXXXXXXXXXX S.S for Variational stateXXXXXXXXXXXXXXXX"<<endl;
        //     _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length,  _MODEL.State_ , _MODEL.PBC);
        //     cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

        //     cout<<"Orbital symmetry of Ansatz for half-filling:"<<endl;
        //     _MODEL.Check_orbital_symmetry(_BASIS, _MODEL.State_);

        STATE_RVB = State_;

    }


}


void MODEL_2_orb_Hubb_chain_KSector::Perform_RVB_Analysis_at_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec ){


    if( (basis.Ndn==basis.Ndn)
            &&
            (basis.Nup== basis.Length
             ))
    {

        if(true){
            Get_overlap_matrix_for_Anzatz_basis(basis,Eig_vec);}
        else{
            Read_Anzatz_basis(basis,Eig_vec);
        }
        cout<<scientific<<setprecision(4)<<endl;
        //cout<<"+++++++++OVERLAP MATRIX+++++++++++++"<<endl;
        //Print_Matrix(overlap_matrix_for_Anzatz_basis);
        //cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
        //    _MODEL.Create_OS_TS_states_by_reading(_BASIS,_LANCZOS.Eig_vec);

        //assert(false);
        cout<<scientific<<setprecision(10)<<endl;
        Variational_state_optimization(basis,Eig_vec);


        Mat_1_doub vec1_ ;
        double_type norm1_temp, overlap_temp;
        double norm1_, overlap_temp2;

        vec1_ = Act_Orbital_Exchange(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|RVB> is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;



        vec1_ = Act_Reflection_about_Central_site(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*(1.0/sqrt(norm1_));
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|RVB> is an eigenstate of Parity Operator [Ref. about plane cutting through center] with eigenvalue = "<< overlap_temp<<endl;


        //     cout<<"XXXXXXXXXXXXX S.S for Variational stateXXXXXXXXXXXXXXXX"<<endl;
        //     _LANCZOS.Measure_two_point_observables(_MODEL.two_point_obs, _MODEL.Two_point_oprts, _BASIS.Length,  _MODEL.State_ , _MODEL.PBC);
        //     cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

        //     cout<<"Orbital symmetry of Ansatz for half-filling:"<<endl;
        //     _MODEL.Check_orbital_symmetry(_BASIS, _MODEL.State_);

        STATE_RVB = State_;

    }
}


void MODEL_2_orb_Hubb_chain_KSector::Perform_Product_state_Analysis_at_2_hole_doped(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec ){

    double_type overlap_temp, norm1_temp;
    double overlap_;

    if( (basis.Ndn==basis.Ndn)
            &&
            (basis.Nup== basis.Length -1
             ))
    {

        Variational_state_pair_Geometry.resize(basis.Length);
        Variational_state_pair_spin_symmetry.resize(basis.Length);
        Variational_state_pair_orbital_symmetry.resize(basis.Length);
        Variational_state_pair_sites.resize(basis.Length);

        for(int i=0;i<basis.Length -1;i++){
            Variational_state_pair_Geometry[i]="AR";
            Variational_state_pair_spin_symmetry[i]="T";
            Variational_state_pair_orbital_symmetry[i]="NOS";
            Variational_state_pair_sites[i].first=i;
            Variational_state_pair_sites[i].second=i;
        }

        Get_Variational_State(basis,basis.Length - 1);

        overlap_temp=dot_product(Eig_vec,State_);
        overlap_ = abs(overlap_temp);
        cout<<"<GS|Triplet_Product_state> = "<<overlap_temp<<endl;


        Mat_1_doub vec1_ ;
        double norm1_, overlap_temp2;

        vec1_ = Act_Orbital_Exchange(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|Triplet_Product_state> is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;


        vec1_ = Act_Reflection_about_Central_site(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|Triplet_Product_state> is an eigenstate of Parity Operator [Ref. about plane cutting through center] with eigenvalue = "<< overlap_temp<<endl;


        STATE_TPS = State_;

    }
}

void MODEL_2_orb_Hubb_chain_KSector::Perform_Product_state_Analysis_at_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec ){

    double_type overlap_temp, norm1_temp;
    double overlap_;

    if( (basis.Ndn==basis.Ndn)
            &&
            (basis.Nup== basis.Length
             ))
    {

        Variational_state_pair_Geometry.resize(basis.Length);
        Variational_state_pair_spin_symmetry.resize(basis.Length);
        Variational_state_pair_orbital_symmetry.resize(basis.Length);
        Variational_state_pair_sites.resize(basis.Length);

        for(int i=0;i<basis.Length;i++){
            Variational_state_pair_Geometry[i]="AR";
            Variational_state_pair_spin_symmetry[i]="T";
            Variational_state_pair_orbital_symmetry[i]="NOS";
            Variational_state_pair_sites[i].first=i;
            Variational_state_pair_sites[i].second=i;
        }

        Get_Variational_State(basis,basis.Length);

        overlap_temp=dot_product(Eig_vec,State_);
        overlap_ = abs(overlap_temp);
        cout<<"<GS|Triplet_Product_state> = "<<overlap_temp<<endl;


        Mat_1_doub vec1_ ;
        double norm1_, overlap_temp2;

        vec1_ = Act_Orbital_Exchange(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|Triplet_Product_state> is an eigenstate of Orbital Symmetry Operator with eigenvalue = "<< overlap_temp<<endl;


        vec1_ = Act_Reflection_about_Central_site(basis, State_);
        norm1_temp = dot_product(vec1_,vec1_);
        norm1_ = abs(norm1_temp);
        for(int i=0;i<vec1_.size();i++){
            vec1_[i] = vec1_[i]*sqrt(1.0/norm1_);
        }
        overlap_temp = dot_product(State_,vec1_);
        overlap_temp2 = abs(overlap_temp)*abs(overlap_temp);
        cout <<"|Triplet_Product_state> is an eigenstate of Parity Operator [Ref. about plane cutting through center] with eigenvalue = "<< overlap_temp<<endl;


        STATE_TPS = State_;

    }
}

void MODEL_2_orb_Hubb_chain_KSector::Read_parameters(BASIS_2_orb_Hubb_chain_KSector &basis, string filename){


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
    string momentum_n_, Momentum_n_ = "Momentum_n = ";

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

            if ((offset = line.find(Momentum_n_, 0)) != string::npos) {
                momentum_n_ = line.substr (offset + Momentum_n_.length());		}

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
    basis.Momentum_n=atoi(momentum_n_.c_str());


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


//#endif

/*
This class includes the Model for which Lanczos is being done
*/
#ifndef USE_COMPLEX
#include "Model_3_orb_Hubbard_chain.h"
#include <stdlib.h>
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

void MODEL_3_orb_Hubb_chain::Add_diagonal_terms(BASIS_3_orb_Hubb_chain &basis){

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
            for(int gamma=0;gamma<3;gamma++){
                for(int gamma_p=gamma+1;gamma_p<3;gamma_p++){
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
            for(int gamma=0;gamma<3;gamma++){
                for(int gamma_p=gamma+1;gamma_p<3;gamma_p++){
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


            //Crystal Field Splitting (CFE):
            for(int gamma=0;gamma<3;gamma++){
                for(int site=0;site<basis.Length;site++){
                    value+=(CFS[gamma])*
                            ( ( bit_value(basis.D_up_basis[i],gamma*basis.Length + site) +
                                bit_value(basis.D_dn_basis[j],gamma*basis.Length + site) )
                              );
                }
            }



            //magnetic Field
            for(int gamma=0;gamma<3;gamma++){
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
        }
    }

}

void MODEL_3_orb_Hubb_chain::Add_non_diagonal_terms(BASIS_3_orb_Hubb_chain &basis){


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


                    } // site
                } //gamma_p
            } //gamma

        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals

}

void MODEL_3_orb_Hubb_chain::Add_connections(BASIS_3_orb_Hubb_chain &basis){

    double Hopp_connection;
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
                                        PBC
                                        )

                                    )
                            { // nearest neighbour
                                if(neigh==1){
                                    Hopp_connection = Hopping_mat_NN[gamma_p][gamma];}
                                else if(neigh==-1){
                                    Hopp_connection = Hopping_mat_NN[gamma][gamma_p];
                                }
                                else{
                                    cout <<" Only nearest neighbout allowed "<<endl;
                                    assert(abs(neigh)==1);
                                }

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

                                    if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopp_connection)*one);
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



                                    if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Hamil.value.push_back(-1.0*sign_FM*(Hopp_connection)*one);
                                        Hamil.rows.push_back((m_new));
                                        Hamil.columns.push_back((m));}

                                } // if up hopping possible


                            }//nearest neighbour
                        } //gamma_p

                    }//site_p

                } //gamma


            } // site
        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals

}


void MODEL_3_orb_Hubb_chain::Read_parameters(BASIS_3_orb_Hubb_chain &basis, string filename){


    string filepath = filename;


    double temp_val;
    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string ucoul, Ucoul = "U = ";
    string jhund, Jhund = "JHund = ";
    string upcoul, Upcoul = "Uprime = ";
    string hmag, Hmag = "H_mag = ";
    string cfs_, CFS_ = "CFS = ";
    string hopp0_, Hopp0_ = "Hopping_mat[0][orb] = ";
    string hopp1_, Hopp1_ = "Hopping_mat[1][orb] = ";
    string hopp2_, Hopp2_ = "Hopping_mat[2][orb] = ";



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

            if ((offset = line.find(Hopp2_, 0)) != string::npos) {
                hopp2_ = line.substr (offset+Hopp2_.length());				}




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

    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }


    stringstream cfs_stream;
    cfs_stream<<cfs_;

    CFS.clear();
    CFS.resize(3);


    for(int n=0;n<3;n++){
        cfs_stream >> temp_val;
        CFS[n]=temp_val;

    }



    Hopping_mat_NN.clear();
    Hopping_mat_NN.resize(3);
    for (int i=0;i<3;i++){
        Hopping_mat_NN[i].resize(3);
    }

    //Hopping_mat_NN[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}


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


void MODEL_3_orb_Hubb_chain::Read_parameters_for_dynamics(string filename){

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



void MODEL_3_orb_Hubb_chain::Initialize_macro_oprs_to_calculate(BASIS_3_orb_Hubb_chain &basis){
    macro_obs.resize(6);

    macro_obs[0]="U_intra";
    macro_obs[1]="U_inter";
    macro_obs[2]="Hunds_z";
    macro_obs[3]="Hunds_pm";
    macro_obs[4]="U_Pair";
    macro_obs[5]="H_KE";

    Macro_oprts.resize(6);



    int T_no_oprs=6;





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


            if(value2!=0){
                Macro_oprts[1].value.push_back(value2);
                Macro_oprts[1].rows.push_back(m);
                Macro_oprts[1].columns.push_back(m);
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

            if(value3!=0){
                Macro_oprts[2].value.push_back(value3);
                Macro_oprts[2].rows.push_back(m);
                Macro_oprts[2].columns.push_back(m);
            }

            //CFS
            value4=0;
            for(int gamma=0;gamma<3;gamma++){
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



                                    if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Macro_oprts[5].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma]));
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



                                    if((Hopping_mat_NN[gamma_p][gamma])!=0){
                                        assert(m_new<m);
                                        Macro_oprts[5].value.push_back(-1.0*sign_FM*(Hopping_mat_NN[gamma_p][gamma]));
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




void MODEL_3_orb_Hubb_chain::Initialize_one_point_to_calculate(BASIS_3_orb_Hubb_chain &basis){
    one_point_obs.resize(7);
    one_point_obs[0]="n_0_up";
    one_point_obs[1]="n_1_up";
    one_point_obs[2]="n_2_up";
    one_point_obs[3]="n_0_dn";
    one_point_obs[4]="n_1_dn";
    one_point_obs[5]="n_2_dn";
    one_point_obs[6]="n_i"; //write now n_k is n_i


    int T_no_oprs=one_point_obs.size();

    One_point_oprts.resize(T_no_oprs);



    int orb;
    int spin;


    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }




    for(int opr_no=0;opr_no<6;opr_no++){


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
        One_point_oprts[6][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
        One_point_oprts[6][site].ncols = One_point_oprts[6][site].nrows;
    }



    Hamiltonian_1_COO n_i;
    n_i.resize(basis.Length);

    Matrix_COO temp;

    for(int site=0;site<basis.Length;site++){
        temp = One_point_oprts[0][site];
        for(int dof=1;dof<6;dof++){
            Sum(temp, One_point_oprts[dof][site], temp, 1.0, 1.0);
        }

        n_i[site]=temp;
        One_point_oprts[6][site]=temp;
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



void MODEL_3_orb_Hubb_chain::Initialize_two_point_to_calculate(BASIS_3_orb_Hubb_chain &basis){
    two_point_obs.resize(3);
    two_point_obs[0]="SzSz";
    two_point_obs[1]="SpSm";
    two_point_obs[2]="SmSp";
    Two_point_oprts.resize(3);


    int T_no_oprs=3;



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



    }


}


void MODEL_3_orb_Hubb_chain::Initialize_Opr_for_Dynamics(BASIS_3_orb_Hubb_chain &basis){

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
                        for (int orb=0;orb<3;orb++){
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
        if(Dyn_opr_string == "J_2"){
            gamma=2;
        }

        if(Dyn_opr_string == "J_0" || Dyn_opr_string == "J_1" || Dyn_opr_string == "J_2"){
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



                            if((Hopping_mat_NN[gamma][gamma])!=0){

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



                            if((Hopping_mat_NN[gamma][gamma])!=0){

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



void MODEL_3_orb_Hubb_chain::Calculate_Local_Obs_for_States_to_Look(LANCZOS & lanczos, BASIS_3_orb_Hubb_chain & basis){

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
#endif

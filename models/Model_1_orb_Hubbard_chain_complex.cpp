/*
This class includes the Model for which Lanczos is being done
*/

#ifdef NOTUSE_COMPLEX
#include "Model_1_orb_Hubbard_chain.h"
#include <stdlib.h>
using namespace std;
#define PI 3.14159265

/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ]
    site----->[012....................(L-1)]


2)  similarly for "down spin" basis

3)  For total
    m=basis.D_dn_basis.size()*i + j;
*/



void MODEL_1_orb_Hubb_chain::Act_Hamil(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert (Vec_in.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
    Vec_out.clear();
    Vec_out.resize(basis.D_up_basis.size()*basis.D_dn_basis.size());
    Act_diagonal_terms(basis, Vec_in, Vec_out);
    Act_non_diagonal_terms(basis, Vec_in, Vec_out);
    Act_connections(basis, Vec_in, Vec_out);

}





void MODEL_1_orb_Hubb_chain::Act_connections(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

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
                for(int site_p=0;site_p<basis.Length ;site_p++){

                    if((Hopping_mat_NN[site_p][site])!=zero)

                    {

                        //---------------Hopping for up electrons-------------------//
                        //there have to be one up electron in site
                        //there have to be no up electron in site_p
                        if(
                                (bit_value(basis.D_up_basis[i], site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i], site_p)==0)
                                )
                        {

                            D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                          - pow(2,site) );


                            //i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            i_new = basis.inverse_Dup[D_up - basis.DupMin_];
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);

                            if(m_new>=m){
                                cout<<" Hopping: "<<site_p<<"  "<< site<<"  "<<Hopping_mat_NN[site_p][site]<<endl;
                            }
                            assert(m_new<m);

                            Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            Vec_out[m] += Vec_in[m_new]*conj(-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one);


                        } // if up hopping possible


                        //---------------Hopping for dn electrons-------------------//
                        //there have to be one dn electron in site
                        //there have to be no dn electron in site_p
                        if(
                                (bit_value(basis.D_dn_basis[j],site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[j],site_p)==0)
                                )
                        {

                            D_dn = (int) (basis.D_dn_basis[j] + pow(2,site_p)
                                          - pow(2,site) );


                            //j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            j_new = basis.inverse_Ddn[D_dn - basis.DdnMin_];
                            i_new = i;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);

                            assert(m_new<m);

                            Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            Vec_out[m] += Vec_in[m_new]*conj(-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one);


                        } // if up hopping possible


                    }//nearest neighbour


                }//site_p




            } // site
        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals

}


void MODEL_1_orb_Hubb_chain::Act_non_diagonal_terms(BASIS_1_orb_Hubb_chain &basis,  Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    for(int type_ind=0;type_ind<three_point_intrs.size();type_ind++){

        int m;
        double_type value;
        if(three_point_intrs[type_ind]=="SzSpSm"){
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            int site1, site2, site3;

            for(int sites_set=0;sites_set<three_point_intrs_sites_set[type_ind].size();sites_set++){
                site1=three_point_intrs_sites_set[type_ind][sites_set][0];
                site2=three_point_intrs_sites_set[type_ind][sites_set][1];
                site3=three_point_intrs_sites_set[type_ind][sites_set][2];

                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //Sz[site1]Sp[site2]*Sm[site3]:
                        //there have to be ony up electron at site2
                        //there have to be only down electron at site

                        assert(site1!=site2);
                        assert(site1!=site3);
                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==1)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                                ((bit_value(basis.D_up_basis[i], site3)==1)
                                 &&
                                 (bit_value(basis.D_dn_basis[j], site3)==0)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = (int) (basis.D_dn_basis[j] + pow(2, site3)
                                          - pow(2, site2) );

                            //i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            //j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);
                            i_new = basis.inverse_Dup[D_up - basis.DupMin_];
                            j_new = basis.inverse_Ddn[D_dn - basis.DdnMin_];

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                            value = sign_FM*(0.5*( bit_value(basis.D_up_basis[i_new], site1) -
                                                   bit_value(basis.D_dn_basis[j_new], site1) ))*three_point_intrs_vals[type_ind][sites_set];

                            //assert(m_new<m);
                            if(abs(value)>0.0000000001){
                                Vec_out[m_new] += Vec_in[m]*value*one;
                            }

                        }
                    }
                }

            }

        }


    }



}

void MODEL_1_orb_Hubb_chain::Act_diagonal_terms(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    assert (Vec_out.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
    assert (Vec_in.size() == Vec_out.size());


    //Remember H[l][m]=<l|H|m>
    int m;
    double value;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;
            //coulomb repulsion:
            value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);



            //magnetic Field
            for(int site=0;site<basis.Length;site++){
                value+=0.5*(H_field[site])*
                        ( ( bit_value(basis.D_up_basis[i],site) -
                            bit_value(basis.D_dn_basis[j],site) )
                          );
            }

            //Onsite_Energy
            for(int site=0;site<basis.Length;site++){
                value+=1.0*(Onsite_Energy[site])*
                        ( ( bit_value(basis.D_up_basis[i],site) +
                            bit_value(basis.D_dn_basis[j],site) )
                          );
                //  cout<<"site = "<<site<<" : "<<Onsite_Energy[site]<<endl;
            }

            //LongRange interactions ninj
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){

                    if(NonLocalInteractions_mat[site_i][site_j]!=0.0){

                        value+=1.0*NonLocalInteractions_mat[site_i][site_j]*(
                                    ( bit_value(basis.D_up_basis[i],site_i) + bit_value(basis.D_dn_basis[j],site_i))*
                                    ( bit_value(basis.D_up_basis[i],site_j) + bit_value(basis.D_dn_basis[j],site_j))
                                    );

                    }
                }
            }




            if(value!=0){
                Vec_out[m] +=value*one*Vec_in[m];
            }
        }
    }

}


void MODEL_1_orb_Hubb_chain::Add_diagonal_terms(BASIS_1_orb_Hubb_chain &basis){

    Hamil.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    Hamil.ncols = Hamil.nrows;



    //Remember H[l][m]=<l|H|m>
    int m;
    double_type value;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=zero;
            //on-site coulomb repulsion:
            value+=one*U*(1.0*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]));





            //magnetic Field
            for(int site=0;site<basis.Length;site++){
                value+=one*0.5*(H_field[site])*
                        ( 1.0*( bit_value(basis.D_up_basis[i],site) -
                                bit_value(basis.D_dn_basis[j],site) )
                          );
            }

            //long-range coulomb repulsion b/w n_sitei and n_sitej
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){
                    if(DenDen_Interaction_mat[site_i][site_j]!=zero){

                        value+=  DenDen_Interaction_mat[site_i][site_j]*
                                (0.5*(  ( bit_value(basis.D_up_basis[i], site_i) +
                                          bit_value(basis.D_dn_basis[j], site_i) )*
                                        ( bit_value(basis.D_up_basis[i], site_j) +
                                          bit_value(basis.D_dn_basis[j], site_j) )
                                        ));




                    }

                }

            }




            if(value!=zero){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m);
                Hamil.columns.push_back(m);
            }
        }
    }

}
void MODEL_1_orb_Hubb_chain::Add_non_diagonal_terms(BASIS_1_orb_Hubb_chain &basis){}
void MODEL_1_orb_Hubb_chain::Add_connections(BASIS_1_orb_Hubb_chain &basis){

    double_type value;
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

            value=zero;


            for(int site=0;site<basis.Length ;site++){

                for(int site_p=0;site_p<basis.Length ;site_p++){


                    if(Hopping_mat[site_p][site]!=zero)

                    {

                        //---------------Hopping for up electrons-------------------//
                        //there have to be one up electron in site
                        //there have to be no up electron in site_p
                        if(
                                (bit_value(basis.D_up_basis[i], site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i], site_p)==0)
                                )
                        {

                            D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                          - pow(2,site) );


                            //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            //i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            i_new = basis.inverse_Dup[D_up - basis.D_up_basis[0]];
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);


                            assert(m_new<m);
                            Hamil.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one);
                            Hamil.rows.push_back((m_new));
                            Hamil.columns.push_back((m));


                        } // if up hopping possible


                        //---------------Hopping for dn electrons-------------------//
                        //there have to be one dn electron in site
                        //there have to be no dn electron in site_p
                        if(
                                (bit_value(basis.D_dn_basis[j],site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[j],site_p)==0)
                                )
                        {

                            D_dn = (int) (basis.D_dn_basis[j] + pow(2,site_p)
                                          - pow(2,site) );


                            //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                            //j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);
                            j_new = basis.inverse_Ddn[D_dn - basis.D_dn_basis[0]];
                            i_new = i;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);

                            assert(m_new<m);
                            Hamil.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one);
                            Hamil.rows.push_back((m_new));
                            Hamil.columns.push_back((m));

                        } // if dn hopping possible


                    }//if hopping finite


                }//site_p




            } // site
        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals

}


void MODEL_1_orb_Hubb_chain::Read_parameters(BASIS_1_orb_Hubb_chain &basis, string filename){



    string filepath = filename;
    string pbc_,PBC_ ="PBC = ";
    string length, Length = "Length = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string ucoul, Ucoul = "U = ";

    string hmag, Hmag = "H_mag = ";

    string hopp_, Hopp_ = "Hopping NN = ";

    string longrangehopping_, LongRangeHopping_ = "LongRangeHopping = ";
    string LongRangeHoppingfile_ = "LongRangeHopping file = ";

    string longrangeinteraction_, LongRangeInteraction_ = "LongRangeInteraction = ";
    string LongRangeInteractionfile_ = "LongRangeInteraction file = ";

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


            if ((offset = line.find(LongRangeInteraction_, 0)) != string::npos) {
                longrangeinteraction_ = line.substr (offset+LongRangeInteraction_.length());  }

            if ((offset = line.find(No_Of_Onepoint_Obs_, 0)) != string::npos) {
                no_of_onepoint_obs_ = line.substr (offset+No_Of_Onepoint_Obs_.length());  }

            if ((offset = line.find(FourPointObservablesSitesFile_, 0)) != string::npos) {
                fourpointobservablessitesfile_ = line.substr (offset+FourPointObservablesSitesFile_.length());  }


            if ((offset = line.find(LongRangeInteractionfile_, 0)) != string::npos) {
                LongRangeInteractionfilepath = line.substr (offset+LongRangeInteractionfile_.length());  }


            if ((offset = line.find(LongRangeHopping_, 0)) != string::npos) {
                longrangehopping_ = line.substr (offset+LongRangeHopping_.length());  }

            if ((offset = line.find(LongRangeHoppingfile_, 0)) != string::npos) {
                LongRangeHoppingfilepath = line.substr (offset+LongRangeHoppingfile_.length());  }

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

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(Hopp_, 0)) != string::npos) {
                hopp_ = line.substr (offset+Hopp_.length());				}



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

    No_of_onepoint_obs=atoi(no_of_onepoint_obs_.c_str());

    U=atof(ucoul.c_str());


    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }


    Hopping_mat_NN.clear();
    Hopping_mat_NN.resize(1);

    Hopping_mat_NN[0].resize(1);


    //Hopping_mat_NN[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}


    stringstream hopp_stream(hopp_);


    hopp_stream >> Hopping_mat_NN[0][0];


    if(longrangeinteraction_=="true"){
        LongRangeInteraction=true;
        Read_matrix_from_file(LongRangeInteractionfilepath,
                              DenDen_Interaction_mat,basis.Length,basis.Length);
    }
    else{
        DenDen_Interaction_mat.resize(basis.Length);
        for(int site_i=0;site_i<basis.Length;site_i++){
            DenDen_Interaction_mat[site_i].resize(basis.Length);
            for(int site_j=0;site_j<basis.Length;site_j++){
                DenDen_Interaction_mat[site_i][site_j]=zero;
            }
        }
    }

    if(longrangehopping_=="true"){
        LongRangeHopping=true;
        Read_matrix_from_file(LongRangeHoppingfilepath,
                              Hopping_mat,basis.Length,basis.Length);
        for(int site=0;site<basis.Length ;site++){
            for(int site_p=0;site_p<basis.Length ;site_p++){
                if(site_p>=site){
                    Hopping_mat[site_p][site]=zero;
                }

            }}
    }
    else{
        Hopping_mat.resize(basis.Length);
        for(int site=0;site<basis.Length ;site++){
            Hopping_mat[site].resize(basis.Length);
        }


        for(int site=0;site<basis.Length ;site++){
            for(int site_p=0;site_p<basis.Length ;site_p++){
                int neigh =site-site_p;
                if(PBC==true){
                    if(site_p==0 && site==(basis.Length -1)){
                        neigh=1;}
                }

                if(neigh==1){
                    Hopping_mat[site_p][site]=one*Hopping_mat_NN[0][0];
                }
                else{
                    Hopping_mat[site_p][site]=zero;
                }

            }
        }



    }




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



}

void MODEL_1_orb_Hubb_chain::Read_parameters_for_dynamics(string filename){

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

void MODEL_1_orb_Hubb_chain::Initialize_one_point_operator_site_specific(string opr_type , Matrix_COO &OPR, int site, BASIS_1_orb_Hubb_chain &basis){

    int orb=0;
    int spin;
    if(opr_type=="n_up"){
        spin=0;
    }
    if(opr_type=="n_dn"){
        spin=1;
    }
    assert(opr_type=="n_up" || opr_type =="n_dn");


    OPR.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    OPR.ncols = OPR.nrows;

    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double_type value;


    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            //n_orb_spin[site]:
            if(spin==0){
                value=one*(1.0*bit_value(basis.D_up_basis[i],orb*basis.Length + site));
            }
            else{
                value=one*(1.0*bit_value(basis.D_dn_basis[j],orb*basis.Length + site));
            }

            if(value!=zero){
                OPR.value.push_back(value);
                OPR.rows.push_back(m);
                OPR.columns.push_back(m);
            }
        }
    }
}


void MODEL_1_orb_Hubb_chain::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR,
                                                                          int site, int site2, BASIS_1_orb_Hubb_chain &basis){

    OPR.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    OPR.ncols = OPR.nrows;


    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double value;


    if(opr_type=="SzSz"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;
                value=0.0;
                value+=(0.25*( ( bit_value(basis.D_up_basis[i], site) -
                                bit_value(basis.D_dn_basis[j], site) )*
                              ( bit_value(basis.D_up_basis[i], site2) -
                                bit_value(basis.D_dn_basis[j], site2) )
                              ));
                if(value!=0.0){
                    OPR.value.push_back(value*one);
                    OPR.rows.push_back(m);
                    OPR.columns.push_back(m);
                }
            }
        }
    }


    if(opr_type=="SpSm"){
        //Remember OPR[l][m]=<l|OPR|m>
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;

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

                    i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                    j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site;
                    lp= site2;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                    //assert(m_new<m);
                    OPR.value.push_back(one*sign_FM);
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
                        OPR.value.push_back(1.0*one);
                        OPR.rows.push_back(m);
                        OPR.columns.push_back(m);
                    }
                }
            }
        }
    }

    if(opr_type=="SmSp"){
        //Remember OPR[l][m]=<l|OPR|m>
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;

                /*if site !=site2
                Sm[site]*Sp[site2]=Sp[site2]*Sm[site]
                */
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

                    i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                    j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site2;
                    lp= site;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                    //assert(m_new<m);
                    OPR.value.push_back(sign_FM*one);
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
                        OPR.value.push_back(1.0*one);
                        OPR.rows.push_back(m);
                        OPR.columns.push_back(m);
                    }
                }
            }
        }
    }




}


void MODEL_1_orb_Hubb_chain::Initialize_one_point_to_calculate(BASIS_1_orb_Hubb_chain &basis){



    One_point_oprts_onsite.resize(No_of_onepoint_obs);

    for(int i=0;i<No_of_onepoint_obs;i++){
        Read_matrix_from_file(One_point_strs[i],
                              One_point_oprts_onsite[i],4,4);
    }



    int T_no_oprs=No_of_onepoint_obs;
    int orb;
    int spin;

    One_point_oprts.resize(T_no_oprs);
    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }


    one_point_obs=One_point_strs;

    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        for(int site=0;site<basis.Length;site++){
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }


        //Remember OPR[l][m]=<l|OPR|m>
        int m;
        double_type value;


        for(int site=0;site<basis.Length;site++){

            One_point_oprts[opr_no][site].value.clear();
            One_point_oprts[opr_no][site].rows.clear();
            One_point_oprts[opr_no][site].columns.clear();


            for (int i=0;i<basis.D_up_basis.size();i++){
                for (int j=0;j<basis.D_dn_basis.size();j++){

                    value=zero;

                    m=basis.D_dn_basis.size()*i + j;

                    if( (bit_value(basis.D_up_basis[i],site)==0)
                            &&
                            (bit_value(basis.D_dn_basis[j],site)==0)

                            )
                    {
                        value += One_point_oprts_onsite[opr_no][0][0];
                    }

                    if( (bit_value(basis.D_up_basis[i],site)==0)
                            &&
                            (bit_value(basis.D_dn_basis[j],site)==1)

                            ){
                        value +=One_point_oprts_onsite[opr_no][1][1];
                    }

                    if( (bit_value(basis.D_up_basis[i],site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[j],site)==0)

                            ){
                        value +=One_point_oprts_onsite[opr_no][2][2];
                    }

                    if( (bit_value(basis.D_up_basis[i],site)==1)
                            &&
                            (bit_value(basis.D_dn_basis[j],site)==1)

                            )
                    {
                        value +=One_point_oprts_onsite[opr_no][3][3];
                    }

                    if(value!=zero){
                        One_point_oprts[opr_no][site].value.push_back(value);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }
            }


        }












    }


}

void MODEL_1_orb_Hubb_chain::Initialize_two_point_to_calculate(BASIS_1_orb_Hubb_chain &basis){
    two_point_obs.resize(3);
    two_point_obs[0]="SzSz";
    two_point_obs[1]="SpSm";
    two_point_obs[2]="cdagger_upc_up";
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



                            //Sp_site[site]*Sm_site[site2]  Hunds coupling:
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

                                //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                i_new = basis.inverse_Dup[D_up - basis.D_up_basis[0]];
                                //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                                j_new = basis.inverse_Ddn[D_dn - basis.D_dn_basis[0]];

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l= site;
                                lp= site2;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                //assert(m_new<m);

                                Two_point_oprts[opr_no][site][site2].value.push_back(one*sign_FM);
                                Two_point_oprts[opr_no][site][site2].rows.push_back(m_new);
                                Two_point_oprts[opr_no][site][site2].columns.push_back(m);


                            }

                            if(site==site2){


                                if(
                                        ((bit_value(basis.D_up_basis[i], site2)==1)
                                         &&
                                         (bit_value(basis.D_dn_basis[j], site2)==0)
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

        if(two_point_obs[opr_no]=="cdagger_upc_up"){


            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){
                    Two_point_oprts[opr_no][site][site2].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
                    Two_point_oprts[opr_no][site][site2].ncols = Two_point_oprts[opr_no][site][site2].nrows;
                }
            }


            //Remember OPR[l][m]=<l|OPR|m>
            int m;


            int D_up,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            for(int site=0;site<basis.Length;site++){
                for(int site2=site;site2<basis.Length;site2++){


                    for (int i=0;i<basis.D_up_basis.size();i++){
                        for (int j=0;j<basis.D_dn_basis.size();j++){
                            m=basis.D_dn_basis.size()*i + j;



                            //cdagger_up[site]*c_up[site2]:
                            //there have to be up electron at site2
                            //there have to be 0 up electron at site

                            if((bit_value(basis.D_up_basis[i], site)==0
                                )
                                    &&
                                    (bit_value(basis.D_up_basis[i], site2)==1
                                     ))
                            {

                                D_up = (int) (basis.D_up_basis[i] - pow(2, site2)
                                              + pow(2, site) );


                                i_new = basis.inverse_Dup[D_up - basis.D_up_basis[0]];
                                j_new = j;

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l= site2;
                                lp= site;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                                sign_FM = pow(-1.0, sign_pow_up);

                                //assert(m_new<m);

                                Two_point_oprts[opr_no][site][site2].value.push_back(one*sign_FM);
                                Two_point_oprts[opr_no][site][site2].rows.push_back(m_new);
                                Two_point_oprts[opr_no][site][site2].columns.push_back(m);


                            }

                            if(site==site2){

                                Two_point_oprts[opr_no][site][site2].value.push_back(one*(1.0*bit_value(basis.D_up_basis[i], site2)) );
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

void MODEL_1_orb_Hubb_chain::Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_chain &basis){


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

                        value = value + 0.5*(bit_value(basis.D_up_basis[i], site)
                                             - bit_value(basis.D_dn_basis[j], site) );


                        if(value!=0){
                            Oprs_local[site].value.push_back(value);
                            Oprs_local[site].rows.push_back(m);
                            Oprs_local[site].columns.push_back(m);
                        }
                    }
                }

            }
            //local operators created ----------------------------------------------

            //In Momentum space for OBC only-----------------------------------------------------

            Matrix_COO temp;
            temp=Oprs_local[0];
            double value1, value2;
            for(int site=0;site<basis.Length-1;site++){

                value2=sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                if(site==0){
                    value1=sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
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

        if(Dyn_opr_string == "J"){
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
                        //there have to be one up electron in site
                        //there have to be no up electron in site
                        if(
                                (bit_value(basis.D_up_basis[i], site)==1)
                                &&
                                (bit_value(basis.D_up_basis[i], site_p)==0)
                                )
                        {

                            D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                          - pow(2, site) );


                            i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);



                            if((Hopping_mat_NN[0][0])!=0){

                                Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[0][0]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[0][0]));
                                Opr.rows.push_back((m));
                                Opr.columns.push_back((m_new));

                            }

                        } // if up hopping possible


                        //---------------Hopping for dn electrons-------------------//
                        //there have to be one dn electron in gamma, site
                        //there have to be no dn electron in gamma, site
                        if(
                                (bit_value(basis.D_dn_basis[j], site)==1)
                                &&
                                (bit_value(basis.D_dn_basis[j], site_p)==0)
                                )
                        {

                            D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                          - pow(2, site) );


                            j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);
                            i_new = i;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);



                            if((Hopping_mat_NN[0][0])!=0){

                                Opr.value.push_back(1.0*sign_FM*(Hopping_mat_NN[0][0]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[0][0]));
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
#endif

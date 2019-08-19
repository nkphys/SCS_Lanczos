/*
This class includes the Model for which Lanczos is being done
*/

#ifdef USE_COMPLEX
#include "Model_1_orb_tJ.h"
#include <stdlib.h>
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

void MODEL_1_orb_tJ::Add_diagonal_terms(BASIS_1_orb_tJ &basis, string run_type){


    Hamil.value.clear();
    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;




    //Remember H[l][m]=<l|H|m>
    int m,j;

    double_type value;
    for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            value=zero;

            //magnetic Field
            for(int site=0;site<basis.Length;site++){
                value+=one*0.5*(H_field[site])*
                        ( 1.0*(bit_value(basis.D_up_basis[i],site) -
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

            //SzSz exchange
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){
                    if(Jzz_Exchange_mat[site_i][site_j]!=zero){

                        value+=  Jzz_Exchange_mat[site_i][site_j]*
                                (0.5*0.25*(  ( bit_value(basis.D_up_basis[i], site_i) -
                                          bit_value(basis.D_dn_basis[j], site_i) )*
                                        ( bit_value(basis.D_up_basis[i], site_j) -
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
void MODEL_1_orb_tJ::Add_non_diagonal_terms(BASIS_1_orb_tJ &basis){


}

void MODEL_1_orb_tJ::Add_connections(BASIS_1_orb_tJ &basis, string run_type){

    string time_evolution_type= "pump_probe_tJ_Ladder";
    double_type phase;
    double_type value;
    int m,j;
    int D_up,D_dn;
    int i_new,j_new;
    int m_new;
    double sign_FM;
    int sign_pow_up, sign_pow_dn;
    int l,lp;
    int nup_temp, ndn_temp;
    int N1_temp,N2_temp;
    complex<double> iota_(0.0,1.0);

    for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;
            value=zero;


            for(int site_p=0;site_p<basis.Length ;site_p++){

                for(int site=site_p+1;site<basis.Length ;site++){


                    if(Hopping_mat[site_p][site]!=zero)

                    {

                        //-----------
                        if(time_evolution_type=="pump_probe_tJ_Ladder"
                                &&
                           run_type=="time_evolution"){
                            if(abs(site-site_p)==2){
                               phase=exp(iota_*Ax);
                            }
                            else if(abs(site-site_p)==1){
                               phase=exp(iota_*Ay);
                            }
                            else{
                                phase=one;
                            }
                        }
                        else{
                            phase=one;
                        }
                        //----------------




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
                            D_dn = basis.D_dn_basis[i];


                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                            N1_temp =0;
                            N2_temp =0;

                            N1_temp += bit_value(D_up,site) +
                                    bit_value(D_dn, site);
                            N2_temp += bit_value(D_up,site_p) +
                                    bit_value(D_dn,site_p);



                            l= site;
                            lp= site_p;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);







                            if(N1_temp <2 && N2_temp <2){
                                assert(m_new<m);
                            Hamil.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one*phase);
                            Hamil.rows.push_back((m_new));
                            Hamil.columns.push_back((m));
                            }


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
                            D_up = basis.D_up_basis[j];


                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                            //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                            N1_temp =0;
                            N2_temp =0;

                            N1_temp += bit_value(D_up,site) +
                                    bit_value(D_dn, site);
                            N2_temp += bit_value(D_up,site_p) +
                                    bit_value(D_dn,site_p);



                            l= site;
                            lp= site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);


                            if(N1_temp <2 && N2_temp <2){
                            assert(m_new<m);
                            Hamil.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one*phase);
                            Hamil.rows.push_back((m_new));
                            Hamil.columns.push_back((m));
                            }

                        } // if dn hopping possible


                    }//if hopping finite







                    if((Jpm_Exchange_mat[site_p][site]!=zero)){
                    //SpSm-connections
                    //Sp_site_p*Sm_site  Exchange coupling:
                    //there have to be ony up electron in site, site
                    //there have to be only down electron in site_p, site
                    if(
                            ((bit_value(basis.D_dn_basis[j],site_p)==1)
                             &&
                             (bit_value(basis.D_up_basis[i],site_p)==0)
                             )
                            &&
                            ((bit_value(basis.D_up_basis[i],site)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j],site)==0)
                             )
                            )
                    {

                        D_up = (int) (basis.D_up_basis[i] - pow(2,site)
                                      + pow(2,site_p) );
                        D_dn = (int) (basis.D_dn_basis[j] + pow(2,site)
                                      - pow(2,site_p) );



                            nup_temp = __builtin_popcount(D_up);
                            ndn_temp = __builtin_popcount(D_dn);
                            m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];


                        l=site_p;
                        lp=site;

                        sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                        sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                        sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                        assert(m_new<m);

                        Hamil.value.push_back(sign_FM*(0.5*(Jpm_Exchange_mat[site_p][site]))*one);
                        Hamil.rows.push_back(m_new);
                        Hamil.columns.push_back(m);

                    } // if SpSm possible

                }

                }//site_p




            } // site

    } // "i" i.e up_decimals

}


void MODEL_1_orb_tJ::Read_parameters(BASIS_1_orb_tJ &basis, string filename){



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


    string longrangeexchangezz_, LongRangeExchangeZZ_ = "LongRangeExchangeZZ = ";
    string LongRangeExchangeZZfile_ = "LongRangeExchangeZZ file = ";

    string longrangeexchangepm_, LongRangeExchangePM_ = "LongRangeExchangePM = ";
    string LongRangeExchangePMfile_ = "LongRangeExchangePM file = ";

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

            if ((offset = line.find(LongRangeExchangeZZ_, 0)) != string::npos) {
                longrangeexchangezz_ = line.substr (offset+LongRangeExchangeZZ_.length());  }

            if ((offset = line.find(LongRangeExchangePM_, 0)) != string::npos) {
                longrangeexchangepm_ = line.substr (offset+LongRangeExchangePM_.length());  }

            if ((offset = line.find(No_Of_Onepoint_Obs_, 0)) != string::npos) {
                no_of_onepoint_obs_ = line.substr (offset+No_Of_Onepoint_Obs_.length());  }

            if ((offset = line.find(FourPointObservablesSitesFile_, 0)) != string::npos) {
                fourpointobservablessitesfile_ = line.substr (offset+FourPointObservablesSitesFile_.length());  }


            if ((offset = line.find(LongRangeInteractionfile_, 0)) != string::npos) {
                LongRangeInteractionfilepath = line.substr (offset+LongRangeInteractionfile_.length());  }

            if ((offset = line.find(LongRangeExchangeZZfile_, 0)) != string::npos) {
                LongRangeExchangeZZfilepath = line.substr (offset+LongRangeExchangeZZfile_.length());  }

            if ((offset = line.find(LongRangeExchangePMfile_, 0)) != string::npos) {
                LongRangeExchangePMfilepath = line.substr (offset+LongRangeExchangePMfile_.length());  }


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
    basis.N_total=basis.Ndn+basis.Nup;

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



    if(longrangeexchangepm_=="true"){
        LongRangeExchangePM_=true;
        Read_matrix_from_file(LongRangeExchangePMfilepath,
                              Jpm_Exchange_mat ,basis.Length,basis.Length);
    }
    else{
        Jpm_Exchange_mat.resize(basis.Length);
        for(int site_i=0;site_i<basis.Length;site_i++){
            Jpm_Exchange_mat[site_i].resize(basis.Length);
            for(int site_j=0;site_j<basis.Length;site_j++){
                Jpm_Exchange_mat[site_i][site_j]=zero;
            }
        }
    }

    if(longrangeexchangezz_=="true"){
        LongRangeExchangeZZ_=true;
        Read_matrix_from_file(LongRangeExchangeZZfilepath,
                              Jzz_Exchange_mat ,basis.Length,basis.Length);
    }
    else{
        Jzz_Exchange_mat.resize(basis.Length);
        for(int site_i=0;site_i<basis.Length;site_i++){
            Jzz_Exchange_mat[site_i].resize(basis.Length);
            for(int site_j=0;site_j<basis.Length;site_j++){
                Jzz_Exchange_mat[site_i][site_j]=zero;
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

void MODEL_1_orb_tJ::Read_parameters_for_dynamics(string filename){

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


void MODEL_1_orb_tJ::Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis){

 /*   int T_no_oprs=2;
    int orb;
    int spin;



   //  0 n
   //  1 Sz

    One_point_oprts.resize(T_no_oprs);
    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }


*/




}

void MODEL_1_orb_tJ::Initialize_one_point_to_calculate_from_file(BASIS_1_orb_tJ &basis){



    One_point_oprts_onsite.resize(No_of_onepoint_obs);

    for(int i=0;i<No_of_onepoint_obs;i++){
        Read_matrix_from_file(One_point_strs[i],
                              One_point_oprts_onsite[i],3,3);
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
            One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size();
            One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
        }


        //Remember OPR[l][m]=<l|OPR|m>
        int m,j;
        double_type value;


        for(int site=0;site<basis.Length;site++){

            One_point_oprts[opr_no][site].value.clear();
            One_point_oprts[opr_no][site].rows.clear();
            One_point_oprts[opr_no][site].columns.clear();


            for (int i=0;i<basis.D_up_basis.size();i++){
               j=i;

                    value=zero;

                    m=i;

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



                    if(value!=zero){
                        One_point_oprts[opr_no][site].value.push_back(value);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }

            }


        }












    }




    H_Total.value.clear();
    H_Total.nrows = basis.D_up_basis.size();
    H_Total.ncols = H_Total.nrows;
    H_Total=Hamil;

    //H_KE Macro observable to measure
    //H_KE
    H_KE.value.clear();
    H_KE.nrows = basis.D_up_basis.size();
    H_KE.ncols = H_KE.nrows;


       string time_evolution_type= "pump_probe_tJ_Ladder";
       double_type phase;
       double_type value;
       int m,j;
       int D_up,D_dn;
       int i_new,j_new;
       int m_new;
       double sign_FM;
       int sign_pow_up, sign_pow_dn;
       int l,lp;
       int nup_temp, ndn_temp;
       int N1_temp,N2_temp;
       complex<double> iota_(0.0,1.0);

       for (int i=0;i<basis.D_up_basis.size();i++){

               m=i;
               j=i;
               value=zero;


               for(int site_p=0;site_p<basis.Length ;site_p++){

                   for(int site=site_p+1;site<basis.Length ;site++){


                       if(Hopping_mat[site_p][site]!=zero)

                       {

                           //-----------
                           if(time_evolution_type=="pump_probe_tJ_Ladder"){
                               if(abs(site-site_p)==2){
                                  phase=exp(iota_*Ax);
                               }
                               else if(abs(site-site_p)==1){
                                  phase=exp(iota_*Ay);
                               }
                               else{
                                   phase=one;
                               }
                           }
                           else{
                               phase=one;
                           }
                           //----------------




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
                               D_dn = basis.D_dn_basis[i];


                               nup_temp = __builtin_popcount(D_up);
                               ndn_temp = __builtin_popcount(D_dn);
                               // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                               //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                               m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                               N1_temp =0;
                               N2_temp =0;

                               N1_temp += bit_value(D_up,site) +
                                       bit_value(D_dn, site);
                               N2_temp += bit_value(D_up,site_p) +
                                       bit_value(D_dn,site_p);



                               l= site;
                               lp= site_p;

                               sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                               sign_FM = pow(-1.0, sign_pow_up);







                               if(N1_temp <2 && N2_temp <2){
                                   assert(m_new<m);
                               H_KE.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one*phase);
                               H_KE.rows.push_back((m_new));
                               H_KE.columns.push_back((m));
                               }


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
                               D_up = basis.D_up_basis[j];


                               nup_temp = __builtin_popcount(D_up);
                               ndn_temp = __builtin_popcount(D_dn);
                               // m_new = Find_commont_int(basis.D_up_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]],
                               //       basis.D_dn_reverse[ndn_temp][D_dn-basis.D_dn_min[ndn_temp]]);
                               m_new = basis.D_updn_reverse[nup_temp][D_up-basis.D_up_min[nup_temp]][D_dn-basis.D_dn_min[ndn_temp]];

                               N1_temp =0;
                               N2_temp =0;

                               N1_temp += bit_value(D_up,site) +
                                       bit_value(D_dn, site);
                               N2_temp += bit_value(D_up,site_p) +
                                       bit_value(D_dn,site_p);



                               l= site;
                               lp= site_p;

                               sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                               sign_FM = pow(-1.0, sign_pow_dn);


                               if(N1_temp <2 && N2_temp <2){
                               assert(m_new<m);
                               H_KE.value.push_back(1.0*sign_FM*((Hopping_mat[site_p][site]))*one*phase);
                               H_KE.rows.push_back((m_new));
                               H_KE.columns.push_back((m));
                               }

                           } // if dn hopping possible


                       }//if hopping finite


                   }//site_p




               } // site

       } // "i" i.e up_decimals





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

void MODEL_1_orb_tJ::Initialize_Opr_for_Dynamics(BASIS_1_orb_tJ &basis){



}
#endif

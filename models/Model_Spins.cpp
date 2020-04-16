/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "Model_Spins.h"
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
            value+=one*(H_field[2])*
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
    int i_new, i_new_temp;
    int m_new;
    int Val_site, Val_site_p;
    int Val_site_new, Val_site_p_new;
    double_type Factor;
    complex<double> iota_(0.0,1.0);

    for (int i=basis.D_min;i<basis.D_max+1;i++){

        value=zero;

        for(int site_p=0;site_p<basis.Length ;site_p++){
            Val_site_p = value_at_pos(i, site_p, basis.BASE);

            for(int site=site_p+1;site<basis.Length ;site++){
                Val_site = value_at_pos(i, site, basis.BASE);


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


                        i_new_temp = Updated_decimal_with_value_at_pos(i, site, basis.BASE, Val_site_new);
                        i_new = Updated_decimal_with_value_at_pos(i_new_temp, site_p, basis.BASE, Val_site_p_new);


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

    } // "i" i.e up_decimals

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

            if ((offset = line.find(LongRangeExchangeZZfile_, 0)) != string::npos) {
                LongRangeExchangeZZfilepath = line.substr (offset+LongRangeExchangeZZfile_.length());  }

            if ((offset = line.find(LongRangeExchangePMfile_, 0)) != string::npos) {
                LongRangeExchangePMfilepath = line.substr (offset+LongRangeExchangePMfile_.length());  }

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



    Read_matrix_from_file(LongRangeExchangePMfilepath,
                          Jpm_Exchange_mat ,basis.Length,basis.Length);

    Read_matrix_from_file(LongRangeExchangeZZfilepath,
                          Jzz_Exchange_mat ,basis.Length,basis.Length);



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
    int i_new;


    if(Dyn_opr_string=="Sz"){
        for (int i=basis.D_min;i<basis.D_max + 1;i++){
            value=zero;
            //Sz[site]
            value=one*
                    ( ( (1.0*value_at_pos(i, site_, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );


            if(value!=zero){
                Dyn_opr.value.push_back(value*one);
                Dyn_opr.rows.push_back(i);
                Dyn_opr.columns.push_back(i);
            }
        }
    }
    else if(Dyn_opr_string=="Sm"){

        for (int i=basis.D_min;i<basis.D_max+1;i++){

            value=zero;
            Val_site = value_at_pos(i, site_, basis.BASE);

            //Sm[site]:
            //site cannot be in -Spin

            if(Val_site != 0)
            {
                Val_site_new = Val_site - 1;
                i_new = Updated_decimal_with_value_at_pos(i, site_, basis.BASE, Val_site_new);

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
    int i_new;


    if(Dyn_opr_string=="Sz"){

        for(int site_=0;site_<basis.Length;site_++){
            for (int i=basis.D_min;i<basis.D_max + 1;i++){
                value=zero;
                //Sz[site]
                value=one*
                        ( ( (1.0*value_at_pos(i, site_, basis.BASE)) - (0.5*basis.TwoTimesSpin)) );

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
            for (int i=basis.D_min;i<basis.D_max+1;i++){

                value=zero;
                Val_site = value_at_pos(i, site_, basis.BASE);

                //Sm[site]:
                //site cannot be in -Spin

                if(Val_site != 0)
                {
                    Val_site_new = Val_site - 1;
                    i_new = Updated_decimal_with_value_at_pos(i, site_, basis.BASE, Val_site_new);

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

#include "Model_1_orb_Hubbard_GC.h"
#include <stdlib.h>
using namespace std;
#define PI 3.14159265
/*convention for basis:

    1)  for "up-spin" basis
                  [_______________________  _  ]
        site----->[012....................(L-1)]


    2)  similarly for "down-spin" basis

    */
template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_diagonal_terms(){

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
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_non_diagonal_terms(){
    //NOTHING FOR THIS MODEL
}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Add_connections(){


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
void MODEL_1_orb_Hubbard_GC<Basis_type>::Read_parameters(string filename){


    string filepath = filename;


    double temp_val;
    string length, Length = "Length = ";
    string ntotal, Ntotal = "N_Total = ";

    string ucoul, Ucoul = "U = ";
    string hmag, Hmag = "H_mag = ";
    string cfs_site_resolved, CFS_SITE_RESOLVED = "CFS_SITE_RESOLVED = ";

    string LongRangeHoppingfile_ = "LongRangeHopping_file = ";

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

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(CFS_SITE_RESOLVED, 0)) != string::npos) {
                cfs_site_resolved = line.substr (offset+CFS_SITE_RESOLVED.length());				}

            if ((offset = line.find(LongRangeHoppingfile_, 0)) != string::npos) {
                LongRangeHoppingfilepath = line.substr (offset+LongRangeHoppingfile_.length());  }

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

    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }
    cout<<"H_field = "<<h<<endl;



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

    //Print_Matrix(Hopping_mat_LongRange);

    //THINK ABOUT IT LATER :)
    /*for(int site=0;site<basis.Length ;site++){
        for(int site_p=0;site_p<basis.Length ;site_p++){
            if(site_p>=site){
                Hopping_mat[site_p][site]=zero;
            }

        }}
        */

    cout<<"READING PARAMETERS"<<endl;


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


    int TOTAL_NO_OBS=4;
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

            cout<<one_point_obs[obs_no]<<"["<<site<<"] = "<<value_.real() << "  "<<value_.imag()<<endl;
        }

        cout<<endl;

        vector< int >().swap( OPR_.columns );
        vector< int >().swap( OPR_.rows );
        vector< double_type >().swap( OPR_.value );
    }

}

template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site){


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
            for(int sigma_p=0;sigma_p<2;sigma_p++){

                //A[sigma][sigma_p]*c^{dagger}(sigma)*c(sigma_p)
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

        OPR.value.push_back(value_diagonal*one);
        OPR.rows.push_back(i);
        OPR.columns.push_back(i);
    }




    //OFF-diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){

                //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma)*c(sigma_p,gamma_p)
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

                            //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);


                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;


                            //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                            //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                            //m_new = basis.D_dn_basis.size()*i_new + j_new;


                            l=site;
                            lp=site;

                            max_up = basis.Length -1;
                            min_up = 0;
                            max_dn = basis.Length -1;
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

                            OPR.value.push_back(sign_FM*AMat[(sigma)][(sigma_p)]*one);
                            OPR.rows.push_back(m_new);
                            OPR.columns.push_back(i);

                        }
                    }

                }


            }

        }
    }


}


template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site, int site_p){

    //sum_{sigma,sigma_p} AMat[sigma][sigma_p]
    //c_{sigma}^{dagger,site}*c_{sigma_p,site_p}
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
            for(int sigma_p=0;sigma_p<2;sigma_p++){

                //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma)*c(sigma_p,gamma_p)
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

        OPR.value.push_back(value_diagonal*one);
        OPR.rows.push_back(i);
        OPR.columns.push_back(i);
    }




    //OFF-diagonal
    for (int i=0;i<basis.D_up_basis.size();i++){

        for(int sigma=0;sigma<2;sigma++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){

                //A[sigma,gamma][sigma_p,gamma_p]*c^{dagger}(sigma,gamma,site)*c(sigma_p,gamma_p,site_p)
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

                            //i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                            //j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                            i_new = Find_int_in_intarray(D_up,basis.Canonical_partition_up[__builtin_popcount(D_up)]);
                            j_new = Find_int_in_intarray(D_dn,basis.Canonical_partition_dn[__builtin_popcount(D_up)]);

                            m_new = (basis.Canonical_partition_dn[__builtin_popcount(D_up)].size()*i_new + j_new) +
                                    basis.Nup_offsets[__builtin_popcount(D_up)].first;


                            //m_new = Find_intpair_in_intarraypair(D_up,D_dn,basis.D_up_basis,
                            //basis.D_dn_basis,__builtin_popcount(D_up),basis.Nup_offsets);

                            //m_new = basis.D_dn_basis.size()*i_new + j_new;


                            l=site;
                            lp=site_p;

                            max_up = basis.Length -1;
                            min_up = 0;
                            max_dn = basis.Length -1;
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

                            OPR.value.push_back(sign_FM*AMat[(sigma)][(sigma_p)]*one);
                            OPR.rows.push_back(m_new);
                            OPR.columns.push_back(i);

                        }
                    }
                }
            }
        }
    }



}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Initialize_Opr_for_Dynamics(LANCZOS &lanczos_GS){

}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_c_on_GS(LANCZOS & lanczos, BASIS_1_orb_Hubbard_GC & basis_Nm1,
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

                    value = sign_FM*lanczos.Eig_vec[i]*value_in;

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

                    value = sign_FM*lanczos.Eig_vec[i]*value_in;

                    State_c_on_GS[i_new] += value;


                }
            }



        }


    }



}



template <typename Basis_type>
void MODEL_1_orb_Hubbard_GC<Basis_type>::Get_cdagger_on_GS(LANCZOS & lanczos, BASIS_1_orb_Hubbard_GC & basis_Np1,
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






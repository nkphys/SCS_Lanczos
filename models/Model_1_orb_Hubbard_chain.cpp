/*
This class includes the Model for which Lanczos is being done
*/

//#ifndef USE_COMPLEX
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

void MODEL_1_orb_Hubb_chain::Add_diagonal_terms(BASIS_1_orb_Hubb_chain &basis){

    Hamil.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    Hamil.ncols = Hamil.nrows;

    double EPS_=0.0000001;

    //Remember H[l][m]=<l|H|m>
    int m;
    double_type value;
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            m=basis.D_dn_basis.size()*i + j;

            value=0;
            //coulomb repulsion:
            //value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);
            for(int site=0;site<basis.Length;site++){
                value+=1.0*(Onsite_U[site])*
                         ( 1.0*( bit_value(basis.D_up_basis[i],site)*
                           bit_value(basis.D_dn_basis[j],site) )
                          );
                //  cout<<"site = "<<site<<" : "<<Onsite_Energy[site]<<endl;
            }


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

             //LongRange Direct Exchange Szi X Szj
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){

                    if(DirectExchange_mat[site_i][site_j]!=0.0){

                        value+=0.25*DirectExchange_mat[site_i][site_j]*(
                                    ( bit_value(basis.D_up_basis[i],site_i) - bit_value(basis.D_dn_basis[j],site_i))*
                                    ( bit_value(basis.D_up_basis[i],site_j) - bit_value(basis.D_dn_basis[j],site_j))
                                    );

                    }
                }
            }




            if(abs(value)>EPS_){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(m);
                Hamil.columns.push_back(m);
            }
        }
    }

}
void MODEL_1_orb_Hubb_chain::Add_non_diagonal_terms(BASIS_1_orb_Hubb_chain &basis){


    bool DirectExchange_term=true;
    bool PairHopping_term=true;
    bool InteractionAssistedHopping_term=true;

    if(DirectExchange_term){
        int m;
        double_type value;
        
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            

                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if(DirectExchange_mat[site2][site3]!=0.0){

                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //Sp[site2]*Sm[site3]:

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

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                            value = sign_FM*0.5*DirectExchange_mat[site2][site3];

                            //assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }

            }
    }
    }
        
    }


       if(PairHopping_term){
        int m;
        double_type value;
        
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            

                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if(PairHopping_mat[site2][site3]!=0.0){

                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //P*[site2]*P[site3]:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==0)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                                ((bit_value(basis.D_up_basis[i], site3)==1)
                                 &&
                                 (bit_value(basis.D_dn_basis[j], site3)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site3)
                                          + pow(2, site2) );

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);


                            value = sign_FM*PairHopping_mat[site2][site3];

                            assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }

            }
    }
    }
        
    }


       if(InteractionAssistedHopping_term){
        int m;
        double_type value;
        
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            

                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if(InteractionAssistedHopping_mat[site2][site3]!=0.0){


//---------------------s=up term1---------------------------------//
                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //n_{site2,up}c_{site2,dn}* c_{site3,dn}:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==0)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==1)
                            )
                                &&
                            ((bit_value(basis.D_dn_basis[j], site3)==1)
                                 ))
                        {

                            D_up = basis.D_up_basis[i];
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site3)
                                          + pow(2, site2) );

                            i_new = i;
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_dn);


                            value = sign_FM*InteractionAssistedHopping_mat[site2][site3];

                            assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }
//--------------s=up term 1 done----------------------//


//---------------------s=up term2---------------------------------//
                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //n_{site3,up}c_{site2,dn}* c_{site3,dn}:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==0)
                            )
                                &&
                            ((bit_value(basis.D_dn_basis[j], site3)==1) &&
                             (bit_value(basis.D_up_basis[i], site3)==1)
                                 ))
                        {

                            D_up = basis.D_up_basis[i];
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site3)
                                          + pow(2, site2) );

                            i_new = i;
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_dn);


                            value = sign_FM*InteractionAssistedHopping_mat[site2][site3];

                            assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }
//--------------s=up term 2 done----------------------//


//--------------s=dn term 1-----------------------//
                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //n_{site2,dn}c_{site2,up}* c_{site3,up}:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==1)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                            ((bit_value(basis.D_up_basis[i], site3)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = basis.D_dn_basis[j];

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);

                            value = sign_FM*InteractionAssistedHopping_mat[site2][site3];

                            assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }

//---------------s=dn term 1 done-----------------//


//--------------s=dn term 2-----------------------//
                for (int i=0;i<basis.D_up_basis.size();i++){
                    for (int j=0;j<basis.D_dn_basis.size();j++){
                        m=basis.D_dn_basis.size()*i + j;

                        //n_{site3,dn}c_{site2,up}* c_{site3,up}:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                            ((bit_value(basis.D_up_basis[i], site3)==1) &&
                             (bit_value(basis.D_dn_basis[j], site3)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = basis.D_dn_basis[j];

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);

                            value = sign_FM*InteractionAssistedHopping_mat[site2][site3];

                            //assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }

//---------------s=dn term 2 done-----------------//



            }
    }
    }
        
    }




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

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                            value = sign_FM*(0.5*( bit_value(basis.D_up_basis[i_new], site1) -
                                                   bit_value(basis.D_dn_basis[j_new], site1) ))*three_point_intrs_vals[type_ind][sites_set];

                            //assert(m_new<m);
                            if(m_new<m){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((m_new));
                                    Hamil.columns.push_back((m));
                                }
                            }
                        }
                    }
                }

            }

        }


    }



}
void MODEL_1_orb_Hubb_chain::Add_connections(BASIS_1_orb_Hubb_chain &basis){

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


                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
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
                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one);
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


                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);
                            i_new = i;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site;
                            lp= site_p;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

                            sign_FM = pow(-1.0, sign_pow_dn);

                            assert(m_new<m);
                            Hamil.value.push_back(-1.0*sign_FM*conjugate(Hopping_mat_NN[site_p][site])*one);
                            Hamil.rows.push_back((m_new));
                            Hamil.columns.push_back((m));

                        } // if up hopping possible


                    }//nearest neighbour


                }//site_p




            } // site
        }// "j" i.e dn_decimals
    } // "i" i.e up_decimals

}





void MODEL_1_orb_Hubb_chain::Act_Hamil(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    assert (Vec_in.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
    Vec_out.clear();
    Vec_out.resize(basis.D_up_basis.size()*basis.D_dn_basis.size());
    Act_diagonal_terms(basis, Vec_in, Vec_out);
    cout<<"Diagonal done"<<endl;
    Act_non_diagonal_terms(basis, Vec_in, Vec_out);
    cout<<"Non diagonal done"<<endl;
    Act_connections(basis, Vec_in, Vec_out);
    cout<<"Connections done"<<endl;

}





void MODEL_1_orb_Hubb_chain::Act_connections(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){



    //Symmetrized, so that parallelization can be done efficiently.
    Mat_2_doub Hopping_Mat;
    Hopping_Mat.resize(basis.Length);
    for(int sitei=0;sitei<basis.Length;sitei++){
        Hopping_Mat[sitei].resize(basis.Length);
        for(int sitej=0;sitej<basis.Length;sitej++){
            if(sitej>=sitei){
                Hopping_Mat[sitei][sitej] = Hopping_mat_NN[sitei][sitej];
            }
            else{
                Hopping_Mat[sitei][sitej] = conjugate(Hopping_mat_NN[sitej][sitei]);
            }
        }
    }

    //    Mat_2_doub Vec_out_temp;
    int no_of_proc;
    no_of_proc=1;

#ifdef _OPENMP
    no_of_proc=min(basis.Length, NProcessors_);
    omp_set_num_threads(no_of_proc);
    cout<<"Connections acting: "<<no_of_proc<<" processors"<<endl;
    //    Vec_out_temp.resize(no_of_proc);
    //    for(int i=0;i<no_of_proc;i++){
    //        Vec_out_temp[i].resize(Vec_in.size());
    //    }
#endif





    //    for (int i=0;i<basis.D_up_basis.size();i++){
    //        for (int j=0;j<basis.D_dn_basis.size();j++){
    //            m=basis.D_dn_basis.size()*i + j;


    for(int site=0;site<basis.Length ;site++){
        for(int site_p=0;site_p<basis.Length ;site_p++){

            if(abs(Hopping_Mat[site_p][site])>0.0000000001)
            {

#ifdef _OPENMP
#pragma omp parallel
                {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
                    for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){

                        double value;
                        int D_up,D_dn;
                        int i_new,j_new;
                        int m_new;
                        double sign_FM;
                        int sign_pow_up, sign_pow_dn;
                        int l,lp;

                        int mytid;
#ifdef _OPENMP
                        mytid = omp_get_thread_num();
#endif

                        int i,j;
                        j = m%basis.D_dn_basis.size();
                        i = int (m/basis.D_dn_basis.size());



                        value=0;


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

                            //if(m_new>=m){
                            //cout<<" Hopping: "<<site_p<<"  "<< site<<"  "<<Hopping_mat_NN[site_p][site]<<endl;
                            //}
                            //assert(m_new<m);

                            //HERE

#ifdef _OPENMP
                            // Vec_out_temp[mytid][m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            //Vec_out_temp[mytid][m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));
                            Vec_out[m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));

#else
                            // Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            Vec_out[m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));
#endif


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

                            //                            assert(m_new<m);

#ifdef _OPENMP
                            // Vec_out_temp[mytid][m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            //Vec_out_temp[mytid][m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(Hopping_Mat[site_p][site])*one);
                            Vec_out[m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]))*one);
#else
                            //Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat_NN[site_p][site])*one;
                            Vec_out[m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]))*one);
#endif


                        } // if dn hopping possible



                    }//m

#ifdef _OPENMP
                }
#endif


            }//Hopping non-zero
        }//site_p

    } // site

    //        }// "j" i.e dn_decimals
    //    } // "i" i.e up_decimals










    //#ifdef _OPENMP
    //#pragma omp parallel for default(shared)
    //    for(int comp=0;comp<Vec_in.size();comp++){
    //        for(int Np=0;Np<no_of_proc;Np++){
    //            Vec_out[comp] += Vec_out_temp[Np][comp];
    //        }
    //    }

    //    for(int Np=0;Np<no_of_proc;Np++){
    //        vector < double_type >().swap(Vec_out_temp[Np]);
    //    }
    //#endif





}


void MODEL_1_orb_Hubb_chain::Act_non_diagonal_terms(BASIS_1_orb_Hubb_chain &basis,  Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){




    bool DirectExchange_term=true;
    bool PairHopping_term=true;
    bool InteractionAssistedHopping_term=true;

    if(DirectExchange_term){
        Mat_2_doub Vec_out_temp;
        int no_of_proc;
        no_of_proc=1;
        int temp_int;
#ifdef _OPENMP
            no_of_proc=NProcessors_;
            omp_set_num_threads(no_of_proc);
            cout<<"Non diagonal acting: "<<no_of_proc<<" processors"<<endl;
#endif


#ifdef _OPENMP
#pragma omp parallel
            {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
        for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){
            int mytid;
            
#ifdef _OPENMP
                    mytid = omp_get_thread_num();
#endif

                    int i,j;
                    j = m%basis.D_dn_basis.size();
                    i = int (m/basis.D_dn_basis.size());

        for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if( (DirectExchange_mat[site2][site3] + DirectExchange_mat[site3][site2]  )!=0.0){

        
            double_type value;
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;
            

                        //Sp[site2]*Sm[site3]:

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

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                            value = sign_FM*0.5*(1.0*(DirectExchange_mat[site2][site3] +  DirectExchange_mat[site3][site2] )  );

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
                                //Above is needed instead of eqn below, because while parallely running 
                                //different "m" running prallely can lead to same m_new and overwriting + operation
                                //Vec_out[m_new] += Vec_in[m]*value*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif

                                }
                          
                        }
            }
    }
    }

                }
    #ifdef _OPENMP
            }
#endif
        
    }

//---------------------------------------------------------------------



       if(PairHopping_term){
        Mat_2_doub Vec_out_temp;
        int no_of_proc;
        no_of_proc=1;
        int temp_int;
#ifdef _OPENMP
            no_of_proc=NProcessors_;
            omp_set_num_threads(no_of_proc);
            cout<<"Non diagonal acting: "<<no_of_proc<<" processors"<<endl;
#endif


#ifdef _OPENMP
#pragma omp parallel
            {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
        for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){
            int mytid;

#ifdef _OPENMP
                    mytid = omp_get_thread_num();
#endif

                    int i,j;
                    j = m%basis.D_dn_basis.size();
                    i = int (m/basis.D_dn_basis.size());


                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if( (PairHopping_mat[site2][site3] +  PairHopping_mat[site3][site2])!=0.0){
                        double_type value;
                        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
                        double sign_FM;
                        

                        //P*[site2]*P[site3]:

                        assert(site2!=site3);

                        if(((bit_value(basis.D_dn_basis[j], site2)==0)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                                ((bit_value(basis.D_up_basis[i], site3)==1)
                                 &&
                                 (bit_value(basis.D_dn_basis[j], site3)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site3)
                                          + pow(2, site2) );

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn);


                            value = sign_FM*(PairHopping_mat[site2][site3]+PairHopping_mat[site3][site2]);

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif
                                }
                            
                        


            }
    }
    }
        }
       
    }
#ifdef _OPENMP
            }
#endif

       }

//-------------------------------------------------------------------
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//-------------------------------------------------------------------


       if(InteractionAssistedHopping_term){
        
         Mat_2_doub Vec_out_temp;
        int no_of_proc;
        no_of_proc=1;
        int temp_int;
#ifdef _OPENMP
            no_of_proc=NProcessors_;
            omp_set_num_threads(no_of_proc);
            cout<<"Non diagonal acting: "<<no_of_proc<<" processors"<<endl;
#endif


#ifdef _OPENMP
#pragma omp parallel
            {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
        for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){
            int mytid;

#ifdef _OPENMP
                    mytid = omp_get_thread_num();
#endif

                    int i,j;
                    j = m%basis.D_dn_basis.size();
                    i = int (m/basis.D_dn_basis.size());



            
            double_type value;
            int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
            double sign_FM;

                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){


            if(  (InteractionAssistedHopping_mat[site2][site3])!=0.0 ){
                        


//---------------------s=up term1---------------------------------//
                

                        //n_{site2,up}c_{site2,dn}* c_{site3,dn}:

                        assert(site2!=site3);
                        if(true){
                        if(((bit_value(basis.D_dn_basis[j], site2)==0)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==1)
                            )
                                &&
                            ((bit_value(basis.D_dn_basis[j], site3)==1)
                                 ))
                        {

                            D_up = basis.D_up_basis[i];
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site3)
                                          + pow(2, site2) );

                            i_new = i;
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_dn);


                            value = sign_FM*(InteractionAssistedHopping_mat[site2][site3] );

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif

                                }
                            //}
                        }
            }
//--------------s=up term 1 done----------------------//


//---------------------s=up term2---------------------------------//
          
                        
                        //n_{site2,up}c_{site3,dn}* c_{site2,dn}:

                        assert(site2!=site3);

                        if(true){
                        if(((bit_value(basis.D_dn_basis[j], site3)==0)
                            )
                                &&
                            ((bit_value(basis.D_dn_basis[j], site2)==1) &&
                             (bit_value(basis.D_up_basis[i], site2)==1)
                                 ))
                        {

                            D_up = basis.D_up_basis[i];
                            D_dn = (int) (basis.D_dn_basis[j] - pow(2, site2)
                                          + pow(2, site3) );

                            i_new = i;
                            j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site2;
                            lp= site3;

                            sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                            sign_FM = pow(-1.0, sign_pow_dn);


                            value = sign_FM*(InteractionAssistedHopping_mat[site2][site3]  );

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif
                                }
                           // }
                        }
            }
             
//--------------s=up term 2 done----------------------//


//--------------s=dn term 1-----------------------//
               
                        

                        //n_{site2,dn}c_{site2,up}* c_{site3,up}:

                        assert(site2!=site3);

                        if(true){
                        if(((bit_value(basis.D_dn_basis[j], site2)==1)
                            &&
                            (bit_value(basis.D_up_basis[i], site2)==0)
                            )
                                &&
                            ((bit_value(basis.D_up_basis[i], site3)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
                                          + pow(2, site2) );
                            D_dn = basis.D_dn_basis[j];

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site3;
                            lp= site2;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);

                            value = sign_FM*(InteractionAssistedHopping_mat[site2][site3] );

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif
                                }
                            //}
                        }
                        }

//---------------s=dn term 1 done-----------------//


//--------------s=dn term 2-----------------------//

                        //n_{site2,dn}c_{site3,up}* c_{site2,up}:

                        assert(site2!=site3);
if(true){
                        if(((bit_value(basis.D_up_basis[i], site3)==0)
                            )
                                &&
                            ((bit_value(basis.D_up_basis[i], site2)==1) &&
                             (bit_value(basis.D_dn_basis[j], site2)==1)
                                 ))
                        {

                            D_up = (int) (basis.D_up_basis[i] - pow(2, site2)
                                          + pow(2, site3) );
                            D_dn = basis.D_dn_basis[j];

                            i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                            j_new = j;

                            m_new = basis.D_dn_basis.size()*i_new + j_new;

                            l= site2;
                            lp= site3;

                            sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

                            sign_FM = pow(-1.0, sign_pow_up);

                            value = sign_FM*(InteractionAssistedHopping_mat[site2][site3] );

                            //assert(m_new<m);
                            //if(m_new<m){
                                if(abs(value)>0.0000000001){
#ifdef _OPENMP
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif
                                }
                            //}
                        }
                        }

//---------------s=dn term 2 done-----------------//



            }
    }
    }
        
    }

#ifdef _OPENMP
            }
#endif
         
       }




//-------------------------------------------------------------------
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//--------------------------------------------------------------------

    for(int type_ind=0;type_ind<three_point_intrs.size();type_ind++){


        if(three_point_intrs[type_ind]=="SzSpSm"){

            Mat_2_doub Vec_out_temp;
            int no_of_proc;
            no_of_proc=1;
            int temp_int= three_point_intrs_sites_set[type_ind].size();

#ifdef _OPENMP
            no_of_proc=min(temp_int, NProcessors_);
            omp_set_num_threads(no_of_proc);
            cout<<"Non diagonal acting: "<<no_of_proc<<" processors"<<endl;
            //            Vec_out_temp.resize(no_of_proc);
            //            for(int i=0;i<no_of_proc;i++){
            //                Vec_out_temp[i].resize(Vec_in.size());
            //            }
#endif



#ifdef _OPENMP
#pragma omp parallel
            {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
                for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){

                    int mytid;
#ifdef _OPENMP
                    mytid = omp_get_thread_num();
#endif

                    int i,j;
                    j = m%basis.D_dn_basis.size();
                    i = int (m/basis.D_dn_basis.size());


                    for(int sites_set=0;sites_set<three_point_intrs_sites_set[type_ind].size();sites_set++){


                        double_type value;
                        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
                        double sign_FM;
                        int site1, site2, site3;

                        site1=three_point_intrs_sites_set[type_ind][sites_set][0];
                        site2=three_point_intrs_sites_set[type_ind][sites_set][1];
                        site3=three_point_intrs_sites_set[type_ind][sites_set][2];



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
                            if(abs(value)>0.00000000){

#ifdef _OPENMP
                                //Vec_out_temp[mytid][m_new] += Vec_in[m]*value*one;
                                Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
#else
                                Vec_out[m_new] += Vec_in[m]*value*one;
#endif


                            }

                        }
                    }
                }

#ifdef _OPENMP
            }
#endif

            //#ifdef _OPENMP
            //#pragma omp parallel for default(shared)
            //            for(int comp=0;comp<Vec_in.size();comp++){
            //                for(int Np=0;Np<no_of_proc;Np++){
            //                    Vec_out[comp] += Vec_out_temp[Np][comp];
            //                }
            //            }

            //            for(int Np=0;Np<no_of_proc;Np++){
            //                vector < double_type >().swap(Vec_out_temp[Np]);
            //            }
            //#endif

        }


    }



}

void MODEL_1_orb_Hubb_chain::Act_diagonal_terms(BASIS_1_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    double EPS_=0.00000001;
    assert (Vec_out.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
    assert (Vec_in.size() == Vec_out.size());


    //Remember H[l][m]=<l|H|m>

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for (int i=0;i<basis.D_up_basis.size();i++){
        for (int j=0;j<basis.D_dn_basis.size();j++){
            int m;
            double_type value;
            m=basis.D_dn_basis.size()*i + j;

            value=0;
            //coulomb repulsion:
            //value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);
            for(int site=0;site<basis.Length;site++){
                value+=1.0*(Onsite_U[site])*
                         ( 1.0*( bit_value(basis.D_up_basis[i],site)*
                           bit_value(basis.D_dn_basis[j],site) )
                          );
                //  cout<<"site = "<<site<<" : "<<Onsite_Energy[site]<<endl;
            }


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

             //LongRange Direct Exchange Szi X Szj
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){

                    if(DirectExchange_mat[site_i][site_j]!=0.0){

                        value+=0.25*DirectExchange_mat[site_i][site_j]*(
                                    ( bit_value(basis.D_up_basis[i],site_i) - bit_value(basis.D_dn_basis[j],site_i))*
                                    ( bit_value(basis.D_up_basis[i],site_j) - bit_value(basis.D_dn_basis[j],site_j))
                                    );

                    }
                }
            }




            if(abs(value)>EPS_){
                Vec_out[m] +=value*one*Vec_in[m];
            }
        }
    }

}

void MODEL_1_orb_Hubb_chain::Read_parameters(BASIS_1_orb_Hubb_chain &basis, string filename){



    string filepath = filename;
    string pbc_,PBC_ ="PBC = ";
    string length_x, Length_X = "Length_X = ";
    string length_y, Length_Y = "Length_Y = ";
    string total_sites_, Total_Sites_ = "Total_Sites = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string ucoul, Ucoul = "U = ";

    string hmag, Hmag = "H_mag = ";

    string hopp_, Hopp_ = "Hopping = ";

    string nn_hopp_, NN_Hopp_ = "Hopping_NN = ";

    string geometry_, Geometry_ = "Geometry = ";

    string file_onsite_u_, File_Onsite_U_ = "File_Onsite_U = ";
    string file_onsite_energies_, File_Onsite_Energies_ = "File_Onsite_Energies = ";
    string file_hopping_connections_, File_Hopping_Connections_ = "File_Hopping_Connections = ";
    string file_nonlocal_int_connections_, File_NonLocal_Int_Connections_ = "File_NonLocal_Int_Connections = ";
    string file_directexchange_int_connections_, File_DirectExchange_Int_Connections_ = "File_Direct_Exchange_Connections = ";
    
    //Interaction_Assisted_Hopping, _Pair_Hopping
    string file_interactionassistedhopping_connections_, File_InteractionAssistedHopping_Connections_ = "File_Interaction_Assisted_Hopping_Connections = ";
    string file_pairhopping_connections_, File_PairHopping_Connections_ = "File_Pair_Hopping_Connections = ";


    string file_three_point_observation_, File_Three_Point_Observation_ = "File_Three_Point_Observation = ";
    string file_three_point_interaction_, File_Three_Point_Interaction_ = "File_Three_Point_Interaction = ";

    string processors_, Processors_ = "Processors = ";

    string read_onsite_energies;
    string read_onsite_u;



    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(Geometry_, 0)) != string::npos) {
                geometry_ = line.substr (offset+Geometry_.length());				}

            if ((offset = line.find(PBC_, 0)) != string::npos) {
                pbc_ = line.substr (offset+PBC_.length());				}

            if ((offset = line.find(Length_X, 0)) != string::npos) {
                length_x = line.substr (offset + Length_X.length());		}

            if ((offset = line.find(Length_Y, 0)) != string::npos) {
                length_y = line.substr (offset + Length_Y.length());		}

            if ((offset = line.find(Total_Sites_, 0)) != string::npos) {
                total_sites_ = line.substr (offset + Total_Sites_.length());		}

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

            if ((offset = line.find(NN_Hopp_, 0)) != string::npos) {
                nn_hopp_ = line.substr (offset+NN_Hopp_.length());				}

            if ((offset = line.find(File_Three_Point_Observation_, 0)) != string::npos) {
                file_three_point_observation_ = line.substr (offset+File_Three_Point_Observation_.length());}

            if ((offset = line.find(File_Three_Point_Interaction_, 0)) != string::npos) {
                file_three_point_interaction_ = line.substr (offset+File_Three_Point_Interaction_.length());}

            if ((offset = line.find(File_Onsite_U_, 0)) != string::npos) {
                file_onsite_u_ = line.substr (offset+File_Onsite_U_.length());          }

            if ((offset = line.find(File_Onsite_Energies_, 0)) != string::npos) {
                file_onsite_energies_ = line.substr (offset+File_Onsite_Energies_.length());				}

            if ((offset = line.find(File_Hopping_Connections_, 0)) != string::npos) {
                file_hopping_connections_ = line.substr (offset+File_Hopping_Connections_.length());				}

            if ((offset = line.find(File_NonLocal_Int_Connections_, 0)) != string::npos) {
                file_nonlocal_int_connections_ = line.substr (offset+File_NonLocal_Int_Connections_.length());				}

            if ((offset = line.find(File_DirectExchange_Int_Connections_, 0)) != string::npos) {
                file_directexchange_int_connections_ = line.substr (offset+File_DirectExchange_Int_Connections_.length());          }

             if ((offset = line.find(File_InteractionAssistedHopping_Connections_, 0)) != string::npos) {
                file_interactionassistedhopping_connections_ = line.substr (offset+File_InteractionAssistedHopping_Connections_.length());          }

            if ((offset = line.find(File_PairHopping_Connections_, 0)) != string::npos) {
                file_pairhopping_connections_ = line.substr (offset+File_PairHopping_Connections_.length());          }

            if ((offset = line.find(Processors_ , 0)) != string::npos) {
                processors_ = line.substr (offset+Processors_ .length());	}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}



    NProcessors_=atoi(processors_.c_str());


    if(pbc_ == "true"){
        PBC =true;
    }
    else{
        PBC=false;
    }

    int Length_X_int, Length_Y_int, Total_Sites_int;
    Length_X_int=atoi(length_x.c_str());
    Length_Y_int=atoi(length_y.c_str());

    Total_Sites_int = atoi(total_sites_.c_str());
    basis.Length=Total_Sites_int;
    basis.Ndn=atoi(ndn.c_str());
    basis.Nup=atoi(nup.c_str());



    // Mat_1_string three_point_oprs;
    // Mat_1_trio_int three_point_oprs_sites_set;



    int N_three_point_oprs;
    int n_sites_set;
    // stringstream _file_three_point_observation_(file_three_point_observation_);
    ifstream inputfile_three_point_observation(file_three_point_observation_.c_str());

    inputfile_three_point_observation>>N_three_point_oprs;
    three_point_oprs.resize(N_three_point_oprs);
    three_point_oprs_sites_set.resize(N_three_point_oprs);

    for(int n=0;n<three_point_oprs.size();n++){
        inputfile_three_point_observation>>three_point_oprs[n];
        inputfile_three_point_observation>> n_sites_set;
        three_point_oprs_sites_set[n].resize(n_sites_set);
        for(int m=0;m<n_sites_set;m++){
            three_point_oprs_sites_set[n][m].resize(3);
            for(int i=0;i<3;i++){
                inputfile_three_point_observation>>three_point_oprs_sites_set[n][m][i];
            }
        }
    }



    int N_three_point_intrs;
    ifstream inputfile_three_point_interaction(file_three_point_interaction_.c_str());

    inputfile_three_point_interaction>>N_three_point_intrs;
    three_point_intrs.resize(N_three_point_intrs);
    three_point_intrs_sites_set.resize(N_three_point_intrs);
    three_point_intrs_vals.resize(N_three_point_intrs);

    for(int n=0;n<three_point_intrs.size();n++){
        inputfile_three_point_interaction>>three_point_intrs[n];
        inputfile_three_point_interaction>> n_sites_set;
        three_point_intrs_sites_set[n].resize(n_sites_set);
        three_point_intrs_vals[n].resize(n_sites_set);
        for(int m=0;m<n_sites_set;m++){
            three_point_intrs_sites_set[n][m].resize(3);
            for(int i=0;i<3;i++){
                inputfile_three_point_interaction>>three_point_intrs_sites_set[n][m][i];
            }
            inputfile_three_point_interaction>>three_point_intrs_vals[n][m];
        }
    }





    stringstream _file_onsite_energies_(file_onsite_energies_);
    _file_onsite_energies_ >> read_onsite_energies;
    Onsite_Energy.resize(basis.Length);
    string filename_Onsite_Energy;
    string line_temp;

    // string temp_x_, temp_y_, temp_site_, Ener_val_ ;
    int temp_x, temp_y, temp_site;
    double Ener_val;

    if(read_onsite_energies == "true"){
        _file_onsite_energies_ >> filename_Onsite_Energy;

        cout<<"reading Onsite energies from '"<<filename_Onsite_Energy<<"'"<<endl;
        ifstream inputfile_Onsite_Energy(filename_Onsite_Energy.c_str());
        //getline(inputfile_Onsite_Energy,line_temp);

        //for(int iy=0;iy<Length_Y_int;iy++){
            //for(int ix=0;ix<Length_X_int;ix++){
        for(int i=0;i<basis.Length;i++){
                //inputfile_Onsite_Energy >> temp_x_ >> temp_y_ >> temp_site_ >> Ener_val_;
                //temp_x = atoi(temp_x_.c_str());
                //temp_y = atoi(temp_y_.c_str());
                //temp_site = atoi(temp_site_.c_str());
                //Ener_val = atof(Ener_val_.c_str());

                inputfile_Onsite_Energy >> temp_site >> Ener_val;
                cout<<temp_site<<"  "<<Ener_val<<endl;
                assert(temp_site==i);
                Onsite_Energy[temp_site]=Ener_val;
            }
        //}

    }
    else{
        for(int i=0;i<basis.Length;i++){
            Onsite_Energy[i]=0.0;
        }
    }


    U=atof(ucoul.c_str());
    stringstream _file_onsite_u_(file_onsite_u_);
    _file_onsite_u_ >> read_onsite_u;
    Onsite_U.resize(basis.Length);
    string filename_Onsite_U;

    // string temp_x_, temp_y_, temp_site_, Ener_val_ ;
    double U_val_temp;

    if(read_onsite_u == "true"){
        _file_onsite_u_ >> filename_Onsite_U;

        cout<<"reading Onsite Coulomb U from '"<<filename_Onsite_U<<"'"<<endl;
        ifstream inputfile_Onsite_U(filename_Onsite_U.c_str());
        //getline(inputfile_Onsite_Energy,line_temp);

        //for(int iy=0;iy<Length_Y_int;iy++){
        //for(int ix=0;ix<Length_X_int;ix++){
        for(int i=0;i<basis.Length;i++){
            //inputfile_Onsite_Energy >> temp_x_ >> temp_y_ >> temp_site_ >> Ener_val_;
            //temp_x = atoi(temp_x_.c_str());
            //temp_y = atoi(temp_y_.c_str());
            //temp_site = atoi(temp_site_.c_str());
            //Ener_val = atof(Ener_val_.c_str());

            inputfile_Onsite_U >> temp_site >> U_val_temp;
            cout<<temp_site<<"  "<<U_val_temp<<endl;
            assert(temp_site==i);
            Onsite_U[temp_site]=U_val_temp;
        }
        //}

    }
    else{
        for(int i=0;i<basis.Length;i++){
            Onsite_U[i]=U;
        }
    }








    double h;
    int temp_Tsites;
    stringstream _sstring_hmag_(hmag);
    //h=atof(hmag.c_str());
    _sstring_hmag_>> temp_Tsites;
    if(temp_Tsites!=basis.Length){
        cout<<"temp_Tsites = "<<temp_Tsites<<endl;
        cout<<"basis.Length = "<<basis.Length<<endl;
    }
    assert(temp_Tsites==basis.Length);
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        _sstring_hmag_>>h;
        H_field[i]=h;
    }

    cout<<"H field used-----------------"<<endl;
    for(int i=0;i<basis.Length;i++){
        cout<<H_field[i]<<"  ";
    }
    cout<<"------------------------"<<endl;


    double hopping_double;
    stringstream hopp_stream(hopp_);
    hopp_stream >> hopping_double;
    cout<<"NN hopping = "<<hopping_double<<endl;


    double nn_hopping_double;
    stringstream nn_hopp_stream(nn_hopp_);
    nn_hopp_stream >> nn_hopping_double;

    Hopping_mat_NN.clear();


    assert(geometry_ == "NN_chain" || geometry_ == "NN_2D_Lattice" || geometry_ == "NNN_2D_Lattice" || geometry_ == "LongRange");

    if(geometry_=="NN_chain"){

        Total_Sites_int=Length_X_int;
        cout<<"1 dimensional chain with nearest neighbour hopping is solved"<<endl;
        assert(Length_Y_int==1);
        Hopping_mat_NN.resize(basis.Length);
        for(int site=0;site<basis.Length;site++){
            Hopping_mat_NN[site].resize(basis.Length);
        }

        Mat_2_int Neighs_;
        Neighs_.resize(basis.Length);
        for(int site=0;site<basis.Length;site++){

            if(PBC){
                Neighs_[site].resize(2);

                // <---- "-x*(1)" direction
                Neighs_[site][0]=(site + basis.Length - 1)%(basis.Length);

                // ----> "+x*(1)" direction
                Neighs_[site][1]=(site + 1)%(basis.Length);
            }
            else{
                if(site==0){
                    Neighs_[site].resize(1);
                    // ----> "+x*(1)" direction
                    Neighs_[site][0]=(site + 1)%(basis.Length);
                }
                else if(site==(basis.Length -1)){
                    Neighs_[site].resize(1);
                    // <---- "-x*(1)" direction
                    Neighs_[site][0]=(site + basis.Length - 1)%(basis.Length);
                }
                else{
                    Neighs_[site].resize(2);
                    // <---- "-x*(1)" direction
                    Neighs_[site][0]=(site + basis.Length - 1)%(basis.Length);

                    // ----> "+x*(1)" direction
                    Neighs_[site][1]=(site + 1)%(basis.Length);
                }
            }
        }

        for(int site=0;site<basis.Length;site++){
            for(int neigh_no=0;neigh_no<Neighs_[site].size();neigh_no++){
                //Only upper diagonal part is created
                if(site > Neighs_[site][neigh_no] ){
                    Hopping_mat_NN[Neighs_[site][neigh_no]][site]=hopping_double;
                }
            }
        }
    }

    if(geometry_=="NN_2D_Lattice"){
        cout<<"2 dimensional (Lx)"<<Length_X_int<<" X (Ly)"<<Length_Y_int<< "Lattice with nearest neighbour hopping is solved"<<endl;


        //SITE LABELING XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//

        int Total_sites;
        Total_Sites_int = Length_X_int*Length_Y_int;
        basis.Length=Total_Sites_int;
        Total_sites = Total_Sites_int;

        vector<int> indx_, indy_;
        Mat_2_int Nc_;
        indx_.clear(); 	indx_.resize(Total_sites);
        indy_.clear();	indy_.resize(Total_sites);
        Nc_.resize(Length_X_int);
        for(int ix=0;ix<Length_X_int;ix++){
            Nc_[ix].resize(Length_Y_int);
        }

        int icount=0;
        for(int j=0;j<Length_Y_int;j++){
            for(int i=0;i<Length_X_int;i++){
                indx_[icount]=i;
                indy_[icount]=j;
                Nc_[i][j]=icount;
                icount++;
            }}
        //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//

        Hopping_mat_NN.resize(Total_sites);
        for(int site=0;site<Total_sites;site++){
            Hopping_mat_NN[site].resize(Total_sites);
        }


        Mat_2_int Neighs_;
        Neighs_.resize(basis.Length);
        for(int iy=0;iy<Length_Y_int;iy++){
            for(int ix=0;ix<Length_X_int;ix++){
                int site = Nc_[ix][iy];
                int x_new, y_new;

                if(PBC){

                    if(Length_X_int >1 && Length_Y_int >1){
                        // <---- "-x*(1)" direction
                        Neighs_[site].resize(1);
                        x_new = (ix + Length_X_int - 1)%(Length_X_int);
                        Neighs_[site][0]=Nc_[x_new][iy];

                        // ----> "+x*(1)" direction
                        Neighs_[site].resize(2);
                        x_new = (ix + 1)%(Length_X_int);
                        Neighs_[site][1]=Nc_[x_new][iy];

                        // <---- "-y*(1)" direction
                        Neighs_[site].resize(3);
                        y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                        Neighs_[site][2]=Nc_[ix][y_new];

                        // ----> "+y*(1)" direction
                        Neighs_[site].resize(4);
                        y_new = (iy + 1)%(Length_Y_int);
                        Neighs_[site][3]=Nc_[ix][y_new];
                    }
                    if(Length_X_int ==1 && Length_Y_int >1){

                        // <---- "-y*(1)" direction
                        Neighs_[site].resize(1);
                        y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                        Neighs_[site][0]=Nc_[ix][y_new];

                        // ----> "+y*(1)" direction
                        Neighs_[site].resize(2);
                        y_new = (iy + 1)%(Length_Y_int);
                        Neighs_[site][1]=Nc_[ix][y_new];
                    }
                    if(Length_X_int >1 && Length_Y_int ==1){
                        // <---- "-x*(1)" direction
                        Neighs_[site].resize(1);
                        x_new = (ix + Length_X_int - 1)%(Length_X_int);
                        Neighs_[site][0]=Nc_[x_new][iy];

                        // ----> "+x*(1)" direction
                        Neighs_[site].resize(2);
                        x_new = (ix + 1)%(Length_X_int);
                        Neighs_[site][1]=Nc_[x_new][iy];
                    }
                }
                else{
                    cout<<"Only PBC is allowed right now"<<endl;
                    assert(PBC);
                }
            }
        }

        for(int site=0;site<basis.Length;site++){
            for(int neigh_no=0;neigh_no<Neighs_[site].size();neigh_no++){
                //Only upper diagonal part is created
                if(site > Neighs_[site][neigh_no] ){
                    Hopping_mat_NN[Neighs_[site][neigh_no]][site]=hopping_double;
                }
            }
        }
    }

    if(geometry_=="NNN_2D_Lattice"){
        cout<<"2 dimensional (Lx)"<<Length_X_int<<" X (Ly)"<<Length_Y_int<< "Lattice with next nearest neighbour hopping is solved"<<endl;


        //SITE LABELING XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//
        int Total_sites;
        Total_Sites_int = Length_X_int*Length_Y_int;
        basis.Length=Total_Sites_int;
        Total_sites = Total_Sites_int;

        vector<int> indx_, indy_;
        Mat_2_int Nc_;
        indx_.clear(); 	indx_.resize(Total_sites);
        indy_.clear();	indy_.resize(Total_sites);
        Nc_.resize(Length_X_int);
        for(int ix=0;ix<Length_X_int;ix++){
            Nc_[ix].resize(Length_Y_int);
        }

        int icount=0;
        for(int j=0;j<Length_Y_int;j++){
            for(int i=0;i<Length_X_int;i++){
                indx_[icount]=i;
                indy_[icount]=j;
                Nc_[i][j]=icount;
                icount++;
            }}
        //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//

        Hopping_mat_NN.resize(Total_sites);
        for(int site=0;site<Total_sites;site++){
            Hopping_mat_NN[site].resize(Total_sites);
        }


        Mat_2_int Neighs_;
        Neighs_.resize(basis.Length);
        for(int iy=0;iy<Length_Y_int;iy++){
            for(int ix=0;ix<Length_X_int;ix++){
                int site = Nc_[ix][iy];
                int x_new, y_new;

                if(Length_X_int >1 && Length_Y_int >1){
                    // <---- "-x*(1)" direction
                    Neighs_[site].resize(1);
                    y_new=iy;
                    x_new = (ix + Length_X_int - 1)%(Length_X_int);
                    Neighs_[site][0]=Nc_[x_new][y_new];

                    // ----> "+x*(1)" direction
                    Neighs_[site].resize(2);
                    y_new=iy;
                    x_new = (ix + 1)%(Length_X_int);
                    Neighs_[site][1]=Nc_[x_new][y_new];

                    // <---- "-y*(1)" direction
                    Neighs_[site].resize(3);
                    x_new=ix;
                    y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                    Neighs_[site][2]=Nc_[x_new][y_new];

                    // ----> "+y*(1)" direction
                    Neighs_[site].resize(4);
                    x_new=ix;
                    y_new = (iy + 1)%(Length_Y_int);
                    Neighs_[site][3]=Nc_[x_new][y_new];

                    // "+x+y*(1)" direction
                    Neighs_[site].resize(5);
                    y_new = (iy + 1)%(Length_Y_int);
                    x_new = (ix + 1)%(Length_X_int);
                    Neighs_[site][4]=Nc_[x_new][y_new];

                    // "+x-y*(1)" direction
                    Neighs_[site].resize(6);
                    x_new = (ix + 1)%(Length_X_int);
                    y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                    Neighs_[site][5]=Nc_[x_new][y_new];

                    // "-x+y*(1)" direction
                    Neighs_[site].resize(7);
                    x_new = (ix + Length_X_int - 1)%(Length_X_int);
                    y_new = (iy + 1)%(Length_Y_int);
                    Neighs_[site][6]=Nc_[x_new][y_new];

                    // "-x-y*(1)" direction
                    Neighs_[site].resize(8);
                    x_new = (ix + Length_X_int - 1)%(Length_X_int);
                    y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                    Neighs_[site][7]=Nc_[x_new][y_new];

                }
                if(Length_X_int ==1 && Length_Y_int >1){

                    // <---- "-y*(1)" direction
                    Neighs_[site].resize(1);
                    y_new = (iy + Length_Y_int - 1)%(Length_Y_int);
                    Neighs_[site][0]=Nc_[ix][y_new];

                    // ----> "+y*(1)" direction
                    Neighs_[site].resize(2);
                    y_new = (iy + 1)%(Length_Y_int);
                    Neighs_[site][1]=Nc_[ix][y_new];
                }
                if(Length_X_int >1 && Length_Y_int ==1){
                    // <---- "-x*(1)" direction
                    Neighs_[site].resize(1);
                    x_new = (ix + Length_X_int - 1)%(Length_X_int);
                    Neighs_[site][0]=Nc_[x_new][iy];

                    // ----> "+x*(1)" direction
                    Neighs_[site].resize(2);
                    x_new = (ix + 1)%(Length_X_int);
                    Neighs_[site][1]=Nc_[x_new][iy];
                }
            }
        }

        int site_neigh;
        int neigh_x, neigh_y;
        int _x, _y;
        for(int site=0;site<basis.Length;site++){
            _x=indx_[site]; _y=indy_[site];
            for(int neigh_no=0;neigh_no<Neighs_[site].size();neigh_no++){
                site_neigh=Neighs_[site][neigh_no];
                neigh_x=indx_[site_neigh]; neigh_y=indy_[site_neigh];
                //Only upper diagonal part is created
                if(site > site_neigh ){
                    if(neigh_no<4){
                        if(PBC){
                            Hopping_mat_NN[site_neigh][site]=hopping_double;
                        }
                        else{
                            if(abs(neigh_x - _x)<2 && abs(neigh_y - _y)<2){
                                Hopping_mat_NN[site_neigh][site]=hopping_double;
                            }
                        }
                    }
                    else{
                        if(PBC){
                            Hopping_mat_NN[site_neigh][site]=nn_hopping_double;
                        }
                        else{
                            if(abs(neigh_x - _x)<2 && abs(neigh_y - _y)<2){
                                Hopping_mat_NN[site_neigh][site]=nn_hopping_double;
                            }

                        }

                    }
                }
            }
        }
    }


    if(geometry_ == "LongRange"){


        cout<<"Length_X, Length_Y, Hopping, Hopping_NN, and PBC from the input file are not used in LongRange, Only Total_Sites is used."<<endl;
        Total_Sites_int = atoi(total_sites_.c_str());
        basis.Length=Total_Sites_int;


        //Hoppings Mat(i,j)ci^{dagger}cj
        Hopping_mat_NN.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            Hopping_mat_NN[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_hopping_connections(file_hopping_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_hopping_connections>>Hopping_mat_NN[site_i][site_j];

                if(site_i>=site_j){
                    if(abs(Hopping_mat_NN[site_i][site_j])>0.000000001){
                        cout<<"ONLY UPPER TRIANGLE IS ALLOWED IN HOPPING"<<endl;
                        assert(false);
                    }
                }

            }
        }


        //Interactions  Mat(i,j)ninj
        NonLocalInteractions_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            NonLocalInteractions_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_nonlocal_int_connections(file_nonlocal_int_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_nonlocal_int_connections>>NonLocalInteractions_mat[site_i][site_j];
                if(site_i>=site_j){
                    if(abs(NonLocalInteractions_mat[site_i][site_j])>0.000000001){
                        cout<<"ONLY UPPER TRIANGLE IS ALLOWED IN NonLocalInteractions"<<endl;
                        assert(false);
                    }
                }
            }
        }

        //Direct Exchanges  Mat(i,j)Si.Sj
        DirectExchange_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            DirectExchange_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_directexchange_int_connections(file_directexchange_int_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_directexchange_int_connections>>DirectExchange_mat[site_i][site_j];
                if(site_i>=site_j){
                    if(abs(DirectExchange_mat[site_i][site_j])>0.000000001){
                        cout<<"ONLY UPPER TRIANGLE IS ALLOWED IN DirectExchange"<<endl;
                        assert(false);
                    }
                }
            }
        }


        //PairHopping  Mat(i,j) Pi* X Pj +h.c
        PairHopping_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            PairHopping_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_pairhopping_connections(file_pairhopping_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_pairhopping_connections>>PairHopping_mat[site_i][site_j];

                if(site_i>=site_j){
                    if(abs(PairHopping_mat[site_i][site_j])>0.000000001){
                        cout<<"ONLY UPPER TRIANGLE IS ALLOWED IN PairHopping"<<endl;
                        assert(false);
                    }
                }

            }
        }


        //InteractionAssistedHopping  Mat(i,j) { [n(i,s) X c(i,sbar)* X c(j,sbar)] + h.c.}
        InteractionAssistedHopping_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            InteractionAssistedHopping_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_interactionassistedhopping_connections(file_interactionassistedhopping_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_interactionassistedhopping_connections>>InteractionAssistedHopping_mat[site_i][site_j];
                cout<<InteractionAssistedHopping_mat[site_i][site_j]<<"  ";
            }
            cout<<endl;
        }



    }


    cout<<"PRINTING HOPPING MATRIX"<<endl;
    Print_Matrix(Hopping_mat_NN);

    cout<<""<<endl;
    cout<<"PRINTING Interaction MATRIX"<<endl;
    Print_Matrix(NonLocalInteractions_mat);
    cout<<"**************************"<<endl;




    /*

    Hopping_mat_NN.resize(1);
    Hopping_mat_NN[0].resize(1);
    //Hopping_mat_NN[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}
    stringstream hopp_stream(hopp_);
    hopp_stream >> Hopping_mat_NN[0][0];

*/





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

double_type MODEL_1_orb_Hubb_chain::Measure_one_point_operator_site_specific(string opr_type , Mat_1_doub &EigVec_, int site, BASIS_1_orb_Hubb_chain &basis){

    int orb=0;
    int spin;
    if(opr_type=="n_up"){
        spin=0;
    }
    if(opr_type=="n_dn"){
        spin=1;
    }
    assert(opr_type=="n_up" || opr_type =="n_dn");


    Mat_1_doub Vec_new;
    Vec_new.resize(EigVec_.size());

    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double_type value;
    double_type value_final;

    value_final=zero;
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
                value_final +=value*EigVec_[m]*conjugate(EigVec_[m]);
            }
        }
    }

    return value_final;
}

void MODEL_1_orb_Hubb_chain::Initialize_one_point_to_calculate(BASIS_1_orb_Hubb_chain &basis){



    one_point_obs.resize(2);
    one_point_obs[0]="n_up";
    one_point_obs[1]="n_dn";
    One_point_oprts.resize(2);



    int T_no_oprs=2;
    int orb;
    int spin;


    for(int i=0;i<T_no_oprs;i++){
        One_point_oprts[i].resize(basis.Length);
    }




    for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


        if(one_point_obs[opr_no]=="n_up" || one_point_obs[opr_no]=="n_dn"){
            orb=0;
        }

        if(one_point_obs[opr_no]=="n_up"){
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
        double_type value;


        for(int site=0;site<basis.Length;site++){

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
                        One_point_oprts[opr_no][site].value.push_back(value);
                        One_point_oprts[opr_no][site].rows.push_back(m);
                        One_point_oprts[opr_no][site].columns.push_back(m);
                    }
                }
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



    if(opr_type=="denden"){
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;
                value=0;
                value+=1.0*( ( bit_value(basis.D_up_basis[i], site) +
                               bit_value(basis.D_dn_basis[j], site) )*
                             ( bit_value(basis.D_up_basis[i], site2) +
                               bit_value(basis.D_dn_basis[j], site2) )
                             );
                if(value!=0){
                    OPR.value.push_back(value*one);
                    OPR.rows.push_back(m);
                    OPR.columns.push_back(m);
                }
            }
        }
    }

    if(opr_type=="SzSz"){
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

//                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
//                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                    i_new = basis.inverse_Dup[D_up - basis.DupMin_];
                   j_new = basis.inverse_Ddn[D_dn - basis.DdnMin_];

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site;
                    lp= site2;

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
                            ((bit_value(basis.D_up_basis[i], site2)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j], site2)==0)
                             )
                            )
                    {
                        OPR.value.push_back(one);
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

                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

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
                        OPR.value.push_back(one);
                        OPR.rows.push_back(m);
                        OPR.columns.push_back(m);
                    }
                }
            }
        }
    }




}








double_type MODEL_1_orb_Hubb_chain::Measure_two_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_, int site, int site2, BASIS_1_orb_Hubb_chain &basis){



    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double value;
    double_type value_final;
    value_final=zero;

    if(opr_type=="denden"){
        value_final=zero;
        for (int i=0;i<basis.D_up_basis.size();i++){
            for (int j=0;j<basis.D_dn_basis.size();j++){
                m=basis.D_dn_basis.size()*i + j;
                value=0;
                value+=1.0*( ( bit_value(basis.D_up_basis[i], site) +
                               bit_value(basis.D_dn_basis[j], site) )*
                             ( bit_value(basis.D_up_basis[i], site2) +
                               bit_value(basis.D_dn_basis[j], site2) )
                             );
                if(value!=0){
                    value_final += value*EigVec_[m]*conjugate(EigVec_[m]);
                }
            }
        }
    }

    if(opr_type=="SzSz"){
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
                    value_final += value*EigVec_[m]*conjugate(EigVec_[m]);
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

                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site;
                    lp= site2;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                    //assert(m_new<m);

                    value_final += sign_FM*EigVec_[m]*conjugate(EigVec_[m_new]);
                }

                if((site==site2)){
                    if(
                            ((bit_value(basis.D_up_basis[i], site2)==1)
                             &&
                             (bit_value(basis.D_dn_basis[j], site2)==0)
                             )
                            )
                    {

                        value_final += one*EigVec_[m]*conjugate(EigVec_[m]);
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

                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site2;
                    lp= site;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                    //assert(m_new<m);
                    value_final += sign_FM*EigVec_[m]*conjugate(EigVec_[m_new]);
                }

                if((site==site2)){
                    if(
                            ((bit_value(basis.D_up_basis[i], site2)==0)
                             &&
                             (bit_value(basis.D_dn_basis[j], site2)==1)
                             )
                            )
                    {
                        value_final += EigVec_[m]*conjugate(EigVec_[m]);
                    }
                }
            }
        }
    }



    return value_final;


}





double_type MODEL_1_orb_Hubb_chain::Measure_three_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_, int site1, int site2, int site3, BASIS_1_orb_Hubb_chain &basis){


    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double value;
    double_type value_final;
    value_final=zero;

    if(opr_type=="SzSpSm"){
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;

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

                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site3;
                    lp= site2;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                    value = sign_FM*(0.5*( bit_value(basis.D_up_basis[i_new], site1) -
                                           bit_value(basis.D_dn_basis[j_new], site1) ));

                    //assert(m_new<m);
                    value_final += value*EigVec_[m]*conjugate(EigVec_[m_new]);

                }


            }
        }
    }



    return value_final;

}


void MODEL_1_orb_Hubb_chain::Initialize_three_point_operator_sites_specific(string opr_type , Matrix_COO &OPR,
                                                                            int site1, int site2, int site3, BASIS_1_orb_Hubb_chain &basis){

    OPR.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    OPR.ncols = OPR.nrows;


    //Remember OPR[l][m]=<l|OPR|m>
    int m;
    double value;


    if(opr_type=="SzSpSm"){
        int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
        double sign_FM;

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

                    i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
                    j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);

                    m_new = basis.D_dn_basis.size()*i_new + j_new;

                    l= site3;
                    lp= site2;

                    sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                    sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                    sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


                    value = sign_FM*(0.5*( bit_value(basis.D_up_basis[i_new], site1) -
                                           bit_value(basis.D_dn_basis[j_new], site1) ));

                    //assert(m_new<m);
                    OPR.value.push_back(value*one);
                    OPR.rows.push_back(m_new);
                    OPR.columns.push_back(m);
                }


            }
        }
    }


}


void MODEL_1_orb_Hubb_chain::Initialize_two_point_to_calculate(BASIS_1_orb_Hubb_chain &basis){
    two_point_obs.resize(2);
    two_point_obs[0]="SzSz";
    two_point_obs[1]="SpSm";
    //two_point_obs[2]="SmSp";
    Two_point_oprts.resize(2);


    int T_no_oprs=2;



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

                                i_new = Find_int_in_intarray(D_up,basis.D_up_basis);
                                j_new = Find_int_in_intarray(D_dn,basis.D_dn_basis);

                                m_new = basis.D_dn_basis.size()*i_new + j_new;

                                l= site;
                                lp= site2;

                                sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
                                sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
                                sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);



                                //assert(m_new<m);

                                Two_point_oprts[opr_no][site][site2].value.push_back(sign_FM*one);
                                Two_point_oprts[opr_no][site][site2].rows.push_back(m_new);
                                Two_point_oprts[opr_no][site][site2].columns.push_back(m);


                            }

                            if((site==site2)){


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
                            Oprs_local[site].value.push_back(value*one);
                            Oprs_local[site].rows.push_back(m);
                            Oprs_local[site].columns.push_back(m);
                        }
                    }
                }

            }
            //local operators created ----------------------------------------------

            //In Momentum space for OBC only-----------------------------------------------------

            Matrix_COO temp;
            vector< int >().swap( Dyn_opr.columns );
            vector< int >().swap( Dyn_opr.rows );
            vector< double_type >().swap( Dyn_opr.value );

            if(!PBC){
                temp=Oprs_local[0];
                double_type value1, value2;
                for(int site=0;site<basis.Length-1;site++){

                    value2=one*sin((site+2)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                    if(site==0){
                        value1=one*sin((site+1)*Dyn_Momentum*PI)*sqrt(2.0/(basis.Length +1));
                        Sum(temp, Oprs_local[site+1], temp, value1, value2);}
                    else{
                        Sum(temp, Oprs_local[site+1], temp, 1.0, value2);
                    }

                }
            }
            else{
                double_type value1, value2;
                temp.value.clear();
                temp.rows.clear();
                temp.columns.clear();
                temp.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size() ;
                temp.ncols = temp.nrows;

                for(int site=0;site<basis.Length;site++){

#ifdef USE_COMPLEX
                    value2=exp(iota_comp*((Dyn_Momentum*PI*site) ))*sqrt(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                    cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif


                    Sum(temp, Oprs_local[site], temp, 1.0, value2);


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



                            if((Hopping_mat_NN[0][0])!=zero){

                                Opr.value.push_back(one*sign_FM*(Hopping_mat_NN[0][0]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*one*sign_FM*(Hopping_mat_NN[0][0]));
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



                            if((Hopping_mat_NN[0][0])!=zero){

                                Opr.value.push_back(one*sign_FM*(Hopping_mat_NN[0][0]));
                                Opr.rows.push_back((m_new));
                                Opr.columns.push_back((m));
                                Opr.value.push_back(-1.0*one*sign_FM*(Hopping_mat_NN[0][0]));
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

void MODEL_1_orb_Hubb_chain::Get_c_on_GS(Mat_1_doub & EigVec_, BASIS_1_orb_Hubb_chain & basis_Nm1, BASIS_1_orb_Hubb_chain & basis,
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
        orb_val = 0;
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


                    if(bit_value(basis.D_up_basis[i], site_val)==1){
                        l = site_val;

                        D_up_new = (int) (basis.D_up_basis[i] - pow(2,site_val) );
                        D_dn_new = basis.D_dn_basis[j];

                        i_up = Find_int_in_intarray_smartly(D_up_new,basis_Nm1.D_up_basis,basis_Nm1.partitions_up,basis_Nm1.Dup_val_at_partitions);
                        i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Nm1.D_dn_basis,basis_Nm1.partitions_dn,basis_Nm1.Ddn_val_at_partitions);

                        m_new = basis_Nm1.D_dn_basis.size()*i_up + i_dn;

                        max_up = basis.Length -1;
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

                    if(bit_value(basis.D_dn_basis[j],site_val)==1){
                        l = site_val;

                        D_dn_new = (int) (basis.D_dn_basis[j] - pow(2, site_val));
                        D_up_new = basis.D_up_basis[i];

                        i_up = Find_int_in_intarray_smartly(D_up_new,basis_Nm1.D_up_basis,basis_Nm1.partitions_up,basis_Nm1.Dup_val_at_partitions);
                        i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Nm1.D_dn_basis,basis_Nm1.partitions_dn,basis_Nm1.Ddn_val_at_partitions);

                        m_new = basis_Nm1.D_dn_basis.size()*i_up + i_dn;

                        max_dn = basis.Length -1;
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



void MODEL_1_orb_Hubb_chain::Get_cdagger_on_GS(Mat_1_doub & EigVec_, BASIS_1_orb_Hubb_chain & basis_Np1, BASIS_1_orb_Hubb_chain & basis,
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


                    if(bit_value(basis.D_up_basis[i],site_val)==0){
                        l = site_val;

                        D_up_new = (int) (basis.D_up_basis[i] + pow(2, site_val) );
                        D_dn_new = basis.D_dn_basis[j];

                        i_up = Find_int_in_intarray_smartly(D_up_new,basis_Np1.D_up_basis,basis_Np1.partitions_up,basis_Np1.Dup_val_at_partitions);
                        i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Np1.D_dn_basis,basis_Np1.partitions_dn,basis_Np1.Ddn_val_at_partitions);

                        m_new = basis_Np1.D_dn_basis.size()*i_up + i_dn;

                        max_up = basis.Length -1;
                        min_up=0;
                        sign_pow_up = one_bits_in_bw(min_up ,l, basis.D_up_basis[i]) ;
                        if(l != min_up){
                            sign_pow_up += bit_value(basis.D_up_basis[i],min_up);
                        }

                        sign_FM = pow(-1.0, sign_pow_up);

                        value = sign_FM*EigVec_[m]*(value_in);

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

                    if(bit_value(basis.D_dn_basis[j],site_val)==0){
                        l = site_val;

                        D_dn_new = (int) (basis.D_dn_basis[j] + pow(2, site_val) );
                        D_up_new = basis.D_up_basis[i];

                        i_up = Find_int_in_intarray_smartly(D_up_new,basis_Np1.D_up_basis,basis_Np1.partitions_up,basis_Np1.Dup_val_at_partitions);
                        i_dn = Find_int_in_intarray_smartly(D_dn_new,basis_Np1.D_dn_basis,basis_Np1.partitions_dn,basis_Np1.Ddn_val_at_partitions);

                        m_new = basis_Np1.D_dn_basis.size()*i_up + i_dn;

                        max_dn = basis.Length -1;
                        min_dn=0;
                        sign_pow_dn = one_bits_in_bw(min_dn ,l, basis.D_dn_basis[j]) ;
                        if(l != min_dn){
                            sign_pow_dn += bit_value(basis.D_dn_basis[j],min_dn);
                        }
                        sign_pow_dn += __builtin_popcount(D_up_new); //jump over all c^{\dagger}_up

                        sign_FM = pow(-1.0, sign_pow_dn);

                        value = sign_FM*EigVec_[m]*(value_in);

                        State_cdagger_on_GS[m_new] += value;
                    }
                }
            }
        }


    }



}

//#endif

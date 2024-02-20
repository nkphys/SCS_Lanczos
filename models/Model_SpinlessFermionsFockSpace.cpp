/*
This class includes the Model for which Lanczos is being done
*/

//#ifndef USE_COMPLEX
#include "Model_SpinlessFermionsFockSpace.h"
#include <stdlib.h>
using namespace std;
#define PI 3.14159265


void MODEL_SpinlessFermionsFockSpace::Add_diagonal_terms(BASIS_SpinlessFermionsFockSpace &basis){

    Hamil.nrows = basis.basis_size;
    Hamil.ncols = Hamil.nrows;


    //Remember H[l][m]=<l|H|m>
    int m;
    double value;
    for (unsigned long long int i=basis.DMin_;i<basis.DMax_+1;i++){
      
            value=0;
            
            //Onsite_Energy
            for(int site=0;site<basis.Length;site++){
                value+=1.0*(Onsite_Energy[site])*
                        (value_at_pos(i, site, basis.Base));
                //  cout<<"site = "<<site<<" : "<<Onsite_Energy[site]<<endl;
            }

            //LongRange interactions ninj
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){

                    if(abs(NonLocalInteractions_mat[site_i][site_j])>0.000000001){

                        value+=1.0*NonLocalInteractions_mat[site_i][site_j]*(
                                    ( value_at_pos(i, site_i, basis.Base))*
                                    ( value_at_pos(i, site_j, basis.Base))
                                    );

                    }
                }
            }


            if(value!=0){
                Hamil.value.push_back(value*one);
                Hamil.rows.push_back(i);
                Hamil.columns.push_back(i);
            }
        
    }

}

void MODEL_SpinlessFermionsFockSpace::Add_non_diagonal_terms(BASIS_SpinlessFermionsFockSpace &basis){


    bool Pairing_term=true;

    if(Pairing_term){
        int m;
        double_type value;
        
            int D_up, D_dn,i_new,i_new_temp,j_new,m_new, l, lp, sign_pow;
            double sign_FM;
            

for (unsigned long long int i=basis.DMin_;i<basis.DMax_+1;i++){
                for(int site2=0;site2<basis.Length;site2++){
                for(int site3=0;site3<basis.Length;site3++){

                    if(NonLocalPairingField_mat[site2][site3]!=0.0){


                        //c[site2]*c[site3]:

                        assert(site2!=site3);

                        if( (value_at_pos(i, site2, basis.Base)==1) &&
                            (value_at_pos(i, site3, basis.Base)==1)
                        )
                        {


                        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site2, basis.Base, 0);
                        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site3, basis.Base, 0);

                           assert(i_new<i);

                            l= site3;
                            lp= site2;

                            sign_pow = one_bits_in_bw(l,lp,i);
                            sign_FM = pow(-1.0, sign_pow);


                            value = sign_FM*NonLocalPairingField_mat[site2][site3];

                            //assert(m_new<m);
                            if(i_new<i){
                                if(abs(value)>0.000000000){
                                    Hamil.value.push_back(value*one);
                                    Hamil.rows.push_back((i_new));
                                    Hamil.columns.push_back((i));
                                }
                            }
                        }
            }
    }
    }
}
        
    }




}


void MODEL_SpinlessFermionsFockSpace::Add_connections(BASIS_SpinlessFermionsFockSpace &basis){

    double value;
    int m;
    int D_up,D_dn;
    int i_new, i_new_temp, j_new;
    int m_new;
    double sign_FM;
    int sign_pow;
    int l,lp;
   
   for (unsigned long long int i=basis.DMin_;i<basis.DMax_+1;i++){
            
            value=0;

            for(int site=0;site<basis.Length ;site++){
                for(int site_p=0;site_p<basis.Length ;site_p++){

                    if((Hopping_mat[site_p][site])!=zero)

                    {

                        //---------------Hopping-------------------//
                        //there have to be oneelectron in site
                        if( (value_at_pos(i, site, basis.Base)==1) &&
                            (value_at_pos(i, site_p, basis.Base)==0)
                        )
                        {

                        i_new_temp = Updated_decimal_with_value_at_pos<unsigned long long int>(i, site, basis.Base, 0);
                        i_new = Updated_decimal_with_value_at_pos<unsigned long long int>(i_new_temp, site_p, basis.Base, 1);

                        

                            l= site;
                            lp= site_p;

                            sign_pow = one_bits_in_bw(l,lp,i);
                            sign_FM = pow(-1.0, sign_pow);

                            if(i_new>=i){
                                cout<<"This Hopping is not required: "<<site_p<<"  "<< site<<"  "<<Hopping_mat[site_p][site]<<endl;
                            }
                            assert(i_new<i);
                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_mat[site_p][site])*one);
                            Hamil.rows.push_back((i_new));
                            Hamil.columns.push_back((i));


                        } // if hopping 


                    }


                }//site_p
            } // site
       
    } 

}





void MODEL_SpinlessFermionsFockSpace::Act_Hamil(BASIS_SpinlessFermionsFockSpace &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

 cout<<"NOT WORKING AT PRESENT"<<endl;
    // assert (Vec_in.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
    // Vec_out.clear();
    // Vec_out.resize(basis.D_up_basis.size()*basis.D_dn_basis.size());
    // Act_diagonal_terms(basis, Vec_in, Vec_out);
    // cout<<"Diagonal done"<<endl;
    // Act_non_diagonal_terms(basis, Vec_in, Vec_out);
    // cout<<"Non diagonal done"<<endl;
    // Act_connections(basis, Vec_in, Vec_out);
    // cout<<"Connections done"<<endl;

}





void MODEL_SpinlessFermionsFockSpace::Act_connections(BASIS_SpinlessFermionsFockSpace &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){



//     //Symmetrized, so that parallelization can be done efficiently.
//     Mat_2_doub Hopping_Mat;
//     Hopping_Mat.resize(basis.Length);
//     for(int sitei=0;sitei<basis.Length;sitei++){
//         Hopping_Mat[sitei].resize(basis.Length);
//         for(int sitej=0;sitej<basis.Length;sitej++){
//             if(sitej>=sitei){
//                 Hopping_Mat[sitei][sitej] = Hopping_mat[sitei][sitej];
//             }
//             else{
//                 Hopping_Mat[sitei][sitej] = conjugate(Hopping_mat[sitej][sitei]);
//             }
//         }
//     }

//     //    Mat_2_doub Vec_out_temp;
//     int no_of_proc;
//     no_of_proc=1;

// #ifdef _OPENMP
//     no_of_proc=min(basis.Length, NProcessors_);
//     omp_set_num_threads(no_of_proc);
//     cout<<"Connections acting: "<<no_of_proc<<" processors"<<endl;
//     //    Vec_out_temp.resize(no_of_proc);
//     //    for(int i=0;i<no_of_proc;i++){
//     //        Vec_out_temp[i].resize(Vec_in.size());
//     //    }
// #endif





//     //    for (int i=0;i<basis.D_up_basis.size();i++){
//     //        for (int j=0;j<basis.D_dn_basis.size();j++){
//     //            m=basis.D_dn_basis.size()*i + j;


//     for(int site=0;site<basis.Length ;site++){
//         for(int site_p=0;site_p<basis.Length ;site_p++){


//             if(abs(Hopping_Mat[site_p][site])>0.0000000001)

//             {



// #ifdef _OPENMP
// #pragma omp parallel
//                 {
// #endif


// #ifdef _OPENMP
// #pragma omp for nowait
// #endif
//                     for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){

//                         double value;
//                         int D_up,D_dn;
//                         int i_new,j_new;
//                         int m_new;
//                         double sign_FM;
//                         int sign_pow_up, sign_pow_dn;
//                         int l,lp;

//                         int mytid;
// #ifdef _OPENMP
//                         mytid = omp_get_thread_num();
// #endif

//                         int i,j;
//                         j = m%basis.D_dn_basis.size();
//                         i = int (m/basis.D_dn_basis.size());



//                         value=0;


//                         //---------------Hopping for up electrons-------------------//
//                         //there have to be one up electron in site
//                         //there have to be no up electron in site_p
//                         if(
//                                 (bit_value(basis.D_up_basis[i], site)==1)
//                                 &&
//                                 (bit_value(basis.D_up_basis[i], site_p)==0)
//                                 )
//                         {

//                             D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
//                                           - pow(2,site) );


//                             //i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
//                             i_new = basis.inverse_Dup[D_up - basis.DupMin_];
//                             j_new = j;

//                             m_new = basis.D_dn_basis.size()*i_new + j_new;

//                             l= site;
//                             lp= site_p;

//                             sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);

//                             sign_FM = pow(-1.0, sign_pow_up);

//                             //if(m_new>=m){
//                             //cout<<" Hopping: "<<site_p<<"  "<< site<<"  "<<Hopping_mat[site_p][site]<<endl;
//                             //}
//                             //assert(m_new<m);

//                             //HERE

// #ifdef _OPENMP
//                             // Vec_out_temp[mytid][m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat[site_p][site])*one;
//                             //Vec_out_temp[mytid][m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));
//                             Vec_out[m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));

// #else
//                             // Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat[site_p][site])*one;
//                             Vec_out[m] += Vec_in[m_new]*-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]));
// #endif


//                         } // if up hopping possible


//                         //---------------Hopping for dn electrons-------------------//
//                         //there have to be one dn electron in site
//                         //there have to be no dn electron in site_p
//                         if(
//                                 (bit_value(basis.D_dn_basis[j],site)==1)
//                                 &&
//                                 (bit_value(basis.D_dn_basis[j],site_p)==0)
//                                 )
//                         {

//                             D_dn = (int) (basis.D_dn_basis[j] + pow(2,site_p)
//                                           - pow(2,site) );


//                             //j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);
//                             j_new = basis.inverse_Ddn[D_dn - basis.DdnMin_];
//                             i_new = i;

//                             m_new = basis.D_dn_basis.size()*i_new + j_new;

//                             l= site;
//                             lp= site_p;

//                             sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);

//                             sign_FM = pow(-1.0, sign_pow_dn);

//                             //                            assert(m_new<m);

// #ifdef _OPENMP
//                             // Vec_out_temp[mytid][m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat[site_p][site])*one;
//                             //Vec_out_temp[mytid][m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(Hopping_Mat[site_p][site])*one);
//                             Vec_out[m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]))*one);
// #else
//                             //Vec_out[m_new] += Vec_in[m]*-1.0*sign_FM*(Hopping_mat[site_p][site])*one;
//                             Vec_out[m] += Vec_in[m_new]*conjugate(-1.0*sign_FM*(conjugate(Hopping_Mat[site_p][site]))*one);
// #endif


//                         } // if dn hopping possible



//                     }//m

// #ifdef _OPENMP
//                 }
// #endif


//             }//Hopping non-zero
//         }//site_p

//     } // site

//     //        }// "j" i.e dn_decimals
//     //    } // "i" i.e up_decimals










//     //#ifdef _OPENMP
//     //#pragma omp parallel for default(shared)
//     //    for(int comp=0;comp<Vec_in.size();comp++){
//     //        for(int Np=0;Np<no_of_proc;Np++){
//     //            Vec_out[comp] += Vec_out_temp[Np][comp];
//     //        }
//     //    }

//     //    for(int Np=0;Np<no_of_proc;Np++){
//     //        vector < double_type >().swap(Vec_out_temp[Np]);
//     //    }
//     //#endif





}


void MODEL_SpinlessFermionsFockSpace::Act_non_diagonal_terms(BASIS_SpinlessFermionsFockSpace &basis,  Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


//     for(int type_ind=0;type_ind<three_point_intrs.size();type_ind++){


//         if(three_point_intrs[type_ind]=="SzSpSm"){

//             Mat_2_doub Vec_out_temp;
//             int no_of_proc;
//             no_of_proc=1;
//             int temp_int= three_point_intrs_sites_set[type_ind].size();

// #ifdef _OPENMP
//             no_of_proc=min(temp_int, NProcessors_);
//             omp_set_num_threads(no_of_proc);
//             cout<<"Non diagonal acting: "<<no_of_proc<<" processors"<<endl;
//             //            Vec_out_temp.resize(no_of_proc);
//             //            for(int i=0;i<no_of_proc;i++){
//             //                Vec_out_temp[i].resize(Vec_in.size());
//             //            }
// #endif



// #ifdef _OPENMP
// #pragma omp parallel
//             {
// #endif


// #ifdef _OPENMP
// #pragma omp for nowait
// #endif
//                 for(int m=0;m<basis.D_up_basis.size()*basis.D_dn_basis.size();m++){

//                     int mytid;
// #ifdef _OPENMP
//                     mytid = omp_get_thread_num();
// #endif

//                     int i,j;
//                     j = m%basis.D_dn_basis.size();
//                     i = int (m/basis.D_dn_basis.size());


//                     for(int sites_set=0;sites_set<three_point_intrs_sites_set[type_ind].size();sites_set++){


//                         double_type value;
//                         int D_up, D_dn,i_new,j_new,m_new, l, lp, sign_pow_up , sign_pow_dn;
//                         double sign_FM;
//                         int site1, site2, site3;

//                         site1=three_point_intrs_sites_set[type_ind][sites_set][0];
//                         site2=three_point_intrs_sites_set[type_ind][sites_set][1];
//                         site3=three_point_intrs_sites_set[type_ind][sites_set][2];



//                         //Sz[site1]Sp[site2]*Sm[site3]:
//                         //there have to be ony up electron at site2
//                         //there have to be only down electron at site

//                         assert(site1!=site2);
//                         assert(site1!=site3);
//                         assert(site2!=site3);

//                         if(((bit_value(basis.D_dn_basis[j], site2)==1)
//                             &&
//                             (bit_value(basis.D_up_basis[i], site2)==0)
//                             )
//                                 &&
//                                 ((bit_value(basis.D_up_basis[i], site3)==1)
//                                  &&
//                                  (bit_value(basis.D_dn_basis[j], site3)==0)
//                                  ))
//                         {

//                             D_up = (int) (basis.D_up_basis[i] - pow(2, site3)
//                                           + pow(2, site2) );
//                             D_dn = (int) (basis.D_dn_basis[j] + pow(2, site3)
//                                           - pow(2, site2) );

//                             //i_new = Find_int_in_intarray_smartly(D_up,basis.D_up_basis,basis.partitions_up,basis.Dup_val_at_partitions);
//                             //j_new = Find_int_in_intarray_smartly(D_dn,basis.D_dn_basis,basis.partitions_dn,basis.Ddn_val_at_partitions);
//                             i_new = basis.inverse_Dup[D_up - basis.DupMin_];
//                             j_new = basis.inverse_Ddn[D_dn - basis.DdnMin_];

//                             m_new = basis.D_dn_basis.size()*i_new + j_new;

//                             l= site3;
//                             lp= site2;

//                             sign_pow_up = one_bits_in_bw(l,lp,basis.D_up_basis[i]);
//                             sign_pow_dn = one_bits_in_bw(l,lp,basis.D_dn_basis[j]);
//                             sign_FM = pow(-1.0, sign_pow_up + sign_pow_dn+1);


//                             value = sign_FM*(0.5*( bit_value(basis.D_up_basis[i_new], site1) -
//                                                    bit_value(basis.D_dn_basis[j_new], site1) ))*three_point_intrs_vals[type_ind][sites_set];

//                             //assert(m_new<m);
//                             if(abs(value)>0.00000000){

// #ifdef _OPENMP
//                                 //Vec_out_temp[mytid][m_new] += Vec_in[m]*value*one;
//                                 Vec_out[m] += Vec_in[m_new]*conjugate(value)*one;
// #else
//                                 Vec_out[m_new] += Vec_in[m]*value*one;
// #endif


//                             }

//                         }
//                     }
//                 }

// #ifdef _OPENMP
//             }
// #endif

//             //#ifdef _OPENMP
//             //#pragma omp parallel for default(shared)
//             //            for(int comp=0;comp<Vec_in.size();comp++){
//             //                for(int Np=0;Np<no_of_proc;Np++){
//             //                    Vec_out[comp] += Vec_out_temp[Np][comp];
//             //                }
//             //            }

//             //            for(int Np=0;Np<no_of_proc;Np++){
//             //                vector < double_type >().swap(Vec_out_temp[Np]);
//             //            }
//             //#endif

//         }


//     }



}

void MODEL_SpinlessFermionsFockSpace::Act_diagonal_terms(BASIS_SpinlessFermionsFockSpace &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


//     assert (Vec_out.size() == basis.D_up_basis.size()*basis.D_dn_basis.size());
//     assert (Vec_in.size() == Vec_out.size());


//     //Remember H[l][m]=<l|H|m>

// #ifdef _OPENMP
// #pragma omp parallel for default(shared)
// #endif
//     for (int i=0;i<basis.D_up_basis.size();i++){
//         for (int j=0;j<basis.D_dn_basis.size();j++){
//             int m;
//             double value;
//             m=basis.D_dn_basis.size()*i + j;

//             value=0;
//             //coulomb repulsion:
//             value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);



//             //magnetic Field
//             for(int site=0;site<basis.Length;site++){
//                 value+=0.5*(H_field[site])*
//                         ( ( bit_value(basis.D_up_basis[i],site) -
//                             bit_value(basis.D_dn_basis[j],site) )
//                           );
//             }

//             //Onsite_Energy
//             for(int site=0;site<basis.Length;site++){
//                 value+=1.0*(Onsite_Energy[site])*
//                         ( ( bit_value(basis.D_up_basis[i],site) +
//                             bit_value(basis.D_dn_basis[j],site) )
//                           );
//                 //  cout<<"site = "<<site<<" : "<<Onsite_Energy[site]<<endl;
//             }

//             //LongRange interactions ninj
//             for(int site_i=0;site_i<basis.Length;site_i++){
//                 for(int site_j=0;site_j<basis.Length;site_j++){

//                     if(NonLocalInteractions_mat[site_i][site_j]!=0.0){

//                         value+=1.0*NonLocalInteractions_mat[site_i][site_j]*(
//                                     ( bit_value(basis.D_up_basis[i],site_i) + bit_value(basis.D_dn_basis[j],site_i))*
//                                     ( bit_value(basis.D_up_basis[i],site_j) + bit_value(basis.D_dn_basis[j],site_j))
//                                     );

//                     }
//                 }
//             }




//             if(value!=0){
//                 Vec_out[m] +=value*one*Vec_in[m];
//             }
//         }
//     }

}

void MODEL_SpinlessFermionsFockSpace::Read_parameters(BASIS_SpinlessFermionsFockSpace &basis, string filename){



    string filepath = filename;
    
    string total_sites_, Total_Sites_ = "Total_Sites = ";
    string geometry_, Geometry_ = "Geometry = ";

    string file_onsite_energies_, File_Onsite_Energies_ = "File_Onsite_Energies = ";
    string file_hopping_connections_, File_Hopping_Connections_ = "File_Hopping_Connections = ";
    string file_nonlocal_int_connections_, File_NonLocal_Int_Connections_ = "File_NonLocal_Int_Connections = ";
    string file_pairing_connections_, File_Pairing_Connections_ = "File_Pairing_Connections = ";
    
    string processors_, Processors_ = "Processors = ";

    string read_onsite_energies;



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

            if ((offset = line.find(Total_Sites_, 0)) != string::npos) {
                total_sites_ = line.substr (offset + Total_Sites_.length());		}

            if ((offset = line.find(File_Onsite_Energies_, 0)) != string::npos) {
                file_onsite_energies_ = line.substr (offset+File_Onsite_Energies_.length());				}

            if ((offset = line.find(File_Hopping_Connections_, 0)) != string::npos) {
                file_hopping_connections_ = line.substr (offset+File_Hopping_Connections_.length());				}

            if ((offset = line.find(File_NonLocal_Int_Connections_, 0)) != string::npos) {
                file_nonlocal_int_connections_ = line.substr (offset+File_NonLocal_Int_Connections_.length());				}

            if ((offset = line.find(File_Pairing_Connections_, 0)) != string::npos) {
                file_pairing_connections_ = line.substr (offset+File_Pairing_Connections_.length());				}

            if ((offset = line.find(Processors_ , 0)) != string::npos) {
                processors_ = line.substr (offset+Processors_ .length());	}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    NProcessors_=atoi(processors_.c_str());


    int  Total_Sites_int;
    Total_Sites_int = atoi(total_sites_.c_str());
    basis.Length=Total_Sites_int;
  

    stringstream _file_onsite_energies_(file_onsite_energies_);
    _file_onsite_energies_ >> read_onsite_energies;
    Onsite_Energy.resize(basis.Length);
    string filename_Onsite_Energy;
    string line_temp;

    // string temp_x_, temp_y_, temp_site_, Ener_val_ ;
    int temp_site;
    double Ener_val;
    if(read_onsite_energies == "true"){
        _file_onsite_energies_ >> filename_Onsite_Energy;

        cout<<"reading Onsite energies from '"<<filename_Onsite_Energy<<"'"<<endl;
        ifstream inputfile_Onsite_Energy(filename_Onsite_Energy.c_str());
        getline(inputfile_Onsite_Energy,line_temp);

        for(int i_=0;i_<Total_Sites_int;i_++){
                inputfile_Onsite_Energy >> temp_site >> Ener_val;
                Onsite_Energy[temp_site]=Ener_val;
        }
    }
    else{
        for(int i=0;i<basis.Length;i++){
            Onsite_Energy[i]=0.0;
        }
    }




    Hopping_mat.clear();

    assert(geometry_ == "LongRange");

    if(geometry_ == "LongRange"){

        //Hoppings Mat(i,j)ci^{dagger}cj
        Hopping_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            Hopping_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_hopping_connections(file_hopping_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_hopping_connections>>Hopping_mat[site_i][site_j];
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
            }
        }

        //Pairingfield  Mat(i,j)[c*jXc*i + ciXcj]
        NonLocalPairingField_mat.resize(Total_Sites_int);
        for(int site_=0;site_<Total_Sites_int;site_++){
            NonLocalPairingField_mat[site_].resize(Total_Sites_int);
        }

        ifstream inputfile_pairing_connections(file_pairing_connections_.c_str());
        for(int site_i=0;site_i<Total_Sites_int;site_i++){
            for(int site_j=0;site_j<Total_Sites_int;site_j++){
                inputfile_pairing_connections>>NonLocalPairingField_mat[site_i][site_j];
            }
        }



    }


    cout<<"PRINTING HOPPING MATRIX"<<endl;
    Print_Matrix(Hopping_mat);

    cout<<""<<endl;
    cout<<"PRINTING Interaction MATRIX"<<endl;
    Print_Matrix(NonLocalInteractions_mat);
    cout<<"**************************"<<endl;


    cout<<""<<endl;
    cout<<"PRINTING PairingField MATRIX"<<endl;
    Print_Matrix(NonLocalPairingField_mat);
    cout<<"**************************"<<endl;


    /*

    Hopping_mat.resize(1);
    Hopping_mat[0].resize(1);
    //Hopping_mat[alpha][beta] comes in front of c^{\dagger}_{alpha\sigma}c_{beta\sigma}
    stringstream hopp_stream(hopp_);
    hopp_stream >> Hopping_mat[0][0];

*/





}

void MODEL_SpinlessFermionsFockSpace::Read_parameters_for_dynamics(string filename){

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

void MODEL_SpinlessFermionsFockSpace::Initialize_one_point_operator_site_specific(string opr_type , Matrix_COO &OPR, int site, BASIS_SpinlessFermionsFockSpace &basis){

    // int orb=0;
    // int spin;
    // if(opr_type=="n_up"){
    //     spin=0;
    // }
    // if(opr_type=="n_dn"){
    //     spin=1;
    // }
    // assert(opr_type=="n_up" || opr_type =="n_dn");


    // OPR.nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    // OPR.ncols = OPR.nrows;

    // //Remember OPR[l][m]=<l|OPR|m>
    // int m;
    // double_type value;


    // for (int i=0;i<basis.D_up_basis.size();i++){
    //     for (int j=0;j<basis.D_dn_basis.size();j++){
    //         m=basis.D_dn_basis.size()*i + j;

    //         //n_orb_spin[site]:
    //         if(spin==0){
    //             value=one*(1.0*bit_value(basis.D_up_basis[i],orb*basis.Length + site));
    //         }
    //         else{
    //             value=one*(1.0*bit_value(basis.D_dn_basis[j],orb*basis.Length + site));
    //         }

    //         if(value!=zero){
    //             OPR.value.push_back(value);
    //             OPR.rows.push_back(m);
    //             OPR.columns.push_back(m);
    //         }
    //     }
    // }
}

double_type MODEL_SpinlessFermionsFockSpace::Measure_one_point_operator_site_specific(string opr_type , Mat_1_doub &EigVec_, int site, BASIS_SpinlessFermionsFockSpace &basis){

    // int orb=0;
    // int spin;
    // if(opr_type=="n_up"){
    //     spin=0;
    // }
    // if(opr_type=="n_dn"){
    //     spin=1;
    // }
    // assert(opr_type=="n_up" || opr_type =="n_dn");


    // Mat_1_doub Vec_new;
    // Vec_new.resize(EigVec_.size());

    // //Remember OPR[l][m]=<l|OPR|m>
    // int m;
    // double_type value;
    // double_type value_final;

    // value_final=zero;
    // for (int i=0;i<basis.D_up_basis.size();i++){
    //     for (int j=0;j<basis.D_dn_basis.size();j++){
    //         m=basis.D_dn_basis.size()*i + j;

    //         //n_orb_spin[site]:
    //         if(spin==0){
    //             value=one*(1.0*bit_value(basis.D_up_basis[i],orb*basis.Length + site));
    //         }
    //         else{
    //             value=one*(1.0*bit_value(basis.D_dn_basis[j],orb*basis.Length + site));
    //         }

    //         if(value!=zero){
    //             value_final +=value*EigVec_[m]*conjugate(EigVec_[m]);
    //         }
    //     }
    // }

    // return value_final;
}

void MODEL_SpinlessFermionsFockSpace::Initialize_one_point_to_calculate(BASIS_SpinlessFermionsFockSpace &basis){



    // one_point_obs.resize(1);
    // one_point_obs[0]="n";
    // One_point_oprts.resize(1);



    // int T_no_oprs=1;
    // int orb;
    // int spin;


    // for(int i=0;i<T_no_oprs;i++){
    //     One_point_oprts[i].resize(basis.Length);
    // }




    // for(int opr_no=0;opr_no<T_no_oprs;opr_no++){


    //     for(int site=0;site<basis.Length;site++){
    //         One_point_oprts[opr_no][site].nrows = basis.D_up_basis.size()*basis.D_dn_basis.size();
    //         One_point_oprts[opr_no][site].ncols = One_point_oprts[opr_no][site].nrows;
    //     }


    //     //Remember OPR[l][m]=<l|OPR|m>
    //     int m;
    //     double_type value;


    //     for(int site=0;site<basis.Length;site++){

    //         for (int i=0;i<basis.D_up_basis.size();i++){
    //             for (int j=0;j<basis.D_dn_basis.size();j++){
    //                 m=basis.D_dn_basis.size()*i + j;

    //                 //n_orb_spin[site]:
    //                 if(spin==0){
    //                     value=one*(1.0*bit_value(basis.D_up_basis[i],orb*basis.Length + site));
    //                 }
    //                 else{
    //                     value=one*(1.0*bit_value(basis.D_dn_basis[j],orb*basis.Length + site));
    //                 }





    //                 if(value!=zero){
    //                     One_point_oprts[opr_no][site].value.push_back(value);
    //                     One_point_oprts[opr_no][site].rows.push_back(m);
    //                     One_point_oprts[opr_no][site].columns.push_back(m);
    //                 }
    //             }
    //         }

    //     }

    // }


}

void MODEL_SpinlessFermionsFockSpace::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR,
                                                                          int site, int site2, BASIS_SpinlessFermionsFockSpace &basis){

 

}








double_type MODEL_SpinlessFermionsFockSpace::Measure_two_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_, int site, int site2, BASIS_SpinlessFermionsFockSpace &basis){



}





double_type MODEL_SpinlessFermionsFockSpace::Measure_three_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_, int site1, int site2, int site3, BASIS_SpinlessFermionsFockSpace &basis){
   

}


void MODEL_SpinlessFermionsFockSpace::Initialize_three_point_operator_sites_specific(string opr_type , Matrix_COO &OPR,
                                                                            int site1, int site2, int site3, BASIS_SpinlessFermionsFockSpace &basis){


}


void MODEL_SpinlessFermionsFockSpace::Initialize_two_point_to_calculate(BASIS_SpinlessFermionsFockSpace &basis){
   


}

void MODEL_SpinlessFermionsFockSpace::Initialize_Opr_for_Dynamics(BASIS_SpinlessFermionsFockSpace &basis){


}

void MODEL_SpinlessFermionsFockSpace::Get_c_on_GS(Mat_1_doub & EigVec_, BASIS_SpinlessFermionsFockSpace & basis_Nm1, BASIS_SpinlessFermionsFockSpace & basis,
                                         Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


}



void MODEL_SpinlessFermionsFockSpace::Get_cdagger_on_GS(Mat_1_doub & EigVec_, BASIS_SpinlessFermionsFockSpace & basis_Np1, BASIS_SpinlessFermionsFockSpace & basis,
                                               Mat_1_trio_int TRIO_VEC, Mat_1_doub values){


}

//#endif

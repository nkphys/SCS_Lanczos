/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "assert.h"
#include "Basis_1_orb_tJ.h"
using namespace std;


void BASIS_1_orb_tJ::Construct_basis(){

    //Mat_1_int temp=decimal_to_binary(2147483647);
    int d_up_max,d_up_min;
    int d_dn_max,d_dn_min;

    //int d_up_old;
    n_orb=1;
    //Calculating min and max decimal_{up,dn}


    cout<<"here"<<endl;

    //putting correct D_'s in the D arrays
    int numup, numdn;


    int N_up_local, N_dn_local , N_local;
    bool allowed=false;

    D_up_basis.clear();
    D_dn_basis.clear();
    Nup_offsets.clear();

    partitions_up.clear();
    Dup_val_at_partitions.clear();

    Nup_offsets.resize(N_total+1);
    D_updn_reverse.resize(N_total+1);



    //   D_up_reverse.resize(N_total+1);
    //   D_dn_reverse.resize(N_total+1);
    D_up_min.resize(N_total+1);D_up_max.resize(N_total+1);
    D_dn_min.resize(N_total+1);D_dn_max.resize(N_total+1);

    int basis_count=0;

    // for(int Nup=max(0, N_total - n_orb*Length);Nup<=min( n_orb*Length, N_total);Nup++){

    if(Nup==Nup){

        cout<<"Restricted basis constructionfor Nup = "<<Nup<<
              "sector is being done"<<endl;
        // Nup_offsets[Nup].first = D_up_basis.size();

        d_up_min=0;
        for (int i=0;i<Nup;i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_min=d_up_min+ temp ;
        }
        D_up_min[Nup]=d_up_min;

        d_up_max=0;
        for (int i=(n_orb*Length - Nup);i<(n_orb*Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_max=d_up_max+ temp;
        }
        D_up_max[Nup]=d_up_max;

        //  cout <<d_up_min<<"   "<<d_up_max<<endl;




        d_dn_min=0;
        for (int i=0;i<Ndn;i++){
            int temp = (int)(pow(2,i)+0.5);
            d_dn_min=d_dn_min+ temp;
        }
        D_dn_min[Ndn]=d_dn_min;


        d_dn_max=0;
        for (int i=(n_orb*Length - Ndn);i<(n_orb*Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_dn_max=d_dn_max+ temp;
        }
        D_dn_max[Ndn]=d_dn_max;





        //        D_updn_reverse[Nup].resize(d_up_max-d_up_min+1);
        //        for(int i=0;i<d_up_max-d_up_min+1;i++){
        //            D_updn_reverse[Nup][i].resize(d_dn_max-d_dn_min+1);
        //        }



        char state_str[500];
        int length_state_str;
        Mat_1_ullint Dec_vec;
        Dec_vec.clear();
        fromDeci(state_str, 2, d_up_max);

        stringstream ss_temp;
        string state_actual_str;
        ss_temp << state_str;
        ss_temp >> state_actual_str;
        cout<<state_actual_str<<"  |  "<< d_up_max<<endl;

        length_state_str = Length; //strlen(state_str);
        findPermutations(state_str, 0, length_state_str, Dec_vec, 2);
        cout <<"No. of up basis = "<<Dec_vec.size()<<endl;







        Mat_1_ullint Dec_vec_dn;
        Dec_vec_dn.clear();
        fromDeci(state_str, 2, d_dn_max);

        stringstream ss_temp_dn;
        string state_actual_str_dn;
        ss_temp_dn << state_str;
        ss_temp_dn >> state_actual_str_dn;
        cout<<state_actual_str_dn<<"  |  "<< d_dn_max<<endl;

        length_state_str = Length; //strlen(state_str);
        findPermutations(state_str, 0, length_state_str, Dec_vec_dn, 2);
        cout <<"No. of dn basis = "<<Dec_vec_dn.size()<<endl;



        int dup_old=-100;
        int dup;
        int ddn;
        partitions_up.push_back(0);
        for(int dup_ind=Dec_vec.size()-1;dup_ind>=0;dup_ind--){

            dup=Dec_vec[dup_ind];
            //
            //numup= __builtin_popcount (dup);
            //assert( numup == Nup );

            //for(int ddn_ind=0;ddn_ind<Dec_vec_dn.size();ddn_ind++){

            assert(Nup==(Length-Ndn));

            ddn = (int)(pow(2,Length)+0.5) - dup -1;

                //ddn = Dec_vec_dn[ddn_ind];
                //numdn= __builtin_popcount (ddn);
                //assert( numdn == Ndn );

                    //check dup and ddn satisfy conditions
                    bool check;
                    allowed=true;
                    for(int site=0;site<Length;site++){
                        N_up_local=0;
                        N_dn_local=0;

                        for(int gamma=0;gamma<n_orb;gamma++){
                            N_up_local+=bit_value(dup,gamma*Length + site);
                            N_dn_local+=bit_value(ddn,gamma*Length + site);
                        }
                        N_local=N_up_local+N_dn_local;
                        // if(N_local == 1){
                        //   cout<<N_local;
                        //}
                        //if(N_local==4 || N_local==5 || N_local==3){
                        if( N_local==0 || N_local==1 ){
                            check=true;
                            //cout<<dup<<","<<ddn<<"allowed"<<endl;
                        }
                        else{
                            check=false;
                            //cout<<dup<<","<<ddn<<"not allowed"<<endl;
                        }
                        allowed = allowed && check;
                    }
                    //temp_pair.first = dup;
                    //temp_pair.second = ddn;
                    if(allowed==true){
                        D_up_basis.push_back(dup);
                        D_dn_basis.push_back(ddn);

                        if(dup!=dup_old){
                            Dup_val_at_partitions.push_back(dup); //if at ith pos

                            if(D_up_basis.size()>1){
                                partitions_up.push_back(D_up_basis.size()-1); //search  p[i]<= and <=p[i+1]-1
                            }
                            dup_old = dup;
                        }

                        //D_updn_reverse[Nup][dup-d_up_min][ddn-d_dn_min]=basis_count;
                        basis_count++;
                    }




            //}


            if(dup%10000==0){
                cout<<"dup ="<<dup<<" done"<<endl;
            }

        }
        partitions_up.push_back(D_up_basis.size());

    }


    /*
        for(int b=0;b<D_up_basis.size();b++){
        cout<<"-----------------------"<<endl;
        cout<<"UP spins:"<<endl;
        print_binary_of_decimal(D_up_basis[b]);
        cout<<"DN spins:"<<endl;
        print_binary_of_decimal(D_dn_basis[b]);
        }
        */

    cout<<"Total basis = "<<basis_count<<endl;

    //assert(false);
    //}
}



void BASIS_1_orb_tJ::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_multi_orb_Hubbard_chain_GC_restricted.h"
using namespace std;

/*
1)  for "up-spin" basis
              [_______________________  _  ] [_______________________  _  ]
    site----->[012....................(L-1)] [012....................(L-1)]
    orbital-->[.........orb - "0"..........] [.........orb - "1"..........]

2)  similarly for "down-spin" basis

4) orb_0 = d_xz, orb_1 = d_yz, orb_2 = d_xy
*/

void BASIS_multi_orb_Hubb_chain_GC_Restricted::Construct_basis(){
    Restricted=true;

    int d_up_max,d_up_min;
    int d_dn_max,d_dn_min;

    
    //Calculating min and max decimal_{up,dn}



    //putting correct D_'s in the D arrays
    int numup, numdn;
    int Ndn;

    int N_up_local, N_dn_local , N_local;
    bool allowed=false;

    D_up_basis.clear();
    D_dn_basis.clear();
    Nup_offsets.clear();

    Nup_offsets.resize(N_total+1);
    D_updn_reverse.resize(N_total+1);


    Mat_1_intpair D_up_D_dn_allowed;
    pair_int temp_pair;

    Canonical_partition_up.resize(N_total+1);
    Canonical_partition_dn.resize(N_total+1);

 //   D_up_reverse.resize(N_total+1);
 //   D_dn_reverse.resize(N_total+1);
    D_up_min.resize(N_total+1);D_up_max.resize(N_total+1);
    D_dn_min.resize(N_total+1);D_dn_max.resize(N_total+1);

    int basis_count=0;

    for(int Nup=max(0, N_total - N_orb*Length);Nup<=min( N_orb*Length, N_total);Nup++){

        Ndn = N_total - Nup;
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
        for (int i=(N_orb*Length - Nup);i<(N_orb*Length);i++){
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
        for (int i=(N_orb*Length - Ndn);i<(N_orb*Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_dn_max=d_dn_max+ temp;
        }
        D_dn_max[Ndn]=d_dn_max;

    //cout <<d_dn_min<<"   "<<d_dn_max<<endl;

     //  D_up_reverse[Nup].resize(d_up_max-d_up_min+1);
     //  D_dn_reverse[Ndn].resize(d_dn_max-d_dn_min+1);
       D_updn_reverse[Nup].resize(d_up_max-d_up_min+1);
       for(int i=0;i<d_up_max-d_up_min+1;i++){
           D_updn_reverse[Nup][i].resize(d_dn_max-d_dn_min+1);
       }


        Canonical_partition_up[Nup].clear();
        Canonical_partition_dn[Nup].clear();


        for(int dup=d_up_min;dup<=d_up_max;dup++){
            numup= __builtin_popcount (dup);
            if( numup == Nup ){

                for(int ddn=d_dn_min;ddn<=d_dn_max;ddn++){
                    numdn= __builtin_popcount (ddn);
                    if( numdn == Ndn ){

                        //check dup and ddn satisfy conditions
                        bool check;
                        allowed=true;
                        for(int site=0;site<Length;site++){
                            N_up_local=0;
                            N_dn_local=0;

                            for(int gamma=0;gamma<3;gamma++){
                        N_up_local+=bit_value(dup,gamma*Length + site);
                        N_dn_local+=bit_value(ddn,gamma*Length + site);
                            }
                            N_local=N_up_local+N_dn_local;
                           // if(N_local == 1){
                             //   cout<<N_local;
                            //}
                            //if(N_local==4 || N_local==5 || N_local==3){
                            if( Is_int_in_array(N_local,Local_occupations_allowed) ){
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
                        //D_up_D_dn_allowed.push_back(temp_pair);
                       // D_up_reverse[Nup][dup-d_up_min].push_back(basis_count);
                       // D_dn_reverse[Ndn][ddn-d_dn_min].push_back(basis_count);
                        D_updn_reverse[Nup][dup-d_up_min][ddn-d_dn_min]=basis_count;
                        basis_count++;
                        }



                    }
                }
            }

        }







/*
        for(int dup=d_up_min;dup<=d_up_max;dup++){

            numup= __builtin_popcount (dup);

            if( numup == Nup ){

                for(int ddn=d_dn_min;ddn<=d_dn_max;ddn++){
                    numdn= __builtin_popcount (ddn);
                    if( numdn == Ndn ){
                        D_up_basis.push_back(dup);
                        D_dn_basis.push_back(ddn);
                    }



                }
            }

        }

        if(D_up_basis.size() ==0){
            Nup_offsets[Nup].second =0;
        }
        else{
            Nup_offsets[Nup].second = D_up_basis.size() - 1;
        }


        */

        cout<<"basis upto now "<<basis_count<<endl;

    }


    cout<<"Size of Restricted Hilbert space size = "<<D_up_basis.size()<<"   "<<basis_count<<endl;



    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_multi_orb_Hubb_chain_GC_Restricted::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

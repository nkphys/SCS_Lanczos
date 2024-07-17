/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_Moire_Kspace.h"
using namespace std;


void BASIS_Moire_Kspace::Construct_basis(){

    ulli d_up_max,d_up_min;
    ulli d_dn_max,d_dn_min;

    Length = Length1*Length2;
    assert(n_orb>0);


    cout<<"Length = "<<Length<<endl;
    cout<<"Nup = "<<Nup<<endl;
    cout<<"Ndn = "<<Ndn<<endl;
    cout<<"K1 = "<<K1_target <<endl;
    cout<<"K2 = "<<K2_target <<endl;

    int K1_up_bare, K2_up_bare;
    int K1_dn_bare, K2_dn_bare;
    int K1_, K2_;
    ulli d_up, d_dn;

    D_up_basis.clear();
    D_dn_basis.clear();

    Mat_1_ullint Dec_vec_up;
    int length_state_str_up;
    Dec_vec_up.clear();
    char state_str_up[n_orb*Length];
    for(int i=0;i<n_orb*Length;i++){
        if(i<Nup){
        state_str_up[i]='1';
        }
        else{
        state_str_up[i]='0';
        }
    }
    state_str_up[Length] = '\0';
    length_state_str_up = strlen(state_str_up);
    findPermutations(state_str_up, 0, length_state_str_up, Dec_vec_up, 2);


    Mat_1_ullint Dec_vec_dn;
    int length_state_str_dn;
    Dec_vec_dn.clear();
    char state_str_dn[n_orb*Length];
    for(int i=0;i<n_orb*Length;i++){
        if(i<Ndn){
        state_str_dn[i]='1';
        }
        else{
        state_str_dn[i]='0';
        }
    }
    state_str_dn[Length] = '\0';
    length_state_str_dn = strlen(state_str_dn);
    findPermutations(state_str_dn, 0, length_state_str_dn, Dec_vec_dn, 2);


    Mat_1_int K1_up_bare_array, K2_up_bare_array;
    Mat_1_int K1_dn_bare_array, K2_dn_bare_array;

    for(int d_up_ind=0;d_up_ind<Dec_vec_up.size();d_up_ind++){
        d_up= Dec_vec_up[d_up_ind];
        K1_up_bare=0; K2_up_bare=0;
        for(int k1_=0;k1_<Length1;k1_++){
        for(int k2_=0;k2_<Length2;k2_++){
        for(int gamma=0;gamma<n_orb;gamma++){
            K1_up_bare += bit_value(d_up,gamma*Length + (k1_+k2_*Length1))*k1_;
            K2_up_bare += bit_value(d_up,gamma*Length + (k1_+k2_*Length1))*k2_;
        }
    }
    }
     K1_up_bare_array.push_back(K1_up_bare);
     K2_up_bare_array.push_back(K2_up_bare);
    }

    for(int d_dn_ind=0;d_dn_ind<Dec_vec_dn.size();d_dn_ind++){
        d_dn= Dec_vec_dn[d_dn_ind];
        K1_dn_bare=0; K2_dn_bare=0;
        for(int k1_=0;k1_<Length1;k1_++){
        for(int k2_=0;k2_<Length2;k2_++){
        for(int gamma=0;gamma<n_orb;gamma++){
            K1_dn_bare += bit_value(d_dn,gamma*Length + (k1_+k2_*Length1))*k1_;
            K2_dn_bare += bit_value(d_dn,gamma*Length + (k1_+k2_*Length1))*k2_;
        }
    }
    }
     K1_dn_bare_array.push_back(K1_dn_bare);
     K2_dn_bare_array.push_back(K2_dn_bare);
    }

    ulli d_up_old;
    bool new_up_set=true;
    for(int d_up_ind=0;d_up_ind<Dec_vec_up.size();d_up_ind++){
        d_up= Dec_vec_up[d_up_ind];
        for(int d_dn_ind=0;d_dn_ind<Dec_vec_dn.size();d_dn_ind++){
            d_dn= Dec_vec_dn[d_dn_ind];
             K1_ = (K1_up_bare_array[d_up_ind] + K1_dn_bare_array[d_dn_ind])%Length1;
             K2_ = (K2_up_bare_array[d_up_ind] + K2_dn_bare_array[d_dn_ind])%Length2;

             if(K1_==K1_target && K2_==K2_target){

             D_up_basis.push_back(d_up);
             D_dn_basis.push_back(d_dn);

             new_up_set = ((d_up !=d_up_old) || (D_up_basis.size()==1) );

             if(new_up_set){
              d_up_old = d_up;
              D_up_basis_range.push_back(d_up);
              D_up_basis_range_min.push_back(D_up_basis.size()-1);
              if(D_up_basis.size()>1){
              D_up_basis_range_max.push_back(D_up_basis.size()-2);
              }

             }

             }

        }
    }


    D_up_basis_range_max.push_back(D_up_basis.size()-1);

    K1_up_bare_array.clear();
    K2_up_bare_array.clear();
    K1_dn_bare_array.clear();
    K2_dn_bare_array.clear();

 cout<<"Basis size : "<<D_up_basis.size()<<endl;
 basis_size = D_up_basis.size();

}


void BASIS_Moire_Kspace::Construct_basis_old(){

    ulli d_up_max,d_up_min;
    ulli d_dn_max,d_dn_min;

    Length = Length1*Length2;

    assert(n_orb>0);

    //Calculating min and max decimal_{up,dn}
    d_up_min=0;
    for (int i=0;i<Nup;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_min=d_up_min+ temp ;
    }

    d_dn_min=0;
    for (int i=0;i<Ndn;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_min=d_dn_min+ temp;
    }

    d_up_max=0;
    for (int i=(n_orb*Length - Nup);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_max=d_up_max+ temp;
    }

    d_dn_max=0;
    for (int i=(n_orb*Length - Ndn);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_max=d_dn_max+ temp;
    }


    //putting correct D_'s in the D arrays
    int num_up, num_dn;
    int K1_up_bare, K2_up_bare;
    int K1_dn_bare, K2_dn_bare;
    int K1_, K2_;
    D_up_basis.clear();
    D_dn_basis.clear();
    for(ulli d_up=d_up_min;d_up<=d_up_max;d_up++){
       for(ulli d_dn=d_dn_min;d_dn<=d_dn_max;d_dn++){
            num_up = __builtin_popcount (d_up);
            num_dn = __builtin_popcount (d_dn);
           if(num_up==Nup && num_dn==Ndn){

               K1_up_bare=0;K1_dn_bare=0;
               K2_up_bare=0;K2_dn_bare=0;
               for(int k1_=0;k1_<Length1;k1_++){
               for(int k2_=0;k2_<Length2;k2_++){

               for(int gamma=0;gamma<n_orb;gamma++){
               K1_up_bare += bit_value(d_up,gamma*Length + (k1_+k2_*Length1))*k1_;
               K2_up_bare += bit_value(d_up,gamma*Length + (k1_+k2_*Length1))*k2_;
               K1_dn_bare += bit_value(d_dn,gamma*Length + (k1_+k2_*Length1))*k1_;
               K2_dn_bare += bit_value(d_dn,gamma*Length + (k1_+k2_*Length1))*k2_;
               }

               }
               }


                K1_ = (K1_up_bare + K1_dn_bare)%Length1;
                K2_ = (K2_up_bare + K2_dn_bare)%Length2;


                if(K1_==K1_target && K2_==K2_target){
                D_up_basis.push_back(d_up);
                D_dn_basis.push_back(d_dn);
                }
           }
       }
    }


    cout<<"Basis size : "<<D_up_basis.size()<<endl;
    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_Moire_Kspace::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

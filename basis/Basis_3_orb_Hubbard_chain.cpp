/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_3_orb_Hubbard_chain.h"
using namespace std;


void BASIS_3_orb_Hubb_chain::Construct_basis(){
    int d_up_max,d_up_min;
    int d_dn_max,d_dn_min;
    n_orb=3;
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
    int num;
    D_up_basis.clear();
    for(int d=d_up_min;d<=d_up_max;d++){
        num= __builtin_popcount (d);
        if(num == Nup){
            D_up_basis.push_back(d);
        }

    }


    D_dn_basis.clear();
    for(int d=d_dn_min;d<=d_dn_max;d++){
        num= __builtin_popcount (d);
        if(num== Ndn){
            D_dn_basis.push_back(d);
        }

    }

    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_3_orb_Hubb_chain::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

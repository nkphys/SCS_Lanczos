/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_3_orb_Hubbard_chain_two_SzSectors.h"
using namespace std;


void BASIS_3_orb_Hubb_chain_two_SzSectors::Construct_basis(){

    Nup_sec1=Nup+1;
    Ndn_sec1=Ndn;

    Nup_sec2=Nup;
    Ndn_sec2=Ndn+1;

    int d_up_max_sec1,d_up_min_sec1;
    int d_dn_max_sec1,d_dn_min_sec1;
    int d_up_max_sec2,d_up_min_sec2;
    int d_dn_max_sec2,d_dn_min_sec2;

    n_orb=3;

    //Calculating min and max decimal_{up,dn} sector-1
    d_up_min_sec1=0;
    for (int i=0;i<Nup_sec1;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_min_sec1=d_up_min_sec1+ temp ;
    }

    d_dn_min_sec1=0;
    for (int i=0;i<Ndn_sec1;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_min_sec1=d_dn_min_sec1+ temp;
    }

    d_up_max_sec1=0;
    for (int i=(n_orb*Length - Nup_sec1);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_max_sec1=d_up_max_sec1+ temp;
    }

    d_dn_max_sec1=0;
    for (int i=(n_orb*Length - Ndn_sec1);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_max_sec1=d_dn_max_sec1+ temp;
    }

    //Calculating min and max decimal_{up,dn} sector-2
    d_up_min_sec2=0;
    for (int i=0;i<Nup_sec2;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_min_sec2=d_up_min_sec2+ temp ;
    }

    d_dn_min_sec2=0;
    for (int i=0;i<Ndn_sec2;i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_min_sec2=d_dn_min_sec2+ temp;
    }

    d_up_max_sec2=0;
    for (int i=(n_orb*Length - Nup_sec2);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_max_sec2=d_up_max_sec2+ temp;
    }

    d_dn_max_sec2=0;
    for (int i=(n_orb*Length - Ndn_sec2);i<(n_orb*Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_max_sec2 =d_dn_max_sec2 + temp;
    }

    //putting correct D_'s in the D arrays
    int num;
    D_up_basis_sec1.clear();
    inverse_Dup_sec1.resize( (d_up_max_sec1-d_up_min_sec1+1) );
    int temp_i=0;

    for(int d=d_up_min_sec1;d<=d_up_max_sec1;d++){
        num= __builtin_popcount (d);
        if(num == Nup_sec1){
            inverse_Dup_sec1[d-d_up_min_sec1]=temp_i;
            D_up_basis_sec1.push_back(d);
            temp_i++;
        }

    }


    D_dn_basis_sec1.clear();
    inverse_Ddn_sec1.resize( (d_dn_max_sec1-d_dn_min_sec1+1) );
    temp_i=0;
    for(int d=d_dn_min_sec1;d<=d_dn_max_sec1;d++){
        num= __builtin_popcount (d);
        if(num== Ndn_sec1){
            inverse_Ddn_sec1[d-d_dn_min_sec1]=temp_i;
            D_dn_basis_sec1.push_back(d);
            temp_i++;
        }

    }



    D_up_basis_sec2.clear();
    inverse_Dup_sec2.resize( (d_up_max_sec2-d_up_min_sec2+1) );
    temp_i=0;

    for(int d=d_up_min_sec2;d<=d_up_max_sec2;d++){
        num= __builtin_popcount (d);
        if(num == Nup_sec2){
            inverse_Dup_sec2[d-d_up_min_sec2]=temp_i;
            D_up_basis_sec2.push_back(d);
            temp_i++;
        }

    }


    D_dn_basis_sec2.clear();
    inverse_Ddn_sec2.resize( (d_dn_max_sec2-d_dn_min_sec2+1) );
    temp_i=0;
    for(int d=d_dn_min_sec2;d<=d_dn_max_sec2;d++){
        num= __builtin_popcount (d);
        if(num== Ndn_sec2){
            inverse_Ddn_sec2[d-d_dn_min_sec2]=temp_i;
            D_dn_basis_sec2.push_back(d);
            temp_i++;
        }

    }


 cout<<"Hilbert space size (sec1) = "<<(D_dn_basis_sec1.size()*D_up_basis_sec1.size())
          <<endl;
 cout<<"Hilbert space size (sec2) = "<<(D_dn_basis_sec2.size()*D_up_basis_sec2.size())
          <<endl;
 cout<<"Hilbert space size = "<<(D_dn_basis_sec1.size()*D_up_basis_sec1.size())
       + (D_dn_basis_sec2.size()*D_up_basis_sec2.size()) <<endl;

 basis_size_sec1=D_dn_basis_sec1.size()*D_up_basis_sec1.size();
 basis_size_sec2=D_dn_basis_sec2.size()*D_up_basis_sec2.size();

    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_3_orb_Hubb_chain_two_SzSectors::clear(){

    D_up_basis_sec1.clear();
    D_dn_basis_sec1.clear();
    D_up_basis_sec2.clear();
    D_dn_basis_sec2.clear();

}

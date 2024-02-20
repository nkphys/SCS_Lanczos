/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_1_orb_Hubbard_chain.h"
using namespace std;


void BASIS_1_orb_Hubb_chain::Construct_basis(){
    int d_up_max,d_up_min;
    int d_dn_max,d_dn_min;

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
    for (int i=(Length - Nup);i<(Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_up_max=d_up_max+ temp;
    }

    d_dn_max=0;
    for (int i=(Length - Ndn);i<(Length);i++){
        int temp = (int)(pow(2,i)+0.5);
        d_dn_max=d_dn_max+ temp;
    }

    DdnMax_ = d_dn_max;
    DdnMin_ = d_dn_min;
    DupMax_ = d_up_max;
    DupMin_ = d_up_min;

    

    //putting correct D_'s in the D arrays
    int num;
    D_up_basis.clear();
    inverse_Dup.resize(d_up_max-d_up_min+1);
    int temp_i=0;
    for(int d=d_up_min;d<=d_up_max;d++){
        num= __builtin_popcount (d);
        if(num == Nup){
            inverse_Dup[d-d_up_min]=temp_i;
            D_up_basis.push_back(d);
            temp_i++;
        }

    }

//cout<<"here"<<endl;

    D_dn_basis.clear();
    inverse_Ddn.resize(d_dn_max-d_dn_min+1);
    temp_i=0;
    for(int d=d_dn_min;d<=d_dn_max;d++){
        num= __builtin_popcount (d);
        if(num== Ndn){
            inverse_Ddn[d-d_dn_min]=temp_i;
            D_dn_basis.push_back(d);
            temp_i++;
        }

    }



        cout<<"Hilbert space size = "<<D_dn_basis.size()*D_up_basis.size()<<endl;

//----------makes finding basis index faster----------------------------------//



    partitions_length_up=int (sqrt (D_up_basis.size()));
    if(partitions_length_up>D_up_basis.size()){
        partitions_length_up=D_up_basis.size();
    }

    partitions_length_dn=int (sqrt (D_dn_basis.size()));
    if(partitions_length_dn>D_dn_basis.size()){
        partitions_length_dn=D_dn_basis.size();
    }



    int temp_size_up = 1 + ((D_up_basis.size())/partitions_length_up );
    int temp_size_dn = 1 + ((D_dn_basis.size())/partitions_length_dn );


    partitions_up.clear();
    Dup_val_at_partitions.clear();

    for (int i=0;i<(temp_size_up-1);i++){
        partitions_up.push_back(i*partitions_length_up);
        Dup_val_at_partitions.push_back(D_up_basis[i*partitions_length_up]);
    }
    partitions_up.push_back(D_up_basis.size()-1);
    Dup_val_at_partitions.push_back(D_up_basis[D_up_basis.size()-1]);



    partitions_dn.clear();
    Ddn_val_at_partitions.clear();

    for (int i=0;i<(temp_size_dn-1);i++){
        partitions_dn.push_back(i*partitions_length_dn);
        Ddn_val_at_partitions.push_back(D_dn_basis[i*partitions_length_dn]);
    }
    partitions_dn.push_back(D_dn_basis.size()-1);
    Ddn_val_at_partitions.push_back(D_dn_basis[D_dn_basis.size()-1]);


//-------------------------------------------------------------------------------//

    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_1_orb_Hubb_chain::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

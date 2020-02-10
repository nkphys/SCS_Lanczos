/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_3_orb_Hubbard_chain_GC.h"
using namespace std;


void BASIS_3_orb_Hubb_chain_GC::Construct_basis(){
    Restricted=false;
    int d_up_max,d_up_min;
    int d_dn_max,d_dn_min;
    n_orb=3;
    //Calculating min and max decimal_{up,dn}



    //putting correct D_'s in the D arrays
    int numup, numdn;
    int Ndn;

    D_up_basis.clear();
    D_dn_basis.clear();
    Nup_offsets.clear();

    Nup_offsets.resize(N_total+1);

    Canonical_partition_up.resize(N_total+1);
    Canonical_partition_dn.resize(N_total+1);


    //cout<<"1"<<endl;
    for(int Nup=(max(N_total - (n_orb*Length),0));Nup<=(min(N_total,n_orb*Length));Nup++){

        //cout<<Nup<<endl;

        Nup_offsets[Nup].first = D_up_basis.size();

        d_up_min=0;
        for (int i=0;i<Nup;i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_min=d_up_min+ temp ;
        }

        d_up_max=0;
        for (int i=(n_orb*Length - Nup);i<(n_orb*Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_max=d_up_max+ temp;
        }

        Ndn = N_total - Nup;
        assert(Ndn<32);
        assert(Nup<32);

        d_dn_min=0;
        for (int i=0;i<Ndn;i++){
            int temp = (int)(pow(2,i)+0.5);
            d_dn_min=d_dn_min+ temp;
        }

        d_dn_max=0;
        for (int i=(n_orb*Length - Ndn);i<(n_orb*Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_dn_max=d_dn_max+ temp;
        }


        Canonical_partition_up[Nup].clear();
        Canonical_partition_dn[Nup].clear();


        for(int dup=d_up_min;dup<=d_up_max;dup++){
            numup= __builtin_popcount (dup);
            if( numup == Nup ){
                Canonical_partition_up[Nup].push_back(dup);
            }
        }



        for(int ddn=d_dn_min;ddn<=d_dn_max;ddn++){
            numdn= __builtin_popcount (ddn);
            if( numdn == Ndn ){
                Canonical_partition_dn[Nup].push_back(ddn);
            }
        }



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

    }


    bool PRINTING_BASIS =false;
    if(PRINTING_BASIS){
        int up_el, dn_el;
        string value_string;
        for(int basis_index=0;basis_index<D_up_basis.size();basis_index++){
            cout<<"basis = "<<basis_index<<" XXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
            for(int orb_no=2;orb_no>=0;orb_no--){
                for(int site=0;site<Length;site++){
                    up_el=bit_value(D_up_basis[basis_index],orb_no*Length + site);
                    dn_el=bit_value(D_dn_basis[basis_index],orb_no*Length + site);
                    if(up_el==1 && dn_el==1){
                        value_string="ud";
                    }
                    else if(up_el==1 && dn_el==0){
                        value_string ="up";
                    }
                    else if(up_el==0 && dn_el==1){
                        value_string="dn";
                    }
                    else{
                        value_string="00";
                    }

                    cout<<value_string<<" ";
                }
                cout<<endl;
            }
            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

        }

    }

    //cout<<"here 2"<<endl;

    //print_decimal_of_binary(2147483647);//2147483647
    //Mat_1_int temp=decimal_to_binary(2147483647);
    //cout<<endl;
}


void BASIS_3_orb_Hubb_chain_GC::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

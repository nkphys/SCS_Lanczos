/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_2_orb_Hubbard_chain_KSector.h"
#include <assert.h>
//#include <bits/stdc++.h>
using namespace std;


/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ] [_______________________  _  ]
    site----->[012....................(L-1)] [012....................(L-1)]
    orbital-->[.........orb - "0"..........] [.........orb - "1"..........]

2)  similarly for "down spin" basis

3)  For total
    m=basis.D_dn_basis.size()*i_up + j_dn;
*/


void BASIS_2_orb_Hubb_chain_KSector::Construct_basis(){

    n_orb=2;
    Write_Basis=false;
    Read_Basis=false;

    if(!Read_Basis){
        Mat_1_int D_up_basis_All;
        Mat_1_int D_dn_basis_All;

        int d_up_basis_temp, d_dn_basis_temp;
        int d_up_basis_temp2, d_dn_basis_temp2;
        int d_up_max,d_up_min;
        int d_dn_max,d_dn_min;


        //Calculating min and max decimal_{up,dn}
        d_up_min=0;
        for (int i=0;i<Nup;i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_min=d_up_min+ temp;
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
        D_up_basis_All.clear();
        for(int d=d_up_min;d<=d_up_max;d++){
            num= __builtin_popcount (d);
            if(num == Nup){
                D_up_basis_All.push_back(d);
            }
        }


        D_dn_basis_All.clear();
        for(int d=d_dn_min;d<=d_dn_max;d++){
            num= __builtin_popcount (d);
            if(num== Ndn){
                D_dn_basis_All.push_back(d);
            }
        }


        cout<<"No of Basis states without targeting Momentum Sector = "<<D_dn_basis_All.size()*D_up_basis_All.size()<<endl;


        //Now Saving only representative states for the groups, for spin=UP sector
        D_Period.clear();
        D_Norm.clear();
        D_up_basis.clear();
        D_dn_basis.clear();
        Dup_Range.clear();

        int Max_Dup;
        Max_Dup = D_up_basis_All[D_up_basis_All.size() - 1];
        Dup_Range.resize(Max_Dup+1); //This size is temporary


        bool Added;
        int state_up_count;


        for(int state_up=0;state_up<D_up_basis_All.size();state_up++){
            state_up_count=0;

            for(int state_dn=0;state_dn<D_dn_basis_All.size();state_dn++){

                //            cout<<"|state_up>("<<D_up_basis_All[state_up]<<") = ";
                //            print_binary_of_decimal(D_up_basis_All[state_up]);
                //            cout<<"|state_dn>("<<D_dn_basis_All[state_dn]<<") = ";
                //            print_binary_of_decimal(D_dn_basis_All[state_dn]);

                Added=false;

                d_up_basis_temp2=D_up_basis_All[state_up];
                d_dn_basis_temp2=D_dn_basis_All[state_dn];

                for(int R=1;R<Length;R++){
                    //Translation on orb-0,spin_up
                    d_up_basis_temp = Act_Translation_assuming_PBC(d_up_basis_temp2,0,Length-1);
                    //Translation on orb-1,spin_up
                    d_up_basis_temp2 = Act_Translation_assuming_PBC(d_up_basis_temp,Length,n_orb*Length-1);

                    //Translation on orb-0,spin_dn
                    d_dn_basis_temp = Act_Translation_assuming_PBC(d_dn_basis_temp2,0,Length-1);
                    //Translation on orb-1,spin_dn
                    d_dn_basis_temp2 = Act_Translation_assuming_PBC(d_dn_basis_temp,Length,n_orb*Length-1);



                    //                cout<<"T^{"<<R<<"}|state_up> ("<<d_up_basis_temp2<<") = ";
                    //                print_binary_of_decimal(d_up_basis_temp2);
                    //                cout<<"T^{"<<R<<"}|state_dn> ("<<d_dn_basis_temp2<<") = ";
                    //                print_binary_of_decimal(d_dn_basis_temp2);

                    if(d_up_basis_temp2 < D_up_basis_All[state_up]){
                        //It means the state "d_up_basis_temp2" has been reported before and
                        //or in other words the representative state of D_up_basis_All[state_up]
                        //is already present
                        Added=true;
                        break;
                    }
                    if ( (d_up_basis_temp2==D_up_basis_All[state_up])
                         &&
                         (d_dn_basis_temp2 < D_dn_basis_All[state_dn])
                         ){
                        //It means the state "d_up_basis_temp2,d_dn_basis_temp2" has been reported before and
                        //or in other words the representative state of D_up_basis_All[state_up]
                        //is already present
                        Added=true;
                        break;
                    }

                    if((d_up_basis_temp2==D_up_basis_All[state_up])
                            &&
                            (d_dn_basis_temp2==D_dn_basis_All[state_dn])
                            ){

                        //Check Compatibility with Momentum targeted
                        if( (Momentum_n%(Length/R)) == 0){
                            D_Period.push_back(R);
                            D_Norm.push_back((Length*Length*1.0)/(R*1.0));
                            D_up_basis.push_back(D_up_basis_All[state_up]);
                            D_dn_basis.push_back(D_dn_basis_All[state_dn]);

                            if(state_up_count==0){
                                Dup_Range[D_up_basis_All[state_up]].first = D_up_basis.size() - 1;
                            }

                            state_up_count++;
                            Added=true;
                            break;
                        }
                        else{
                            Added=true;
                            break;
                        }
                    }
                }

                if(Added==false){
                    D_Period.push_back(Length);
                    D_Norm.push_back((Length*Length*1.0)/(Length*1.0));
                    D_up_basis.push_back(D_up_basis_All[state_up]);
                    D_dn_basis.push_back(D_dn_basis_All[state_dn]);

                    if(state_up_count==0){
                        Dup_Range[D_up_basis_All[state_up]].first = D_up_basis.size() - 1;
                    }
                    state_up_count++;
                }
            }

            if(state_up_count==0){
                Dup_Range[D_up_basis_All[state_up]].first = -1;
                Dup_Range[D_up_basis_All[state_up]].second = -1;
            }
            if(state_up_count!=0){
                Dup_Range[D_up_basis_All[state_up]].second = D_up_basis.size() - 1;
            }

            if((state_up%500)==0){
                cout<<"state_up = "<<state_up<<" done"<<endl;
            }
        }


        cout<<"No of Basis states in Momentum Sector ="<<Momentum_n<<" is "<<D_Period.size()<<endl;

        Max_Dup = D_up_basis[D_up_basis.size() -1];
        Dup_Range.resize(Max_Dup+1);


    }
    else{
        string first_line;
        string line;
        int index_n, temp_int;
        double temp_double;
        pair_int temp_pair_int;
        assert(Read_Basis);
        int check_Length;
        int check_N_up;
        int check_N_dn;
        int check_n_orb;
        int check_Momentum_n;
        string file_in_D_up_down_Period_Norm_basis = "IN_D_up_down_Period_Norm_basis.txt";
        ifstream infile_file_in_D_up_down_Period_Norm_basis(file_in_D_up_down_Period_Norm_basis.c_str());

        string file_in_Dup_Range = "IN_Dup_Range.txt";
        ifstream infile_file_in_Dup_Range(file_in_Dup_Range.c_str());


        getline(infile_file_in_D_up_down_Period_Norm_basis,first_line);
        stringstream first_line_ss;
        first_line_ss<<first_line;
        first_line_ss>>check_Length>>check_N_up>>check_N_dn>>check_n_orb>>check_Momentum_n;
        assert(check_Length==Length);
        assert(check_N_up==Nup);
        assert(check_N_dn==Ndn);
        assert(check_n_orb==n_orb);
        assert(check_Momentum_n==Momentum_n);

        getline(infile_file_in_Dup_Range,first_line);
        stringstream first_line_ss_2;
        first_line_ss_2<<first_line;
        first_line_ss_2>>check_Length>>check_N_up>>check_N_dn>>check_n_orb>>check_Momentum_n;
        assert(check_Length==Length);
        assert(check_N_up==Nup);
        assert(check_N_dn==Ndn);
        assert(check_n_orb==n_orb);
        assert(check_Momentum_n==Momentum_n);


        cout<<"Reading Basis"<<endl;
        D_up_basis.clear();
        D_dn_basis.clear();
        D_Period.clear();
        D_Norm.clear();
        while (getline(infile_file_in_D_up_down_Period_Norm_basis, line))
        {
            stringstream line_ss(line);
            line_ss>>temp_int;D_up_basis.push_back(temp_int);
            line_ss>>temp_int;D_dn_basis.push_back(temp_int);
            line_ss>>temp_int;D_Period.push_back(temp_int);
            line_ss>>temp_double;D_Norm.push_back(temp_double);
        }
        cout<<"Basis read"<<endl;

        Dup_Range.clear();
        while (getline(infile_file_in_Dup_Range, line))
        {
            stringstream line_ss(line);
            line_ss>>index_n;
            line_ss>>temp_int;temp_pair_int.first=temp_int;
            line_ss>>temp_int;temp_pair_int.second=temp_int;
            Dup_Range.push_back(temp_pair_int);

        }


    }


    /*

     Mat_1_int checks;
     checks.clear();
     checks.push_back(928);checks.push_back(876);checks.push_back(822);checks.push_back(770);
     checks.push_back(636);checks.push_back(567);checks.push_back(498);checks.push_back(429);
     for(int i=0;i<checks.size();i++){
         cout<<"-----------------------------"<<endl;
         cout << "FOR Nup("<<checks[i]<<"):----------------"<<endl;
         print_binary_of_decimal(D_up_basis[checks[i]]);
         cout << "FOR Ndn("<<checks[i]<<"):----------------"<<endl;
         print_binary_of_decimal(D_dn_basis[checks[i]]);
         cout<<"-----------------------------"<<endl;
         cout<<"-----------------------------"<<endl;
         cout<<"-----------------------------"<<endl;
     }

    cout <<"END"<<endl;
    */


    if(Write_Basis){

        string file_out_D_up_down_Period_Norm_basis = "OUT_D_up_down_Period_Norm_basis.txt";
        ofstream outfile_file_out_D_up_down_Period_Norm_basis(file_out_D_up_down_Period_Norm_basis.c_str());

        string file_out_Dup_Range = "OUT_Dup_Range.txt";
        ofstream outfile_file_out_Dup_Range(file_out_Dup_Range.c_str());


        outfile_file_out_D_up_down_Period_Norm_basis<<Length<<"  "<<Nup<<"  "<<Ndn<<"  "<<n_orb<<"  "<<Momentum_n<<endl;
        outfile_file_out_Dup_Range<<Length<<"  "<<Nup<<"  "<<Ndn<<"  "<<n_orb<<"  "<<Momentum_n<<endl;

        for(int n=0;n<D_up_basis.size();n++){
            outfile_file_out_D_up_down_Period_Norm_basis<<D_up_basis[n]<<"   "<<
                                                          D_dn_basis[n]<<"   "<<
                                                          D_Period[n]<<"   "<<
                                                          D_Norm[n]<<"   "<<endl;
        }

        for(int n=0;n<Dup_Range.size();n++){
            outfile_file_out_Dup_Range<<n<<"    "<<Dup_Range[n].first<<"    "<<Dup_Range[n].second<<endl;
        }

    }


}


void BASIS_2_orb_Hubb_chain_KSector::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

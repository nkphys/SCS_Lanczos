/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_1_orb_Hubbard_2D_KSector.h"
#include <assert.h>
#define PI 3.14159265358979323846
//#include <bits/stdc++.h>
using namespace std;


/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ]
    site----->[012....................(L-1)]

2)  similarly for "down spin" basis

3)  For total
    m=basis.D_dn_basis.size()*i_up + j_dn;
*/


void BASIS_1_orb_Hubb_2D_KSector::Construct_basis(){


    double epsilon=0.000001;
    int n_int;
    double n_double, n_new_double;

    if(!Read_Basis){
        Mat_1_int D_up_basis_All;
        Mat_1_int D_dn_basis_All;

        int d_up_basis_temp, d_dn_basis_temp;
        int d_up_basis_temp_Xtrans, d_dn_basis_temp_Xtrans;
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
        for (int i=(Length - Nup);i<(Length);i++){
            int temp = (int)(pow(2,i)+0.5);
            d_up_max=d_up_max+ temp;
        }

        d_dn_max=0;
        for (int i=(Length - Ndn);i<(Length);i++){
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
        Dx_Period.clear();
        Dy_Period.clear();
        D_Norm.clear();
        D_up_basis.clear();
        D_dn_basis.clear();
        Dup_Range.clear();

        int Max_Dup;
        Max_Dup = D_up_basis_All[D_up_basis_All.size() - 1];
        Dup_Range.resize(Max_Dup+1); //This size is temporary

        int Ry_min, Ry_max;
        int Rx_min, Rx_max;
        int m_bar;
        int Ry_, Rx_;
        int sign_pow_min_R, sign_up_min_R, sign_dn_min_R;
        double sign_min_R;
        int sign_up_min_R_orb0, sign_dn_min_R_orb0;
        int sign_up_min_R_Xtrans, sign_dn_min_R_Xtrans;
        complex<double> gamma_;

        bool Already_Added;
        bool Addition_Allowed;
        bool new_and_allowed_Tx, new_and_allowed_Ty;
        int state_up_count;


        for(int state_up=0;state_up<D_up_basis_All.size();state_up++){
            state_up_count=0;

            for(int state_dn=0;state_dn<D_dn_basis_All.size();state_dn++){

                //            cout<<"|state_up>("<<D_up_basis_All[state_up]<<") = ";
                //            print_binary_of_decimal(D_up_basis_All[state_up]);
                //            cout<<"|state_dn>("<<D_dn_basis_All[state_dn]<<") = ";
                //            print_binary_of_decimal(D_dn_basis_All[state_dn]);

                d_up_basis_temp=D_up_basis_All[state_up];
                d_dn_basis_temp=D_dn_basis_All[state_dn];

                Already_Added=false;
                Addition_Allowed=false;
                sign_up_min_R=0;
                sign_dn_min_R=0;
                for(int Rx=1;Rx<=Lx;Rx++){

                    if(Lx>1){
                        for(int iy=0;iy<Ly;iy++){

                            sign_dn_min_R_orb0 = one_bits_in_bw(iy*Lx, (iy+1)*Lx - 1, d_dn_basis_temp) +
                                    1*bit_value(d_dn_basis_temp,iy*Lx);
                            if(bit_value(d_dn_basis_temp, (iy+1)*Lx - 1 )==1){
                                sign_dn_min_R += 1*sign_dn_min_R_orb0;
                            }

                            sign_up_min_R_orb0 = one_bits_in_bw(iy*Lx, (iy+1)*Lx - 1, d_up_basis_temp) +
                                    1*bit_value(d_up_basis_temp,iy*Lx);
                            if(bit_value(d_up_basis_temp, (iy+1)*Lx - 1 )==1){
                                sign_up_min_R += 1*sign_up_min_R_orb0;
                            }


                            //Translation on spin_up
                            d_up_basis_temp = Act_Translation_2D_alongX_assuming_PBC(d_up_basis_temp, Lx ,Ly, iy);
                            //Translation on spin_dn
                            d_dn_basis_temp = Act_Translation_2D_alongX_assuming_PBC(d_dn_basis_temp, Lx ,Ly, iy);
                        }
                    }


                    sign_up_min_R_Xtrans = sign_up_min_R;
                    sign_dn_min_R_Xtrans = sign_dn_min_R;

                    d_up_basis_temp_Xtrans=d_up_basis_temp;
                    d_dn_basis_temp_Xtrans=d_dn_basis_temp;

                    for(int Ry=1;Ry<=Ly;Ry++){

                        if(Ly>1){
                            for(int ix=0;ix<Lx;ix++){

                                for(int iy_=0;iy_<Ly-1;iy_++){
                                    sign_dn_min_R_orb0 = one_bits_in_bw(ix + iy_*Lx, ix + (iy_+1)*Lx, d_dn_basis_temp);
                                    sign_dn_min_R_orb0 = sign_dn_min_R_orb0*bit_value(d_dn_basis_temp, ix + iy_*Lx);

                                    sign_dn_min_R += sign_dn_min_R_orb0;
                                }
                                sign_dn_min_R_orb0 = one_bits_in_bw(ix, ix + (Ly-1)*Lx, d_dn_basis_temp) +
                                        bit_value(d_dn_basis_temp, ix) ;
                                sign_dn_min_R_orb0 = sign_dn_min_R_orb0*bit_value(d_dn_basis_temp, ix + (Ly-1)*Lx);
                                sign_dn_min_R += sign_dn_min_R_orb0;

                                for(int iy_=0;iy_<Ly-1;iy_++){
                                    sign_up_min_R_orb0 = one_bits_in_bw(ix + iy_*Lx, ix + (iy_+1)*Lx, d_up_basis_temp);
                                    sign_up_min_R_orb0 = sign_up_min_R_orb0*bit_value(d_up_basis_temp, ix + iy_*Lx);

                                    sign_up_min_R += sign_up_min_R_orb0;
                                }
                                sign_up_min_R_orb0 = one_bits_in_bw(ix, ix + (Ly-1)*Lx, d_up_basis_temp) +
                                        bit_value(d_up_basis_temp, ix) ;
                                sign_up_min_R_orb0 = sign_up_min_R_orb0*bit_value(d_up_basis_temp, ix + (Ly-1)*Lx);
                                sign_up_min_R += sign_up_min_R_orb0;


                                //Translation on spin_up
                                d_up_basis_temp = Act_Translation_2D_alongY_assuming_PBC(d_up_basis_temp, Lx ,Ly, ix);
                                //Translation on spin_dn
                                d_dn_basis_temp = Act_Translation_2D_alongY_assuming_PBC(d_dn_basis_temp, Lx ,Ly, ix);
                            }
                        }


                        if(d_up_basis_temp < D_up_basis_All[state_up]){
                            //It means the state "d_up_basis_temp" has been reported before and
                            //or in other words the representative state of D_up_basis_All[state_up]
                            //is already present
                            Already_Added=true;
                            break;
                        }
                        if ( (d_up_basis_temp==D_up_basis_All[state_up])
                             &&
                             (d_dn_basis_temp < D_dn_basis_All[state_dn])
                             ){
                            //It means the state "d_up_basis_temp,d_dn_basis_temp" has been reported before and
                            //or in other words the representative state of D_up_basis_All[state_up]
                            //is already present
                            Already_Added=true;
                            break;
                        }

                        if(   (  (d_up_basis_temp==D_up_basis_All[state_up])
                                &&
                                (d_dn_basis_temp==D_dn_basis_All[state_dn])  )
                              &&
                              !Addition_Allowed

                                ){

                            Addition_Allowed =true;
                            Rx_=Rx;
                            Ry_=Ry;
                            sign_pow_min_R = sign_up_min_R + sign_dn_min_R;

                        }
                    }
                    d_up_basis_temp=d_up_basis_temp_Xtrans;
                    d_dn_basis_temp=d_dn_basis_temp_Xtrans;
                    sign_up_min_R = sign_up_min_R_Xtrans;
                    sign_dn_min_R = sign_dn_min_R_Xtrans;

                }


                if( (!Already_Added) && Addition_Allowed ){

                    for(int m_=1;m_<=max(Lx,Ly);m_++){
                        if( ((m_*Rx_)%Lx)==0  &&   ((m_*Ry_)%Ly)==0    ){
                            m_bar=m_;
                            break;
                        }
                    }

                    sign_min_R = pow(-1.0, 1.0*sign_pow_min_R);

                    gamma_=zero_comp;
                    for(int m=0;m<m_bar;m++){
                        gamma_ += exp(-1.0*iota_comp*( ((2.0*PI*Momentum_nx*Rx_*m)/(Lx*1.0)) +
                                                       ((2.0*PI*Momentum_ny*Ry_*m)/(Ly*1.0))
                                                       ))* ( pow(1.0*sign_min_R, 1.0*m) );
                    }

                    if(abs(gamma_)>=epsilon)
                    {
                        Dx_Period.push_back(Rx_);
                        Dy_Period.push_back(Ry_);
                        D_Norm.push_back((Lx*Ly*abs(gamma_)*abs(gamma_)*1.0)/(m_bar*1.0));
                        D_up_basis.push_back(D_up_basis_All[state_up]);
                        D_dn_basis.push_back(D_dn_basis_All[state_dn]);

                        if(state_up_count==0){
                            Dup_Range[D_up_basis_All[state_up]].first = D_up_basis.size() - 1;
                        }
                        state_up_count++;
                    }

                }

            }

            if(state_up_count==0){
                Dup_Range[D_up_basis_All[state_up]].first = -1;
                Dup_Range[D_up_basis_All[state_up]].second = -1;
            }
            if(state_up_count!=0){
                Dup_Range[D_up_basis_All[state_up]].second = D_up_basis.size() - 1;
            }

            if((state_up%200)==0){
                cout<<"state_up = "<<state_up<<" done"<<endl;
            }
        }


        cout<<"No of Basis states in Momentum Sector = ("<<Momentum_nx<<", "<< Momentum_ny<<") is "<<Dx_Period.size()<<endl;

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
        int check_Lx;
        int check_Ly;
        int check_N_up;
        int check_N_dn;
        int check_Momentum_nx;
        int check_Momentum_ny;
        string file_in_D_up_down_Period_Norm_basis =  file_read_basis + "_D_up_down_Period_Norm.txt";
        ifstream infile_file_in_D_up_down_Period_Norm_basis(file_in_D_up_down_Period_Norm_basis.c_str());

        string file_in_Dup_Range = file_read_basis + "_Dup_Range.txt";
        ifstream infile_file_in_Dup_Range(file_in_Dup_Range.c_str());


        getline(infile_file_in_D_up_down_Period_Norm_basis,first_line);
        stringstream first_line_ss;
        first_line_ss<<first_line;
        first_line_ss>>check_Lx>>check_Ly>>check_N_up>>check_N_dn>>check_Momentum_nx>>check_Momentum_ny;
        assert(check_Lx==Lx);
        assert(check_Ly==Ly);
        assert(check_N_up==Nup);
        assert(check_N_dn==Ndn);
        assert(check_Momentum_nx==Momentum_nx);
        assert(check_Momentum_ny==Momentum_ny);

        getline(infile_file_in_Dup_Range,first_line);
        stringstream first_line_ss_2;
        first_line_ss_2<<first_line;
        first_line_ss_2>>check_Lx>>check_Ly>>check_N_up>>check_N_dn>>check_Momentum_nx>>check_Momentum_ny;
        assert(check_Lx==Lx);
        assert(check_Ly==Ly);
        assert(check_N_up==Nup);
        assert(check_N_dn==Ndn);
        assert(check_Momentum_nx==Momentum_nx);
        assert(check_Momentum_ny==Momentum_ny);


        cout<<"Reading Basis"<<endl;
        D_up_basis.clear();
        D_dn_basis.clear();
        Dx_Period.clear();
        Dy_Period.clear();
        D_Norm.clear();
        while (getline(infile_file_in_D_up_down_Period_Norm_basis, line))
        {
            stringstream line_ss(line);
            line_ss>>temp_int;D_up_basis.push_back(temp_int);
            line_ss>>temp_int;D_dn_basis.push_back(temp_int);
            line_ss>>temp_int;Dx_Period.push_back(temp_int);
            line_ss>>temp_int;Dy_Period.push_back(temp_int);
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


        cout<<"No of Basis states in Momentum Sector = ("<<Momentum_nx<<", "<< Momentum_ny<<") is "<<Dx_Period.size()<<endl;

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

        string file_out_D_up_down_Period_Norm_basis = file_write_basis + "_D_up_down_Period_Norm_basis.txt";
        ofstream outfile_file_out_D_up_down_Period_Norm_basis(file_out_D_up_down_Period_Norm_basis.c_str());

        string file_out_Dup_Range = file_write_basis + "_Dup_Range.txt";
        ofstream outfile_file_out_Dup_Range(file_out_Dup_Range.c_str());


        outfile_file_out_D_up_down_Period_Norm_basis<<Lx<<"  "<<Ly<<"  "<<Nup<<"  "<<Ndn<<"  "<<Momentum_nx<<"  "<<Momentum_ny<<endl;
        outfile_file_out_Dup_Range<<Lx<<"  "<<Ly<<"  "<<Nup<<"  "<<Ndn<<"  "<<Momentum_nx<<"  "<<Momentum_ny<<endl;

        for(int n=0;n<D_up_basis.size();n++){
            outfile_file_out_D_up_down_Period_Norm_basis<<D_up_basis[n]<<"   "<<
                                                          D_dn_basis[n]<<"   "<<
                                                          Dx_Period[n]<<"   "<<
                                                          Dy_Period[n]<<"   "<<
                                                          D_Norm[n]<<"   "<<endl;
        }

        for(int n=0;n<Dup_Range.size();n++){
            outfile_file_out_Dup_Range<<n<<"    "<<Dup_Range[n].first<<"    "<<Dup_Range[n].second<<endl;
        }

    }


}


void BASIS_1_orb_Hubb_2D_KSector::clear(){
    D_up_basis.clear();
    D_dn_basis.clear();
}

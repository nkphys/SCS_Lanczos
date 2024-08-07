/*
This class includes the Model for which Lanczos is being done
*/
//#ifndef USE_COMPLEX
#include "Model_1_orb_Hubbard_2D_KSector.h"
#include <stdlib.h>
#include <string>
using namespace std;
#define PI 3.14159265358979323846
#ifdef _OPENMP
#include <omp.h>
#endif
/*convention for basis:

1)  for "up-spin" basis
              [_______________________  _  ]
    site----->[012....................(L-1)]

2)  similarly for "down spin" basis

3)  For total
    m=basis.D_dn_basis.size()*i + j;
*/

void MODEL_1_orb_Hubb_2D_KSector::Act_Hamil(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){

    Vec_out.clear();
    Vec_out.resize(Vec_in.size());

    Act_diagonal_terms(basis, Vec_in, Vec_out);
    Act_connections(basis, Vec_in, Vec_out);

}


void MODEL_1_orb_Hubb_2D_KSector::Act_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());

    int no_of_proc;
    no_of_proc=1;

#ifdef _OPENMP
    int temp_int;
    temp_int = basis.D_up_basis.size();
    no_of_proc=min(temp_int, NProcessors_);
    omp_set_num_threads(no_of_proc);
    cout<<"Diagonal terms acting: "<<no_of_proc<<" processors"<<endl;
#endif



#ifdef _OPENMP
#pragma omp parallel
    {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (int i=0;i<basis.D_up_basis.size();i++){
            //Remember H[l][m]=<l|H|m>
            int m,j;
            double value;
            m=i;
            j=i;

            value=0;
            //on-site coulomb repulsion  ni_up ni_dn:
            value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);


            //longrange coulomb interaction ni.nj
            for(int site_i=0;site_i<basis.Length;site_i++){
                for(int site_j=0;site_j<basis.Length;site_j++){
                    value += NonLocalInteractions_mat[site_i][site_j]*( bit_value(basis.D_up_basis[i], site_i) +
                                                                        bit_value(basis.D_dn_basis[j], site_i) )*
                            ( bit_value(basis.D_up_basis[i], site_j) +
                              bit_value(basis.D_dn_basis[j], site_j) );
                }
            }

            //magnetic Field
            for(int site=0;site<basis.Length;site++){
                value+=0.5*(H_field[site])*
                        ( ( bit_value(basis.D_up_basis[i], site) -
                            bit_value(basis.D_dn_basis[j], site) )
                          );
            }


            if(value!=0){
                Vec_out[m] += Vec_in[m]*value*one;
            }


        }

#ifdef _OPENMP
    }
#endif


    cout<<"Done Hamiltonian construction: Diagonal"<<endl;

}



void MODEL_1_orb_Hubb_2D_KSector::Act_connections(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out){


    int no_of_proc;
    no_of_proc=1;

#ifdef _OPENMP
    double begin_time_connections, end_time_connections;
    begin_time_connections = omp_get_wtime();
#endif



#ifdef _OPENMP
    int temp_int;
    temp_int = basis.D_up_basis.size();
    no_of_proc=min(temp_int, NProcessors_);
    omp_set_num_threads(no_of_proc);
    cout<<"Connections acting: "<<no_of_proc<<" processors"<<endl;
    //    Vec_out_temp.resize(no_of_proc);
    //    for(int i=0;i<no_of_proc;i++){
    //        Vec_out_temp[i].resize(Vec_in.size());
    //    }
#endif



#ifdef _OPENMP
#pragma omp parallel
    {
#endif


#ifdef _OPENMP
#pragma omp for nowait
#endif

        for (int i=0;i<basis.D_up_basis.size();i++){

            double value;
            int m,j;
            int D_up,D_dn;
            int i_new;
            int m_new;
            double sign_FM;
            int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
            int sign_pow_up, sign_pow_dn;
            int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
            int l,lp;
            int range_min, range_max;
            bool row_found_;
            double_type phase_;
            int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
            complex<double> iota_(0.0,1.0);
            int D_up_temp ,D_dn_temp;
            int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
            bool repeating_rows;
            int row_counter;
            int check_min, check_max;
            int site, site_p;
            int Lxm1, Lym1;
            Lxm1=basis.Lx-1;
            Lym1=basis.Ly-1;



            /******TO REMOVE*********************
        cout<<"up("<<i<<") : ";
        print_binary_of_decimal(basis.D_up_basis[i]);
        cout<<"dn("<<i<<") : ";
        print_binary_of_decimal(basis.D_dn_basis[i]);
        cout<<endl;
        ************************************/

            m=i;
            j=i;

            value=0;

            row_counter=0;
            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);

                    for(int ix_p=0;ix_p<basis.Lx ;ix_p++){
                        for(int iy_p=0;iy_p<basis.Ly ;iy_p++){
                            site_p=ix_p + (iy_p*basis.Lx);

                            if(abs(Hopping_mat_NN[site_p][site])>0.00000001)// && Dis_y==0)
                            { // LongRange
                                Hopping_NN=Hopping_mat_NN[site_p][site];

                                //---------------Hopping for up electrons-------------------//
                                //there have to be one up electron in site
                                //there have to be no up electron in site_p
                                if(
                                        (bit_value(basis.D_up_basis[i],site)==1)
                                        &&
                                        (bit_value(basis.D_up_basis[i],site_p)==0)

                                        )
                                {

                                    sign_pow_up=0;
                                    sign_pow_dn=0;

                                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                                  - pow(2,site) );
                                    D_dn = basis.D_dn_basis[m] ;

                                    D_up_temp=D_up;
                                    D_dn_temp=D_dn;
                                    row_found_=false;

                                    for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                        if(inv_trnsltns_x>0 && basis.Lx>1){

                                            for(int iy_=0;iy_<basis.Ly;iy_++){

                                                //Inv Translation on spin_dn
                                                sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                        1*bit_value(D_dn_temp,iy_*basis.Lx);
                                                if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                    sign_pow_dn += 1*sign_pow_dn_orb0;
                                                }

                                                D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                                //Inv Translation on spin_up
                                                sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                        1*bit_value(D_up_temp,iy_*basis.Lx);
                                                if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                    sign_pow_up += 1*sign_pow_up_orb0;
                                                }

                                                D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                            }
                                        }
                                        D_dn_temp_Xtrans = D_dn_temp;
                                        D_up_temp_Xtrans = D_up_temp;
                                        sign_pow_up_Xtrans = sign_pow_up;
                                        sign_pow_dn_Xtrans = sign_pow_dn;
                                        for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                            if(inv_trnsltns_y>0 && basis.Ly>1){

                                                for(int ix_=0;ix_<basis.Lx;ix_++){

                                                    //Inv Translation on spin_dn
                                                    for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                        sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                        sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                        sign_pow_dn += sign_pow_dn_orb0;
                                                    }
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                            bit_value(D_dn_temp, ix_) ;
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                    sign_pow_dn += sign_pow_dn_orb0;

                                                    D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                    //Inv Translation on spin_up
                                                    for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                        sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                        sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                        sign_pow_up += sign_pow_up_orb0;
                                                    }
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                            bit_value(D_up_temp, ix_) ;
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                    sign_pow_up += sign_pow_up_orb0;

                                                    D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                                }
                                            }


                                            if(D_up_temp>=basis.Dup_Range.size())
                                            {
                                                row_found_=false;
                                            }
                                            else
                                            {
                                                assert(D_up_temp<basis.Dup_Range.size());
                                                range_min=basis.Dup_Range[D_up_temp].first;
                                                range_max=basis.Dup_Range[D_up_temp].second;
                                                if(range_min==-1)
                                                {
                                                    row_found_=false;
                                                    assert(range_max==-1);
                                                }
                                                else
                                                {
                                                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                    if(i_new==-1){
                                                        row_found_=false;
                                                    }
                                                    else{
                                                        row_found_=true;
                                                        Inv_Trnsltns_x_=inv_trnsltns_x;
                                                        Inv_Trnsltns_y_=inv_trnsltns_y;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        if(row_found_){break;}

                                        D_up_temp=D_up_temp_Xtrans ;
                                        D_dn_temp=D_dn_temp_Xtrans ;
                                        sign_pow_up = sign_pow_up_Xtrans;
                                        sign_pow_dn = sign_pow_dn_Xtrans;
                                    }


                                    if(row_found_==true){
                                        m_new = i_new;

#ifdef USE_COMPLEX
                                        phase_=exp(-1.0*iota_*( ((2.0*PI*(Inv_Trnsltns_x_)*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*(Inv_Trnsltns_y_)*basis.Momentum_ny)/(basis.Ly))   )
                                                   )
                                                *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                        if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                            cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                        }
                                        assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                        phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                        sign_pow_up += one_bits_in_bw(site,site_p,basis.D_up_basis[i]);

                                        sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);

                                        //                                    if(m_new<=m)
                                        //                                    {
                                        //                                        repeating_rows=false;
                                        //                                        check_min=Hamil.rows.size()-1;
                                        //                                        check_max=(Hamil.rows.size()-1)-row_counter;
                                        //                                        // cout<<check_min<<endl;
                                        //                                        // cout<<check_max<<endl;
                                        //                                        for(int check_=check_min;check_>check_max;check_--){
                                        //                                            if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                        //                                                Hamil.value[check_] +=-1.0*sign_FM*(Hopping_NN)*one*phase_;
                                        //                                                repeating_rows=true;
                                        //                                                break;
                                        //                                            }
                                        //                                        }

                                        //                                        if(!repeating_rows){
                                        //                                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_NN)*one*phase_);
                                        //                                            Hamil.rows.push_back(m_new);
                                        //                                            Hamil.columns.push_back(m);
                                        //                                            row_counter++;
                                        //                                        }


                                        //                                      Vec_out[m_new] +=  -1.0*sign_FM*(Hopping_NN)*one*phase_*Vec_in[m];
                                        Vec_out[m] +=  -1.0*sign_FM*(conjugate(Hopping_NN))*one*phase_*Vec_in[m_new];

                                        //                                    }
                                    }

                                } // if up hopping possible


                                //---------------Hopping for dn electrons-------------------//
                                //there have to be one dn electron in gamma, site
                                //there have to be no dn electron in gamma_p, site_p
                                if(
                                        (bit_value(basis.D_dn_basis[j], site)==1)
                                        &&
                                        (bit_value(basis.D_dn_basis[j], site_p)==0)
                                        )
                                {

                                    sign_pow_dn=0;
                                    sign_pow_up=0;

                                    D_up = basis.D_up_basis[m];
                                    D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                                  - pow(2, site) );

                                    D_up_temp=D_up;
                                    D_dn_temp=D_dn;

                                    row_found_=false;

                                    for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                        if(inv_trnsltns_x>0 && basis.Lx>1){

                                            for(int iy_=0;iy_<basis.Ly;iy_++){

                                                //Translation on spin_dn
                                                sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                        1*bit_value(D_dn_temp,iy_*basis.Lx);
                                                if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                    sign_pow_dn += 1*sign_pow_dn_orb0;
                                                }

                                                D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                                //Translation on spin_up
                                                sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                        1*bit_value(D_up_temp,iy_*basis.Lx);
                                                if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                    sign_pow_up += 1*sign_pow_up_orb0;
                                                }

                                                D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                            }
                                        }
                                        D_dn_temp_Xtrans = D_dn_temp;
                                        D_up_temp_Xtrans = D_up_temp;
                                        sign_pow_up_Xtrans = sign_pow_up;
                                        sign_pow_dn_Xtrans = sign_pow_dn;
                                        for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){

                                            if(inv_trnsltns_y>0 && basis.Ly>1){

                                                for(int ix_=0;ix_<basis.Lx;ix_++){

                                                    //Translation on spin_dn
                                                    for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                        sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                        sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                        sign_pow_dn += sign_pow_dn_orb0;
                                                    }
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                            bit_value(D_dn_temp, ix_) ;
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                    sign_pow_dn += sign_pow_dn_orb0;

                                                    D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                    //Translation on spin_up
                                                    for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                        sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                        sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                        sign_pow_up += sign_pow_up_orb0;
                                                    }
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                            bit_value(D_up_temp, ix_) ;
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                    sign_pow_up += sign_pow_up_orb0;

                                                    D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                                }
                                            }




                                            if(D_up_temp>=basis.Dup_Range.size())
                                            {
                                                row_found_=false;
                                            }
                                            else
                                            {
                                                assert(D_up_temp<basis.Dup_Range.size());
                                                range_min=basis.Dup_Range[D_up_temp].first;
                                                range_max=basis.Dup_Range[D_up_temp].second;
                                                if(range_min==-1)
                                                {
                                                    row_found_=false;
                                                    assert(range_max==-1);
                                                }
                                                else
                                                {
                                                    i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                    if(i_new==-1){
                                                        row_found_=false;
                                                    }
                                                    else{
                                                        row_found_=true;
                                                        Inv_Trnsltns_x_=inv_trnsltns_x;
                                                        Inv_Trnsltns_y_=inv_trnsltns_y;
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                        if(row_found_){break;}
                                        D_up_temp=D_up_temp_Xtrans ;
                                        D_dn_temp=D_dn_temp_Xtrans ;
                                        sign_pow_up = sign_pow_up_Xtrans;
                                        sign_pow_dn = sign_pow_dn_Xtrans;
                                    }


                                    if(row_found_==true){
                                        m_new = i_new;

#ifdef USE_COMPLEX
                                        phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                                   )
                                                *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));

#endif
#ifndef USE_COMPLEX
                                        if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                            cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                        }
                                        assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                        phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));

#endif

                                        sign_pow_dn += one_bits_in_bw(site,site_p,basis.D_dn_basis[j]);

                                        sign_FM = pow(-1.0, sign_pow_dn+sign_pow_up);


                                        //                                    if(m_new<=m)
                                        //                                    {
                                        //                                        repeating_rows=false;
                                        //                                        check_min=Hamil.rows.size()-1;
                                        //                                        check_max=(Hamil.rows.size()-1)-row_counter;
                                        //                                        for(int check_=check_min;check_>check_max;check_--){
                                        //                                            if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                        //                                                Hamil.value[check_] +=-1.0*sign_FM*(Hopping_NN)*one*phase_;
                                        //                                                repeating_rows=true;
                                        //                                                break;
                                        //                                            }
                                        //                                        }

                                        //                                        if(!repeating_rows){
                                        //                                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_NN)*one*phase_);
                                        //                                            Hamil.rows.push_back(m_new);
                                        //                                            Hamil.columns.push_back(m);
                                        //                                            row_counter++;
                                        //                                        }

                                        //                                        Vec_out[m_new] += -1.0*sign_FM*(Hopping_NN)*one*phase_*Vec_in[m];
                                        Vec_out[m] += -1.0*sign_FM*(conjugate(Hopping_NN))*one*phase_*Vec_in[m_new];

                                        //                                    }

                                    }
                                } // if dn hopping possible


                            }//nearest neighbour

                        }  //iy_p
                    } //ix_p

                } //iy
            } // ix

            if(m%10 ==1){
                //cout<<"Connection: done "<<m<<" basis"<<endl;
            }

        } // "i" i.e up_decimals


#ifdef _OPENMP
    }
#endif

    cout<<"Done Hamiltonian construction: Connections"<<endl;

#ifdef _OPENMP
    end_time_connections = omp_get_wtime();
    cout<<"Time for adding connections to Hamil [using OpenMP] = "<<double(end_time_connections - begin_time_connections)<<endl;
#endif

}




void MODEL_1_orb_Hubb_2D_KSector::Add_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis){


    cout<<"Started Hamiltonian construction: Diagonal"<<endl;
    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());

    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;


    //Remember H[l][m]=<l|H|m>
    int m,j;
    double value;
    for (int i=0;i<basis.D_up_basis.size();i++){
        m=i;
        j=i;

        value=0;
        //on-site coulomb repulsion  ni_up ni_dn:
        value+=U*countCommonBits(basis.D_up_basis[i],basis.D_dn_basis[j]);


        //longrange coulomb interaction ni.nj
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                value += NonLocalInteractions_mat[site_i][site_j]*( bit_value(basis.D_up_basis[i], site_i) +
                                                                    bit_value(basis.D_dn_basis[j], site_i) )*
                        ( bit_value(basis.D_up_basis[i], site_j) +
                          bit_value(basis.D_dn_basis[j], site_j) );
            }
        }

        //magnetic Field
        for(int site=0;site<basis.Length;site++){
            value+=0.5*(H_field[site])*
                    ( ( bit_value(basis.D_up_basis[i], site) -
                        bit_value(basis.D_dn_basis[j], site) )
                      );
        }


        if(value!=0){
            Hamil.value.push_back(value*one);
            Hamil.rows.push_back(m);
            Hamil.columns.push_back(m);
        }


        if(m%1000 ==1){
            //cout<<"done "<<m<<" basis"<<endl;
        }

    }


    cout<<"Done Hamiltonian construction: Diagonal"<<endl;

}

void MODEL_1_orb_Hubb_2D_KSector::Add_non_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis){

    //    cout<<"Started Hamiltonian construction: Non Diagonal"<<endl;
    //    cout<<"Done Hamiltonian construction: Non Diagonal"<<endl;

}

void MODEL_1_orb_Hubb_2D_KSector::Add_connections(BASIS_1_orb_Hubb_2D_KSector &basis){


#ifdef _OPENMP
    double begin_time_connections, end_time_connections;
    begin_time_connections = omp_get_wtime();
#endif


    cout<<"Started Hamiltonian construction: Connections"<<endl;
    assert(basis.D_up_basis.size()==basis.D_dn_basis.size());



    Hamil.nrows = basis.D_up_basis.size();
    Hamil.ncols = Hamil.nrows;

    double value;
    int m,j;
    int D_up,D_dn;
    int i_new;
    int m_new;
    double sign_FM;
    int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
    int sign_pow_up, sign_pow_dn;
    int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
    int l,lp;
    int range_min, range_max;
    bool row_found_;
    double_type phase_;
    int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
    complex<double> iota_(0.0,1.0);
    int D_up_temp ,D_dn_temp;
    int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
    bool repeating_rows;
    int row_counter;
    int check_min, check_max;
    int site, site_p;
    int Lxm1, Lym1;
    Lxm1=basis.Lx-1;
    Lym1=basis.Ly-1;

    for (int i=0;i<basis.D_up_basis.size();i++){

        /******TO REMOVE*********************
        cout<<"up("<<i<<") : ";
        print_binary_of_decimal(basis.D_up_basis[i]);
        cout<<"dn("<<i<<") : ";
        print_binary_of_decimal(basis.D_dn_basis[i]);
        cout<<endl;
        ************************************/

        m=i;
        j=i;

        value=0;

        row_counter=0;
        for(int ix=0;ix<basis.Lx ;ix++){
            for(int iy=0;iy<basis.Ly ;iy++){
                site=ix + (iy*basis.Lx);

                for(int ix_p=0;ix_p<basis.Lx ;ix_p++){
                    for(int iy_p=0;iy_p<basis.Ly ;iy_p++){
                        site_p=ix_p + (iy_p*basis.Lx);

                        if(abs(Hopping_mat_NN[site_p][site])!=0.0)// && Dis_y==0)
                        { // LongRange
                            Hopping_NN=Hopping_mat_NN[site_p][site];

                            //---------------Hopping for up electrons-------------------//
                            //there have to be one up electron in site
                            //there have to be no up electron in site_p
                            if(
                                    (bit_value(basis.D_up_basis[i],site)==1)
                                    &&
                                    (bit_value(basis.D_up_basis[i],site_p)==0)

                                    )
                            {

                                sign_pow_up=0;
                                sign_pow_dn=0;

                                D_up = (int) (basis.D_up_basis[i] + pow(2, site_p)
                                              - pow(2,site) );
                                D_dn = basis.D_dn_basis[m] ;

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Inv Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Inv Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Inv Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Inv Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }


                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if(row_found_){break;}

                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*(Inv_Trnsltns_x_)*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*(Inv_Trnsltns_y_)*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                    sign_pow_up += one_bits_in_bw(site,site_p,basis.D_up_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);

                                    if(m_new<=m)
                                    {
                                        repeating_rows=false;
                                        check_min=Hamil.rows.size()-1;
                                        check_max=(Hamil.rows.size()-1)-row_counter;
                                        // cout<<check_min<<endl;
                                        // cout<<check_max<<endl;
                                        for(int check_=check_min;check_>check_max;check_--){
                                            if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                                Hamil.value[check_] +=-1.0*sign_FM*(Hopping_NN)*one*phase_;
                                                repeating_rows=true;
                                                break;
                                            }
                                        }

                                        if(!repeating_rows){
                                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_NN)*one*phase_);
                                            Hamil.rows.push_back(m_new);
                                            Hamil.columns.push_back(m);
                                            row_counter++;
                                        }
                                    }
                                }

                            } // if up hopping possible


                            //---------------Hopping for dn electrons-------------------//
                            //there have to be one dn electron in gamma, site
                            //there have to be no dn electron in gamma_p, site_p
                            if(
                                    (bit_value(basis.D_dn_basis[j], site)==1)
                                    &&
                                    (bit_value(basis.D_dn_basis[j], site_p)==0)
                                    )
                            {

                                sign_pow_dn=0;
                                sign_pow_up=0;

                                D_up = basis.D_up_basis[m];
                                D_dn = (int) (basis.D_dn_basis[j] + pow(2, site_p)
                                              - pow(2, site) );

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;

                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){

                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }




                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }

                                        }
                                    }
                                    if(row_found_){break;}
                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));

#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));

#endif

                                    sign_pow_dn += one_bits_in_bw(site,site_p,basis.D_dn_basis[j]);

                                    sign_FM = pow(-1.0, sign_pow_dn+sign_pow_up);


                                    if(m_new<=m)
                                    {
                                        repeating_rows=false;
                                        check_min=Hamil.rows.size()-1;
                                        check_max=(Hamil.rows.size()-1)-row_counter;
                                        for(int check_=check_min;check_>check_max;check_--){
                                            if(Hamil.rows[check_]==m_new && Hamil.columns[check_]==m){
                                                Hamil.value[check_] +=-1.0*sign_FM*(Hopping_NN)*one*phase_;
                                                repeating_rows=true;
                                                break;
                                            }
                                        }

                                        if(!repeating_rows){
                                            Hamil.value.push_back(-1.0*sign_FM*(Hopping_NN)*one*phase_);
                                            Hamil.rows.push_back(m_new);
                                            Hamil.columns.push_back(m);
                                            row_counter++;
                                        }

                                    }

                                }
                            } // if dn hopping possible


                        }//nearest neighbour

                    }  //iy_p
                } //ix_p

            } //iy
        } // ix

        if(m%10 ==1){
            //cout<<"Connection: done "<<m<<" basis"<<endl;
        }

    } // "i" i.e up_decimals


    cout<<"Done Hamiltonian construction: Connections"<<endl;

#ifdef _OPENMP
    end_time_connections = omp_get_wtime();
    cout<<"Time for adding connections to Hamil [using OpenMP] = "<<double(end_time_connections - begin_time_connections)<<endl;
#endif

}

void MODEL_1_orb_Hubb_2D_KSector::Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_2D_KSector &basis ,BASIS_1_orb_Hubb_2D_KSector & basis_Kminusq){

#ifdef _OPENMP
    double begin_time_dynopr, end_time_dynopr;
    begin_time_dynopr = omp_get_wtime();
#endif

    vector< int >().swap( Dyn_opr.columns );
    vector< int >().swap( Dyn_opr.rows );
    vector< double_type >().swap( Dyn_opr.value );

    Dyn_opr.value.clear();Dyn_opr.rows.clear();
    Dyn_opr.columns.clear();
    Dyn_opr.ncols = basis.D_up_basis.size();
    Dyn_opr.nrows = basis_Kminusq.D_up_basis.size();

    int ix, iy, site;
    double value;
    double_type value2, szq_val, value_final;
    int j, D_up_temp, D_dn_temp, i_new;
    bool row_found_;
    int range_min, range_max;

    for (int i=0;i<basis.D_up_basis.size();i++){
        D_up_temp = basis.D_up_basis[i];
        D_dn_temp = basis.D_dn_basis[i];
        //for (int j=0;j<basis_Kminusq.D_up_basis.size();j++){
        //            if( (basis_Kminusq.D_up_basis[j]==basis.D_up_basis[i])
        //                    &&
        //                    (basis_Kminusq.D_dn_basis[j]==basis.D_dn_basis[i])
        //                    ){

        //---searching same up,dn basis in Kmq sector-----------
        if(D_up_temp>=basis_Kminusq.Dup_Range.size())
        {
            row_found_=false;
        }
        else
        {
            assert(D_up_temp<basis_Kminusq.Dup_Range.size());
            range_min=basis_Kminusq.Dup_Range[D_up_temp].first;
            range_max=basis_Kminusq.Dup_Range[D_up_temp].second;
            if(range_min==-1)
            {
                row_found_=false;
                assert(range_max==-1);
            }
            else
            {
                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis_Kminusq.D_dn_basis, range_min, range_max);
                if(i_new==-1){
                    row_found_=false;
                }
                else{
                    row_found_=true;
                }
            }
        }

        j = i_new;
        //------------------------------------------------------



        if(row_found_){
            assert(basis.Dm_bar[i]==basis_Kminusq.Dm_bar[j]);
            szq_val=zero;
            for(ix=0;ix<basis.Lx;ix++){
                for(iy=0;iy<basis.Ly;iy++){

                    site = ix + iy*(basis.Lx);

                    value=0.5*(1.0)*
                            ( ( bit_value(basis.D_up_basis[i], site) -
                                bit_value(basis.D_dn_basis[i], site) )
                              );

#ifdef USE_COMPLEX
                    value2=exp(iota_comp*(((2.0/(1.0*basis.Lx))*Dyn_Momentum_x*PI*ix) + ((2.0/(1.0*basis.Ly))*Dyn_Momentum_y*PI*iy)))*sqrt(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                    cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif

                    szq_val += value2*value;

                }
            }

            if(szq_val !=zero){

                value_final=(conjugate(basis_Kminusq.Dgamma[j])*basis.Dgamma[i])*
                        szq_val*( (1.0*basis.Lx*basis.Ly)/(basis.Dm_bar[i])  )*
                        (1.0/sqrt(basis_Kminusq.D_Norm[j]*basis.D_Norm[i]));

                Dyn_opr.columns.push_back(i);
                Dyn_opr.rows.push_back(j);
                Dyn_opr.value.push_back(value_final);
            }


        }


        //} //if both have same Representative state
        //}
    }


#ifdef _OPENMP
    end_time_dynopr = omp_get_wtime();
    cout<<"Time for Creating DYN_OPR [using OpenMP] = "<<double(end_time_dynopr - begin_time_dynopr)<<endl;
#endif

}

/*
void MODEL_1_orb_Hubb_2D_KSector::Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_2D_KSector &basis){


    vector< int >().swap( Dyn_opr.columns );
    vector< int >().swap( Dyn_opr.rows );
    vector< double_type >().swap( Dyn_opr.value );
    Dyn_opr.ncols = Hamil.ncols;
    Dyn_opr.nrows = Hamil.nrows;




    //    Hamiltonian_2_COO Oprs_local;
    //    Oprs_local.resize(basis.Lx);
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        Oprs_local[ix].resize(basis.Ly);
    //    }
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        for(int iy=0;iy<basis.Ly;iy++){
    //            Oprs_local[ix][iy].nrows = basis.D_up_basis.size() ;
    //            Oprs_local[ix][iy].ncols = Oprs_local[ix][iy].nrows;
    //        }
    //    }

    if(Dyn_opr_string == "Sz"){


        //Remember O[l][m]=<l|O|m>
        int m,j, site;
        double value;
        double_type value2;
        Matrix_COO Oprs_local;
        Oprs_local.nrows = basis.D_up_basis.size() ;
        Oprs_local.ncols = Oprs_local.nrows;
        Matrix_COO temp;
        temp.nrows = basis.D_up_basis.size() ;
        temp.ncols = temp.nrows;

        for(int ix=0;ix<basis.Lx;ix++){
            for(int iy=0;iy<basis.Ly;iy++){
                Oprs_local.value.clear();
                Oprs_local.rows.clear();
                Oprs_local.columns.clear();

                for (int i=0;i<basis.D_up_basis.size();i++){
                    m=i;
                    j=i;
                    site = ix + iy*(basis.Lx);

                    value=0.5*(1.0)*
                            ( ( bit_value(basis.D_up_basis[i], site) -
                                bit_value(basis.D_dn_basis[j], site) )
                              );

                    if(value!=0){
                        Oprs_local.value.push_back(value*one);
                        Oprs_local.rows.push_back(m);
                        Oprs_local.columns.push_back(m);
                    }
                }


#ifdef USE_COMPLEX
                value2=exp(iota_comp*((Dyn_Momentum_x*PI*ix) + (Dyn_Momentum_y*PI*iy)))*sqrt(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif



                Sum(temp, Oprs_local, temp, 1.0, value2);


                vector< int >().swap( Oprs_local.columns );
                vector< int >().swap( Oprs_local.rows );
                vector< double_type >().swap( Oprs_local.value );


            }
        }

        Dyn_opr = temp;
        vector< int >().swap( temp.columns );
        vector< int >().swap( temp.rows );
        vector< double_type >().swap( temp.value );



        //Now multiplying with sum_{lx,ly}exp(iota*k.l) where l leads to distinct classes
        Mat_1_int lx_; //saves l for creating distinct states of class
        Mat_1_int ly_;
        Mat_1_int D_up_states,D_dn_states;
        int D_up, D_up_temp, D_up_temp_Xtrans;
        int D_dn, D_dn_temp, D_dn_temp_Xtrans;
        int distinct_states;
        bool distinct_, distinct_up, distinct_dn;
        int index;
        double_type coeff;

        for (int row=0;row<Dyn_opr.rows.size();row++){

            index=Dyn_opr.rows[row];
            D_up =basis.D_up_basis[index];
            D_dn =basis.D_dn_basis[index];

            D_up_temp=D_up;
            D_dn_temp=D_dn;
            distinct_states=0;
            lx_.clear();
            ly_.clear();
            D_up_states.clear();
            D_dn_states.clear();
            coeff=zero;
            for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                if(inv_trnsltns_x>0 && basis.Lx>1){

                    for(int iy_=0;iy_<basis.Ly;iy_++){
                        D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);
                        D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);
                    }
                }
                D_dn_temp_Xtrans = D_dn_temp;
                D_up_temp_Xtrans = D_up_temp;

                for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){

                    if(inv_trnsltns_y>0 && basis.Ly>1){

                        for(int ix_=0;ix_<basis.Lx;ix_++){
                            D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);
                            D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);
                        }

                    }

                    if(distinct_states==0)
                    {
                        lx_.push_back(inv_trnsltns_x);
                        ly_.push_back(inv_trnsltns_y);
                        D_up_states.push_back(D_up_temp);
                        D_dn_states.push_back(D_dn_temp);
                        distinct_states++;
                    }
                    else{
                        distinct_up = !(Is_int_in_array(D_up_temp, D_up_states));
                        distinct_dn = !(Is_int_in_array(D_dn_temp, D_dn_states));
                        distinct_ = (distinct_up && distinct_dn);
                        if(distinct_){
                            lx_.push_back(inv_trnsltns_x);
                            ly_.push_back(inv_trnsltns_y);
                            D_up_states.push_back(D_up_temp);
                            D_dn_states.push_back(D_dn_temp);
                            distinct_states++;
                        }
                    }
                }

                D_up_temp=D_up_temp_Xtrans ;
                D_dn_temp=D_dn_temp_Xtrans ;

            }


            for(int l=0;l<lx_.size();l++){
#ifdef USE_COMPLEX
                coeff +=exp(iota_comp*((Dyn_Momentum_x*PI*lx_[l]) + (Dyn_Momentum_y*PI*ly_[l])));
#endif
#ifndef USE_COMPLEX
                cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif
            }

            Dyn_opr.value[row] = Dyn_opr.value[row]*(coeff)*(1.0/(distinct_states*1.0));



        } //i basis


    }
    else{
        cout<<"Only Sz operator is implemented for dynamics"<<endl;
        assert(false);
    }




}

*/
void MODEL_1_orb_Hubb_2D_KSector::Getting_Local_Sz_Opr(BASIS_1_orb_Hubb_2D_KSector &basis, Matrix_COO & OPR_LOCAL, int site){

    vector< int >().swap( OPR_LOCAL.columns );
    vector< int >().swap( OPR_LOCAL.rows );
    vector< double_type >().swap(OPR_LOCAL.value );
    OPR_LOCAL.ncols = Hamil.ncols;
    OPR_LOCAL.nrows = Hamil.nrows;

    double value;

    assert(false);


}

void MODEL_1_orb_Hubb_2D_KSector::Initialize_Opr_for_Structure_factor(BASIS_1_orb_Hubb_2D_KSector &basis){


    vector< int >().swap( OPR_SF.columns );
    vector< int >().swap( OPR_SF.rows );
    vector< double_type >().swap(OPR_SF.value );
    OPR_SF.ncols = basis.D_up_basis.size();
    OPR_SF.nrows = basis.D_up_basis.size();




    //    Hamiltonian_2_COO Oprs_local;
    //    Oprs_local.resize(basis.Lx);
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        Oprs_local[ix].resize(basis.Ly);
    //    }
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        for(int iy=0;iy<basis.Ly;iy++){
    //            Oprs_local[ix][iy].nrows = basis.D_up_basis.size() ;
    //            Oprs_local[ix][iy].ncols = Oprs_local[ix][iy].nrows;
    //        }
    //    }

    if(Dyn_opr_string == "Sz"){


        //Remember O[l][m]=<l|O|m>
        int m,j, site_i, site_j;
        double value;
        double_type value2;
        Matrix_COO Oprs_ij;
        Oprs_ij.nrows = basis.D_up_basis.size() ;
        Oprs_ij.ncols = Oprs_ij.nrows;
        Matrix_COO temp;
        temp.nrows = basis.D_up_basis.size() ;
        temp.ncols = temp.nrows;

        for(int jx=0;jx<basis.Lx;jx++){
            for(int jy=0;jy<basis.Ly;jy++){
                site_j = jx + jy*(basis.Lx);


                for(int ix=0;ix<basis.Lx;ix++){
                    for(int iy=0;iy<basis.Ly;iy++){
                        site_i = ix + iy*(basis.Lx);

                        //if(site_i==site_j){


                        Oprs_ij.value.clear();
                        Oprs_ij.rows.clear();
                        Oprs_ij.columns.clear();

                        for (int i=0;i<basis.D_up_basis.size();i++){

                            m=i;
                            j=i;

                            value=0.25*(1.0)*
                                    ( ( bit_value(basis.D_up_basis[i], site_i) -
                                        bit_value(basis.D_dn_basis[j], site_i) )*
                                      ( bit_value(basis.D_up_basis[i], site_j) -
                                        bit_value(basis.D_dn_basis[j], site_j) )
                                      );

                            if(value!=0){
                                Oprs_ij.value.push_back(value*one);
                                Oprs_ij.rows.push_back(m);
                                Oprs_ij.columns.push_back(m);
                            }
                        }


#ifdef USE_COMPLEX
                        value2=exp(iota_comp*( ((2.0/(1.0*basis.Lx))*Dyn_Momentum_x*PI*(ix-jx)) + ((2.0/(1.0*basis.Ly))*Dyn_Momentum_y*PI*(iy-jy))  ))*(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                        cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif


                        Sum(temp, Oprs_ij, temp, 1.0, value2);


                        vector< int >().swap( Oprs_ij.columns );
                        vector< int >().swap( Oprs_ij.rows );
                        vector< double_type >().swap( Oprs_ij.value );


                        //}
                    }
                }

            }
        }


        OPR_SF = temp;
        vector< int >().swap( temp.columns );
        vector< int >().swap( temp.rows );
        vector< double_type >().swap( temp.value );


    }

}




double_type MODEL_1_orb_Hubb_2D_KSector::Measure_Opr_for_Structure_factor(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &EigVec_){


    double_type value_final;

    if(Dyn_opr_string == "Sz"){


        //Remember O[l][m]=<l|O|m>

        value_final =zero;
        int m,j, site_i, site_j;
        double value;
        double_type value2;
        Matrix_COO Oprs_ij;
        Oprs_ij.nrows = basis.D_up_basis.size() ;
        Oprs_ij.ncols = Oprs_ij.nrows;
        Matrix_COO temp;
        temp.nrows = basis.D_up_basis.size() ;
        temp.ncols = temp.nrows;



        //if(site_i==site_j){


        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            for(int jx=0;jx<basis.Lx;jx++){
                for(int jy=0;jy<basis.Ly;jy++){
                    site_j = jx + jy*(basis.Lx);


                    for(int ix=0;ix<basis.Lx;ix++){
                        for(int iy=0;iy<basis.Ly;iy++){
                            site_i = ix + iy*(basis.Lx);

                            value=0.25*(1.0)*
                                    ( ( bit_value(basis.D_up_basis[i], site_i) -
                                        bit_value(basis.D_dn_basis[j], site_i) )*
                                      ( bit_value(basis.D_up_basis[i], site_j) -
                                        bit_value(basis.D_dn_basis[j], site_j) )
                                      );





#ifdef USE_COMPLEX
                            value2=exp(iota_comp*( ((2.0/(1.0*basis.Lx))*Dyn_Momentum_x*PI*(ix-jx)) + ((2.0/(1.0*basis.Ly))*Dyn_Momentum_y*PI*(iy-jy))  ))*(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                            cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif

                            if(value!=0){

                                //                                Oprs_ij.value.push_back(value*one);
                                //                                Oprs_ij.rows.push_back(m);
                                //                                Oprs_ij.columns.push_back(m);
                                value_final += conjugate(EigVec_[m])*value*one*value2*EigVec_[m];

                            }
                        }

                        //}
                    }
                }

            }
        }

    }


    return value_final;
}


double_type MODEL_1_orb_Hubb_2D_KSector::Measure_Opr_for_LocalNupNdn(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &EigVec_){


    double_type value_final;
    value_final=zero;


    for (int i=0;i<basis.D_up_basis.size();i++){
        int m,j, site_i;
        double value;

        for(int ix=0;ix<basis.Lx;ix++){
            for(int iy=0;iy<basis.Ly;iy++){
                site_i = ix + iy*(basis.Lx);

                m=i;
                j=i;

                value=(1.0)*
                        (bit_value(basis.D_up_basis[i], site_i))*
                        (bit_value(basis.D_dn_basis[i], site_i));

                if(value!=0){
                    value_final += EigVec_[m]*conjugate(EigVec_[m])*value*one*(1.0/(basis.Length));
                }
            }
        }

    }

    return value_final;


}


void MODEL_1_orb_Hubb_2D_KSector::Initialize_Oprs_for_meausurement(BASIS_1_orb_Hubb_2D_KSector &basis){

    //    Mat_1_string obs_string;
    //    Hamiltonian_2_COO Oprts_array;

    obs_string.clear();
    obs_string.push_back("SzSz_00"); //q1=(0,0) i.e. FM
    obs_string.push_back("SzSz_pipi"); //q1=(1,1) i.e. AFM

    obs_string.push_back("Local_nup_ndn");


    Oprts_array.clear();
    Oprts_array.resize(3);



    Dyn_Momentum_x=0;Dyn_Momentum_y=0;
    //***************Sq***********************************//
    vector< int >().swap( OPR_SF.columns );
    vector< int >().swap( OPR_SF.rows );
    vector< double_type >().swap(OPR_SF.value );
    OPR_SF.ncols = basis.D_up_basis.size();
    OPR_SF.nrows = basis.D_up_basis.size();




    //    Hamiltonian_2_COO Oprs_local;
    //    Oprs_local.resize(basis.Lx);
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        Oprs_local[ix].resize(basis.Ly);
    //    }
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        for(int iy=0;iy<basis.Ly;iy++){
    //            Oprs_local[ix][iy].nrows = basis.D_up_basis.size() ;
    //            Oprs_local[ix][iy].ncols = Oprs_local[ix][iy].nrows;
    //        }
    //    }


    //Remember O[l][m]=<l|O|m>
    int m,j, site_i, site_j;
    double value;
    double_type value2;
    Matrix_COO Oprs_ij;
    Oprs_ij.nrows = basis.D_up_basis.size() ;
    Oprs_ij.ncols = Oprs_ij.nrows;
    Matrix_COO temp;
    temp.nrows = basis.D_up_basis.size() ;
    temp.ncols = temp.nrows;

    for(int jx=0;jx<basis.Lx;jx++){
        for(int jy=0;jy<basis.Ly;jy++){
            site_j = jx + jy*(basis.Lx);


            for(int ix=0;ix<basis.Lx;ix++){
                for(int iy=0;iy<basis.Ly;iy++){
                    site_i = ix + iy*(basis.Lx);

                    //if(site_i==site_j){


                    Oprs_ij.value.clear();
                    Oprs_ij.rows.clear();
                    Oprs_ij.columns.clear();

                    for (int i=0;i<basis.D_up_basis.size();i++){

                        m=i;
                        j=i;

                        value=0.25*(1.0)*
                                ( ( bit_value(basis.D_up_basis[i], site_i) -
                                    bit_value(basis.D_dn_basis[j], site_i) )*
                                  ( bit_value(basis.D_up_basis[i], site_j) -
                                    bit_value(basis.D_dn_basis[j], site_j) )
                                  );

                        if(value!=0){
                            Oprs_ij.value.push_back(value*one);
                            Oprs_ij.rows.push_back(m);
                            Oprs_ij.columns.push_back(m);
                        }
                    }


#ifdef USE_COMPLEX
                    value2=exp(iota_comp*( ((2.0/(1.0*basis.Lx))*Dyn_Momentum_x*PI*(ix-jx)) + ((2.0/(1.0*basis.Ly))*Dyn_Momentum_y*PI*(iy-jy))  ))*(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                    cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif

                    Sum(temp, Oprs_ij, temp, 1.0*one, value2);

                    vector< int >().swap( Oprs_ij.columns );
                    vector< int >().swap( Oprs_ij.rows );
                    vector< double_type >().swap( Oprs_ij.value );

                    //}
                }
            }
        }
    }

    OPR_SF = temp;
    vector< int >().swap( temp.columns );
    vector< int >().swap( temp.rows );
    vector< double_type >().swap( temp.value );


    Oprts_array[0]=OPR_SF;

    vector< int >().swap( OPR_SF.columns );
    vector< int >().swap( OPR_SF.rows );
    vector< double_type >().swap(OPR_SF.value );

    //****************************************************//



    Dyn_Momentum_x=basis.Lx/2;
    Dyn_Momentum_y=basis.Ly/2;
    //***************Sq***********************************//
    vector< int >().swap( OPR_SF.columns );
    vector< int >().swap( OPR_SF.rows );
    vector< double_type >().swap(OPR_SF.value );
    OPR_SF.ncols = Hamil.ncols;
    OPR_SF.nrows = Hamil.nrows;




    //    Hamiltonian_2_COO Oprs_local;
    //    Oprs_local.resize(basis.Lx);
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        Oprs_local[ix].resize(basis.Ly);
    //    }
    //    for(int ix=0;ix<basis.Lx;ix++){
    //        for(int iy=0;iy<basis.Ly;iy++){
    //            Oprs_local[ix][iy].nrows = basis.D_up_basis.size() ;
    //            Oprs_local[ix][iy].ncols = Oprs_local[ix][iy].nrows;
    //        }
    //    }


    //Remember O[l][m]=<l|O|m>

    Oprs_ij.nrows = basis.D_up_basis.size() ;
    Oprs_ij.ncols = Oprs_ij.nrows;

    temp.nrows = basis.D_up_basis.size() ;
    temp.ncols = temp.nrows;

    for(int jx=0;jx<basis.Lx;jx++){
        for(int jy=0;jy<basis.Ly;jy++){
            site_j = jx + jy*(basis.Lx);


            for(int ix=0;ix<basis.Lx;ix++){
                for(int iy=0;iy<basis.Ly;iy++){
                    site_i = ix + iy*(basis.Lx);

                    Oprs_ij.value.clear();
                    Oprs_ij.rows.clear();
                    Oprs_ij.columns.clear();

                    for (int i=0;i<basis.D_up_basis.size();i++){

                        m=i;
                        j=i;

                        value=0.25*(1.0)*
                                ( ( bit_value(basis.D_up_basis[i], site_i) -
                                    bit_value(basis.D_dn_basis[j], site_i) )*
                                  ( bit_value(basis.D_up_basis[i], site_j) -
                                    bit_value(basis.D_dn_basis[j], site_j) )
                                  );

                        if(value!=0){
                            Oprs_ij.value.push_back(value*one);
                            Oprs_ij.rows.push_back(m);
                            Oprs_ij.columns.push_back(m);
                        }
                    }


#ifdef USE_COMPLEX
                    value2=exp(iota_comp*( ((2.0/(1.0*basis.Lx))*Dyn_Momentum_x*PI*(ix-jx)) + ((2.0/(1.0*basis.Ly))*Dyn_Momentum_y*PI*(iy-jy))  ))*(1.0/(basis.Length));
#endif
#ifndef USE_COMPLEX
                    cout<<"For PBC=true and Dynamics=true, compile with USE_COMPLEX"<<endl;
#endif

                    Sum(temp, Oprs_ij, temp, 1.0*one, value2);

                    vector< int >().swap( Oprs_ij.columns );
                    vector< int >().swap( Oprs_ij.rows );
                    vector< double_type >().swap( Oprs_ij.value );

                }
            }
        }
    }

    OPR_SF = temp;
    vector< int >().swap( temp.columns );
    vector< int >().swap( temp.rows );
    vector< double_type >().swap( temp.value );


    Oprts_array[1]=OPR_SF;

    vector< int >().swap( OPR_SF.columns );
    vector< int >().swap( OPR_SF.rows );
    vector< double_type >().swap(OPR_SF.value );

    //****************************************************//



    //******Local NupNdn**********************************//

    vector< int >().swap( Oprts_array[2].columns );
    vector< int >().swap( Oprts_array[2].rows );
    vector< double_type >().swap(Oprts_array[2].value );
    Oprts_array[2].ncols = Hamil.ncols;
    Oprts_array[2].nrows = Hamil.nrows;


    //Remember O[l][m]=<l|O|m>
    Matrix_COO Oprs_i;
    Oprs_i.nrows = basis.D_up_basis.size() ;
    Oprs_i.ncols = Oprs_i.nrows;

    temp.columns.clear();temp.rows.clear();temp.value.clear();
    temp.nrows = basis.D_up_basis.size();
    temp.ncols = temp.nrows;


    for(int ix=0;ix<basis.Lx;ix++){
        for(int iy=0;iy<basis.Ly;iy++){
            site_i = ix + iy*(basis.Lx);

            Oprs_i.value.clear();
            Oprs_i.rows.clear();
            Oprs_i.columns.clear();

            for (int i=0;i<basis.D_up_basis.size();i++){

                m=i;
                j=i;

                value=(1.0)*
                        (bit_value(basis.D_up_basis[i], site_i))*
                        (bit_value(basis.D_dn_basis[i], site_i));

                if(value!=0){
                    Oprs_i.value.push_back(value*one);
                    Oprs_i.rows.push_back(m);
                    Oprs_i.columns.push_back(m);
                }
            }


            value2=(1.0/(basis.Length));


            Sum(temp, Oprs_i, temp, 1.0*one, value2);

            vector< int >().swap( Oprs_i.columns );
            vector< int >().swap( Oprs_i.rows );
            vector< double_type >().swap( Oprs_i.value );

        }
    }

    Oprts_array[2] = temp;
    vector< int >().swap( temp.columns );
    vector< int >().swap( temp.rows );
    vector< double_type >().swap( temp.value );


    //****************************************************//


}

void MODEL_1_orb_Hubb_2D_KSector::Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site_x, int site_y, BASIS_1_orb_Hubb_2D_KSector &basis){


    OPR_.nrows = basis.D_up_basis.size();
    OPR_.ncols = OPR_.nrows;


    if(opr_type=="SzSz"){
        double value;
        int m,j;
        int site, site_p;
        int ix_p,iy_p;

        OPR_.value.clear();
        OPR_.rows.clear();
        OPR_.columns.clear();
        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;
            value=0;

            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);

                    ix_p = (ix + site_x)%basis.Lx;
                    iy_p = (iy + site_y)%basis.Ly;
                    site_p=ix_p + (iy_p*basis.Lx);

                    //---------------Sz[site]Sz[site_p]-------------------//

                    value +=0.25*( bit_value(basis.D_up_basis[i], site) - bit_value(basis.D_dn_basis[m], site) )*
                            ( bit_value(basis.D_up_basis[i], site_p) - bit_value(basis.D_dn_basis[m], site_p) );




                } //iy
            } // ix

            OPR_.value.push_back(one*value*(1.0/(basis.Ly*basis.Lx)));
            OPR_.rows.push_back(m);
            OPR_.columns.push_back(m);

        } // "i" i.e up_decimals

    }




    if(opr_type=="SpSm"){

        double value;
        int m,j;
        int D_up,D_dn;
        int i_new;
        int m_new;
        double sign_FM;
        int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
        int sign_pow_up, sign_pow_dn;
        int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
        int l,lp;
        int range_min, range_max;
        bool row_found_;
        double_type phase_;
        int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
        complex<double> iota_(0.0,1.0);
        int D_up_temp ,D_dn_temp;
        int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
        bool repeating_rows;
        int row_counter;
        int check_min, check_max;
        int site, site_p;
        int Lxm1, Lym1;
        int ix_p,iy_p;
        Lxm1=basis.Lx-1;
        Lym1=basis.Ly-1;

        OPR_.value.clear();
        OPR_.rows.clear();
        OPR_.columns.clear();
        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            value=0;
            row_counter=0;


            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);

                    ix_p = (ix + site_x)%basis.Lx;
                    iy_p = (iy + site_y)%basis.Ly;
                    site_p=ix_p + (iy_p*basis.Lx);



                    //---------------Sp[site]Sm[site_p]-------------------//


                    //---------------Hopping for up electrons-------------------//
                    //there have to be only one electron in site with down spin
                    //there have to be only one electron in site_p with up spin
                    bool SpSm_allowed;

                    if(site!=site_p){
                        SpSm_allowed = (bit_value(basis.D_dn_basis[i],site)==1) && (bit_value(basis.D_up_basis[i],site)==0)
                                && (bit_value(basis.D_up_basis[i],site_p)==1) && (bit_value(basis.D_dn_basis[i],site_p)==0);
                    }
                    else{
                        SpSm_allowed = (bit_value(basis.D_up_basis[i],site)==1) && (bit_value(basis.D_dn_basis[i],site)==0);
                    }

                    if(SpSm_allowed)
                    {

                        sign_pow_up=0;
                        sign_pow_dn=0;

                        D_up = (int) (basis.D_up_basis[i] - pow(2, site_p)
                                      + pow(2,site) );
                        D_dn = (int) (basis.D_dn_basis[i] + pow(2, site_p)
                                      - pow(2,site) ) ;

                        D_up_temp=D_up;
                        D_dn_temp=D_dn;
                        row_found_=false;

                        for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                            if(inv_trnsltns_x>0 && basis.Lx>1){

                                for(int iy_=0;iy_<basis.Ly;iy_++){

                                    //Inv Translation on spin_dn
                                    sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                            1*bit_value(D_dn_temp,iy_*basis.Lx);
                                    if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                        sign_pow_dn += 1*sign_pow_dn_orb0;
                                    }

                                    D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                    //Inv Translation on spin_up
                                    sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                            1*bit_value(D_up_temp,iy_*basis.Lx);
                                    if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                        sign_pow_up += 1*sign_pow_up_orb0;
                                    }

                                    D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                }
                            }
                            D_dn_temp_Xtrans = D_dn_temp;
                            D_up_temp_Xtrans = D_up_temp;
                            sign_pow_up_Xtrans = sign_pow_up;
                            sign_pow_dn_Xtrans = sign_pow_dn;
                            for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                if(inv_trnsltns_y>0 && basis.Ly>1){

                                    for(int ix_=0;ix_<basis.Lx;ix_++){

                                        //Inv Translation on spin_dn
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                            sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_dn += sign_pow_dn_orb0;
                                        }
                                        sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                bit_value(D_dn_temp, ix_) ;
                                        sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_dn += sign_pow_dn_orb0;

                                        D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                        //Inv Translation on spin_up
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                            sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_up += sign_pow_up_orb0;
                                        }
                                        sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                bit_value(D_up_temp, ix_) ;
                                        sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_up += sign_pow_up_orb0;

                                        D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                    }
                                }


                                if(D_up_temp>=basis.Dup_Range.size())
                                {
                                    row_found_=false;
                                }
                                else
                                {
                                    assert(D_up_temp<basis.Dup_Range.size());
                                    range_min=basis.Dup_Range[D_up_temp].first;
                                    range_max=basis.Dup_Range[D_up_temp].second;
                                    if(range_min==-1)
                                    {
                                        row_found_=false;
                                        assert(range_max==-1);
                                    }
                                    else
                                    {
                                        i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                        if(i_new==-1){
                                            row_found_=false;
                                        }
                                        else{
                                            row_found_=true;
                                            Inv_Trnsltns_x_=inv_trnsltns_x;
                                            Inv_Trnsltns_y_=inv_trnsltns_y;
                                            break;
                                        }
                                    }
                                }
                            }
                            if(row_found_){break;}

                            D_up_temp=D_up_temp_Xtrans ;
                            D_dn_temp=D_dn_temp_Xtrans ;
                            sign_pow_up = sign_pow_up_Xtrans;
                            sign_pow_dn = sign_pow_dn_Xtrans;
                        }


                        if(row_found_==true){
                            m_new = i_new;

#ifdef USE_COMPLEX
                            phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                       )
                                    *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                            if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                            }
                            assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                            phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                            sign_pow_up += one_bits_in_bw(site,site_p,basis.D_up_basis[i]);
                            sign_pow_dn += one_bits_in_bw(site,site_p,basis.D_dn_basis[i]); //may be +1 is needed here
                            if(site!=site_p){
                                sign_pow_dn +=1;
                            }

                            sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                            repeating_rows=false;
                            check_min=OPR_.rows.size()-1;
                            check_max=(OPR_.rows.size()-1)-row_counter;
                            // cout<<check_min<<endl;
                            // cout<<check_max<<endl;
                            for(int check_=check_min;check_>check_max;check_--){
                                if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                    OPR_.value[check_] +=1.0*sign_FM*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                    repeating_rows=true;
                                    break;
                                }
                            }

                            if(!repeating_rows){
                                OPR_.value.push_back(1.0*sign_FM*one*phase_*(1.0/(basis.Lx*basis.Ly)));
                                OPR_.rows.push_back(m_new);
                                OPR_.columns.push_back(m);
                                row_counter++;
                            }

                        }

                    }

                } //iy
            } // ix



        } // "i" i.e up_decimals

    }


}



double_type MODEL_1_orb_Hubb_2D_KSector::Measure_two_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_ , int site_x, int site_y, BASIS_1_orb_Hubb_2D_KSector &basis){


    double_type value_final;
    value_final=zero;


    if(opr_type=="SzSz"){
        double value;
        int m,j;
        int site, site_p;
        int ix_p,iy_p;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;
            value=0;

            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);

                    ix_p = (ix + site_x)%basis.Lx;
                    iy_p = (iy + site_y)%basis.Ly;
                    site_p=ix_p + (iy_p*basis.Lx);

                    //---------------Sz[site]Sz[site_p]-------------------//

                    value +=0.25*( bit_value(basis.D_up_basis[i], site) - bit_value(basis.D_dn_basis[m], site) )*
                            ( bit_value(basis.D_up_basis[i], site_p) - bit_value(basis.D_dn_basis[m], site_p) );




                } //iy
            } // ix

            value_final += EigVec_[m]*conjugate(EigVec_[m])*one*value*(1.0/(basis.Ly*basis.Lx));

        } // "i" i.e up_decimals

    }




    if(opr_type=="SpSm"){

        double value;
        int m,j;
        int D_up,D_dn;
        int i_new;
        int m_new;
        double sign_FM;
        int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
        int sign_pow_up, sign_pow_dn;
        int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
        int l,lp;
        int range_min, range_max;
        bool row_found_;
        double_type phase_;
        int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
        complex<double> iota_(0.0,1.0);
        int D_up_temp ,D_dn_temp;
        int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
        bool repeating_rows;
        int row_counter;
        int check_min, check_max;
        int site, site_p;
        int Lxm1, Lym1;
        int ix_p,iy_p;
        Lxm1=basis.Lx-1;
        Lym1=basis.Ly-1;

        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            value=0;
            row_counter=0;


            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);

                    ix_p = (ix + site_x)%basis.Lx;
                    iy_p = (iy + site_y)%basis.Ly;
                    site_p=ix_p + (iy_p*basis.Lx);



                    //---------------Sp[site]Sm[site_p]-------------------//


                    //---------------Hopping for up electrons-------------------//
                    //there have to be only one electron in site with down spin
                    //there have to be only one electron in site_p with up spin
                    bool SpSm_allowed;

                    if(site!=site_p){
                        SpSm_allowed = (bit_value(basis.D_dn_basis[i],site)==1) && (bit_value(basis.D_up_basis[i],site)==0)
                                && (bit_value(basis.D_up_basis[i],site_p)==1) && (bit_value(basis.D_dn_basis[i],site_p)==0);
                    }
                    else{
                        SpSm_allowed = (bit_value(basis.D_up_basis[i],site)==1) && (bit_value(basis.D_dn_basis[i],site)==0);
                    }

                    if(SpSm_allowed)
                    {

                        sign_pow_up=0;
                        sign_pow_dn=0;

                        D_up = (int) (basis.D_up_basis[i] - pow(2, site_p)
                                      + pow(2,site) );
                        D_dn = (int) (basis.D_dn_basis[i] + pow(2, site_p)
                                      - pow(2,site) ) ;

                        D_up_temp=D_up;
                        D_dn_temp=D_dn;
                        row_found_=false;

                        for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                            if(inv_trnsltns_x>0 && basis.Lx>1){

                                for(int iy_=0;iy_<basis.Ly;iy_++){

                                    //Inv Translation on spin_dn
                                    sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                            1*bit_value(D_dn_temp,iy_*basis.Lx);
                                    if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                        sign_pow_dn += 1*sign_pow_dn_orb0;
                                    }

                                    D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                    //Inv Translation on spin_up
                                    sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                            1*bit_value(D_up_temp,iy_*basis.Lx);
                                    if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                        sign_pow_up += 1*sign_pow_up_orb0;
                                    }

                                    D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                }
                            }
                            D_dn_temp_Xtrans = D_dn_temp;
                            D_up_temp_Xtrans = D_up_temp;
                            sign_pow_up_Xtrans = sign_pow_up;
                            sign_pow_dn_Xtrans = sign_pow_dn;
                            for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                if(inv_trnsltns_y>0 && basis.Ly>1){

                                    for(int ix_=0;ix_<basis.Lx;ix_++){

                                        //Inv Translation on spin_dn
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                            sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_dn += sign_pow_dn_orb0;
                                        }
                                        sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                bit_value(D_dn_temp, ix_) ;
                                        sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_dn += sign_pow_dn_orb0;

                                        D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                        //Inv Translation on spin_up
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                            sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_up += sign_pow_up_orb0;
                                        }
                                        sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                bit_value(D_up_temp, ix_) ;
                                        sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_up += sign_pow_up_orb0;

                                        D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                    }
                                }


                                if(D_up_temp>=basis.Dup_Range.size())
                                {
                                    row_found_=false;
                                }
                                else
                                {
                                    assert(D_up_temp<basis.Dup_Range.size());
                                    range_min=basis.Dup_Range[D_up_temp].first;
                                    range_max=basis.Dup_Range[D_up_temp].second;
                                    if(range_min==-1)
                                    {
                                        row_found_=false;
                                        assert(range_max==-1);
                                    }
                                    else
                                    {
                                        i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                        if(i_new==-1){
                                            row_found_=false;
                                        }
                                        else{
                                            row_found_=true;
                                            Inv_Trnsltns_x_=inv_trnsltns_x;
                                            Inv_Trnsltns_y_=inv_trnsltns_y;
                                            break;
                                        }
                                    }
                                }
                            }
                            if(row_found_){break;}

                            D_up_temp=D_up_temp_Xtrans ;
                            D_dn_temp=D_dn_temp_Xtrans ;
                            sign_pow_up = sign_pow_up_Xtrans;
                            sign_pow_dn = sign_pow_dn_Xtrans;
                        }


                        if(row_found_==true){
                            m_new = i_new;

#ifdef USE_COMPLEX
                            phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                       )
                                    *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                            if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                            }
                            assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                            phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                            sign_pow_up += one_bits_in_bw(site,site_p,basis.D_up_basis[i]);
                            sign_pow_dn += one_bits_in_bw(site,site_p,basis.D_dn_basis[i]); //may be +1 is needed here
                            if(site!=site_p){
                                sign_pow_dn +=1;
                            }

                            sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


//                            repeating_rows=false;
//                            check_min=OPR_.rows.size()-1;
//                            check_max=(OPR_.rows.size()-1)-row_counter;
//                            // cout<<check_min<<endl;
//                            // cout<<check_max<<endl;
//                            for(int check_=check_min;check_>check_max;check_--){
//                                if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
//                                    OPR_.value[check_] +=1.0*sign_FM*one*phase_*(1.0/(basis.Lx*basis.Ly));
//                                    repeating_rows=true;
//                                    break;
//                                }
//                            }

//                            if(!repeating_rows){
//                                OPR_.value.push_back(1.0*sign_FM*one*phase_*(1.0/(basis.Lx*basis.Ly)));
//                                OPR_.rows.push_back(m_new);
//                                OPR_.columns.push_back(m);
//                                row_counter++;
//                            }

                            value_final += 1.0*sign_FM*one*phase_*(1.0/(basis.Lx*basis.Ly))*conjugate(EigVec_[m_new])*EigVec_[m];

                        }

                    }

                } //iy
            } // ix



        } // "i" i.e up_decimals

    }


    return value_final;


}



void MODEL_1_orb_Hubb_2D_KSector::Initialize_three_point_operator_sites_specific(string opr_type, Matrix_COO &OPR_, int sitejx, int sitejy, int sitelx, int sitely, BASIS_1_orb_Hubb_2D_KSector &basis){

    OPR_.nrows = basis.D_up_basis.size();
    OPR_.ncols = OPR_.nrows;


    if(opr_type=="SzSpSm"){

        assert ( ((sitejx*sitejx) + (sitejy*sitejy)) !=0 );
        assert ( ((sitelx*sitelx) + (sitely*sitely)) !=0 );

        assert (  !((sitejx==sitelx) && (sitejy==sitely))    );
        double value;
        int m,j;
        int D_up,D_dn;
        int i_new;
        int m_new;
        double sign_FM;
        int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
        int sign_pow_up, sign_pow_dn;
        int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
        int l,lp;
        int range_min, range_max;
        bool row_found_;
        double_type phase_;
        double Sz_value;
        int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
        complex<double> iota_(0.0,1.0);
        int D_up_temp ,D_dn_temp;
        int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
        bool repeating_rows;
        int row_counter;
        int check_min, check_max;
        int site, site_j, site_l;
        int Lxm1, Lym1;
        int ix_j,iy_j;
        int ix_l,iy_l;
        Lxm1=basis.Lx-1;
        Lym1=basis.Ly-1;

        OPR_.value.clear();
        OPR_.rows.clear();
        OPR_.columns.clear();
        for (int i=0;i<basis.D_up_basis.size();i++){

            m=i;
            j=i;

            value=0;
            row_counter=0;


            for(int ix=0;ix<basis.Lx ;ix++){
                for(int iy=0;iy<basis.Ly ;iy++){
                    site=ix + (iy*basis.Lx);
                    Sz_value = 0.5*(bit_value(basis.D_up_basis[i], site) - bit_value(basis.D_dn_basis[j], site) );

                    ix_j = (ix + sitejx)%basis.Lx;
                    iy_j = (iy + sitejy)%basis.Ly;
                    site_j=ix_j + (iy_j*basis.Lx);

                    ix_l = (ix + sitelx)%basis.Lx;
                    iy_l = (iy + sitely)%basis.Ly;
                    site_l=ix_l + (iy_l*basis.Lx);



                    //---------------Sz[site]Sp[site_j]Sm[site_l]-------------------//


                    //---------------Hopping for up electrons-------------------//
                    //there have to be only one electron in site_l with up spin
                    //there have to be only one electron in site_j with down spin
                    bool SpSm_allowed;

                    SpSm_allowed = (bit_value(basis.D_dn_basis[i],site_j)==1) && (bit_value(basis.D_up_basis[i],site_j)==0)
                            && (bit_value(basis.D_up_basis[i],site_l)==1) && (bit_value(basis.D_dn_basis[i],site_l)==0);

                    if(SpSm_allowed)
                    {

                        sign_pow_up=0;
                        sign_pow_dn=0;

                        D_up = (int) (basis.D_up_basis[i] - pow(2, site_l)
                                      + pow(2,site_j) );
                        D_dn = (int) (basis.D_dn_basis[i] + pow(2, site_l)
                                      - pow(2,site_j) ) ;

                        D_up_temp=D_up;
                        D_dn_temp=D_dn;
                        row_found_=false;

                        for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                            if(inv_trnsltns_x>0 && basis.Lx>1){

                                for(int iy_=0;iy_<basis.Ly;iy_++){

                                    //Inv Translation on spin_dn
                                    sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                            1*bit_value(D_dn_temp,iy_*basis.Lx);
                                    if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                        sign_pow_dn += 1*sign_pow_dn_orb0;
                                    }

                                    D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                    //Inv Translation on spin_up
                                    sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                            1*bit_value(D_up_temp,iy_*basis.Lx);
                                    if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                        sign_pow_up += 1*sign_pow_up_orb0;
                                    }

                                    D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                }
                            }
                            D_dn_temp_Xtrans = D_dn_temp;
                            D_up_temp_Xtrans = D_up_temp;
                            sign_pow_up_Xtrans = sign_pow_up;
                            sign_pow_dn_Xtrans = sign_pow_dn;
                            for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                if(inv_trnsltns_y>0 && basis.Ly>1){

                                    for(int ix_=0;ix_<basis.Lx;ix_++){

                                        //Inv Translation on spin_dn
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                            sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_dn += sign_pow_dn_orb0;
                                        }
                                        sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                bit_value(D_dn_temp, ix_) ;
                                        sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_dn += sign_pow_dn_orb0;

                                        D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                        //Inv Translation on spin_up
                                        for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                            sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                            sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                            sign_pow_up += sign_pow_up_orb0;
                                        }
                                        sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                bit_value(D_up_temp, ix_) ;
                                        sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                        sign_pow_up += sign_pow_up_orb0;

                                        D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                    }
                                }


                                if(D_up_temp>=basis.Dup_Range.size())
                                {
                                    row_found_=false;
                                }
                                else
                                {
                                    assert(D_up_temp<basis.Dup_Range.size());
                                    range_min=basis.Dup_Range[D_up_temp].first;
                                    range_max=basis.Dup_Range[D_up_temp].second;
                                    if(range_min==-1)
                                    {
                                        row_found_=false;
                                        assert(range_max==-1);
                                    }
                                    else
                                    {
                                        i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                        if(i_new==-1){
                                            row_found_=false;
                                        }
                                        else{
                                            row_found_=true;
                                            Inv_Trnsltns_x_=inv_trnsltns_x;
                                            Inv_Trnsltns_y_=inv_trnsltns_y;
                                            break;
                                        }
                                    }
                                }
                            }
                            if(row_found_){break;}

                            D_up_temp=D_up_temp_Xtrans ;
                            D_dn_temp=D_dn_temp_Xtrans ;
                            sign_pow_up = sign_pow_up_Xtrans;
                            sign_pow_dn = sign_pow_dn_Xtrans;
                        }


                        if(row_found_==true){
                            m_new = i_new;

#ifdef USE_COMPLEX
                            phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                       )
                                    *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                            if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                            }
                            assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                            phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                            sign_pow_up += one_bits_in_bw(site_j,site_l,basis.D_up_basis[i]);
                            sign_pow_dn += one_bits_in_bw(site_j,site_l,basis.D_dn_basis[i]); //may be +1 is needed here

                            sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                            repeating_rows=false;
                            check_min=OPR_.rows.size()-1;
                            check_max=(OPR_.rows.size()-1)-row_counter;
                            // cout<<check_min<<endl;
                            // cout<<check_max<<endl;
                            for(int check_=check_min;check_>check_max;check_--){
                                if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                    OPR_.value[check_] +=1.0*sign_FM*Sz_value*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                    repeating_rows=true;
                                    break;
                                }
                            }

                            if(!repeating_rows){
                                OPR_.value.push_back(1.0*sign_FM*one*Sz_value*phase_*(1.0/(basis.Lx*basis.Ly)));
                                OPR_.rows.push_back(m_new);
                                OPR_.columns.push_back(m);
                                row_counter++;
                            }

                        }

                    }

                } //iy
            } // ix



        } // "i" i.e up_decimals

    }


}

void MODEL_1_orb_Hubb_2D_KSector::Initialize_chiral_corr_operator_sites_specific(string opr_type, Matrix_COO &OPR_, int lx_rel, int ly_rel , BASIS_1_orb_Hubb_2D_KSector &basis){

    OPR_.nrows = basis.D_up_basis.size();
    OPR_.ncols = OPR_.nrows;

    int jx_rel, jy_rel, kx_rel, ky_rel, mx_rel, my_rel, nx_rel, ny_rel;

    //Delta1_Z = Sz[i] S+[i+ (-1,1)] S-[i+(1,1)] + ---
    //Delta2_Z = Sz[l] S+[l+(1,1)] S-[l+(1,0)] + ---
    if(opr_type=="Delta1_Z_Delta2_Z"){

        //Delta1
        jx_rel=-1;jy_rel=1;
        kx_rel=1;ky_rel=1;

        //Delta2
        mx_rel=1;my_rel=1;
        nx_rel=1;ny_rel=0;

    }
    if(opr_type=="Delta1_Z_Delta1_Z"){
        //Delta1
        jx_rel=-1;jy_rel=1;
        kx_rel=1;ky_rel=1;

        //Delta1
        mx_rel=-1;my_rel=1;
        nx_rel=1;ny_rel=1;
    }



    //assert ( ((sitejx*sitejx) + (sitejy*sitejy)) !=0 );
    //assert ( ((sitelx*sitelx) + (sitely*sitely)) !=0 );
    //assert (  !((sitejx==sitelx) && (sitejy==sitely))    );

    double value;
    int m,j;
    int D_up,D_dn;
    int i_new;
    int m_new;
    double sign_FM;
    int sign_pow_dn_orb0, sign_pow_dn_orb1, sign_pow_up_orb0, sign_pow_up_orb1;
    int sign_pow_up, sign_pow_dn;
    int sign_pow_up_Xtrans, sign_pow_dn_Xtrans;
    int l,lp;
    int range_min, range_max;
    bool row_found_;
    double_type phase_;
    double Sz_value, Sz_sitel, Sz_sitei;
    int Inv_Trnsltns_x_, Inv_Trnsltns_y_ ;
    complex<double> iota_(0.0,1.0);
    int D_up_temp ,D_dn_temp;
    int D_up_temp_Xtrans ,D_dn_temp_Xtrans;
    bool repeating_rows;
    int row_counter;
    int check_min, check_max;
    int site_i, site_j, site_k, site_l, site_m, site_n;
    int Lxm1, Lym1;
    int sitejx, sitejy, sitekx, siteky, sitelx, sitely, sitemx, sitemy, sitenx, siteny;
    Lxm1=basis.Lx-1;
    Lym1=basis.Ly-1;
    bool allowed;
    double Value_offset;


    OPR_.value.clear();
    OPR_.rows.clear();
    OPR_.columns.clear();



    for (int i=0;i<basis.D_up_basis.size();i++){

        m=i;
        j=i;

        value=0;
        row_counter=0;


        for(int ix=0;ix<basis.Lx ;ix++){
            for(int iy=0;iy<basis.Ly ;iy++){

                site_i=ix + (iy*basis.Lx);

                sitejx = (ix + jx_rel)%basis.Lx;
                sitejy = (iy + jy_rel)%basis.Ly;
                site_j=sitejx + (sitejy*basis.Lx);

                sitekx = (ix + kx_rel)%basis.Lx;
                siteky = (iy + ky_rel)%basis.Ly;
                site_k=sitekx + (siteky*basis.Lx);


                sitelx = (ix + lx_rel)%basis.Lx;
                sitely = (iy + ly_rel)%basis.Ly;
                site_l=sitelx + (sitely*basis.Lx);

                sitemx = (sitelx + mx_rel)%basis.Lx;
                sitemy = (sitely + my_rel)%basis.Ly;
                site_m=sitemx + (sitemy*basis.Lx);

                sitenx = (sitelx + nx_rel)%basis.Lx;
                siteny = (sitely + ny_rel)%basis.Ly;
                site_n=sitenx + (siteny*basis.Lx);


                //---------------Sz[site_i]Sp[site_j]Sm[site_k]Sz[site_l]Sp[site_m]Sm[site_n]-------------------//

                sign_pow_up=0;
                sign_pow_dn=0;

                //there have to be only one electron in site_l with up spin
                //there have to be only one electron in site_j with down spin

                //For Sm[site_n] : site_n must be in |up> state
                if((bit_value(basis.D_dn_basis[i],site_n)==0) && (bit_value(basis.D_up_basis[i],site_n)==1))
                {
                    D_up = (int) (basis.D_up_basis[i] - pow(2, site_n));
                    D_dn = (int) (basis.D_dn_basis[i] + pow(2, site_n));


                    //For Sp[site_m] : site_m must be in |dn> state
                    if((bit_value(D_dn,site_m)==1) && (bit_value(D_up,site_m)==0))
                    {
                        D_up += (int) (D_up + pow(2, site_m));
                        D_dn += (int) (D_dn - pow(2, site_m));

                        Sz_sitel = 0.5*(bit_value(D_up, site_l) - bit_value(D_dn, site_l) );

                        //For Sm[site_k] : site_k must be in |up> state
                        if((bit_value(D_dn,site_k)==0) && (bit_value(D_up,site_k)==1))
                        {
                            D_up += (int) (D_up - pow(2, site_k));
                            D_dn += (int) (D_dn + pow(2, site_k));

                            //For Sp[site_j] : site_j must be in |dn> state
                            if((bit_value(D_dn,site_j)==1) && (bit_value(D_up,site_j)==0))
                            {
                                D_up += (int) (D_up + pow(2, site_j));
                                D_dn += (int) (D_dn - pow(2, site_j));

                                Sz_sitei = 0.5*(bit_value(D_up, site_i) - bit_value(D_dn, site_i) );

                                Sz_value = Sz_sitei*Sz_sitel;

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Inv Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Inv Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Inv Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Inv Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }


                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if(row_found_){break;}

                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                    //these Fermionic signs does not appear because "S" operator has 2 fermionic operators for same site
                                    //sign_pow_up += one_bits_in_bw(site_j,site_l,basis.D_up_basis[i]);
                                    //sign_pow_dn += one_bits_in_bw(site_j,site_l,basis.D_dn_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                                    repeating_rows=false;
                                    check_min=OPR_.rows.size()-1;
                                    check_max=(OPR_.rows.size()-1)-row_counter;
                                    // cout<<check_min<<endl;
                                    // cout<<check_max<<endl;
                                    for(int check_=check_min;check_>check_max;check_--){
                                        if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                            OPR_.value[check_] +=(-4.0)*sign_FM*Sz_value*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                            repeating_rows=true;
                                            break;
                                        }
                                    }

                                    if(!repeating_rows){
                                        OPR_.value.push_back((-4.0)*sign_FM*one*Sz_value*phase_*(1.0/(basis.Lx*basis.Ly)));
                                        OPR_.rows.push_back(m_new);
                                        OPR_.columns.push_back(m);
                                        row_counter++;
                                    }

                                }
                            }
                        }
                    }
                }





                //---------------(-1.0)*Sz[site_i]Sm[site_j]Sp[site_k]Sz[site_l]Sp[site_m]Sm[site_n]-------------------//

                sign_pow_up=0;
                sign_pow_dn=0;
                Value_offset=-1.0;

                //there have to be only one electron in site_l with up spin
                //there have to be only one electron in site_j with down spin

                //For Sm[site_n] : site_n must be in |up> state
                if((bit_value(basis.D_dn_basis[i],site_n)==0) && (bit_value(basis.D_up_basis[i],site_n)==1))
                {
                    D_up = (int) (basis.D_up_basis[i] - pow(2, site_n));
                    D_dn = (int) (basis.D_dn_basis[i] + pow(2, site_n));


                    //For Sp[site_m] : site_m must be in |dn> state
                    if((bit_value(D_dn,site_m)==1) && (bit_value(D_up,site_m)==0))
                    {
                        D_up += (int) (D_up + pow(2, site_m));
                        D_dn += (int) (D_dn - pow(2, site_m));

                        Sz_sitel = 0.5*(bit_value(D_up, site_l) - bit_value(D_dn, site_l) );

                        //For Sp[site_k] : site_k must be in |dn> state
                        if((bit_value(D_dn,site_k)==1) && (bit_value(D_up,site_k)==0))
                        {
                            D_up += (int) (D_up + pow(2, site_k));
                            D_dn += (int) (D_dn - pow(2, site_k));

                            //For Sm[site_j] : site_j must be in |up> state
                            if((bit_value(D_dn,site_j)==0) && (bit_value(D_up,site_j)==1))
                            {
                                D_up += (int) (D_up - pow(2, site_j));
                                D_dn += (int) (D_dn + pow(2, site_j));

                                Sz_sitei = 0.5*(bit_value(D_up, site_i) - bit_value(D_dn, site_i) );

                                Sz_value = Sz_sitei*Sz_sitel;

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Inv Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Inv Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Inv Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Inv Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }


                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if(row_found_){break;}

                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                    //these Fermionic signs does not appear because "S" operator has 2 fermionic operators for same site
                                    //sign_pow_up += one_bits_in_bw(site_j,site_l,basis.D_up_basis[i]);
                                    //sign_pow_dn += one_bits_in_bw(site_j,site_l,basis.D_dn_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                                    repeating_rows=false;
                                    check_min=OPR_.rows.size()-1;
                                    check_max=(OPR_.rows.size()-1)-row_counter;
                                    // cout<<check_min<<endl;
                                    // cout<<check_max<<endl;
                                    for(int check_=check_min;check_>check_max;check_--){
                                        if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                            OPR_.value[check_] +=(-4.0*Value_offset)*sign_FM*Sz_value*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                            repeating_rows=true;
                                            break;
                                        }
                                    }

                                    if(!repeating_rows){
                                        OPR_.value.push_back((-4.0*Value_offset)*sign_FM*one*Sz_value*phase_*(1.0/(basis.Lx*basis.Ly)));
                                        OPR_.rows.push_back(m_new);
                                        OPR_.columns.push_back(m);
                                        row_counter++;
                                    }

                                }
                            }
                        }
                    }
                }







                //---------------(-1.0)Sz[site_i]Sp[site_j]Sm[site_k]Sz[site_l]Sm[site_m]Sp[site_n]-------------------//

                sign_pow_up=0;
                sign_pow_dn=0;
                Value_offset=-1.0;

                //there have to be only one electron in site_l with up spin
                //there have to be only one electron in site_j with down spin

                //For Sp[site_n] : site_n must be in |dn> state
                if((bit_value(basis.D_dn_basis[i],site_n)==1) && (bit_value(basis.D_up_basis[i],site_n)==0))
                {
                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_n));
                    D_dn = (int) (basis.D_dn_basis[i] - pow(2, site_n));


                    //For Sm[site_m] : site_m must be in |up> state
                    if((bit_value(D_dn,site_m)==0) && (bit_value(D_up,site_m)==1))
                    {
                        D_up += (int) (D_up - pow(2, site_m));
                        D_dn += (int) (D_dn + pow(2, site_m));

                        Sz_sitel = 0.5*(bit_value(D_up, site_l) - bit_value(D_dn, site_l) );

                        //For Sm[site_k] : site_k must be in |up> state
                        if((bit_value(D_dn,site_k)==0) && (bit_value(D_up,site_k)==1))
                        {
                            D_up += (int) (D_up - pow(2, site_k));
                            D_dn += (int) (D_dn + pow(2, site_k));

                            //For Sp[site_j] : site_j must be in |dn> state
                            if((bit_value(D_dn,site_j)==1) && (bit_value(D_up,site_j)==0))
                            {
                                D_up += (int) (D_up + pow(2, site_j));
                                D_dn += (int) (D_dn - pow(2, site_j));

                                Sz_sitei = 0.5*(bit_value(D_up, site_i) - bit_value(D_dn, site_i) );

                                Sz_value = Sz_sitei*Sz_sitel;

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Inv Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Inv Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Inv Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Inv Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }


                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if(row_found_){break;}

                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                    //these Fermionic signs does not appear because "S" operator has 2 fermionic operators for same site
                                    //sign_pow_up += one_bits_in_bw(site_j,site_l,basis.D_up_basis[i]);
                                    //sign_pow_dn += one_bits_in_bw(site_j,site_l,basis.D_dn_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                                    repeating_rows=false;
                                    check_min=OPR_.rows.size()-1;
                                    check_max=(OPR_.rows.size()-1)-row_counter;
                                    // cout<<check_min<<endl;
                                    // cout<<check_max<<endl;
                                    for(int check_=check_min;check_>check_max;check_--){
                                        if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                            OPR_.value[check_] +=(-4.0*Value_offset)*sign_FM*Sz_value*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                            repeating_rows=true;
                                            break;
                                        }
                                    }

                                    if(!repeating_rows){
                                        OPR_.value.push_back((-4.0*Value_offset)*sign_FM*one*Sz_value*phase_*(1.0/(basis.Lx*basis.Ly)));
                                        OPR_.rows.push_back(m_new);
                                        OPR_.columns.push_back(m);
                                        row_counter++;
                                    }

                                }
                            }
                        }
                    }
                }






                //---------------(1.0)Sz[site_i]Sm[site_j]Sp[site_k]Sz[site_l]Sm[site_m]Sp[site_n]-------------------//

                sign_pow_up=0;
                sign_pow_dn=0;
                Value_offset=1.0;

                //there have to be only one electron in site_l with up spin
                //there have to be only one electron in site_j with down spin

                //For Sp[site_n] : site_n must be in |dn> state
                if((bit_value(basis.D_dn_basis[i],site_n)==1) && (bit_value(basis.D_up_basis[i],site_n)==0))
                {
                    D_up = (int) (basis.D_up_basis[i] + pow(2, site_n));
                    D_dn = (int) (basis.D_dn_basis[i] - pow(2, site_n));


                    //For Sm[site_m] : site_m must be in |up> state
                    if((bit_value(D_dn,site_m)==0) && (bit_value(D_up,site_m)==1))
                    {
                        D_up += (int) (D_up - pow(2, site_m));
                        D_dn += (int) (D_dn + pow(2, site_m));

                        Sz_sitel = 0.5*(bit_value(D_up, site_l) - bit_value(D_dn, site_l) );

                        //For Sp[site_k] : site_k must be in |dn> state
                        if((bit_value(D_dn,site_k)==1) && (bit_value(D_up,site_k)==0))
                        {
                            D_up += (int) (D_up + pow(2, site_k));
                            D_dn += (int) (D_dn - pow(2, site_k));

                            //For Sm[site_j] : site_j must be in |up> state
                            if((bit_value(D_dn,site_j)==0) && (bit_value(D_up,site_j)==1))
                            {
                                D_up += (int) (D_up - pow(2, site_j));
                                D_dn += (int) (D_dn + pow(2, site_j));

                                Sz_sitei = 0.5*(bit_value(D_up, site_i) - bit_value(D_dn, site_i) );

                                Sz_value = Sz_sitei*Sz_sitel;

                                D_up_temp=D_up;
                                D_dn_temp=D_dn;
                                row_found_=false;

                                for(int inv_trnsltns_x=0;inv_trnsltns_x<basis.Lx;inv_trnsltns_x++){
                                    if(inv_trnsltns_x>0 && basis.Lx>1){

                                        for(int iy_=0;iy_<basis.Ly;iy_++){

                                            //Inv Translation on spin_dn
                                            sign_pow_dn_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_dn_temp) +
                                                    1*bit_value(D_dn_temp,iy_*basis.Lx);
                                            if(bit_value(D_dn_temp, (iy_+1)*basis.Lx - 1 )==1){
                                                sign_pow_dn += 1*sign_pow_dn_orb0;
                                            }

                                            D_dn_temp = Act_Translation_2D_alongX_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, iy_);

                                            //Inv Translation on spin_up
                                            sign_pow_up_orb0 = one_bits_in_bw(iy_*basis.Lx, (iy_+1)*basis.Lx - 1, D_up_temp) +
                                                    1*bit_value(D_up_temp,iy_*basis.Lx);
                                            if(bit_value(D_up_temp, (iy_+1)*basis.Lx - 1)==1){
                                                sign_pow_up += 1*sign_pow_up_orb0;
                                            }

                                            D_up_temp = Act_Translation_2D_alongX_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, iy_);

                                        }
                                    }
                                    D_dn_temp_Xtrans = D_dn_temp;
                                    D_up_temp_Xtrans = D_up_temp;
                                    sign_pow_up_Xtrans = sign_pow_up;
                                    sign_pow_dn_Xtrans = sign_pow_dn;
                                    for(int inv_trnsltns_y=0;inv_trnsltns_y<basis.Ly;inv_trnsltns_y++){


                                        if(inv_trnsltns_y>0 && basis.Ly>1){

                                            for(int ix_=0;ix_<basis.Lx;ix_++){

                                                //Inv Translation on spin_dn
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_dn_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_dn_temp);
                                                    sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_dn += sign_pow_dn_orb0;
                                                }
                                                sign_pow_dn_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_dn_temp) +
                                                        bit_value(D_dn_temp, ix_) ;
                                                sign_pow_dn_orb0 = sign_pow_dn_orb0*bit_value(D_dn_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_dn += sign_pow_dn_orb0;

                                                D_dn_temp = Act_Translation_2D_alongY_assuming_PBC(D_dn_temp,basis.Lx, basis.Ly, ix_);



                                                //Inv Translation on spin_up
                                                for(int iy_=0;iy_<basis.Ly-1;iy_++){
                                                    sign_pow_up_orb0 = one_bits_in_bw(ix_ + iy_*basis.Lx, ix_ + (iy_+1)*basis.Lx, D_up_temp);
                                                    sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + iy_*basis.Lx);

                                                    sign_pow_up += sign_pow_up_orb0;
                                                }
                                                sign_pow_up_orb0 = one_bits_in_bw(ix_, ix_ + (basis.Ly-1)*basis.Lx, D_up_temp) +
                                                        bit_value(D_up_temp, ix_) ;
                                                sign_pow_up_orb0 = sign_pow_up_orb0*bit_value(D_up_temp, ix_ + (basis.Ly-1)*basis.Lx);
                                                sign_pow_up += sign_pow_up_orb0;

                                                D_up_temp = Act_Translation_2D_alongY_assuming_PBC(D_up_temp,basis.Lx, basis.Ly, ix_);

                                            }
                                        }


                                        if(D_up_temp>=basis.Dup_Range.size())
                                        {
                                            row_found_=false;
                                        }
                                        else
                                        {
                                            assert(D_up_temp<basis.Dup_Range.size());
                                            range_min=basis.Dup_Range[D_up_temp].first;
                                            range_max=basis.Dup_Range[D_up_temp].second;
                                            if(range_min==-1)
                                            {
                                                row_found_=false;
                                                assert(range_max==-1);
                                            }
                                            else
                                            {
                                                i_new = Find_int_in_part_of_intarray(D_dn_temp, basis.D_dn_basis, range_min, range_max);
                                                if(i_new==-1){
                                                    row_found_=false;
                                                }
                                                else{
                                                    row_found_=true;
                                                    Inv_Trnsltns_x_=inv_trnsltns_x;
                                                    Inv_Trnsltns_y_=inv_trnsltns_y;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if(row_found_){break;}

                                    D_up_temp=D_up_temp_Xtrans ;
                                    D_dn_temp=D_dn_temp_Xtrans ;
                                    sign_pow_up = sign_pow_up_Xtrans;
                                    sign_pow_dn = sign_pow_dn_Xtrans;
                                }


                                if(row_found_==true){
                                    m_new = i_new;

#ifdef USE_COMPLEX
                                    phase_=exp(-1.0*iota_*( ((2.0*PI*Inv_Trnsltns_x_*basis.Momentum_nx)/(basis.Lx)) + ((2.0*PI*Inv_Trnsltns_y_*basis.Momentum_ny)/(basis.Ly))   )
                                               )
                                            *sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif
#ifndef USE_COMPLEX
                                    if( !(basis.Momentum_nx==0 && basis.Momentum_ny==0) ){
                                        cout<<"ONLY Kx=0,Ky=0 is allowed in real space calculations"<<endl;
                                    }
                                    assert(basis.Momentum_nx==0 && basis.Momentum_ny==0);
                                    phase_=one*sqrt((1.0*basis.D_Norm[m_new])/(1.0*basis.D_Norm[m]));
#endif

                                    //these Fermionic signs does not appear because "S" operator has 2 fermionic operators for same site
                                    //sign_pow_up += one_bits_in_bw(site_j,site_l,basis.D_up_basis[i]);
                                    //sign_pow_dn += one_bits_in_bw(site_j,site_l,basis.D_dn_basis[i]);

                                    sign_FM = pow(-1.0, sign_pow_up+sign_pow_dn);


                                    repeating_rows=false;
                                    check_min=OPR_.rows.size()-1;
                                    check_max=(OPR_.rows.size()-1)-row_counter;
                                    // cout<<check_min<<endl;
                                    // cout<<check_max<<endl;
                                    for(int check_=check_min;check_>check_max;check_--){
                                        if(OPR_.rows[check_]==m_new && OPR_.columns[check_]==m){
                                            OPR_.value[check_] +=(-4.0*Value_offset)*sign_FM*Sz_value*one*phase_*(1.0/(basis.Lx*basis.Ly));
                                            repeating_rows=true;
                                            break;
                                        }
                                    }

                                    if(!repeating_rows){
                                        OPR_.value.push_back((-4.0*Value_offset)*sign_FM*one*Sz_value*phase_*(1.0/(basis.Lx*basis.Ly)));
                                        OPR_.rows.push_back(m_new);
                                        OPR_.columns.push_back(m);
                                        row_counter++;
                                    }

                                }
                            }
                        }
                    }
                }








            } //iy
        } // ix



    } // "i" i.e up_decimals




}

void MODEL_1_orb_Hubb_2D_KSector::Read_parameters(BASIS_1_orb_Hubb_2D_KSector &basis, string filename){


    string filepath = filename;

    bool read_onsite_energies=false;
    string Geometry = "Geometry = ";
    string pbc_,PBC_ ="PBC = ";
    string lx_, Lx_ = "Lx = ";
    string ly_, Ly_ = "Ly = ";
    string ndn, Ndn = "Ndown = ";
    string nup, Nup = "Nup = ";
    string ucoul, Ucoul = "U = ";
    string hmag, Hmag = "H_mag = ";
    string hopp_, Hopp_ = "Hopping_NN = ";
    string momentum_nx_, Momentum_nx_ = "Momentum_nx = ";
    string momentum_ny_, Momentum_ny_ = "Momentum_ny = ";
    string read_basis_bool_, Read_Basis_BOOL_ = "Read_Basis_bool = ";
    string basis_read_file_, Basis_Read_File_ = "Read_File_Basis = ";
    string write_basis_bool_, Write_Basis_BOOL_ = "Write_Basis_bool = ";
    string basis_write_file_, Basis_Write_File_ = "Write_File_Basis = ";

    string file_onsite_energies_, File_Onsite_Energies_ = "File_Onsite_Energies = ";
    string file_hopping_connections_, File_Hopping_Connections_ = "File_Hopping_Connections = ";
    string file_nonlocal_int_connections_, File_NonLocal_Int_Connections_ = "File_NonLocal_Int_Connections = ";

    string file_three_point_observation_, File_Three_Point_Observation_ = "File_Three_Point_Observation = ";

    string saving_hamil_, Saving_Hamil_ = "Saving_Hamiltonian_and_Oprs = ";
    string processors_, Processors_ = "Processors = ";



    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);


            if ((offset = line.find(Geometry, 0)) != string::npos) {
                geometry_ = line.substr (offset+Geometry.length());				}

            if ((offset = line.find(Read_Basis_BOOL_, 0)) != string::npos) {
                read_basis_bool_ = line.substr (offset+Read_Basis_BOOL_.length());				}


            if ((offset = line.find(Basis_Read_File_, 0)) != string::npos) {
                basis_read_file_ = line.substr (offset+Basis_Read_File_.length());				}


            if ((offset = line.find(Write_Basis_BOOL_, 0)) != string::npos) {
                write_basis_bool_ = line.substr (offset+Write_Basis_BOOL_.length());				}

            if ((offset = line.find(Basis_Write_File_, 0)) != string::npos) {
                basis_write_file_ = line.substr (offset+Basis_Write_File_.length());				}

            if ((offset = line.find(PBC_, 0)) != string::npos) {
                pbc_ = line.substr (offset+PBC_.length());				}

            if ((offset = line.find(Lx_, 0)) != string::npos) {
                lx_ = line.substr (offset + Lx_.length());		}

            if ((offset = line.find(Ly_, 0)) != string::npos) {
                ly_ = line.substr (offset + Ly_.length());		}

            if ((offset = line.find(Momentum_nx_, 0)) != string::npos) {
                momentum_nx_ = line.substr (offset + Momentum_nx_.length());		}

            if ((offset = line.find(Momentum_ny_, 0)) != string::npos) {
                momentum_ny_ = line.substr (offset + Momentum_ny_.length());		}

            if ((offset = line.find(Ndn, 0)) != string::npos) {
                ndn = line.substr (offset + Ndn.length());		}

            if ((offset = line.find(Nup, 0)) != string::npos) {
                nup= line.substr (offset + Nup.length());		}

            if ((offset = line.find(Ucoul, 0)) != string::npos) {
                ucoul= line.substr (offset + Ucoul.length());		}

            if ((offset = line.find(Hmag, 0)) != string::npos) {
                hmag = line.substr (offset + Hmag.length());		}

            if ((offset = line.find(Hopp_, 0)) != string::npos) {
                hopp_ = line.substr (offset+Hopp_.length());				}

            if ((offset = line.find(File_Three_Point_Observation_, 0)) != string::npos) {
                file_three_point_observation_ = line.substr (offset+File_Three_Point_Observation_.length());}

            if ((offset = line.find(File_Onsite_Energies_, 0)) != string::npos) {
                file_onsite_energies_ = line.substr (offset+File_Onsite_Energies_.length());				}

            if ((offset = line.find(File_Hopping_Connections_, 0)) != string::npos) {
                file_hopping_connections_ = line.substr (offset+File_Hopping_Connections_.length());				}

            if ((offset = line.find(File_NonLocal_Int_Connections_, 0)) != string::npos) {
                file_nonlocal_int_connections_ = line.substr (offset+File_NonLocal_Int_Connections_.length());				}

            if ((offset = line.find(Saving_Hamil_, 0)) != string::npos) {
                saving_hamil_ = line.substr (offset + Saving_Hamil_.length());		}

            if ((offset = line.find(Processors_ , 0)) != string::npos) {
                processors_ = line.substr (offset+Processors_ .length());	}


        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}



    int N_three_point_oprs;
    int n_sites_set;
    // stringstream _file_three_point_observation_(file_three_point_observation_);
    ifstream inputfile_three_point_observation(file_three_point_observation_.c_str());

    inputfile_three_point_observation>>N_three_point_oprs;
    three_point_oprs.resize(N_three_point_oprs);
    three_point_oprs_sites_set.resize(N_three_point_oprs);

    for(int n=0;n<three_point_oprs.size();n++){
        inputfile_three_point_observation>>three_point_oprs[n];
        inputfile_three_point_observation>> n_sites_set;
        three_point_oprs_sites_set[n].resize(n_sites_set);
        for(int m=0;m<n_sites_set;m++){
            three_point_oprs_sites_set[n][m].resize(3);
            for(int i=0;i<3;i++){
                inputfile_three_point_observation>>three_point_oprs_sites_set[n][m][i];
            }
        }
    }



    NProcessors_=atoi(processors_.c_str());

    if(pbc_ == "true"){
        PBC =true;
    }
    else{
        PBC=false;
    }

    if(read_basis_bool_=="true"){
        basis.Read_Basis=true;
    }
    else{
        basis.Read_Basis=false;
    }
    basis.file_read_basis = basis_read_file_;

    if(write_basis_bool_=="true"){
        basis.Write_Basis=true;
    }
    else{
        basis.Write_Basis=false;
    }
    basis.file_write_basis = basis_write_file_;



    if(saving_hamil_=="true"){
        Saving_Hamil=true;
    }
    else{
        Saving_Hamil=false;
    }


    basis.Lx=atoi(lx_.c_str());
    basis.Ly=atoi(ly_.c_str());
    basis.Length = basis.Lx*basis.Ly;
    basis.Ndn=atoi(ndn.c_str());
    basis.Nup=atoi(nup.c_str());
    basis.Momentum_nx=atoi(momentum_nx_.c_str());
    basis.Momentum_ny=atoi(momentum_ny_.c_str());


    U=atof(ucoul.c_str());
    Hopping_NN=one*atof(hopp_.c_str());


    double h;
    h=atof(hmag.c_str());
    H_field.resize(basis.Length);
    for(int i=0;i<basis.Length;i++){
        H_field[i]=h;
    }


    stringstream _file_onsite_energies_(file_onsite_energies_);
    _file_onsite_energies_ >> read_onsite_energies;
    Onsite_Energy.resize(basis.Length);
    string filename_Onsite_Energy;
    string line_temp;

    // string temp_x_, temp_y_, temp_site_, Ener_val_ ;
    int temp_x, temp_y, temp_site;
    double Ener_val;
    if(read_onsite_energies){
        _file_onsite_energies_ >> filename_Onsite_Energy;

        ifstream inputfile_Onsite_Energy(filename_Onsite_Energy.c_str());
        getline(inputfile_Onsite_Energy,line_temp);

        for(int iy=0;iy<basis.Ly;iy++){
            for(int ix=0;ix<basis.Lx;ix++){
                inputfile_Onsite_Energy >> temp_x >> temp_y >> temp_site >> Ener_val;
                assert(temp_x==ix);
                assert(temp_y==iy);
                assert(temp_site==iy*basis.Lx + ix);
                Onsite_Energy[temp_site]=Ener_val;
            }
        }

    }
    else{
        cout<< "Onsite energies are 0's"<<endl;
        for(int i=0;i<basis.Length;i++){
            Onsite_Energy[i]=0.0;
        }
    }




    if(geometry_=="LongRange"){

        cout<<"HOPPINGS AND INTERACTIONS MATRICES MUST HAVE TRANSLATIONAL SYMMETRY OF a(x) + a(y), where a is lattice constant"<<endl;

        //Hoppings Mat(i,j)ci^{dagger}cj
        Hopping_mat_NN.resize(basis.Length);
        for(int site_=0;site_<basis.Length;site_++){
            Hopping_mat_NN[site_].resize(basis.Length);
        }

        ifstream inputfile_hopping_connections(file_hopping_connections_.c_str());
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                inputfile_hopping_connections>>Hopping_mat_NN[site_i][site_j];
            }
        }

        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=site_i+1;site_j<basis.Length;site_j++){
                if(abs(Hopping_mat_NN[site_j][site_i]-Hopping_mat_NN[site_i][site_j]) >0.000001){
                    Hopping_mat_NN[site_j][site_i]=Hopping_mat_NN[site_i][site_j];
                }
            }
        }




        //Interactions  Mat(i,j)ninj
        NonLocalInteractions_mat.resize(basis.Length);
        for(int site_=0;site_<basis.Length;site_++){
            NonLocalInteractions_mat[site_].resize(basis.Length);
        }

        ifstream inputfile_nonlocal_int_connections(file_nonlocal_int_connections_.c_str());
        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                inputfile_nonlocal_int_connections>>NonLocalInteractions_mat[site_i][site_j];

                if(site_i==site_j){
                    assert(NonLocalInteractions_mat[site_i][site_j]==0.0);
                }
            }
        }

    }
    else{
        cout<<"Inbuild Hubbard model connections and interaction is used"<<endl;


        NonLocalInteractions_mat.resize(basis.Length);
        for(int site_=0;site_<basis.Length;site_++){
            NonLocalInteractions_mat[site_].resize(basis.Length);
        }

        for(int site_i=0;site_i<basis.Length;site_i++){
            for(int site_j=0;site_j<basis.Length;site_j++){
                NonLocalInteractions_mat[site_i][site_j]=0.0;
            }
        }



        //Hoppings Mat(i,j)ci^{dagger}cj
        Hopping_mat_NN.resize(basis.Length);
        for(int site_=0;site_<basis.Length;site_++){
            Hopping_mat_NN[site_].resize(basis.Length);
        }


        int site, site_p, Dis_x, Dis_y, Dis_sqr;
        for(int ix=0;ix<basis.Lx ;ix++){
            for(int iy=0;iy<basis.Ly ;iy++){
                site=ix + (iy*basis.Lx);

                for(int ix_p=0;ix_p<basis.Lx ;ix_p++){
                    for(int iy_p=0;iy_p<basis.Ly ;iy_p++){
                        site_p=ix_p + (iy_p*basis.Lx);


                        Dis_x = abs(ix_p - ix);
                        Dis_y = abs(iy_p - iy);
                        if(Dis_x==basis.Lx-1){
                            Dis_x = 1;
                        }
                        if(Dis_y==basis.Ly-1 && (basis.Ly!=1)){
                            Dis_y = 1;
                        }

                        Dis_sqr = Dis_x*Dis_x + Dis_y*Dis_y;

                        if(Dis_sqr==1){
                            Hopping_mat_NN[site][site_p]=Hopping_NN;

                        }
                    }
                }
            }
        }



    }


    cout<<"PRINTING HOPPING MATRIX"<<endl;
    Print_Matrix(Hopping_mat_NN);

    cout<<""<<endl;
    cout<<"PRINTING Interaction MATRIX"<<endl;
    Print_Matrix(NonLocalInteractions_mat);
    cout<<"**************************"<<endl;




}



void MODEL_1_orb_Hubb_2D_KSector::Read_parameters_for_dynamics(string filename){

    string dyn_momentum_x_, Dyn_Momentum_x_ = "kx_for_dynamics = ";
    string dyn_momentum_y_, Dyn_Momentum_y_ = "ky_for_dynamics = ";
    string dyn_momentum_resolved_, Dyn_Momentum_Resolved_ = "Momentum_resolved = ";
    string Dyn_opr_string_  = "Opr_for_Dynamics = ";
    string basis_read_file_Kminusq_, Basis_Read_File_Kminusq_ = "Read_File_Basis_Kminusq = ";



    int offset;
    string line;
    ifstream inputfile(filename.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(Basis_Read_File_Kminusq_, 0)) != string::npos) {
                basis_read_file_Kminusq_ = line.substr (offset+Basis_Read_File_Kminusq_.length());	}

            if ((offset = line.find(Dyn_Momentum_Resolved_, 0)) != string::npos) {
                dyn_momentum_resolved_ = line.substr (offset + Dyn_Momentum_Resolved_.length());		}

            if ((offset = line.find(Dyn_Momentum_x_, 0)) != string::npos) {
                dyn_momentum_x_ = line.substr (offset + Dyn_Momentum_x_.length());		}

            if ((offset = line.find(Dyn_Momentum_y_, 0)) != string::npos) {
                dyn_momentum_y_ = line.substr (offset + Dyn_Momentum_y_.length());		}

            if ((offset = line.find(Dyn_opr_string_, 0)) != string::npos) {
                Dyn_opr_string = line.substr (offset + Dyn_opr_string_.length());		}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the Model class."<<endl;}


    file_read_basis_Kminusq = basis_read_file_Kminusq_;

    Dyn_Momentum_x=atoi(dyn_momentum_x_.c_str());
    Dyn_Momentum_y=atoi(dyn_momentum_y_.c_str());

    if(dyn_momentum_resolved_=="true"){
        Dyn_Momentum_Resolved=true;
    }
    else{
        Dyn_Momentum_Resolved=false;
    }

}
//#endif

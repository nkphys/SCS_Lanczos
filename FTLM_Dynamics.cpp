#include "FTLM_Dynamics.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX



void FTLM_DYNAMICS::Perform_FTLM(string inp_filename, Matrix_COO& OPR_){


    LANCZOS Lanczos1_;LANCZOS Lanczos2_;

    double offset_E;

    Lanczos1_.Read_Lanczos_parameters(inp_filename);
    Lanczos1_.Dynamics_performed=false;
    Lanczos2_.Read_Lanczos_parameters(inp_filename);
    Lanczos2_.Dynamics_performed=false;


    Lanczos1_.Save_the_Seed=true;

    Lanczos2_.Get_SeedVec_from_another_routine=true;
    Lanczos2_.Get_insitu_FTLM_overlaps=true;


    if(!Lanczos1_.Get_Full_Spectrum){
        cout<<"Get Full Spectrum must be true"<<endl;
        assert(Lanczos1_.Get_Full_Spectrum);
    }
    if(!Lanczos1_.need_few_eig_vecs){
        cout<<"need_few_eig_vecs must be true"<<endl;
        assert(Lanczos1_.need_few_eig_vecs);
    }

    if(!Lanczos2_.Get_Full_Spectrum){
        cout<<"Get Full Spectrum must be true"<<endl;
        assert(Lanczos2_.Get_Full_Spectrum);
    }
    if(!Lanczos2_.need_few_eig_vecs){
        cout<<"need_few_eig_vecs must be true"<<endl;
        assert(Lanczos2_.need_few_eig_vecs);
    }




    Temperature = 10.0; //CHANGE THIS
    M_=Lanczos1_.M_FTLM;
    Total_Random_States=Lanczos1_.Total_Random_States_for_FTLM;
    Boltzman_const = 1.0;
    Beta = 1.0/(Boltzman_const*Temperature);


    //Save Eigenvectors only for Lanczos-1, Lanczos-2.
    // Lanczos-2 will delete the vectors with in-situ overlap calculations
    Lanczos1_.states_to_look.resize(M_);
    Lanczos2_.states_to_look.resize(M_);
    for(int i=0;i<M_;i++){
        Lanczos1_.states_to_look[i]=i;
        Lanczos2_.states_to_look[i]=i;
    }
    Lanczos2_.Mat_elements.resize(M_);
    for(int m=0;m<M_;m++){
        Lanczos2_.Mat_elements[m].resize(M_);
    }


    int Lanc_steps;

    double Conf_Partition_Func, Conf_Partition_Func_Sqr;
    double Partition_Func_Std_Dev;

    double omega;
    int omega_points;
    omega_points = (int) ((Lanczos1_.omega_max - Lanczos1_.omega_min)/(Lanczos1_.d_omega));
    Mat_1_doub Conf_Cw;
    Conf_Cw.resize(omega_points);



    Sum_Partition_Func=0.0;
    Sum_Partition_Func_Sqr=0.0;
    Sum_Cw.resize(omega_points);
    for(int point=0;point<omega_points;point++){
        Sum_Cw[point]=zero;
    }

    for(int run_no=0;run_no<Total_Random_States;run_no++){




        //LANCZOS RUN-1   ---------------------------------//
        Lanczos1_.Random_seed_value += run_no;
        cout<<"-------LANCZOS-1 PERFORMED FOR CONFIGURATION NO. "<<run_no<<" with seed = |"<<Lanczos1_.Random_seed_value<<">";
        cout<<"------------------"<<endl;

        Lanczos1_.Perform_LANCZOS(Hamil);
        Lanc_steps = Lanczos1_.Evals_Tri_all.size();
        if(run_no==0){
            offset_E = Lanczos1_.Evals_Tri_all[Lanc_steps-1][0];
        }
        Evals1 = Lanczos1_.Evals_Tri_all[Lanc_steps-1];

        //Creating Vecs for Opr|Eig_{i}(run-1)> and seed for Lanczos-2
        Lanczos2_.Vecs_FTLM.resize(M_);
        for(int i=0;i<M_;i++){
            Matrix_COO_vector_multiplication("FULL", OPR_, Lanczos1_.Eig_vecs[i], Lanczos2_.Vecs_FTLM[i]);
            vector < double_type >().swap(Lanczos1_.Eig_vecs[i]);
        }
        Matrix_COO_vector_multiplication("FULL", OPR_, Lanczos1_.Saved_Seed, Lanczos2_.Seed_used);


        cout<<"-------LANCZOS-2 PERFORMED FOR CONFIGURATION NO. "<<run_no<<" with seed = OPR|"<<Lanczos1_.Random_seed_value<<">";
        cout<<"------------------"<<endl;
        Lanczos2_.Perform_LANCZOS(Hamil);
        Evals2 = Lanczos2_.Evals_Tri_all[Lanc_steps-1];



        Conf_Partition_Func = 0.0;
        for(int j=0;j<M_;j++){
            Conf_Partition_Func += exp(-Beta*(Evals1[j] - offset_E))
                    *abs(Lanczos1_.red_eig_vecs[j][0])*abs(Lanczos1_.red_eig_vecs[j][0]);
        }

        Conf_Partition_Func_Sqr = Conf_Partition_Func*Conf_Partition_Func;
        Sum_Partition_Func += Conf_Partition_Func;
        Sum_Partition_Func_Sqr += Conf_Partition_Func_Sqr;
        Partition_Func_Std_Dev = (Sum_Partition_Func_Sqr*Hamil.ncols)/(run_no+1) - (pow((Sum_Partition_Func/(run_no+1)),2.0)*(Hamil.ncols));



        for(int point=0;point<omega_points;point++){
            omega = Lanczos1_.omega_min + (point*Lanczos1_.d_omega);

            Conf_Cw[point]=zero;
            for(int i=0;i<M_;i++){
                for(int j=0;j<M_;j++){
#ifdef USE_COMPLEX
                    Conf_Cw[point] += Lanczos1_.red_eig_vecs[i][0]*conj(Lanczos2_.red_eig_vecs[j][0])*
                            Lanczos2_.Mat_elements[i][j]*
                            Lorentzian(Lanczos1_.eta, (omega - Evals2[j] + Evals1[i]) );
#endif
#ifndef USE_COMPLEX
                    Conf_Cw[point] += Lanczos1_.red_eig_vecs[i][0]*(Lanczos2_.red_eig_vecs[j][0])*
                            Lanczos2_.Mat_elements[i][j]*
                            Lorentzian(Lanczos1_.eta, (omega - Evals2[j] + Evals1[i]) );
#endif
                }
            }
            Sum_Cw[point] += Conf_Cw[point];
        }




        if(run_no%10 ==0){

            char run_no_char[50];
            sprintf(run_no_char,"%d", run_no);
            string out_filerun = "FTLM_Dynamics_out_run" + string(run_no_char) + ".txt";
            ofstream filerun_out(out_filerun.c_str());
            filerun_out<<"#omega    Cw[averaged upto this conf] "<<endl;

            for(int point=0;point<omega_points;point++){
                omega = Lanczos1_.omega_min + (point*Lanczos1_.d_omega);

#ifdef USE_COMPLEX
                filerun_out<<omega<<"    "<<Sum_Cw[point].real()/Sum_Partition_Func<<"   "<<Sum_Cw[point].imag()/Sum_Partition_Func<<endl;
#endif
#ifndef USE_COMPLEX
                filerun_out<<omega<<"    "<<Sum_Cw[point]/Sum_Partition_Func<<endl;
#endif
            }

        }




        cout<<"Run_no = "<<run_no<<"   "<<"completed"<<endl;
        cout<<"Z*Runs/Nstates = "<< Sum_Partition_Func<<endl;

        Lanczos1_.Clear();
        Lanczos2_.Clear();
    }




}

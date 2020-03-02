#include "FTLM_Static.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX



void FTLM_STATIC::Perform_FTLM(string inp_filename){


    LANCZOS Lanczos_;

    double offset_E;

    Lanczos_.Read_Lanczos_parameters(inp_filename);

    Lanczos_.Dynamics_performed=false;

    if(!Lanczos_.Get_Full_Spectrum){
    cout<<"Get Full Spectrum must be true"<<endl;
    assert(Lanczos_.Get_Full_Spectrum);
    }
    if(!Lanczos_.need_few_eig_vecs){
    cout<<"need_few_eig_vecs must be true"<<endl;
    assert(Lanczos_.need_few_eig_vecs);
    }


    Temperature = Lanczos_.Temperature_FTLM;
    M_=Lanczos_.M_FTLM;
    Total_Random_States=Lanczos_.Total_Random_States_for_FTLM;


    Lanczos_.states_to_look.resize(M_);
    for(int i=0;i<M_;i++){
    Lanczos_.states_to_look[i]=i;
    }

    Boltzman_const = 1.0;
    Beta = 1.0/(Boltzman_const*Temperature);



    double Conf_Partition_Func, Conf_Partition_Func_Sqr;
    double Partition_Func_Std_Dev;

    double Conf_Hamil, Conf_Hamil_Sqr;
    double Hamil_Std_Dev_statistical;

    double Conf_Hamil2, Conf_Hamil2_Sqr;
    double Hamil2_Std_Dev_statistical;



    int Lanc_steps;

    Sum_Partition_Func=0.0;
    Sum_Partition_Func_Sqr=0.0;
    Sum_Hamil=0.0;
    Sum_Hamil_Sqr=0.0;
    Sum_Hamil2=0.0;
    Sum_Hamil2_Sqr=0.0;
    for(int run_no=0;run_no<Total_Random_States;run_no++){

    Lanczos_.Random_seed_value += run_no;


    cout<<"-------LANCZOS PREFORMED FOR CONFIGURATION NO. "<<run_no<<" with random seed = "<<Lanczos_.Random_seed_value;
    cout<<"------------------"<<endl;

    Lanczos_.Perform_LANCZOS(Hamil);
    Lanc_steps = Lanczos_.Evals_Tri_all.size();

    if(run_no==0){
        offset_E = Lanczos_.Evals_Tri_all[Lanc_steps-1][0];
    }



    Conf_Partition_Func = 0.0;
    Conf_Hamil =0.0;
    Conf_Hamil2 =0.0;
    for(int j=0;j<M_;j++){

        Conf_Hamil2 += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j]-offset_E))*
                        Lanczos_.Evals_Tri_all[Lanc_steps-1][j]*Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

        Conf_Hamil += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E ))*
                       Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

        Conf_Partition_Func += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E))
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);
    }

    Conf_Partition_Func_Sqr = Conf_Partition_Func*Conf_Partition_Func;
    Conf_Hamil_Sqr = Conf_Hamil*Conf_Hamil;
    Conf_Hamil2_Sqr = Conf_Hamil2*Conf_Hamil2;

    Sum_Partition_Func += Conf_Partition_Func;
    Sum_Partition_Func_Sqr += Conf_Partition_Func_Sqr;

    Sum_Hamil += Conf_Hamil;
    Sum_Hamil_Sqr += Conf_Hamil_Sqr;

    Sum_Hamil2 += Conf_Hamil2;
    Sum_Hamil2_Sqr += Conf_Hamil2_Sqr;


    Partition_Func_Std_Dev = (Sum_Partition_Func_Sqr*Hamil.ncols)/(run_no+1) - (pow((Sum_Partition_Func/(run_no+1)),2.0)*(Hamil.ncols));
    Hamil_Std_Dev_statistical = (Sum_Hamil_Sqr*Hamil.ncols)/(run_no+1) - (pow((Sum_Hamil/(run_no+1)),2.0)*(Hamil.ncols));
    Hamil2_Std_Dev_statistical = (Sum_Hamil2_Sqr*Hamil.ncols)/(run_no+1) - (pow((Sum_Hamil2/(run_no+1)),2.0)*(Hamil.ncols));


    Quantum_Avg_Hamil = Sum_Hamil / Sum_Partition_Func;
    Quantum_Avg_Hamil2 = Sum_Hamil2 / Sum_Partition_Func;

    cout<<"Run_no = "<<run_no<<"  "<<Quantum_Avg_Hamil<<"  "<<Quantum_Avg_Hamil2<<"   "<<Partition_Func_Std_Dev<<endl;

    Lanczos_.Clear();
    }




}

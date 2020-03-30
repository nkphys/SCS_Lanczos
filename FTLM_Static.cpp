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
    //Lanczos_.save_all_Krylov_space_vecs=true;

    Temperature_min = Lanczos_.Temprature_min_FTLM;
    Temperature_max = Lanczos_.Temprature_max_FTLM;
    delta_Temperature = Lanczos_.delta_Temperature_FTLM;
    assert(Temperature_max >= Temperature_min);
    assert(delta_Temperature != 0.0);

    N_Temperature_points = int((Temperature_max - Temperature_min)/(delta_Temperature) + 0.5 );

    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout <<"No. of Temperature points = "<<N_Temperature_points<<", min = "<<Temperature_min<<", max="<<Temperature_max<<endl;
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;



    M_=min(Hamil.nrows, Lanczos_.M_FTLM);
    Lanczos_.M_FTLM=M_;

    Total_Random_States=Lanczos_.Total_Random_States_for_FTLM;


    int M_temp = M_;
    Lanczos_.states_to_look.resize(M_temp);
    for(int i=0;i<M_temp;i++){
    Lanczos_.states_to_look[i]=i;
    }


    Lanczos_.Eig_vecs_required=false; //For : Only <H>, <H^2>

    Boltzman_const = 1.0;


    int Lanc_steps;

    Mat_1_real Conf_Partition_Func, Conf_Partition_Func_Sqr;
    Mat_1_real Partition_Func_Std_Dev;
    Mat_1_real Conf_Hamil, Conf_Hamil_Sqr;
    Mat_1_real Hamil_Std_Dev_statistical;
    Mat_1_real Conf_Hamil2, Conf_Hamil2_Sqr;
    Mat_1_real Hamil2_Std_Dev_statistical;

    Quantum_Avg_Hamil.resize(N_Temperature_points); Quantum_Avg_Hamil2.resize(N_Temperature_points);
    Conf_Partition_Func.resize(N_Temperature_points); Conf_Partition_Func_Sqr.resize(N_Temperature_points);
    Partition_Func_Std_Dev.resize(N_Temperature_points);
    Conf_Hamil.resize(N_Temperature_points); Conf_Hamil_Sqr.resize(N_Temperature_points);
    Hamil_Std_Dev_statistical.resize(N_Temperature_points);
    Conf_Hamil2.resize(N_Temperature_points); Conf_Hamil2_Sqr.resize(N_Temperature_points);
    Hamil2_Std_Dev_statistical.resize(N_Temperature_points);
    Sum_Partition_Func.resize(N_Temperature_points);Sum_Partition_Func_Sqr.resize(N_Temperature_points);
    Sum_Hamil.resize(N_Temperature_points);Sum_Hamil_Sqr.resize(N_Temperature_points);
    Sum_Hamil2.resize(N_Temperature_points);Sum_Hamil2_Sqr.resize(N_Temperature_points);

    for(int t=0;t<N_Temperature_points;t++){
    Sum_Partition_Func[t]=0.0;Sum_Partition_Func_Sqr[t]=0.0;
    Sum_Hamil[t]=0.0; Sum_Hamil_Sqr[t]=0.0;
    Sum_Hamil2[t]=0.0;  Sum_Hamil2_Sqr[t]=0.0;

    Conf_Partition_Func_Sqr[t]=0.0;
    Partition_Func_Std_Dev[t]=0.0;
    Conf_Hamil_Sqr[t]=0.0;
    Hamil_Std_Dev_statistical[t]=0.0;
    Conf_Hamil2_Sqr[t]=0.0;
    Hamil2_Std_Dev_statistical[t]=0.0;

    }

    for(int run_no=0;run_no<Total_Random_States;run_no++){

    Lanczos_.Random_seed_value += run_no;

    cout<<"-------LANCZOS PREFORMED FOR CONFIGURATION NO. "<<run_no<<" with random seed = "<<Lanczos_.Random_seed_value;
    cout<<"------------------"<<endl;

    Lanczos_.Perform_LANCZOS(Hamil);
    Lanc_steps = Lanczos_.Evals_Tri_all.size();
    if(run_no==0){
        offset_E = Lanczos_.Evals_Tri_all[Lanc_steps-1][0];
    }



    for(int Temp_point=0;Temp_point<N_Temperature_points;Temp_point++){
    Temperature = Temperature_min + Temp_point*(delta_Temperature);
    Beta = 1.0/(Boltzman_const*Temperature);


    Conf_Partition_Func[Temp_point] = 0.0;
    Conf_Hamil[Temp_point] =0.0;
    Conf_Hamil2[Temp_point] =0.0;

    for(int j=0;j<M_;j++){

        Conf_Hamil2[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j]-offset_E))*
                        Lanczos_.Evals_Tri_all[Lanc_steps-1][j]*Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

        Conf_Hamil[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E ))*
                       Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

        Conf_Partition_Func[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E))
                              *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);
    }

    Conf_Partition_Func_Sqr[Temp_point] = Conf_Partition_Func[Temp_point]*Conf_Partition_Func[Temp_point];
    Conf_Hamil_Sqr[Temp_point] = Conf_Hamil[Temp_point]*Conf_Hamil[Temp_point];
    Conf_Hamil2_Sqr[Temp_point] = Conf_Hamil2[Temp_point]*Conf_Hamil2[Temp_point];

    Sum_Partition_Func[Temp_point] += Conf_Partition_Func[Temp_point];
    Sum_Partition_Func_Sqr[Temp_point] += Conf_Partition_Func_Sqr[Temp_point];

    Sum_Hamil[Temp_point] += Conf_Hamil[Temp_point];
    Sum_Hamil_Sqr[Temp_point] += Conf_Hamil_Sqr[Temp_point];

    Sum_Hamil2[Temp_point] += Conf_Hamil2[Temp_point];
    Sum_Hamil2_Sqr[Temp_point] += Conf_Hamil2_Sqr[Temp_point];


    Partition_Func_Std_Dev[Temp_point] = (Sum_Partition_Func_Sqr[Temp_point]*Hamil.ncols)/(run_no+1) - (pow((Sum_Partition_Func[Temp_point]/(run_no+1)),2.0)*(Hamil.ncols));
    Hamil_Std_Dev_statistical[Temp_point] = (Sum_Hamil_Sqr[Temp_point]*Hamil.ncols)/(run_no+1) - (pow((Sum_Hamil[Temp_point]/(run_no+1)),2.0)*(Hamil.ncols));
    Hamil2_Std_Dev_statistical[Temp_point] = (Sum_Hamil2_Sqr[Temp_point]*Hamil.ncols)/(run_no+1) - (pow((Sum_Hamil2[Temp_point]/(run_no+1)),2.0)*(Hamil.ncols));


    Quantum_Avg_Hamil[Temp_point] = Sum_Hamil[Temp_point] / Sum_Partition_Func[Temp_point];
    Quantum_Avg_Hamil2[Temp_point] = Sum_Hamil2[Temp_point] / Sum_Partition_Func[Temp_point];

    cout<<"Run_no = "<<run_no<<"  "<<"Temperature = "<<Temperature<<"   "<<Quantum_Avg_Hamil[Temp_point]<<"  "<<Quantum_Avg_Hamil2[Temp_point]<<"   "<<Partition_Func_Std_Dev[Temp_point]<<endl;

    }


    Lanczos_.Clear();


    }




}

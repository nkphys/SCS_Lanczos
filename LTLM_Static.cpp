#include "LTLM_Static.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX



void LTLM_STATIC::Perform_LTLM(string inp_filename, Hamiltonian_1_COO& OPR_){


    LANCZOS Lanczos_;

    int No_of_oprts;
    No_of_oprts=OPR_.size();

    double offset_E;
    Mat_1_doub Opr_val;
    Opr_val.resize(No_of_oprts);

    Mat_2_doub Vec_Temp;
    Vec_Temp.resize(No_of_oprts);

    Sum_Opr_val.resize(No_of_oprts);

    Lanczos_.Read_Lanczos_parameters(inp_filename);
    Lanczos_.Save_the_Seed=true;

    Lanczos_.Dynamics_performed=false;

    if(!Lanczos_.Get_Full_Spectrum){
        cout<<"Get Full Spectrum must be true"<<endl;
        assert(Lanczos_.Get_Full_Spectrum);
    }
    if(!Lanczos_.need_few_eig_vecs){
        cout<<"need_few_eig_vecs must be true"<<endl;
        assert(Lanczos_.need_few_eig_vecs);
    }
    Lanczos_.save_all_Krylov_space_vecs=true;

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


    Lanczos_.Eig_vecs_required=true; //For False if Only <H>, <H^2> required

    Boltzman_const = 1.0;


    int Lanc_steps;

    Mat_1_doub Conf_Partition_Func;
    Mat_1_doub Conf_Hamil;
    Mat_1_doub Conf_Hamil2;
    Mat_2_doub Conf_Opr_val;Conf_Opr_val.resize(No_of_oprts);

    Quantum_Avg_Hamil.resize(N_Temperature_points);
    Quantum_Avg_Hamil2.resize(N_Temperature_points);
    Conf_Partition_Func.resize(N_Temperature_points);
    Conf_Hamil.resize(N_Temperature_points);
    Conf_Hamil2.resize(N_Temperature_points);
    for(int n=0;n<No_of_oprts;n++){
    Conf_Opr_val[n].resize(N_Temperature_points);
    }

    Sum_Partition_Func.resize(N_Temperature_points);
    Sum_Hamil.resize(N_Temperature_points);
    Sum_Hamil2.resize(N_Temperature_points);
    for(int n=0;n<No_of_oprts;n++){
    Sum_Opr_val[n].resize(N_Temperature_points);
    }

    for(int t=0;t<N_Temperature_points;t++){
        Sum_Partition_Func[t]=zero;
        Sum_Hamil[t]=zero;
        Sum_Hamil2[t]=zero;
        for(int n=0;n<No_of_oprts;n++){
        Sum_Opr_val[n][t]=zero;
        }

    }


    offset_E = Lanczos_.Energy_Offset_FTLM;
    for(int run_no=0;run_no<Total_Random_States;run_no++){

        Lanczos_.Random_seed_value += run_no+10;

        cout<<"-------LANCZOS PREFORMED FOR CONFIGURATION NO. "<<run_no<<" with random seed = "<<Lanczos_.Random_seed_value;
        cout<<"------------------"<<endl;

        Lanczos_.Perform_LANCZOS(Hamil);
        Lanc_steps = Lanczos_.Evals_Tri_all.size();

        for(int Temp_point=0;Temp_point<N_Temperature_points;Temp_point++){
            Conf_Partition_Func[Temp_point] = zero;
            Conf_Hamil[Temp_point] =zero;
            Conf_Hamil2[Temp_point] =zero;
            for(int n=0;n<No_of_oprts;n++){
            Conf_Opr_val[n][Temp_point] =zero;
            }
        }


        for(int l=0;l<M_;l++){

            for(int n=0;n<No_of_oprts;n++){
            Matrix_COO_vector_multiplication("FULL", OPR_[n], Lanczos_.Eig_vecs[l], Vec_Temp[n]);
            }

            for(int j=0;j<M_;j++){

                for(int n=0;n<No_of_oprts;n++){
                Opr_val[n] = dot_product(Vec_Temp[n], Lanczos_.Eig_vecs[j]);
                }


                for(int Temp_point=0;Temp_point<N_Temperature_points;Temp_point++){
                    Temperature = Temperature_min + Temp_point*(delta_Temperature);
                    Beta = 1.0/(Boltzman_const*Temperature);


                    for(int n=0;n<No_of_oprts;n++){
                    Conf_Opr_val[n][Temp_point] += exp(-Beta*( (Lanczos_.Evals_Tri_all[Lanc_steps-1][j] + Lanczos_.Evals_Tri_all[Lanc_steps-1][l])*0.5 - offset_E ))*
                            Opr_val[n]
                            *conj(Lanczos_.red_eig_vecs[l][0])*(Lanczos_.red_eig_vecs[j][0]);
                    }

                    Conf_Hamil2[Temp_point] += exp(-Beta*( (Lanczos_.Evals_Tri_all[Lanc_steps-1][j] + Lanczos_.Evals_Tri_all[Lanc_steps-1][l])*0.5 - offset_E ))*
                            Lanczos_.Evals_Tri_all[Lanc_steps-1][j]*Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                            *conj(Lanczos_.red_eig_vecs[l][0])*(Lanczos_.red_eig_vecs[j][0]);

//                    Conf_Hamil[Temp_point] += one*exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E ))*
//                            Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
//                            *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

                    Conf_Hamil[Temp_point] += exp(-Beta*( (Lanczos_.Evals_Tri_all[Lanc_steps-1][j] + Lanczos_.Evals_Tri_all[Lanc_steps-1][l])*0.5 - offset_E ))*
                            Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                            *conj(Lanczos_.red_eig_vecs[l][0])*(Lanczos_.red_eig_vecs[j][0]);

                    Conf_Partition_Func[Temp_point] += exp(-Beta*( (Lanczos_.Evals_Tri_all[Lanc_steps-1][j] + Lanczos_.Evals_Tri_all[Lanc_steps-1][l])*0.5 - offset_E ))
                            *conj(Lanczos_.red_eig_vecs[l][0])*(Lanczos_.red_eig_vecs[j][0]);

                }
            }
        }

        for(int Temp_point=0;Temp_point<N_Temperature_points;Temp_point++){

            Temperature = Temperature_min + Temp_point*(delta_Temperature);

            Conf_Hamil2[Temp_point]=Conf_Hamil2[Temp_point]/(1.0);
            Conf_Hamil[Temp_point]=Conf_Hamil[Temp_point]/(1.0);
            Conf_Partition_Func[Temp_point]=Conf_Partition_Func[Temp_point]/(1.0);

            Sum_Partition_Func[Temp_point] += Conf_Partition_Func[Temp_point];

            for(int n=0;n<No_of_oprts;n++){
            Sum_Opr_val[n][Temp_point] +=Conf_Opr_val[n][Temp_point];
            }

            Sum_Hamil[Temp_point] += Conf_Hamil[Temp_point];

            Sum_Hamil2[Temp_point] += Conf_Hamil2[Temp_point];

            Quantum_Avg_Hamil[Temp_point] = Sum_Hamil[Temp_point] / Sum_Partition_Func[Temp_point];
            Quantum_Avg_Hamil2[Temp_point] = Sum_Hamil2[Temp_point] / Sum_Partition_Func[Temp_point];

            cout<<"Run_no = "<<run_no<<"  "<<"Temperature = "<<Temperature<<"   "<<Sum_Hamil[Temp_point].real()<<"   "<<Sum_Hamil[Temp_point].imag()<<"  "<<Sum_Hamil2[Temp_point].real()<<"   "<<Sum_Hamil2[Temp_point].imag()<<"   "<<Sum_Partition_Func[Temp_point].real()<<"   "<<Sum_Partition_Func[Temp_point].imag()<<"   ";

            for(int n=0;n<No_of_oprts;n++){
              cout<<Sum_Opr_val[n][Temp_point].real()<<"   "<<Sum_Opr_val[n][Temp_point].imag()<<"   ";
            }
            cout<<endl;
        }


    }


    Lanczos_.Clear();


}






#include "FTLM_Static.h"
#include <stdlib.h>

#define PI 3.14159265
using namespace std;
//#define USE_COMPLEX
//#ifdef USE_COMPLEX



#ifndef FTLM_STATIC_functions
#define FTLM_STATIC_functions

template <typename Basis_type, typename Model_type>
void FTLM_STATIC<Basis_type, Model_type>::Perform_FTLM(string inp_filename, Hamiltonian_1_COO& OPR_){


    LANCZOS<Basis_type, Model_type> Lanczos_(basis,model);

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

    Mat_1_real Conf_Partition_Func;
    Mat_1_real Conf_Hamil;
    Mat_2_doub Conf_Opr_val;Conf_Opr_val.resize(No_of_oprts);
    Mat_1_real Conf_Hamil2;

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
        Sum_Partition_Func[t]=0.0;
        Sum_Hamil[t]=0.0;
        Sum_Hamil2[t]=0.0;
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
            Conf_Partition_Func[Temp_point] = 0.0;
            Conf_Hamil[Temp_point] =0.0;
            Conf_Hamil2[Temp_point] =0.0;
            for(int n=0;n<No_of_oprts;n++){
                Conf_Opr_val[n][Temp_point] =zero;
            }
        }


        for(int n=0;n<No_of_oprts;n++){
            Matrix_COO_vector_multiplication("FULL", OPR_[n], Lanczos_.Saved_Seed, Vec_Temp[n]);
        }

        for(int Temp_point=0;Temp_point<N_Temperature_points;Temp_point++){
            Temperature = Temperature_min + Temp_point*(delta_Temperature);
            Beta = 1.0/(Boltzman_const*Temperature);


            for(int j=0;j<M_;j++){

                for(int n=0;n<No_of_oprts;n++){
                    Opr_val[n] = dot_product(Vec_Temp[n], Lanczos_.Eig_vecs[j]);

                    Conf_Opr_val[n][Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E ))*
                            Opr_val[n]
                            *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);
                }

                Conf_Hamil2[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j]-offset_E))*
                        Lanczos_.Evals_Tri_all[Lanc_steps-1][j]*Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                        *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

                Conf_Hamil[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E ))*
                        Lanczos_.Evals_Tri_all[Lanc_steps-1][j]
                        *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);

                Conf_Partition_Func[Temp_point] += exp(-Beta*(Lanczos_.Evals_Tri_all[Lanc_steps-1][j] - offset_E))
                        *abs(Lanczos_.red_eig_vecs[j][0])*abs(Lanczos_.red_eig_vecs[j][0]);
            }



            Sum_Partition_Func[Temp_point] += Conf_Partition_Func[Temp_point];

            for(int n=0;n<No_of_oprts;n++){
                Sum_Opr_val[n][Temp_point] +=Conf_Opr_val[n][Temp_point];
            }

            Sum_Hamil[Temp_point] += Conf_Hamil[Temp_point];

            Sum_Hamil2[Temp_point] += Conf_Hamil2[Temp_point];

            Quantum_Avg_Hamil[Temp_point] = Sum_Hamil[Temp_point] / Sum_Partition_Func[Temp_point];
            Quantum_Avg_Hamil2[Temp_point] = Sum_Hamil2[Temp_point] / Sum_Partition_Func[Temp_point];

            //cout<<"Run_no = "<<run_no<<"  "<<"Temperature = "<<Temperature<<"   "<<Sum_Hamil[Temp_point]<<"  "<<Sum_Hamil2[Temp_point]<<"   "<<Sum_Opr_val[Temp_point].real()<<"   "<<Sum_Opr_val[Temp_point].imag()<<"   "<<Sum_Partition_Func[Temp_point]<<endl;

            cout<<"Run_no = "<<run_no<<"  "<<"Temperature = "<<Temperature<<"   "<<Sum_Hamil[Temp_point]<<"  "<<Sum_Hamil2[Temp_point]<<"   "<<Sum_Partition_Func[Temp_point]<<"   ";

            for(int n=0;n<No_of_oprts;n++){
#ifdef USE_COMPLEX
                cout<<Sum_Opr_val[n][Temp_point].real()<<"   "<<Sum_Opr_val[n][Temp_point].imag()<<"   ";
#endif
#ifndef USE_COMPLEX
                cout<<Sum_Opr_val[n][Temp_point]<<"   ";
#endif

            }
            cout<<endl;

        }

        Lanczos_.Clear();

    }




}

#endif

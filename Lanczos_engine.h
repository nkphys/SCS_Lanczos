#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>

#ifdef USE_COMPLEX
#include "functions_complex.h"
#else
#include "functions_real.h"
#endif


using namespace std;

#ifndef Lanczos_engine
#define Lanczos_engine

class LANCZOS{

public:
    int max_steps;
    double Lanc_Error;
    Mat_2_doub Krylov_space_vecs;
    bool save_all_Krylov_space_vecs;
    string file_eval_spectrum, overlap_out_file;
    Mat_2_real Evals_Tri_all;
    bool Get_Full_Spectrum;
    int Lanc_iter_done;
    bool need_few_eig_vecs,get_overlap_bw_eigstates_and_K_vecs,get_overlap_with_basis;
    bool lanczos_reorthogonalization;
    int few_;
    double GS_energy;
    Mat_2_doub red_eig_vecs;
    Mat_1_doub Eig_vec;
    Mat_2_doub Eig_vecs;
    Mat_2_pair_realint Overlaps;
    Mat_2_pair_Doubint Overlaps_bare;
    Mat_1_doub Dynamics_seed;
    bool Dynamics_performed;
    double check_omega;
    double eta;
    double omega_max;
    double d_omega;
    double omega_min;
    string file_dynamics_out;
    Mat_1_int states_to_look;
    bool calculate_local_obs_for_states_to_look;
    int no_basis_to_check;
    string file_Loc_obs_in_basis_of_states;
    double Numerator_Dynamics;
    bool Read_the_seed;
    string seed_file_in;
    bool Check_Ghosts;
    int Random_seed_value;
    Mat_2_doub One_point_observables_values;
    Mat_1_real B2 ,Norms;
    Mat_1_doub A, red_eig_vec;

    int time_evolution_order;
    Mat_1_doub Evolving_State;
    double omega_sign;

    int Total_Random_States_for_FTLM;
    double delta_Temperature_FTLM, Temprature_min_FTLM, Temprature_max_FTLM ;
    int M_FTLM;

    bool Get_insitu_FTLM_overlaps;
    bool Save_the_Seed;
    bool Get_SeedVec_from_another_routine;
    Mat_2_doub Vecs_FTLM;
    Mat_1_doub Saved_Seed;
    Mat_1_doub Seed_used;
    Mat_2_doub Mat_elements;





void Perform_LANCZOS(Matrix_COO &Hamil);
void Get_Dynamics_seed(Matrix_COO &Opr);
void Get_Dynamics_seed(Mat_1_doub vec_);
void Write_full_spectrum();
void Measure_macro_observables(Mat_1_string macro_obs, Hamiltonian_1_COO &Macro_oprts, int state_no);
void Measure_one_point_observables(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no);
double_type Measure_observable(Matrix_COO &OPR_, int state_no);
void Measure_two_point_observables(Mat_1_string two_point_obs, Hamiltonian_3_COO & Two_point_oprts, int T_sites, Mat_1_doub vec_, bool PBC_check);
void Measure_two_point_observables(Mat_1_string two_point_obs, Hamiltonian_3_COO &Two_point_oprts, int T_sites, int state_no, bool PBC_check);
void Measure_two_point_observables_smartly(Mat_1_string one_point_obs,Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no, string _model_);
void Measure_four_point_observables(Hamiltonian_3_COO &Two_point_oprts, Mat_1_tetra_int sites_set, int state_no);
void Measure_KE(Matrix_COO &Macro_oprt, int state_no);
void Measure_Total_Energy(Matrix_COO &Macro_oprt, int state_no);
void Read_Lanczos_parameters(string filename);
void Get_Dynamics_seed(string filename);
void Time_evolution_type1(Matrix_COO &Hamil, double dt_, double &Energy_);
void Clear();

};

#endif

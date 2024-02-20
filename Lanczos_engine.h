#include <iostream>
#include <math.h>  //fabs(double x) =|x|
#include <algorithm>
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
#include "models/Model_1_orb_Hubbard_chain.h"
#include "basis/Basis_1_orb_Hubbard_chain.h"
#include "models/Model_3_orb_Hubbard_chain.h"
#include "models/Model_2_orb_Hubbard_chain.h"
#include "models/Model_2_orb_Hubbard_chain_KSector.h"
#include "models/Model_3_orb_Hubbard_chain_two_SzSectors.h"
#include "basis/Basis_3_orb_Hubbard_chain.h"
#include "basis/Basis_2_orb_Hubbard_chain.h"
#include "basis/Basis_2_orb_Hubbard_chain_KSector.h"
#include "basis/Basis_1_orb_Hubbard_2D_KSector.h"
#include "models/Model_1_orb_Hubbard_2D_KSector.h"
#include "basis/Basis_3_orb_Hubbard_chain_two_SzSectors.h"
#include "models/Model_1_orb_Hubbard_chain.h"
#include "models/Model_1_orb_tJ.h"
#include "basis/Basis_1_orb_Hubbard_chain.h"
#include "basis/Basis_1_orb_Hubbard_GC.h"
#include "basis/Basis_1_orb_tJ.h"
#include "basis/Basis_Spins.h"
#include "models/Model_Spins.h"
#include "basis/Basis_KondoModel.h"
#include "models/Model_KondoModel.h"
#include "basis/Basis_Spins_Target_Sz.h"
#include "models/Model_Spins_Target_Sz.h"
#include "models/Model_3_orb_Hubbard_chain_GC.h"
#include "models/Model_multi_orb_Hubbard_chain_GC.h"
#include "models/Model_1_orb_Hubbard_GC.h"
#include "basis/Basis_3_orb_Hubbard_chain_GC.h"
#include "basis/Basis_3_orb_Hubbard_chain_GC_restricted.h"
#include "basis/Basis_multi_orb_Hubbard_chain_GC.h"
#include "basis/Basis_multi_orb_Hubbard_chain_GC_restricted.h"
#include "basis/Basis_SpinlessFermionsFockSpace.h"
#include "models/Model_SpinlessFermionsFockSpace.h"

//Remember "cpp" files for templated class over basis need to be included in this code
#include "models/Model_3_orb_Hubbard_chain_GC.cpp"
#include "models/Model_multi_orb_Hubbard_chain_GC.cpp"
#include "models/Model_1_orb_Hubbard_GC.cpp"


#ifdef USE_COMPLEX
#include "functions_complex.h"
#else
#include "functions_real.h"
#endif


using namespace std;

#ifndef Lanczos_engine
#define Lanczos_engine

template <typename Basis_type, typename Model_type>
class LANCZOS{

public:
        LANCZOS(Basis_type& Basis_, Model_type& Model_)
            :basis(Basis_), model(Model_)
        {

        }
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
    double Energy_Offset_FTLM;
    double Temperature_LTLM_Dynamics;

    bool Get_insitu_FTLM_overlaps;
    bool Save_the_Seed;
    bool Get_SeedVec_from_another_routine;
    Mat_2_doub Vecs_FTLM;
    Mat_1_doub Saved_Seed;
    Mat_1_doub Seed_used;
    Mat_2_doub Mat_elements;


    Mat_2_doub Mat_for_TimeEvo;
    int M_TimeEvo;
    Mat_1_doub Vec_new_TimeEvo;
    double dt_TimeEvo;
    bool TimeEvoPerformed;

    bool Eig_vecs_required;

    bool Saving_Hamil;

    Basis_type& basis;
    Model_type& model;






void Perform_LANCZOS(Matrix_COO &Hamil);
void Get_Dynamics_seed(Matrix_COO &Opr);
void Get_Dynamics_seed(Mat_1_doub vec_);
void Write_full_spectrum();
void Measure_macro_observables(Mat_1_string macro_obs, Hamiltonian_1_COO &Macro_oprts, int state_no);
void Measure_one_point_observables(Mat_1_string one_point_obs, Hamiltonian_2_COO &One_point_oprts, int T_sites, int state_no);
double_type Measure_observable(Matrix_COO &OPR_, int state_no);
double_type Measure_observable(Matrix_COO &OPR_, int state_no, string mult_type);
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

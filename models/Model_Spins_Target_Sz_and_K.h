/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "../basis/Basis_Spins_Target_Sz_and_K.h"
#include "../functions_real.h"
#include "../functions_complex.h"
using namespace std;

#ifndef Model_Spins_Target_Sz_and_K
#define Model_Spins_Target_Sz_and_K

class MODEL_Spins_Target_Sz_and_K{

public:

    Mat_1_real H_field;
    Mat_1_real Dz_anisotropy;
    bool PBC;
    Matrix_COO Hamil;
    int no_of_proc;

    Mat_1_string one_point_obs;
    Hamiltonian_2_COO One_point_oprts;
    Mat_1_string two_point_obs;
    Hamiltonian_3_COO Two_point_oprts;

    int No_of_onepoint_obs;
    Mat_3_doub One_point_oprts_onsite;
    Mat_1_string One_point_strs;
    Mat_1_tetra_int Four_point_sites_set;

    Matrix_COO Dyn_opr;
    string Dyn_opr_string;
    Mat_1_int Dyn_opr_int;
    Mat_1_doub Dyn_opr_coeffs;
    double Dyn_Momentum;
    bool Dyn_Momentum_Resolved;


    string LongRangeExchangeZZfilepath;
    string LongRangeExchangePMfilepath;

    string LongRangeJ1zz_filepath, LongRangeJ1pm_filepath, LongRangeJ2_filepath, LongRangeJ3_filepath  ;

    Mat_2_doub Jzz_Exchange_mat;
    Mat_2_doub Jpm_Exchange_mat;
    Mat_2_doub J1zz_mat, J1pm_mat;
    Mat_4_doub J2_mat;
    Mat_4_doub J3_mat;


    int MultisectionSearch_int;

    int No_of_symmetry_trans;


void Read_parameters(BASIS_Spins_Target_Sz_and_K &basis, string filename);
void Add_diagonal_terms(BASIS_Spins_Target_Sz_and_K &basis);
void Add_non_diagonal_terms(BASIS_Spins_Target_Sz_and_K &basis);
void Add_connections_strictly2point(BASIS_Spins_Target_Sz_and_K &basis);
void Read_parameters_for_dynamics(BASIS_Spins_Target_Sz_and_K &basis, string filename);
void Initialize_Opr_for_Dynamics(BASIS_Spins_Target_Sz_and_K &basis);
void Initialize_State_for_Dynamics(BASIS_Spins_Target_Sz_and_K &basis_new, BASIS_Spins_Target_Sz_and_K &basis, Mat_1_doub &Vec_, Mat_1_doub &Vec_Dyn);
void Initialize_Seed_for_Dynamics(BASIS_Spins_Target_Sz_and_K &basis, BASIS_Spins_Target_Sz_and_K &basis_dyn, Mat_1_doub &VecDyn, double & Num_val,  Mat_1_doub &VecGS);
void Initialize_one_point_to_calculate_from_file(BASIS_Spins_Target_Sz_and_K &basis);
void Act_SziSzj(int &site_i, int &site_j, int &m, Mat_1_int &m_out_array, Mat_1_doub &Coeff_out_Array, BASIS_Spins_Target_Sz_and_K &basis);
void Act_SiSj(int &site_i, int &site_j, int &m, Mat_1_int &m_out_array, Mat_1_doub &Coeff_out_Array, BASIS_Spins_Target_Sz_and_K &basis);
void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site1, int site2, BASIS_Spins_Target_Sz_and_K &basis);
void Direct_product_of_localJ_Ansatz(Mat_1_pair_realInt JJzBasis_state, Mat_1_real &Coefficients, Mat_1_int &m_basis, BASIS_Spins_Target_Sz_and_K &basis);
void Create_JJz_Trial_States(Mat_2_pair_realint &JJzBasis_states, Mat_2_real &Coefficients, Mat_2_int &m_basis, BASIS_Spins_Target_Sz_and_K &basis);
void Overlap_of_JJzBasis_with_State(Mat_1_doub &Vec_ , Mat_1_doub &Overlaps_ , Mat_1_int &sorted_indices,  Mat_2_pair_realint &JJzBasis_states, Mat_2_real &Coefficients, Mat_2_int &m_basis, BASIS_Spins_Target_Sz_and_K &basis);
void Act_Hamil(BASIS_Spins_Target_Sz_and_K &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
void Act_translational_opr(BASIS_Spins_Target_Sz_and_K &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
double_type Get_SzSz(int site1, int site2, Mat_1_doub &Vec_, BASIS_Spins_Target_Sz_and_K &basis);
double_type Get_OprDyn_static(Mat_1_doub &VecGS, BASIS_Spins_Target_Sz_and_K &basis);
double_type Get_Sz(int site1, Mat_1_doub &Vec_, BASIS_Spins_Target_Sz_and_K &basis);
//void Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, BASIS_1_orb_tJ &basis);

};

#endif

#endif



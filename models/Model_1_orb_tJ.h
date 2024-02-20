/*
This class includes the Model for which Lanczos is being done
*/

//#ifndef USE_COMPLEX
#include "../basis/Basis_1_orb_tJ.h"
#include "../functions_real.h"
#include "../functions_complex.h"
using namespace std;

#ifndef Model_1_orb_tJ
#define Model_1_orb_tJ

class MODEL_1_orb_tJ{

public:
    double U;

    Mat_2_real Hopping_mat_NN;//Nearest neighbour hopping matrix
    Mat_1_real H_field;
    bool PBC;
    Matrix_COO Hamil;
    Matrix_COO H_KE;
    Matrix_COO H_Total;

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
    double Dyn_Momentum;
    bool Dyn_Momentum_Resolved;

    bool LongRangeHopping;
    bool LongRangeInteraction;
    bool LongRangeExchangeZZ;
    bool LongRangeExchangePM;
    string LongRangeHoppingfilepath;
    string LongRangeInteractionfilepath;
    string LongRangeExchangeZZfilepath;
    string LongRangeExchangePMfilepath;

    Mat_2_doub DenDen_Interaction_mat;
    Mat_2_doub Jzz_Exchange_mat;
    Mat_2_doub Jpm_Exchange_mat;
    Mat_2_doub Hopping_mat;
    double Dz_anisotropy;

    Mat_1_doub State_Sm_on_GS;


    string RingExchange_filepath;
    Mat_4_doub RingExchange_mat;


    //--------------Time evolution---------------//
    double Ax, Ay;

    //-----------------------------------------//



void Read_parameters(BASIS_1_orb_tJ &basis, string filename);
void Read_parameters_for_dynamics(string filename);
void Add_diagonal_terms(BASIS_1_orb_tJ &basis, string run_type);
void Add_non_diagonal_terms(BASIS_1_orb_tJ &basis);
void Add_connections(BASIS_1_orb_tJ &basis, string run_type);
void Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis);
void Initialize_one_point_to_calculate_from_file(BASIS_1_orb_tJ &basis);
void Initialize_two_point_to_calculate(BASIS_1_orb_tJ &basis);
void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, BASIS_1_orb_tJ &basis);
void Initialize_Opr_for_Dynamics(BASIS_1_orb_tJ &basis);
void Get_Sm_on_GS(Mat_1_doub GS_Sz, BASIS_1_orb_tJ & _BASIS_Szm1, BASIS_1_orb_tJ & _BASIS ,int site);
void Act_Hamil(BASIS_1_orb_tJ &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
void Act_four_spin_opr(BASIS_1_orb_tJ &basis, Mat_1_string &oprs, Mat_1_int & sites, int m, int &m_new, double &sign_FM_, double_type &val_);

};

#endif


//#endif



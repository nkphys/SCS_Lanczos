/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "../basis/Basis_Spins.h"
#include "../functions_real.h"
#include "../functions_complex.h"
using namespace std;

#ifndef Model_Spins
#define Model_Spins

class MODEL_Spins{

public:

    Mat_1_real Hx_field;
    Mat_1_real Hy_field;
    Mat_1_real Hz_field;
    Mat_1_real D_anisotropy;

    Mat_1_real Hx_field_const;
    Mat_1_real Hy_field_const;
    Mat_1_real Hz_field_const;
    Mat_1_real D_anisotropy_const;

    bool PBC;
    Matrix_COO Hamil;

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


    string LongRangeExchangeZZfilepath;
    string LongRangeExchangePMfilepath;
    string MagFieldVectorfilepath;

    Mat_2_doub Jzz_Exchange_mat;
    Mat_2_doub Jpm_Exchange_mat;

    Mat_2_doub Jzz_Exchange_mat_const;
    Mat_2_doub Jpm_Exchange_mat_const;

    string Extenstion_to_FilePaths="";


void Read_parameters(BASIS_Spins &basis, string filename);
void MeasureEnergy(BASIS_Spins &basis, Mat_1_doub &Vec, double_type &Energy_);
void MeasureTwoPointOprs(BASIS_Spins &basis, Mat_1_doub &Vec);
void MeasureLocalOprs(BASIS_Spins &basis, Mat_1_doub &Vec);
void Act_Operator(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub &Vec_out, string opr_str, int opr_site);
void Update_Hamiltonian_Params(BASIS_Spins &basis, double Hx_factor,double Hz_factor, double Jpm_factor, double Jzz_factor);
void Add_diagonal_terms(BASIS_Spins &basis);
void Add_non_diagonal_terms(BASIS_Spins &basis);
void Add_connections(BASIS_Spins &basis);
void Read_parameters_for_dynamics(string filename);
void Initialize_Opr_for_Dynamics(BASIS_Spins &basis, int site);
void Initialize_Opr_for_Dynamics(BASIS_Spins &basis);
void Initialize_one_point_to_calculate_from_file(BASIS_Spins &basis);
void Act_Hamil(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
void Act_diagonal_terms(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
void Act_connections(BASIS_Spins &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);


//void Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, BASIS_1_orb_tJ &basis);

};

#endif

#endif



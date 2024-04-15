/*
This class includes the Model for which Lanczos is being done
*/

#ifndef HIDDEN
#include "../basis/Basis_KondoModel.h"
#include "../functions_real.h"
#include "../functions_complex.h"
using namespace std;

#ifndef Model_KondoModel
#define Model_KondoModel

class MODEL_KondoModel{

public:

    int Length;
    int Target_Total_Ne;
    double Target_Total_Sz;

    Mat_2_real J_LSpins_mat;
    Mat_2_real Hopping_mat;
    Mat_2_real KondoExchange_mat;
    double LocalKondo;

    Mat_1_triad_int KondoHoppings_sites;
    Mat_1_real KondoHoppings;

    Matrix_COO Hamil;

    double EPS_;


void Read_parameters(BASIS_KondoModel &basis, string filename);
void Add_diagonal_terms(BASIS_KondoModel &basis);
void Add_LocalSpin_couplings(BASIS_KondoModel &basis);
void Add_FermionHopping(BASIS_KondoModel &basis);
void Add_Kondocouplings(BASIS_KondoModel &basis);

void Act_Hamil(BASIS_KondoModel &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);

void Initializer_opr_LocalizedSpin_SzSz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_LocalizedSpin_SpSm_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_LocalizedSpin_SmSp_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);

void Initializer_opr_Fermionsz_LocalizedSz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_Fermionsminus_LocalizedSplus_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_Fermionsplus_LocalizedSminus_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);



void Initializer_opr_Fermions_LocalOprs(string opr_type, BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i);
void Initializer_opr_Fermions_szsz_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_Fermions_spsm_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);
void Initializer_opr_Fermions_smsp_Correlation(BASIS_KondoModel &basis, Matrix_COO &OPR, int site_i, int site_j);

//void Initialize_one_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_to_calculate(BASIS_1_orb_tJ &basis);
//void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, BASIS_1_orb_tJ &basis);

};

#endif

#endif



/*
This class includes the Model for which Lanczos is being done
*/


#include "../basis/Basis_1_orb_Hubbard_chain.h"
#include "../functions_real.h"
#include "../functions_complex.h"
#include "../Lanczos_engine.h"
using namespace std;

#ifndef Model_1_orb_Hubb_chain
#define Model_1_orb_Hubb_chain

class MODEL_1_orb_Hubb_chain{

public:
    double U;

    Mat_2_doub Hopping_mat_NN;//LongRange hopping matrix(_NN is misnomer)
    Mat_2_real NonLocalInteractions_mat;//LongRange interactions

    Mat_1_real H_field;
    Mat_1_real Onsite_Energy;
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

    Mat_1_string three_point_oprs;
    Mat_3_int three_point_oprs_sites_set;

    Mat_1_string three_point_intrs;
    Mat_3_int three_point_intrs_sites_set;
    Mat_2_doub three_point_intrs_vals;

    Matrix_COO Dyn_opr;
    string Dyn_opr_string;
    double Dyn_Momentum;
    bool Dyn_Momentum_Resolved;

    bool LongRangeHopping;
    bool LongRangeInteraction;
    string LongRangeHoppingfilepath;
    string LongRangeInteractionfilepath;

    Mat_2_doub Hopping_mat;

    Mat_1_doub State_c_on_GS;
    Mat_1_doub State_cdagger_on_GS;



void Read_parameters(BASIS_1_orb_Hubb_chain &basis, string filename);
void Read_parameters_for_dynamics(string filename);
void Add_diagonal_terms(BASIS_1_orb_Hubb_chain &basis);
void Add_non_diagonal_terms(BASIS_1_orb_Hubb_chain &basis);
void Add_connections(BASIS_1_orb_Hubb_chain &basis);
void Initialize_one_point_to_calculate(BASIS_1_orb_Hubb_chain &basis);
void Initialize_two_point_to_calculate(BASIS_1_orb_Hubb_chain &basis);
void Initialize_one_point_operator_site_specific(string opr_type , Matrix_COO &OPR, int site, BASIS_1_orb_Hubb_chain &basis);
void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, BASIS_1_orb_Hubb_chain &basis);
void Initialize_three_point_operator_sites_specific(string opr_type , Matrix_COO &OPR, int site1, int site2, int site3, BASIS_1_orb_Hubb_chain &basis);
void Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_chain &basis);
void Get_cdagger_on_GS(LANCZOS & lanczos, BASIS_1_orb_Hubb_chain & basis_Np1, BASIS_1_orb_Hubb_chain & basis, Mat_1_trio_int TRIO_VEC, Mat_1_doub values);
void Get_c_on_GS(LANCZOS & lanczos, BASIS_1_orb_Hubb_chain & basis_Nm1, BASIS_1_orb_Hubb_chain & basis, Mat_1_trio_int TRIO_VEC, Mat_1_doub values);


};

#endif

//#endif



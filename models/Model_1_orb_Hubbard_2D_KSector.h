/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_1_orb_Hubbard_2D_KSector.h"
#include "../functions_real.h"
#include "../functions_complex.h"

#ifndef Model_1_orb_Hubb_2D_KSector
#define Model_1_orb_Hubb_2D_KSector

class MODEL_1_orb_Hubb_2D_KSector{

public:
    double U;
    double_type Hopping_NN;//Nearest neighbour hopping
    Mat_1_real H_field;

    string geometry_;
    Mat_2_doub Hopping_mat_NN;//LongRange hopping matrix(_NN is misnomer)
    Mat_2_real NonLocalInteractions_mat;//LongRange interactions
    Mat_1_real Onsite_Energy;

    Matrix_COO Hamil;
    bool Calculate_observables_onepoint;
    Mat_1_string one_point_obs;
    Hamiltonian_2_COO One_point_oprts;

    Mat_1_string obs_string;
    Hamiltonian_1_COO Oprts_array;

    Mat_1_string macro_obs;
    Hamiltonian_1_COO Macro_oprts;
    Mat_1_string two_point_obs;
    Hamiltonian_3_COO Two_point_oprts;
    Matrix_COO Dyn_opr;
    Matrix_COO OPR_SF;
    string Dyn_opr_string;
    int Dyn_Momentum_x, Dyn_Momentum_y;
    bool Dyn_Momentum_Resolved;
    bool PBC;
    Mat_1_real Momentum_values;
    Mat_2_doub overlap_matrix_for_Anzatz_basis;


    Mat_1_string three_point_oprs;
    Mat_3_int three_point_oprs_sites_set;

    Mat_2_doub State_Szq_on_GS;
    string file_read_basis_Kminusq;

    bool Saving_Hamil;
     int NProcessors_;


    void Read_parameters(BASIS_1_orb_Hubb_2D_KSector &basis, string filename);
    void Add_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Add_non_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Add_connections(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Read_parameters_for_dynamics(string filename);
    void Getting_Local_Sz_Opr(BASIS_1_orb_Hubb_2D_KSector &basis, Matrix_COO & OPR_LOCAL, int site);
    //void Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Initialize_Opr_for_Dynamics(BASIS_1_orb_Hubb_2D_KSector &basis ,BASIS_1_orb_Hubb_2D_KSector & basis_Kminusq);
    void Initialize_Opr_for_Structure_factor(BASIS_1_orb_Hubb_2D_KSector &basis);

    double_type Measure_Opr_for_Structure_factor(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &EigVec_);
    double_type Measure_Opr_for_LocalNupNdn(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &EigVec_);

    void Initialize_Oprs_for_meausurement(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Calculate_two_point_observables(Mat_1_doub &Vec_);
    void Initialize_two_point_operator_sites_specific(string opr_type , Matrix_COO &OPR_ , int site_x, int site_y, BASIS_1_orb_Hubb_2D_KSector &basis);
    double_type Measure_two_point_operator_sites_specific(string opr_type , Mat_1_doub &EigVec_ , int site_x, int site_y, BASIS_1_orb_Hubb_2D_KSector &basis);
    void Initialize_three_point_operator_sites_specific(string opr_type  , Matrix_COO &OPR_, int sitejx, int sitejy, int sitelx, int sitely, BASIS_1_orb_Hubb_2D_KSector &basis);
    void Initialize_chiral_corr_operator_sites_specific(string opr_type, Matrix_COO &OPR_, int lx_rel, int ly_rel , BASIS_1_orb_Hubb_2D_KSector &basis);
    void Act_Hamil(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);

    void Act_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);
    void Act_connections(BASIS_1_orb_Hubb_2D_KSector &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);


};

#endif

/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_2_orb_Hubbard_chain_KSector.h"

#include "../Lanczos_engine.h"


#ifndef Model_2_orb_Hubb_chain_KSector
#define Model_2_orb_Hubb_chain_KSector

class MODEL_2_orb_Hubb_chain_KSector{

public:
    double U;
    double U_p;
    double J_H;
    double Dz_Anisotropy;
    Mat_2_real Hopping_mat_NN;//Nearest neighbour hopping matrix
    Mat_1_real CFS;
    Mat_1_real H_field;

    Matrix_COO Hamil;
    bool Calculate_observables_onepoint;
    Mat_1_string one_point_obs;
    Hamiltonian_2_COO One_point_oprts;
    Mat_1_string macro_obs;
    Hamiltonian_1_COO Macro_oprts;
    Mat_1_string two_point_obs;
    Hamiltonian_3_COO Two_point_oprts;
    Matrix_COO Dyn_opr;
    string Dyn_opr_string;
    double Dyn_Momentum;
    bool Dyn_Momentum_Resolved;
    bool PBC;
    Mat_1_real Momentum_values;
    Mat_2_doub overlap_matrix_for_Anzatz_basis;


    //[pair no.][Pair_Geometry]
    Mat_1_string Variational_state_pair_Geometry;

    //[pair no.][Spin symmetry]
    Mat_1_string Variational_state_pair_spin_symmetry;

    //[pair no.][orbital symmetry]
    Mat_1_string Variational_state_pair_orbital_symmetry;

    //[pair no.].first=[sitei];[pair no.].second=[sitej]
    Mat_1_intpair  Variational_state_pair_sites;

    //VB states with same overlap
    Mat_2_int Degenerate_states;
    Mat_1_doub Distinct_overlaps;


    Mat_1_doub State_;
    Mat_2_doub BASIS_STATES_ANSATZ;
    Mat_2_doub BASIS_STATES_ANSATZ_2HOLES;
    Mat_2_doub STATES_OS_TS;
    Mat_2_doub STATES_OS_TS_DIAGONAL_HOLES;
    Mat_1_string variational_state_contruction_strings;
    Mat_1_string variational_2holes_classes_contruction_strings;
    Mat_1_doub Ansatz_Basis_Overlap_with_GS;
    Mat_1_doub Ansatz_2holes_Basis_Overlap_with_GS;
    Mat_2_doub OVERLAP_MATRIX_2holes_Basis;
    Matrix_COO Pair_Annihilation;

    Mat_1_doub STATE_RVB;
    Mat_1_doub STATE_TPS; //Triplet_Product_state


    void Read_parameters(BASIS_2_orb_Hubb_chain_KSector &basis, string filename);
    void Add_diagonal_terms(BASIS_2_orb_Hubb_chain_KSector &basis);
    void Add_non_diagonal_terms(BASIS_2_orb_Hubb_chain_KSector &basis);
    void Add_connections(BASIS_2_orb_Hubb_chain_KSector &basis);

    void Initialize_two_point_operator_sites_orbital_specific(string type , Matrix_COO &OPR_, int site1, int gamma1, int site2, int gamma2, BASIS_2_orb_Hubb_chain_KSector &basis);

    Mat_1_doub Act_Orbital_Exchange(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub vec_);

    Mat_1_doub Act_Reflection_about_Central_site(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub vec_);

    string Bond_type_int_tag_to_string(int int_temp);
    void Create_variational_state_strings_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis);
    void Create_variational_state_strings_2_holes(BASIS_2_orb_Hubb_chain_KSector &basis);
    void Perform_RVB_Analysis_at_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec);
    void Perform_RVB_Analysis_at_2_hole_doped(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec);
    void Perform_Product_state_Analysis_at_half_filling(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec);

    void Perform_Product_state_Analysis_at_2_hole_doped(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &Eig_vec);

    void Get_State_by_acting_creation_operators(BASIS_2_orb_Hubb_chain_KSector &basis,
                                                Mat_1_trio_int MAT_,
                                                Mat_1_doub &State_temp);
    void Get_Variational_State(BASIS_2_orb_Hubb_chain_KSector &basis, int no_of_pairs);
    void Read_Mat_2_trio(Mat_2_trio_int &MAT_TEMP, Mat_1_doub &VALUES_TEMP,
                                                 int pair_no);
    void Read_Anzatz_basis(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_);
    void Get_overlap_matrix_for_Anzatz_basis(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub &GS_);
    void Variational_state_optimization(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_);
    void Variational_state_optimization_old(BASIS_2_orb_Hubb_chain_KSector &basis, Mat_1_doub GS_);

};

#endif

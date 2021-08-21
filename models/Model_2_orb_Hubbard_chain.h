/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_2_orb_Hubbard_chain.h"
#include "../functions_real.h"
#include "../functions_complex.h"



#ifndef Model_2_orb_Hubb_chain
#define Model_2_orb_Hubb_chain

class MODEL_2_orb_Hubb_chain{

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



    void Create_STATES_OS_TS_DIAGONAL_HOLES(BASIS_2_orb_Hubb_chain &basis,
                                            BASIS_2_orb_Hubb_chain &basis_nm2,
                                            Mat_2_doub &STATES_OS_TS_);
    void Optimize_Anzatz_basis_2Holes(BASIS_2_orb_Hubb_chain &basis_nm2, Mat_1_doub &OPT_VEC, string read_overlaps);
    void Optimize_overlap_Nm2_Variational_State(BASIS_2_orb_Hubb_chain &basis,
                                                BASIS_2_orb_Hubb_chain &basis_nm2,
                                                Mat_2_doub &STATES_OS_TS_,
                                                Mat_1_doub &Eig_vec,
                                                Mat_1_doub &vec_final);
    void Create_OS_TS_states_by_reading(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
    void Create_states_with_hole_hopping(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub Alpha_, Mat_1_doub &vec);
    double Get_Holes_Projected_state_probability(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec);
    Mat_1_doub Get_2Holes_Projected_state(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec);
    void Get_Holes_Projected_state(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int Hole_positions, Mat_1_doub &vec);

    void Check_orbital_symmetry(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_);
    Mat_1_doub Act_Translational_operator(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_);
    Mat_1_doub Act_Orbital_Exchange(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub vec_);
    void Get_Pair_Operator_Matrix(BASIS_2_orb_Hubb_chain &basis, BASIS_2_orb_Hubb_chain &basis_nm2,
                                  trio_int trio_0,trio_int trio_1, double value_out_ );
    void Variational_state_optimization(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
    void Get_overlap_matrix_for_Anzatz_basis(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
    void Create_Anzatz_basis_2Holes(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
   // void Create_Anzatz_basis_2Holes_new(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
    void Read_Anzatz_basis(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub GS_);
    void Get_State_by_acting_creation_operators(BASIS_2_orb_Hubb_chain &basis, Mat_1_trio_int MAT_,
                                                Mat_1_doub &State_temp);
    void Discard_double_occupancies(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub &vec);
    void Read_parameters(BASIS_2_orb_Hubb_chain &basis, string filename);
    void Read_parameters_for_variational_state(BASIS_2_orb_Hubb_chain &basis, string filename);
    void Read_parameters_for_dynamics(string filename);
    void Add_diagonal_terms(BASIS_2_orb_Hubb_chain &basis);
    void Add_non_diagonal_terms(BASIS_2_orb_Hubb_chain &basis);
    void Add_connections(BASIS_2_orb_Hubb_chain &basis);


    void Initialize_macro_oprs_to_calculate(BASIS_2_orb_Hubb_chain &basis);
    void Initialize_one_point_to_calculate(BASIS_2_orb_Hubb_chain &basis);
    void Initialize_two_point_operator_sites_specific(string type , Matrix_COO &OPR_, int site, int site2, BASIS_2_orb_Hubb_chain &basis);
    void Initialize_two_point_to_calculate(BASIS_2_orb_Hubb_chain &basis);
    void Initialize_Opr_for_Dynamics(BASIS_2_orb_Hubb_chain &basis);
    void Calculate_Local_Obs_for_States_to_Look(bool calculate_local_obs_for_states_to_look,
                                                Mat_1_int & states_to_look,
                                                string file_Loc_obs_in_basis_of_states,
                                                int no_basis_to_check,
                                                Mat_2_pair_realint &Overlaps,
                                                BASIS_2_orb_Hubb_chain & basis);
    void Get_Variational_State(BASIS_2_orb_Hubb_chain &basis, int no_of_pairs);
    void Read_Mat_2_trio(Mat_2_trio_int &MAT_TEMP, Mat_1_doub &VALUES_TEMP,int pair_no);
    void Act_Hamil(BASIS_2_orb_Hubb_chain &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);


};

#endif

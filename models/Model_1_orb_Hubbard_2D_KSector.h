/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_1_orb_Hubbard_2D_KSector.h"
#include "../Lanczos_engine.h"


#ifndef Model_1_orb_Hubb_2D_KSector
#define Model_1_orb_Hubb_2D_KSector

class MODEL_1_orb_Hubb_2D_KSector{

public:
    double U;
    double Hopping_NN;//Nearest neighbour hopping matrix
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



    void Read_parameters(BASIS_1_orb_Hubb_2D_KSector &basis, string filename);
    void Add_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Add_non_diagonal_terms(BASIS_1_orb_Hubb_2D_KSector &basis);
    void Add_connections(BASIS_1_orb_Hubb_2D_KSector &basis);

};

#endif

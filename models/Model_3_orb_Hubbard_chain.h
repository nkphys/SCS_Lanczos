/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_3_orb_Hubbard_chain.h"

#include "../Lanczos_engine.h"


#ifndef Model_3_orb_Hubb_chain
#define Model_3_orb_Hubb_chain

class MODEL_3_orb_Hubb_chain{

public:
    double U;
    double U_p;
    double J_H;
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




void Read_parameters(BASIS_3_orb_Hubb_chain &basis, string filename);
void Read_parameters_for_dynamics(string filename);
void Add_diagonal_terms(BASIS_3_orb_Hubb_chain &basis);
void Add_non_diagonal_terms(BASIS_3_orb_Hubb_chain &basis);
void Add_connections(BASIS_3_orb_Hubb_chain &basis);

void Initialize_macro_oprs_to_calculate(BASIS_3_orb_Hubb_chain &basis);
void Initialize_one_point_to_calculate(BASIS_3_orb_Hubb_chain &basis);
void Initialize_two_point_to_calculate(BASIS_3_orb_Hubb_chain &basis);
void Initialize_Opr_for_Dynamics(BASIS_3_orb_Hubb_chain &basis);
void Calculate_Local_Obs_for_States_to_Look(LANCZOS & lanczos, BASIS_3_orb_Hubb_chain & basis);
};

#endif

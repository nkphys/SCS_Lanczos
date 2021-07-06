/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_3_orb_Hubbard_chain_GC.h"
#include "../basis/Basis_3_orb_Hubbard_chain_GC_restricted.h"
#include "../functions_real.h"
#include "../functions_complex.h"

#define PI 3.14159265

#ifndef Model_3_orb_Hubb_chain_GC_class
#define Model_3_orb_Hubb_chain_GC_class

template <typename Basis_type>
class MODEL_3_orb_Hubb_chain_GC{

public:
    MODEL_3_orb_Hubb_chain_GC(Basis_type& Basis_)
        : basis(Basis_)
    {
    }

    double U;
    double U_p;
    double J_H;
    Mat_2_real Hopping_mat_NN;//Nearest neighbour hopping matrix
    bool USE_LONG_RANGE_HOPPINGS;
    string LongRangeHoppingfilepath;
    Mat_2_doub Hopping_mat_LongRange;
    Mat_1_real CFS;
    Mat_2_real CFS_SITE_RESOLVED_;
    bool CFS_SITE_RESOLVED_bool;
    Mat_1_real H_field;
    double lambda_SOC;

    Mat_1_string macro_obs;
    Hamiltonian_1_COO Macro_oprts;
    Matrix_COO Hamil;
    bool Calculate_observables_onepoint;
    Mat_1_string one_point_obs;
    Hamiltonian_2_COO One_point_oprts;
    Mat_1_string two_point_obs;
    Hamiltonian_3_COO Two_point_oprts;
    Matrix_COO Dyn_opr;
    string Dyn_opr_string;
    double Dyn_Momentum;
    bool Dyn_Momentum_Resolved;
    bool PBC,TBC;
    Mat_1_doub State_c_on_GS;
    Mat_1_doub State_cdagger_on_GS;

    enum {num_Hopping, num_CFS, num_InterOrbRepulsion, num_IntraOrbRepulsion
                   , num_Hunds_zz, num_Hunds_pm_mp, num_PairHopping};

    void Read_parameters(string filename);
    void Read_parameters_for_dynamics(string filename);
    void Add_diagonal_terms();
    void Add_non_diagonal_terms();
    void Add_connections();
    void Add_Spin_Orbit_Coupling();
    void Initialize_one_point_to_calculate();
    void Initialize_two_point_to_calculate();
    void Initialize_Opr_for_Dynamics();
    void Calculate_Local_Obs_for_States_to_Look(bool calculate_local_obs_for_states_to_look,
                                                Mat_1_int & states_to_look,
                                                string file_Loc_obs_in_basis_of_states,
                                                int no_basis_to_check,
                                                Mat_2_pair_realint &Overlaps);
    void Get_c_on_GS(Mat_1_doub &EigVec_, BASIS_3_orb_Hubb_chain_GC & basis_Nm1,
                     Mat_1_trio_int TRIO_VEC, Mat_1_doub values);
    void Get_cdagger_on_GS(Mat_1_doub &EigVec_, BASIS_3_orb_Hubb_chain_GC & basis_Np1,
                     Mat_1_trio_int TRIO_VEC, Mat_1_doub values);
    void Get_Delta_Matrix(Mat_2_doub &EigVecs_);
    void Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site);
    void Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site, int site_p);
    void Get_ExcitonCoherence_Length(Mat_1_doub &vector_used);
    void Act_Hamil(Basis_type &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);


private:
    Basis_type& basis;

};


/*
void Read_parameters(string filename);
void Read_parameters_for_dynamics(string filename);
void Add_diagonal_terms();
void Add_non_diagonal_terms();
void Add_connections();
void Add_Spin_Orbit_Coupling();
void Initialize_one_point_to_calculate();
void Initialize_two_point_to_calculate();
void Initialize_Opr_for_Dynamics(LANCZOS &lanczos_GS);
void Calculate_Local_Obs_for_States_to_Look(LANCZOS & lanczos);
*/





#endif

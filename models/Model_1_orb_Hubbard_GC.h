/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_1_orb_Hubbard_GC.h"
#include "../functions_real.h"
#include "../functions_complex.h"
#define PI 3.14159265

#ifndef Model_1_orb_Hubbard_GC_Class
#define Model_1_orb_Hubbard_GC_Class

template <typename Basis_type>
class MODEL_1_orb_Hubbard_GC{

public:
    MODEL_1_orb_Hubbard_GC(Basis_type& Basis_)
        : basis(Basis_)
    {
    }

    double U;
    Mat_1_real H_field;

    //2Lx2L matrix for L sites, "2" for spin, it includes hopping and SOC both.
    //index=L*spin + site; spin \in {0=UP,1=DOWN}, and L=No of sites
    Mat_2_doub Hopping_mat_LongRange;
    Mat_1_real CFS;

    string LongRangeHoppingfilepath;
    bool CFS_SITE_RESOLVED_bool;

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
    Mat_1_doub State_c_on_GS;
    Mat_1_doub State_cdagger_on_GS;

    void Read_parameters(string filename);
    void Read_parameters_for_dynamics(string filename);
    void Add_diagonal_terms();
    void Add_non_diagonal_terms();
    void Add_connections();
    //    void Initialize_one_point_to_calculate();
    //    void Initialize_two_point_to_calculate();
    void Initialize_Opr_for_Dynamics();
    //    void Calculate_Local_Obs_for_States_to_Look(LANCZOS & lanczos);
    void Get_c_on_GS(Mat_1_doub & EigVec_, BASIS_1_orb_Hubbard_GC & basis_Nm1,
                     Mat_1_trio_int TRIO_VEC, Mat_1_doub values);
    void Get_cdagger_on_GS(Mat_1_doub & EigVec_, BASIS_1_orb_Hubbard_GC & basis_Np1,
                           Mat_1_trio_int TRIO_VEC, Mat_1_doub values);
    void Calculate_one_point_observables(Mat_1_doub &Vec_);
    void Calculate_two_point_observables(Mat_1_doub &Vec_);
    void Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site);
    void Get_CdaggerC_type_Opr(Mat_2_doub AMat, Matrix_COO &OPR, int site, int site_p);

    void Act_Hamil(BASIS_1_orb_Hubbard_GC &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);


private:
    Basis_type& basis;

};



#endif

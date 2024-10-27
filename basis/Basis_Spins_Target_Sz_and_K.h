/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"
#include "../functions_real.h"
#include "../functions_complex.h"

using namespace std;

#ifndef Basis_Spins_Target_Sz_and_K
#define Basis_Spins_Target_Sz_and_K
class BASIS_Spins_Target_Sz_and_K{

public:
    int Length;
    int TwoTimesSpin;
    double Target_Total_Sz;
    int Target_Total_Value;
    double SPIN;
    int BASE;
    unsigned long long int basis_size;
    Mat_1_ullint Partitions_Dec;
    Mat_1_intpair Partitions_pos;


    Mat_1_ullint D_basis;
    unsigned long long int D_min, D_max;


    bool read_basis, write_basis;
    string read_basis_file, write_basis_file;


    //Using Lin tables
    int Length_A, Length_B;

    Mat_1_int Partition_Length;
    Mat_1_int Dec_to_Index_part0_, Dec_to_Index_part1_;
    Mat_1_ullint Index_to_Dec_part0_, Index_to_Dec_part1_;
    Mat_1_ullint MainIndex_to_Dec_part0_, MainIndex_to_Dec_part1_;

    Mat_1_ullint MainIndex_to_Dec_partA_, MainIndex_to_Dec_partB_;
    Mat_1_int DecB_to_Index_partB_;
    Mat_1_int DecB_to_Sections;
    Mat_2_ullint Index_partA_to_DecA_in_givenSection;
    Mat_2_int R_vals_basis;

    Mat_1_real Norm_vals_basis;

    //Symmetries: Translations, may be reflection etc.
    int No_of_symms_trans;
    Mat_1_string TransformationFiles;
    Mat_3_int Transformation_matrices;
    Mat_3_int InverseTransformation_matrices;
    Mat_1_int Transformation_Cardinalities;
    Mat_1_int EigenvaluesTargeted;

bool CheckState(ulli d_tempA, ulli d_tempB, double_type & Gamma_, int & Repetitions_n);
void Check_Commutation_and_Inversion_bw_LatticeSymmetry_Oprs();
void Construct_basis();
void clear();
};

#endif

/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include <assert.h>

using namespace std;

#ifndef Basis_multi_orb_Hubb_chain_GC
#define Basis_multi_orb_Hubb_chain_GC

class BASIS_multi_orb_Hubb_chain_GC{

public:
    int Length;
    int N_total;
    int N_orb;
    int No_of_sys_copies;
    double Boundary_phase_exponent;
    int n_orb;
    Mat_1_int D_up_basis;
    Mat_1_int D_dn_basis;
    Mat_1_int D_up_min, D_up_max;
    Mat_1_int D_dn_min, D_dn_max;
    Mat_3_int D_up_reverse;
    Mat_3_int D_dn_reverse;
    Mat_3_int D_updn_reverse;
    Mat_1_intpair Nup_offsets;
    Mat_2_int Canonical_partition_up;
    Mat_2_int Canonical_partition_dn;
    bool Restricted;
    Mat_1_int Local_occupations_allowed;


void Construct_basis();
void Print_state_in_file_with_basis(Mat_1_doub vec, string filename);
void clear();
};

#endif

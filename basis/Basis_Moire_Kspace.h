/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"


using namespace std;


#ifndef Basis_Moire_Kspace
#define Basis_Moire_Kspace

class BASIS_Moire_Kspace{

public:
    int Length1, Length2, Length;
    int Nup;
    int Ndn;
    int K1_target, K2_target;
    int n_orb;
    Mat_1_ullint D_up_basis;
    Mat_1_ullint D_dn_basis;
    Mat_1_ullint D_up_basis_range;
    Mat_1_int D_up_basis_range_min;
    Mat_1_int D_up_basis_range_max;

    int basis_size;

void Construct_basis_old();
void Construct_basis();
void clear();
};

#endif

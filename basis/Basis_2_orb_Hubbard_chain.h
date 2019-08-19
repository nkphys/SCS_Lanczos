/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"

using namespace std;

#ifndef Basis_2_orb_Hubb_chain
#define Basis_2_orb_Hubb_chain

class BASIS_2_orb_Hubb_chain{

public:
    int Length;
    int Nup;
    int Ndn;
    int n_orb;
    Mat_1_int D_up_basis;
    Mat_1_int D_dn_basis;



void Construct_basis();
void clear();
};

#endif

/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"

using namespace std;

#ifndef Basis_1_orb_Hubb_chain
#define Basis_1_orb_Hubb_chain

class BASIS_1_orb_Hubb_chain{

public:
    int Length;
    int Nup;
    int Ndn;

    Mat_1_int D_up_basis;
    Mat_1_int D_dn_basis;

    Mat_1_int inverse_Dup;
    Mat_1_int inverse_Ddn;

    Mat_1_int partitions_up;
    Mat_1_int Dup_val_at_partitions;

    Mat_1_int partitions_dn;
    Mat_1_int Ddn_val_at_partitions;

    int partitions_length_up;
    int partitions_length_dn;





void Construct_basis();
void clear();
};

#endif

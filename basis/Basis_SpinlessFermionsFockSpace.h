/*
This class includes the Basis used for Model for which Lanczos is being done
*/

#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"

using namespace std;

#ifndef Basis_SpinlessFermionsFockSpace
#define Basis_SpinlessFermionsFockSpace

class BASIS_SpinlessFermionsFockSpace{

public:
    int Length;
    unsigned long long int basis_size;

    Mat_1_int D_basis;

    Mat_1_int inverse_D;

    Mat_1_int partitions_;
    Mat_1_int D_val_at_partitions;

    int partitions_length;

    unsigned long long int DMax_, DMin_;
    int Base;

void Construct_basis();
void clear();
};

#endif

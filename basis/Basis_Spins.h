/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"

using namespace std;

#ifndef Basis_Spins
#define Basis_Spins
class BASIS_Spins{

public:
    int Length;
    int TwoTimesSpin;
    double SPIN;
    int BASE;
    int basis_size;

    Mat_1_int D_basis;
    Mat_1_int inverse_D;
    unsigned long long int D_min, D_max;


void Construct_basis();
void clear();
};

#endif

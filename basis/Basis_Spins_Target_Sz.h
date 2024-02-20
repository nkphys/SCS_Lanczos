/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"

using namespace std;

#ifndef Basis_Spins_Target_Sz
#define Basis_Spins_Target_Sz
class BASIS_Spins_Target_Sz{

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


void Construct_basis();
void clear();
};

#endif

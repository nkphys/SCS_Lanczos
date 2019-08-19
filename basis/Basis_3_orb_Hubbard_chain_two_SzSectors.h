/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"

using namespace std;

#ifndef Basis_3_orb_Hubb_chain_two_SzSectors
#define Basis_3_orb_Hubb_chain_two_SzSectors

class BASIS_3_orb_Hubb_chain_two_SzSectors{

public:
    int Length;
    int Nup_sec1;
    int Ndn_sec1;
    int Nup_sec2;
    int Ndn_sec2;
    int Nup;
    int Ndn;
    int n_orb;
    Mat_1_int D_up_basis_sec1;
    Mat_1_int D_dn_basis_sec1;
    int basis_size_sec1;
    Mat_1_int D_up_basis_sec2;
    Mat_1_int D_dn_basis_sec2;
    int basis_size_sec2;


    Mat_1_int inverse_Dup_sec1;
    Mat_1_int inverse_Ddn_sec1;
    Mat_1_int inverse_Dup_sec2;
    Mat_1_int inverse_Ddn_sec2;



void Construct_basis();
void clear();
};

#endif

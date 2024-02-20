/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"
#include "../Base_to_Decimal.h"

using namespace std;

#ifndef Basis_KondoModel
#define Basis_KondoModel
class BASIS_KondoModel{

public:
    int Length;

    double Target_Total_Sz;
    int Target_Total_Ne;

    int Target_Total_SzWithOffset;


    Mat_2_int Dec_LocalizedSpins;
    Mat_2_int InverseDec_LocalizedSpins;

    Mat_2_intpair Dec_dnup_Fermions;
    Mat_2_int InverseDec_Fermions;

    Mat_1_int Dmax_dn, Dmax_up;

    Mat_1_intpair SzWithOffsetAllowed;
    Mat_1_int Inverse_SzWithOffsetAllowed;
    int No_of_SzSets;

    unsigned long long int basis_size;

//    Mat_2_ullint D_up_basis,D_dn_basis;
//    Mat_2_int inverse_Dup, inverse_Ddn;
//    Mat_2_int D_basis_LocalizedS;
//    unsigned long long int D_min, D_max;


    bool read_basis, write_basis;
    string read_basis_file, write_basis_file;


int Concatenate_Spins_and_Fermions(int i_fermion, int i_Lspins, int SzSetNo);
void Split_to_Spins_and_Fermions(int &i_fermion, int &i_Lspins, int SzSetNo, int i_local);
ulli Get_basis_ind(int i_local, int SzSetNo);
void Get_SzSec_and_localindex(ulli m, int &i_local, int &SzSetNo);
void Get_SzSectors_Allowed();
void Construct_basis();
int Get_max(Mat_1_int A_vec);
int Get_max_first_type(Mat_1_intpair A_vec);
int Get_max_second_type(Mat_1_intpair A_vec);
void clear();
};

#endif

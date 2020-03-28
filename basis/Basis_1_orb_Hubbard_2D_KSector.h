/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include "../tensor_type.h"
#include "../binary_decimal.h"

using namespace std;

#ifndef Basis_1_orb_Hubb_2D_KSector
#define Basis_1_orb_Hubb_2D_KSector

class BASIS_1_orb_Hubb_2D_KSector{

public:
    int Length;
    int Lx, Ly;
    int Nup;
    int Ndn;
    int Momentum_n;
    int Momentum_nx, Momentum_ny;

    bool Write_Basis;
    bool Read_Basis;

    //arrays below save only representative states of the group
    Mat_1_int D_up_basis;
    Mat_1_int D_dn_basis;

    //arrays below save norm of the |a(K)> states
    Mat_1_real D_Norm;

    //arrays below save Periodicity of the representative states
    Mat_1_int Dx_Period;
    Mat_1_int Dy_Period;


    //arrays below save range of the basis index int which Dup is present
    Mat_1_intpair Dup_Range;



void Construct_basis();
void clear();
};

#endif

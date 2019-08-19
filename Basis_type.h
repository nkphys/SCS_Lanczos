/*
This class chooses the Model class
*/


#include "basis/Basis_3_orb_Hubbard_chain.h"
#include "basis/Basis_2_orb_Hubbard_chain.h"
#include "basis/Basis_1_orb_Hubbard_chain.h"


using namespace std;

#ifndef Basis_
#define Basis_


class BASIS_{

public:
    BASIS_(string model_name_string):
    model_name_string_(model_name_string)
    {
    }

    string model_name_string_;

    int Length;
    int Nup;
    int Ndn;

    Mat_1_int D_up_basis;
    Mat_1_int D_dn_basis;




template <class BASIS_TYPE>
void Construct_basis(){

    BASIS_TYPE _Basis;

//NOW CONTINUE

}

};

#endif


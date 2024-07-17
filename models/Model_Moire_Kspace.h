/*
This class includes the Model for which Lanczos is being done
*/

#include "../basis/Basis_Moire_Kspace.h"
#include "../functions_real.h"
#include "../functions_complex.h"
#ifndef Model_Moire_Kspace
#define Model_Moire_Kspace

class MODEL_Moire_Kspace{

public:
    Mat_3_doub W_int;
    Mat_3_doub Ek_dispersion;
    int n_orb;
    int Length1, Length2, Length;
    Mat_1_string DispersionFilepaths;
    string InteractionFilepath;

    Matrix_COO Hamil;

void Read_parameters(BASIS_Moire_Kspace &basis, string filename);
void Read_Interations(BASIS_Moire_Kspace &basis);
void Read_Dispersion(BASIS_Moire_Kspace &basis);
void Add_Interaction_terms(BASIS_Moire_Kspace &basis);
void Add_Dispersion_term(BASIS_Moire_Kspace &basis);

void Act_Hamil(BASIS_Moire_Kspace &basis, Mat_1_doub &Vec_in, Mat_1_doub& Vec_out);

void Initialize_Sk_opr(BASIS_Moire_Kspace &basis, Matrix_COO &OPR_, int q_1, int q_2, string alpha, string beta);
};

#endif

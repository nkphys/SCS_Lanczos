#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include "Lanczos_engine.h"

#ifdef USE_COMPLEX
#include "functions_complex.h"
#else
#include "functions_real.h"
#endif


using namespace std;

#ifndef LTLM_STATIC_engine
#define LTLM_STATIC_engine

class LTLM_STATIC{

public:

    //NOTATION : STATISTICAL AVG is average over random initial seeds
    //           Qauntum Avg is just expectation values.
    int Total_Random_States;
    double Temperature;
    double Temperature_min, Temperature_max, delta_Temperature;
    int N_Temperature_points;
    int M_;
    double Beta;
    double Boltzman_const;

    Mat_2_doub Sum_Opr_val;
    Mat_1_doub Sum_Hamil, Sum_Hamil_Sqr;
    Mat_1_doub Sum_Hamil2, Sum_Hamil2_Sqr;
    Mat_1_doub Sum_Partition_Func, Sum_Partition_Func_Sqr;
    Mat_1_doub Quantum_Avg_Hamil, Quantum_Avg_Hamil2;

    Matrix_COO Hamil;


    void Perform_LTLM(string inp_filename, Hamiltonian_1_COO& OPR_);

};

#endif



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

#ifndef FTLM_DYNAMICS_engine
#define FTLM_DYNAMICS_engine

class FTLM_DYNAMICS{

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


    Mat_2_doub Sum_Cw;
    Mat_1_real Sum_Partition_Func;

    Matrix_COO Hamil;

    Mat_1_real Evals1, Evals2;


    void Perform_FTLM(string inp_filename, Matrix_COO& OPR_);

};

#endif



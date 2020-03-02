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

#ifndef FTLM_STATIC_engine
#define FTLM_STATIC_engine

class FTLM_STATIC{

public:

    //NOTATION : STATISTICAL AVG is average over random initial seeds
    //           Qauntum Avg is just expectation values.
    int Total_Random_States;
    double Temperature;
    int M_;
    double Beta;
    double Boltzman_const;

    double Sum_Hamil, Sum_Hamil_Sqr;
    double Sum_Hamil2, Sum_Hamil2_Sqr;
    double Sum_Partition_Func, Sum_Partition_Func_Sqr;

    double Quantum_Avg_Hamil, Quantum_Avg_Hamil2;

    Matrix_COO Hamil;


    void Perform_FTLM(string inp_filename);

};

#endif



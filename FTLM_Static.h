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

template <typename Basis_type, typename Model_type>
class FTLM_STATIC{


public:
    FTLM_STATIC(Basis_type& Basis_, Model_type& Model_)
        :basis(Basis_), model(Model_)
    {

    }



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
    Mat_1_real Sum_Hamil, Sum_Hamil_Sqr;
    Mat_1_real Sum_Hamil2, Sum_Hamil2_Sqr;
    Mat_1_real Sum_Partition_Func, Sum_Partition_Func_Sqr;
    Mat_1_real Quantum_Avg_Hamil, Quantum_Avg_Hamil2;

    Matrix_COO Hamil;

    Basis_type& basis;
    Model_type& model;


    void Perform_FTLM(string inp_filename, Hamiltonian_1_COO& OPR_);

};

#endif



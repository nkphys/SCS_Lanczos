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

#ifndef LTLM_DYNAMICS_class
#define LTLM_DYNAMICS_class

template <typename Basis_type, typename Model_type>
class LTLM_DYNAMICS{

public:
    LTLM_DYNAMICS(Basis_type& Basis_, Model_type& Model_)
        :basis(Basis_), model(Model_)
    {

    }


    //NOTATION : STATISTICAL AVG is average over random initial seeds
    //           Qauntum Avg is just expectation values.
    int Total_Random_States;
    double Temperature;
    int M_;
    double Beta;
    double Boltzman_const;


    Mat_1_doub Sum_Cw;
    double Sum_Partition_Func;

    Matrix_COO Hamil;
    Matrix_COO Hamil_intrmd;


    Mat_1_real Evals1, Evals2;

    void Perform_LTLM(string inp_filename, Matrix_COO& OPR_);


private :
    Basis_type& basis;
    Model_type& model;


};

#endif



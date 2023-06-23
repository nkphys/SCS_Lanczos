#include <iostream> //for cin and cout
#include <math.h>   // for pow
#include <stdlib.h> //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include "reading_input.h"
#include "Lanczos_engine.h"

#ifdef USE_COMPLEX
#include "functions_complex.h"
#else
#include "functions_real.h"
#endif

using namespace std;

#ifndef TimeDependentLanczos_class
#define TimeDependentLanczos_class

// template <typename Basis_type, typename Model_type>
class TD_Lanczos
{

public:
    double Time_max;
    double Time_slice_width;
    int No_TimeSlices;

    string State_T0_construction_route;
    string inp_filename;
    int M_;

    string EvolutionType;
    string ConstModelName_str;

    void Constructing_State_T0_route();
    void Constructing_State_T0_via_Hamiltonian();
    void Constructing_State_T0_via_GivenAnsatz();

    void Perform_TD_Lanczos();

    void reading_input();

private:
};

#endif

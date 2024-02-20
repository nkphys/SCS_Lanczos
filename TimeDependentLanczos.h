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
    double Time_otoc;
    double Time_slice_width, dt_;
    int No_TimeSlices;

    int KrylovStepsForTimeEvo;
    int M_EigStates;

    int SiteV, SiteW;

    Mat_1_real Time_, Time_bare;
    Mat_1_real Gamma_, Gamma_bare;
    Mat_1_real Js_, Js_bare;
    bool Use_Scheduler;
    string Scheduler_File;


    string State_T0_construction_route;
    string inp_filename;
    int M_;

    string EvolutionType;
    string ModelName_str;


    Mat_1_doub Psi0;

    void Constructing_InitialState();
    void Constructing_State_T0_via_HamiltonianGS();
    void Constructing_State_T0_via_GivenAnsatz();

    void Perform_TD_Lanczos();
    void Perform_TimeEvolution(Mat_1_doub &Psi_initial, Mat_1_doub &Psi_final, int time_slice_init, int time_slice_final);
    void Calculate_OTOC(string opr_type, int siteV, int siteW);
    void reading_input();
    void Create_Scheduler();

    void Read_Scheduler();


private:
};

#endif

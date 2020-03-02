#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
using namespace std;

//void reading_input(string filename, double & J1_p, double & J1_m, double & J1_z, double & J2_p, double & J2_m, double & J_z, double & H_mag, double & T_Sz, int & m_infinite);

void reading_input(string filename, string & Model_Name, bool &Do_Dynamics, bool& Static_Finite_Temp, bool &Perform_RIXS, bool &Restricted_Basis);


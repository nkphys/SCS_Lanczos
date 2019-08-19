/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_Spins.h"
using namespace std;


void BASIS_Spins::Construct_basis(){

    //NOTES:
    /*
    B is Base. We choose B=2S+1.
    D = \sum_{i=0}^{L-1} V_{i} B^{i}, where V_{i} \in {0,1,2,...,B}
    */
    //---------------------------------//

    SPIN = ((1.0*TwoTimesSpin)/2.0);
    BASE = TwoTimesSpin + 1; //(2S+1)
    D_min = 0;
    D_max = pow(BASE,Length) - 1;
    basis_size = D_max + 1;


    cout<<"Total basis = "<< basis_size<<endl;

}



void BASIS_Spins::clear(){
    D_basis.clear();
}

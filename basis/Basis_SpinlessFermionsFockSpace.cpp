/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_SpinlessFermionsFockSpace.h"
using namespace std;


void BASIS_SpinlessFermionsFockSpace::Construct_basis(){
    unsigned long long int d_max,d_min;

    //Calculating min and max decimal
    d_min=0; //No fermions
    d_max=(int)(pow(2,Length)+0.5)-1;; //Full filling i.e. N_fermions = Length

    DMax_ = d_max;
    DMin_ = d_min;
    basis_size = d_max +1;
    Base=2;

    //putting correct D_'s in the D arrays
    
    cout<<"Hilbert space size = "<<basis_size<<endl;

//----------makes finding basis index faster----------------------------------//


}


void BASIS_SpinlessFermionsFockSpace::clear(){

}

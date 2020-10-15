/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_Spins_Target_Sz.h"
using namespace std;


void BASIS_Spins_Target_Sz::Construct_basis(){

    //NOTES:
    /*
    B is Base. We choose B=2S+1.
    D = \sum_{i=0}^{L-1} V_{i} B^{i}, where V_{i} \in {0,1,2,...,B-1}
    */
    //---------------------------------//

    bool Basis_direclty_in_TotalSz=true;

    SPIN = ((1.0*TwoTimesSpin)/2.0);
    BASE = TwoTimesSpin + 1; //(2S+1)
    D_min = 0;
    D_max = pow(BASE,Length) - 1;
    Target_Total_Value = int(Target_Total_Sz + (SPIN*Length));

    double check;
    check=(Target_Total_Value*1.0) - (Target_Total_Sz + (SPIN*Length));
    if(abs(check)>0.0001){
        cout<<"Some issue in targetting Total Sz, not mapping to integer total value"<<endl;
        assert(false);
    }


    int Value;
    D_basis.clear();
    Mat_1_int Value_local;
    Value_local.resize(Length);

    unsigned long long int d;



    if(Basis_direclty_in_TotalSz){
        cout<<"Basis are constructed directly in the targetted Sz sector"<<endl;

        int BASE_new = Length+1;
        int Length_new = BASE;
        int index_max = pow(BASE_new,Length_new) - 1;
        int inputNum = index_max;
        Mat_1_int Nvec;
        Nvec.resize(Length_new);
        for(int l=0;l<Length_new;l++){
            Nvec[l]=0;
        }
        int sum_N, sum_Nl, l_;
        int n_min, n_max;
        unsigned long long int Dec_temp;
        Mat_1_ullint Dec_vec;
        char state_str[500];
        int length_state_str;
        int old_size;

        /*
     NOTE following conditions must be satisfied:
     (1) \sum_{i=0}^{Length-1} n_{i}.i = Total_Value
     (2) \sum_{j=0}^{BASE-1} n_{i} = Length
     n_{i} is Number os sites with value = "i"

    We define:
    index= \sum_{i=0}^{Length_new-1}n_{i}(BASE_new)^{i}
    Basically, we see it as system with "BASE" number of sites, where each site can have n_{i} value
    */

        for(int index=0;index<=index_max;index++){

            inputNum=index;
            l_=0;
            sum_N=0;
            sum_Nl=0;

            Nvec.clear();
            Nvec.resize(Length_new);
            for(int l=0;l<Length_new;l++){
                Nvec[l]=0;
            }
            while (inputNum > 0)
            {
                assert(l_<Length_new);
                Nvec[l_]= inputNum%BASE_new;
                inputNum = inputNum/BASE_new;
                sum_N +=Nvec[l_];
                sum_Nl +=Nvec[l_]*l_;
                l_++;
            }


            if((sum_N==Length) && (sum_Nl==Target_Total_Value)){

                Dec_temp=0;
                n_min=0;
                n_max=0;
                for(int i_=0;i_<BASE;i_++){
                    n_max +=Nvec[i_];
                    for(int site=n_min;site<n_max;site++){
                      Dec_temp +=i_*pow(BASE,site);
                    }
                    n_min +=Nvec[i_];
                }

                //Now all permutations of Dec_temp are required
                Dec_vec.clear();
                fromDeci(state_str, BASE, Dec_temp);
                length_state_str = strlen(state_str);
                findPermutations(state_str, 0, length_state_str, Dec_vec, BASE);

                old_size = D_basis.size();
                D_basis.resize(D_basis.size() + Dec_vec.size());
                for(int x=0;x<Dec_vec.size();x++){
                    D_basis[old_size + x] = Dec_vec[x];
                }

                cout<<"basis under construction : "<<D_basis.size()<<endl;
                vector < unsigned long long int >().swap(Dec_vec);

                //cout<<"here 1"<<endl;

            }



        }


    }
    else{
        for(d=D_min;d<=D_max;d++){

            //Value=0;
            //        for(int site=0;site<Length;site++){
            //        Value += value_at_pos(d, site, BASE);
            //        }

            Value =Sum_of_Values(d, BASE);

            if(Value==Target_Total_Value){
                D_basis.push_back(d);

                // if((D_basis.size()%100000)==0){
                if((d%1000000)==0){
                    cout<<"Basis creation in progress : "<<((d*1.0)/(1.0*D_max))*100.0<< "% checked, d="<<d<<endl;
                }

            }
        }
    }



    basis_size = (unsigned long long int)D_basis.size();


    cout<<endl;
    cout<<"-----------------------------------------------"<<endl;
    cout<<"Total basis = "<< D_max + (unsigned long long int)1<<endl;
    cout<<"In the Sz = "<<Target_Total_Sz<<" sector : "<<basis_size<<endl;
    cout<<"-----------------------------------------------"<<endl<<endl;


}



void BASIS_Spins_Target_Sz::clear(){
    D_basis.clear();
}

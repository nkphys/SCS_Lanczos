/*
This class includes the Basis used for Model for which Lanczos is being done
*/
#include <iostream>
#include <math.h>
#include "Basis_Spins_Target_Sz_and_K.h"
using namespace std;
#define PI 3.14159265

void BASIS_Spins_Target_Sz_and_K::Check_Commutation_and_Inversion_bw_LatticeSymmetry_Oprs(){

    cout<<"No. of Transformations : "<<No_of_symms_trans<<endl;

    for(int i_opr=0;i_opr<No_of_symms_trans;i_opr++){
        cout<<"----------Transformation no. "<<i_opr<<"----------"<<endl;
        for(int i=0;i<Length;i++){
            for(int j=0;j<Length;j++){
                cout<<    Transformation_matrices[i_opr][i][j]<<" ";
                }
        cout<<endl;
        }
        cout<<"-------------------------------------------"<<endl;
    }

    int CommutationMeasure;
    int AB_ij, BA_ij;
    cout<<"------Commutation matrix (b/w lattice transformation given)------"<<endl;
    for(int i_opr=0;i_opr<No_of_symms_trans;i_opr++){
        for(int j_opr=0;j_opr<No_of_symms_trans;j_opr++){

            CommutationMeasure=0;
            for(int i=0;i<Length;i++){
                for(int j=0;j<Length;j++){
                    AB_ij=0;BA_ij=0;
                    for(int k=0;k<Length;k++){
                    AB_ij += Transformation_matrices[i_opr][i][k]*
                             Transformation_matrices[j_opr][k][j];
                    BA_ij += Transformation_matrices[j_opr][i][k]*
                             Transformation_matrices[i_opr][k][j];
                    }
                    CommutationMeasure += abs(AB_ij - BA_ij);
                }
            }
            cout<<CommutationMeasure<<" ";
        }
        cout<<endl;
    }

    cout<<"--------------------------------------"<<endl;



    cout<<"---------Checking Inversion---------------"<<endl;
    for(int i_opr=0;i_opr<No_of_symms_trans;i_opr++){

        cout<<"-------For Transformation no. "<<i_opr<<"-----------"<<endl;
        for(int i=0;i<Length;i++){
            for(int j=0;j<Length;j++){
                AB_ij=0;

                for(int k=0;k<Length;k++){
                AB_ij +=Transformation_matrices[i_opr][i][k]*InverseTransformation_matrices[i_opr][k][j];
                }

                cout<<AB_ij<<" ";
            }
            cout<<endl;
        }

        cout<<"------------------------------------------"<<endl;
    }


}

void BASIS_Spins_Target_Sz_and_K::Construct_basis(){

    //NOTES:
    /*
    B is Base. We choose B=2S+1.
    D = \sum_{i=0}^{L-1} V_{i} B^{i}, where V_{i} \in {0,1,2,...,B-1}
    */
    //---------------------------------//

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


    //Creating bi-partition A and B
    //[A][B]
    //[A]=[ (0) (1) .... (Lp-1) ], [B]=[ (Lp) (Lp + 1) .... (L-1) ]

    Mat_1_ullint Dec_vec;
    Mat_1_int Vec_temp;
    ulli Dec_temp;
    int No_of_partitions=2;
    char state_str[500];
    int length_state_str;
    int Lp = Length/2;
    Length_A = Lp;
    Length_B = Length - Lp;

    Norm_vals_basis.clear();
    R_vals_basis.clear();
    int ValueB_, ValueA_;
    ulli d_tempAB;
    Mat_1_int VecOut;
    Mat_1_int R_vals;
    R_vals.resize(No_of_symms_trans);


    double Prod_cardnlts=1;
    for(int tn_type=0;tn_type<No_of_symms_trans;tn_type++){
        Prod_cardnlts = Prod_cardnlts*Transformation_Cardinalities[tn_type];
    }



    //int basis_count=0;
    ulli Dmin_partA_ = 0;
    ulli Dmax_partA_ = pow(BASE,Length_A) - 1;


    int min_n_up_in_partB = max(0, Target_Total_Value - Length_A);
    ulli Dmin_partB_ = pow(BASE, min_n_up_in_partB)-1;
    ulli Dmax_partB_ = pow(BASE,Length_B) - 1;

    MainIndex_to_Dec_partA_.clear();
    MainIndex_to_Dec_partB_.clear();
    DecB_to_Index_partB_.resize(Dmax_partB_+1);
    DecB_to_Sections.resize(Dmax_partB_+1);
    for(int i=0;i<DecB_to_Index_partB_.size();i++){
        DecB_to_Index_partB_[i]=-100;
        DecB_to_Sections[i]=-100;
    }

    int basis_count=0;
    //int indexB_count=0;
    int Sections_count=-1;


    int indexA_, indexB_;

    ulli d_tempB_old=Dmin_partB_;
    int indexA_old;
    indexA_old=-1;
    indexB_=0;

    for(ulli d_tempB=Dmin_partB_;d_tempB<=Dmax_partB_;d_tempB++){
    //ValueB_ =Sum_of_Values(d_tempB, BASE);

        ValueB_ =Sum_of_Values(d_tempB, BASE);
        if(ValueB_ >= min_n_up_in_partB ){
            Dec_vec.clear();
            Vec_temp.resize(Length_A);
            ValueA_ = Target_Total_Value - ValueB_;
            assert(ValueA_<=Vec_temp.size());
            for(int i=0;i<Vec_temp.size();i++){
                if(i<ValueA_){
                    Vec_temp[i]=1;
                }
                else{
                    Vec_temp[i]=0;
                }
            }

            fromVecint_to_Deci(Vec_temp, BASE, Dec_temp, Length_A);
            fromDeci(state_str, BASE, Dec_temp, Length_A);
            length_state_str = strlen(state_str);
            findPermutations(state_str, 0, length_state_str, Dec_vec, BASE);




    indexA_=0;

    for(int d_tempA_index=0;d_tempA_index<Dec_vec.size();d_tempA_index++){
    ulli d_tempA = Dec_vec[d_tempA_index];
    assert(d_tempA<=Dmax_partA_);

    //for(ulli d_tempA=Dmin_partA_;d_tempA<=Dmax_partA_;d_tempA++){
    //ValueA_ =Sum_of_Values(d_tempA, BASE);


    //Here
    double_type Gamma_;
    int Repetitions_n;
    bool new_allowed_representative = CheckState(d_tempA, d_tempB, Gamma_, Repetitions_n);


    if(new_allowed_representative){

    basis_count +=1;
    MainIndex_to_Dec_partB_.push_back(d_tempB);
    MainIndex_to_Dec_partA_.push_back(d_tempA);

    double norm_ = abs(Gamma_)*abs(Gamma_)*(Prod_cardnlts/(1.0*Repetitions_n));
    Norm_vals_basis.push_back(norm_);

    //R_vals_basis.push_back(R_vals);

    //DecB_to_Sections[d_tempB] = Sections_count;
    //DecB_to_Index_partB_[d_tempB]=indexB_count;

    if(d_tempB!=d_tempB_old ||
        (MainIndex_to_Dec_partB_.size()==1) ){
        //Index_to_Dec_part1_.push_back(d_temp1);
        //cout<<"HERE 1"<<endl;
        DecB_to_Index_partB_[d_tempB]=indexA_old+indexB_+1;
        indexB_ = DecB_to_Index_partB_[d_tempB];
        d_tempB_old = d_tempB;

        Sections_count +=1;
        DecB_to_Sections[d_tempB] = Sections_count;

    }
    //Index_to_Dec_part0_.push_back(d_temp0);

    Index_partA_to_DecA_in_givenSection.resize(Sections_count+1);
    Index_partA_to_DecA_in_givenSection[Sections_count].push_back(d_tempA);
   // Dec_to_Index_partA_[d_tempA]=indexA_;

    assert((indexA_ + indexB_+1) == MainIndex_to_Dec_partB_.size());
    indexA_old = indexA_;
    indexA_ +=1;


/*
    basis_count +=1;
    d_tempAB = d_tempA + pow(BASE,Length_A)*d_tempB;
    fromDeci_to_Vecint(VecOut, BASE, d_tempAB, Length);
    cout<<d_tempAB<<" : ";
    for(int site=0;site<Length;site++){
        cout<<VecOut[site]<<" ";
    }
    cout<<" [ ";
    for(int tn=0;tn<R_vals.size();tn++){
        cout<<R_vals[tn]<<" ";
    }
    cout<<"]"<<endl;
*/

    }

    }

    //indexB_count +=1;

    if(d_tempB%(Dmax_partB_/1) == 0){
    cout<<"Basis fraction done :"<<((1.0*d_tempB)/(1.0*Dmax_partB_))<<endl;
    }
    }
    }

  //  R_vals_basis[4][1]=1;
    cout<<"No. of Basis = "<<MainIndex_to_Dec_partA_.size()<<endl;
    basis_size = MainIndex_to_Dec_partA_.size();

}


bool BASIS_Spins_Target_Sz_and_K::CheckState(ulli d_tempA, ulli d_tempB, double_type & Gamma_, int & Repetitions_n){


    // R_vals.clear();
    // R_vals.resize(No_of_symms_trans);
    // for(int tn=0;tn<R_vals.size();tn++){
    //     R_vals[tn]=Transformation_Cardinalities[tn];
    // }


    //-------------------------

    //------------------------

    Repetitions_n=0;

    //int distinct_vecs;

    Mat_1_ullint d_collected;

    Mat_1_int zeros_set;
    zeros_set.resize(No_of_symms_trans);
    for(int i=0;i<No_of_symms_trans;i++){
        zeros_set[i]=0;
    }
    int Total_no_elements,Total_no_elements_expected;

    bool new_allowed_representative_check=true;
    ulli d_tempAB;
    ulli d_AB_out;
    int ValueB_ =Sum_of_Values(d_tempB, BASE);
    int ValueA_ =Sum_of_Values(d_tempA, BASE);
    Mat_1_int Vec_AB;

    if( (ValueB_ + ValueA_) != Target_Total_Value ){
        new_allowed_representative_check=false;
        }
    else{
        d_tempAB = d_tempA + pow(2,Length_A)*d_tempB;
        //All_vecs.push_back(d_tempAB);

        fromDeci_to_Vecint(Vec_AB, BASE, d_tempAB, Length);
        Mat_1_int Vec_in, Vec_out;
        Vec_in = Vec_AB;

        Mat_2_int IndicesSet;
        Total_no_elements=0;
        Total_no_elements_expected=1;
        for(int i=0;i<Transformation_Cardinalities.size();i++){
            Total_no_elements_expected *=Transformation_Cardinalities[i];
        }

        DirectProduct_IndicesSet_by_recursion(zeros_set, Transformation_Cardinalities, Total_no_elements , IndicesSet);
        assert(Total_no_elements_expected == Total_no_elements);
        assert(IndicesSet.size()==Total_no_elements_expected);


        Mat_1_int Diff;Diff.resize(No_of_symms_trans);
        Mat_1_int Indices_old;
        for(int i=0;i<No_of_symms_trans;i++){
            Indices_old.push_back(0);
        }

        d_collected.clear();
        int counter=0;
        Gamma_=0.0;
        for(int n=0;n<Total_no_elements;n++){


        for(int i=0;i<No_of_symms_trans;i++){
        Diff[i]=IndicesSet[n][i] - Indices_old[i];
        }

        for(int tn_type=(No_of_symms_trans-1);tn_type>=0;tn_type--){

            if(Diff[tn_type]>0){
            for(int m=0;m<Diff[tn_type];m++){
            Matrix_vector_multiplication(Transformation_matrices[tn_type], Vec_in, Vec_out);
            Vec_in=Vec_out;
            }
            }
            if(Diff[tn_type]<0){
            for(int m=0;m<abs(Diff[tn_type]);m++){
            Matrix_vector_multiplication(InverseTransformation_matrices[tn_type], Vec_in, Vec_out);
            Vec_in=Vec_out;
            }
            }
            if(Diff[tn_type]==0){
            Vec_out=Vec_in;
            }

        }


        //d_AB_new
        fromVecint_to_Deci(Vec_out, BASE, d_AB_out, Length);

        //All_vecs.push_back(d_AB_out);


        if(d_AB_out<d_tempAB){
            new_allowed_representative_check=false;
            break;
        }

        if(d_AB_out==d_tempAB){
            Repetitions_n=Repetitions_n+1;
            double_type gamma_local=1.0;
            complex<double> exp_char;
            for(int tn_type=0;tn_type<No_of_symms_trans;tn_type++){
                exp_char = exp(iota_comp*2.0*PI*(1.0*EigenvaluesTargeted[tn_type]*IndicesSet[n][tn_type])*(1.0/(1.0*Transformation_Cardinalities[tn_type])));

#ifdef USE_COMPLEX
            gamma_local = gamma_local*exp_char;
#endif

#ifndef USE_COMPLEX
            gamma_local = gamma_local*exp_char.real();
#endif
            }

            Gamma_ += gamma_local;
            // int no_of_nonzeros=0;
            // int tn_interest;
            // for(int tn_type=0;tn_type<No_of_symms_trans;tn_type++){
            //     if(IndicesSet[n][tn_type]!=0){
            //         no_of_nonzeros +=1;
            //         tn_interest=tn_type;
            //     }
            // }

            // if(no_of_nonzeros==1){
            //     if(IndicesSet[n][tn_interest]<R_vals[tn_interest]){
            //     R_vals[tn_interest]=IndicesSet[n][tn_interest];
            //     }}

            //break; //on first repetition
        }

        // bool found_check=false;
        // int trash_counter;
        // if(d_collected.size()>=1){
        // trash_counter = Find_int_in_intarray(found_check, d_AB_out, d_collected);
        // }
        // if(found_check){
        //     break;
        // }
        // d_collected.push_back(d_AB_out);



        for(int i=0;i<No_of_symms_trans;i++){
            Indices_old[i]=IndicesSet[n][i];
        }
        Vec_in=Vec_out;
        }

        //Matrix_vector_multiplication(Transformation_matrices[i_opr], Vec_in, Vec_out);


         // for(int tn_type=0;tn_type<No_of_symms_trans;tn_type++){
         // int mod_tn = EigenvaluesTargeted[tn_type]%
         //            (Transformation_Cardinalities[tn_type]/R_vals[tn_type]);
         // if(mod_tn!=0){
         //     new_allowed_representative_check=false;
         // }
         // }

        if(abs(Gamma_)<0.0000001){
            new_allowed_representative_check=false;
        }

        }



        // if(new_allowed_representative_check){
        // distinct_vecs = Get_distinct_int_numbers(All_vecs);
        // }



        return new_allowed_representative_check;
}

void BASIS_Spins_Target_Sz_and_K::clear(){
    D_basis.clear();
}


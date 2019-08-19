#include "tensor_type.h"
#include <assert.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include "mkl_spblas.h"
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
//#include <bits/stdc++.h>


int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2,
                                 int num1_sector, Mat_1_intpair &sectors_offset);
int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2);
int Find_int_in_intarray(int num, Mat_1_int &array);
void Print_Matrix_COO(Matrix_COO &A);
void Print_vector_in_file(Mat_1_doub vec, string filename);
void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG);
void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO &A,Mat_1_doub &u,Mat_1_doub &v);
double dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2);
void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_1_doub & vecG, int lanc_iter);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG,
                       Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_2_doub & vecG,int lanc_iter ,int few_);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_2_doub & vecG,int lanc_iter ,int few_, Mat_1_int states_to_look);
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2);
void Calculate_recursive_GF(Mat_1_doub A, Mat_1_doub B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy);
bool comp_greater(double i, double j);
bool comp_greater_pair_double_int(pair_double_int i, pair_double_int j);
static bool sort_using_greater_than(double u, double v);
string NumberToString ( int Number );
bool present_before(Mat_1_int nup_2, Mat_1_int ndn_2, Mat_2_int nup_2_group, Mat_2_int ndn_2_group, int &pos);


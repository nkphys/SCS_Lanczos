#ifdef USE_COMPLEX
#include "tensor_type.h"
#include <assert.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <math.h>  //fabs(double x) =|x|
#include <stdlib.h>  //for div(q,n).rem(quot),rand
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <stdio.h>
//#include "mkl_spblas.h"
//#include <mkl_types.h>
//#include <mkl_cblas.h>
//#include <mkl_lapacke.h>
//#include <bits/stdc++.h>
#include "Matrix.h"

void Perform_SVD(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);
int partition(Mat_1_int &a, int s, int e);
void quicksort(Mat_1_int &a, int s, int e);
void Remove_repetitions(Mat_1_ullint & index_array, Mat_1_doub & val_array, Mat_1_ullint & index_new_array, Mat_1_doub & val_new_array);
void Remove_repetitions(Mat_1_int & index_array, Mat_1_doub & val_array, Mat_1_int & index_new_array, Mat_1_doub & val_new_array);
complex<double> conjugate(complex<double> val);
void Normalize_vec(Mat_1_doub &vec_in);
void value_multiply_vector(complex<double> value, Mat_1_doub &vec_in);
void Direct_product_of_Mat_2_trio_int(Mat_2_trio_int MAT1_, Mat_1_doub values1_,
                                      Mat_2_trio_int MAT2_, Mat_1_doub values2_,
                                      Mat_2_trio_int &MAT_RESULT_, Mat_1_doub &values_result_);
void reading_input_dos_trio(string inp_filename, Mat_1_trio_int &TRIO_VEC, Mat_1_doub &values_ );
bool Is_int_in_array(int num, Mat_1_int array);
int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2,
                                 int num1_sector, Mat_1_intpair &sectors_offset);
int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2);
int Find_int_in_intarray(int num, Mat_1_int &array);
int Find_int_in_intarray(ulli num, Mat_1_ullint &array);
int Find_int_in_intarray(bool &found, ulli num, Mat_1_ullint &array);
int Find_int_in_intarray_using_bisection_(bool &found_, ulli num, Mat_1_ullint &array);
int Find_int_in_intarray_using_bisection(ulli num, Mat_1_ullint &array);
int Find_int_in_intarray_smartly(int num,Mat_1_int &array,Mat_1_int &partition_indices,Mat_1_int &vals_at_partitions);
int Find_int_in_part_of_intarray(int num, Mat_1_int &array, int min_i, int max_i);
int Find_int_in_part_of_intarray(ulli num, Mat_1_ullint &array, int min_i, int max_i);
void Print_Matrix_COO(Matrix_COO &A);
void Print_Matrix_COO_format(Matrix_COO &A);
void Print_Matrix_COO_in_file(Matrix_COO &A, string filename);
void Print_Matrix(Mat_2_doub &A);
void Print_Matrix(Mat_2_real &A);
void Print_vector_in_file(Mat_1_doub vec, string filename);
void Sort_vector_in_decreasing_order_in_file(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old);
void Sort_vector_in_decreasing_order(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old, int n_basis);
void Sort_vector_in_decreasing_order(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old, Mat_1_int &Index_new, int n_basis);
void Print_file_in_vector(Mat_1_doub &vec, string filename , int rows);
void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG);
void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_1_doub & vecG);
void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_2_doub & vecs);
void Matrix_vector_multiplication(Mat_2_int &A, Mat_1_int &u, Mat_1_int &v);
void DirectProduct_IndicesSet_by_recursion(Mat_1_int &TempState, Mat_1_int &VecSizes, int & Total_no_elements ,Mat_2_int &IndicesSet);
void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO &A,Mat_1_doub &u,Mat_1_doub &v);
Matrix_COO Dagger(Matrix_COO &A);
int minSwaps(vector<int> &arr, int n);
complex<double> dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2);
double Norm(Mat_1_doub &vec1);
void Subtract( Mat_1_doub &temp1, double_type x, Mat_1_doub &temp2);
void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_1_doub & vecG, int lanc_iter);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG,
                       Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_2_doub & vecG,int lanc_iter ,int few_);
void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_2_doub & vecG,int lanc_iter ,int few_, Mat_1_int states_to_look);
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2);
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, complex<double> value1, complex<double> value2);
void Calculate_recursive_GF(Mat_1_doub A, Mat_1_real B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy);
bool comp_greater(double i, double j);
bool comp_greater_pair_double_int(pair_real_int i, pair_real_int j);
static bool sort_using_greater_than(double u, double v);
string NumberToString ( int Number );
bool present_before(Mat_1_int nup_2, Mat_1_int ndn_2, Mat_2_int nup_2_group, Mat_2_int ndn_2_group, int &pos);

complex<double> divide(complex<double> z1, complex<double> z2);
void Read_matrix_from_file(string filepath,
                         Mat_2_doub &Mat, int row, int column);
complex<double> reading_pair(string pair_str);
Matrix_COO Identity_COO(int rows_no, int cols_no);
int Find_commont_int(Mat_1_int Vec1, Mat_1_int Vec2);
double Lorentzian(double eta, double x);
#endif

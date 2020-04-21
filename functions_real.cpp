#include "functions_real.h"
#define PI_ 3.14159265


//extern "C" void zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
//                         std::complex<double> *,int *, double *, int *);
extern "C" void dsyev_(char * , char * , int * , double * , int *, double *, double *, int *, int *);

#ifndef USE_COMPLEX


double conjugate(double val){
    return (val);
}

double Lorentzian(double eta, double x){
    double val;
    val = (1.0/PI_)*( (eta/2.0) / (  (x*x) + ((eta*eta)/4)   ) );
    return val;
}

double reading_pair(string double_no_str){

    double value;
    value = atof(double_no_str.c_str());
    return value;
}

void Normalize_vec(Mat_1_doub &vec_in){

    double val;

    val = dot_product(vec_in, vec_in);
    val = sqrt(val);
    for(int i=0;i<vec_in.size();i++){
        vec_in[i] = vec_in[i]*(1.0/val);
    }

}

void value_multiply_vector(double value, Mat_1_doub &vec_in){

    for(int i=0;i<vec_in.size();i++){
    vec_in[i] = vec_in[i]*value;
    }
}

void reading_input_dos_trio(string inp_filename, Mat_1_trio_int &TRIO_VEC, Mat_1_doub &values_ ){


    TRIO_VEC.clear();
    values_.clear();

    string filepath = inp_filename;

    string Dynamics_dos_operators;
    string Dynamics_Dos_Operators_ = "Dynamics_dos_operators = ";
    string no_of_fermionic_oprs_string;
    string temp;
    double temp_double;
    int no_of_fermionic_oprs;


    int offset;
    string line;
    ifstream inputfile(filepath.c_str());


    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            getline(inputfile,line);

            if ((offset = line.find(Dynamics_Dos_Operators_, 0)) != string::npos) {
                Dynamics_dos_operators = line.substr (offset+Dynamics_Dos_Operators_.length());	}

        }
        inputfile.close();
    }
    else
    {cout<<"Unable to open input file while in the functions_real.cpp"<<endl;}

    stringstream Dynamics_dos_operators_(Dynamics_dos_operators);
    Dynamics_dos_operators_ >> no_of_fermionic_oprs_string;

    no_of_fermionic_oprs=atoi(no_of_fermionic_oprs_string.c_str());

    TRIO_VEC.resize(no_of_fermionic_oprs);
    values_.resize(no_of_fermionic_oprs);

    for(int i=0;i<no_of_fermionic_oprs;i++){

        //real part
        Dynamics_dos_operators_ >> temp;
        values_[i]=atof(temp.c_str());

        //orbital
        Dynamics_dos_operators_ >> temp;
        TRIO_VEC[i].orb_=atoi(temp.c_str());

        //spin
        Dynamics_dos_operators_ >> temp;
        TRIO_VEC[i].spin_=atoi(temp.c_str());

        //site
        Dynamics_dos_operators_ >> temp;
        TRIO_VEC[i].site_=atoi(temp.c_str());

    }


}

void swap(int &n1, int &n2)
{
    int temp=n1;
    n1=n2;
    n2=temp;
}

void Print_Matrix(Mat_2_doub &A){

    for(int i=0;i<A.size();i++){
        for(int j=0;j<A.size();j++){
            cout<<A[i][j]<<"  ";
        }
        cout<<endl;
    }


}

void Read_matrix_from_file(string filepath, Mat_2_doub &Mat, int row, int column){

    ifstream infile(filepath.c_str());

    if(!infile.is_open()){
        cout<<filepath<<" not present"<<endl;
        }

    Mat.resize(row);
    for(int i=0;i<row;i++){
        Mat[i].resize(column);
    }

    double tmp_doub;

    for (int i =0;i<row;i++){
        for (int j =0;j<column;j++){

            infile>>tmp_doub;
            Mat[i][j]=tmp_doub;
        }
    }

}

Matrix_COO Dagger(Matrix_COO &A){

    Matrix_COO B;
    B.nrows=A.nrows;
    B.ncols = A.ncols;
    int i=0;
    B.value.resize(A.value.size());
    B.rows.resize(A.rows.size());
    B.columns.resize(A.columns.size());
    while(i<A.value.size()){

        B.value[i] = A.value[i];
        B.rows[i] = A.columns[i];
        B.columns[i] = A.rows[i];

        i=i+1;
    }
    return B;
}

void Direct_product_of_Mat_2_trio_int(Mat_2_trio_int MAT1_, Mat_1_doub values1_,
                                      Mat_2_trio_int MAT2_, Mat_1_doub values2_,
                                      Mat_2_trio_int &MAT_RESULT_, Mat_1_doub &values_result_){

    MAT_RESULT_.clear();
    values_result_.clear();
    MAT_RESULT_.resize(MAT1_.size()*MAT2_.size());
    values_result_.resize(MAT1_.size()*MAT2_.size());
    for(int i=0;i<MAT1_.size();i++){
        for(int j=0;j<MAT2_.size();j++){

            for(int n=0;n<MAT1_[i].size();n++){
                MAT_RESULT_[(i*MAT2_.size()) + j].push_back(MAT1_[i][n]);
            }
            for(int n=0;n<MAT2_[j].size();n++){
                MAT_RESULT_[(i*MAT2_.size()) + j].push_back(MAT2_[j][n]);
            }

            assert(MAT_RESULT_[(i*MAT2_.size()) + j].size() == MAT1_[i].size()+MAT2_[j].size());

            values_result_[(i*MAT2_.size()) + j] = values1_[i]*values2_[j];

        }

    }


}

int minSwaps(vector<int> &arr, int n)
{
    // Create an array of pairs where first
    // element is array element and second element
    // is position of first element
    pair<int, int> arrPos[n];
    for (int i = 0; i < n; i++)
    {
        arrPos[i].first = arr[i];
        arrPos[i].second = i;
    }

    // Sort the array by array element values to
    // get right position of every element as second
    // element of pair.
    sort(arrPos, arrPos + n);

    // To keep track of visited elements. Initialize
    // all elements as not visited or false.
    vector<bool> vis(n, false);

    // Initialize result
    int ans = 0;

    // Traverse array elements
    for (int i = 0; i < n; i++)
    {
        // already swapped and corrected or
        // already present at correct pos
        if (vis[i] || arrPos[i].second == i)
            continue;

        // find out the number of  node in
        // this cycle and add in ans
        int cycle_size = 0;
        int j = i;
        while (!vis[j])
        {
            vis[j] = 1;

            // move to next node
            j = arrPos[j].second;
            cycle_size++;
        }

        // Update answer by adding current cycle.
        if(cycle_size > 0)
        {
            ans += (cycle_size - 1);
        }
    }

    // Return result
    for(int i=0;i<n;i++){
        arr[i]=arrPos[i].first;
    }
    return ans;
}


bool Is_int_in_array(int num, Mat_1_int array){

    bool check=false;
    for (int i = 0; i<array.size(); i++){
        if(array[i]==num){
            check=true;
            break;
        }
    }

    return check;
}

int Find_commont_int(Mat_1_int Vec1, Mat_1_int Vec2){

    int temp;
    bool found =false;
    for (int i = 0; i<Vec1.size(); i++){
        if(!found){
            for(int j = 0; j<Vec2.size(); j++){

                if (Vec1[i] == Vec2[j]){
                    temp = Vec1[i];
                    found=true;
                }

            }
        }
        else{
            break;
        }
    }

    return temp;

}

void Print_vector_in_file(Mat_1_doub vec, string filename){

    ofstream outfile(filename.c_str());

    for(int j=0;j<vec.size();j++){
        outfile<<vec[j]<<endl;
    }

}

void Print_file_in_vector(Mat_1_doub &vec, string filename, int rows){

    cout<<"Vector is being read from : "<<filename<<endl;
    vec.clear();
    ifstream infile(filename.c_str());
    //string line;
    double temp_double;


    for(int j=0;j<rows;j++){
        infile>>temp_double;
        vec.push_back(temp_double);
        //cout<<j<<"  "<<Kvector_n[j]<<endl;
    }

    assert(vec.size()>0);

}


bool present_before(Mat_1_int nup_2, Mat_1_int ndn_2, Mat_2_int nup_2_group, Mat_2_int ndn_2_group, int &pos){


    bool check;
    check=false;


    if(nup_2_group.size()==0){
        pos=0;
        check=false;
    }

    else{
        bool c_up, c_dn, c_both;
        c_up=false;c_dn=false;

        for(int i=0;i<nup_2_group.size();i++){

            if(nup_2_group[i]==nup_2){
                c_up=true;
            }
            else{
                c_up=false;
            }
            if(ndn_2_group[i]==ndn_2){
                c_dn=true;
            }
            else{
                c_dn=false;
            }


            c_both= (c_up && c_dn) ;
            if(c_both==true){
                check=true;
                pos=i;
            }



        }
    }
    return check;


}



string NumberToString ( int Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

static bool sort_using_greater_than(double u, double v)
{
    return u > v;
}

bool comp_greater(double i, double j){
    return i > j;
}


bool comp_greater_pair_double_int(pair_real_int i, pair_real_int j){
    return i.first > j.first;
}


int Find_int_in_part_of_intarray(int num, Mat_1_int &array, int min_i, int max_i){

    int pos;
    bool not_found=true;

    for(int ind=min_i;ind<=max_i;ind++){
        if(num==array[ind]){
            not_found=false;
            pos=ind;
            break;
        }
    }

    if(not_found){
        pos=-1;
    }
    return pos;
}

int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2,
                                 int num1_sector, Mat_1_intpair &sectors_offset){


    int pos;
    bool not_found=true;
    int ind, ind_max;

    assert(num1_sector<sectors_offset.size());
    ind=sectors_offset[num1_sector].first;
    ind_max=sectors_offset[num1_sector].second;
    while( (not_found) && (ind<=ind_max) ){

        assert(ind<array1.size());
        if(num1==array1[ind]  && num2==array2[ind] ){
            pos=ind;
            not_found = false;
        }
        ind++;
    }

    assert(!not_found);
    return pos;



}


int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2){



    int pos;
    bool not_found=true;
    int ind=0;
    while(not_found){
        assert(ind<array1.size());
        if(num1==array1[ind]  && num2==array2[ind] ){
            pos=ind;
            not_found = false;
        }
        ind++;
    }

    return pos;


}

double sign_of_double(double val){

    double val_return;
    if(val==0){
        val_return =0.0;
    }
    else if(val>0.0){
        val_return =1.0;
    }
    else{
        assert(val<0.0);
        val_return=-1.0;
    }

    return val_return;
}

int Find_int_in_intarray(int num, Mat_1_int &array){

    int pos;
    bool not_found=true;
    int ind=0;
    while(not_found){
        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }
        ind++;
    }

    return pos;
}

int Find_int_in_intarray_smartly(int num, Mat_1_int &array, Mat_1_int &partition_indices, Mat_1_int &vals_at_partitions){


    int min_ind;
    for(int i=0;i<(vals_at_partitions.size()-1);i++){
        if(num>=vals_at_partitions[i] && num<=vals_at_partitions[i+1]){
            min_ind=partition_indices[i];
            break;
        }
    }


    int pos;
    bool not_found=true;
    int ind=min_ind;
    while(not_found){
        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }
        ind++;
    }

    return pos;

}

void Print_Matrix_COO(Matrix_COO &A){

    Mat_2_doub B;
    B.resize(A.nrows);
    for(int i=0;i<B.size();i++){
        B[i].resize(A.nrows);
        for(int j=0;j<A.nrows;j++){
            B[i][j]=0.0;
        }
    }

    for(int i=0;i<A.value.size();i++){
        B[A.rows[i]][A.columns[i]]=A.value[i];
    }

    cout<<"--------------------PRINTING THE MATRIX:--------------------"<<endl;
    for(int i=0;i<B.size();i++){
        for(int j=0;j<B.size();j++){

            cout<<B[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<"-------------------------------------------------------------"<<endl;


}

double dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2){
    //This dot_product is parallelized always, create another one with _PARALLELIZE_AT_MATRICES_LEVEL
    double temp;
    double temp1=0;
    assert(vec1.size()==vec2.size());
    // bool Parallelize_dot_product;
    // Parallelize_dot_product=true;

    // if(!Parallelize_dot_product){goto skiploop_143;}
    //#pragma omp parallel for default(shared) reduction(+:temp1)
    //skiploop_143:
    for(int i=0;i<vec1.size();i++){
        temp1 = temp1 + (vec1[i])*(vec2[i]);
    }

    temp = temp1;

    return temp;

}

double Norm(Mat_1_doub &vec1){
    double temp;
    double temp1=zero;

    for(int i=0;i<vec1.size();i++){
        temp1 = temp1 + (vec1[i])*((vec1[i])); //<vec1|vec1>
    }

    temp = temp1;
    return temp;

}

void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO & A,Mat_1_doub &u,Mat_1_doub &v){

    v.clear();
    v.resize(A.nrows);


    for (int i=0;i<v.size();i++){
        v[i]=0;}

    if(COO_type =="U"){
        for (int n=0;n<A.value.size();n++){

            if(A.rows[n]==A.columns[n]){
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
            }
            else{
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
                v[A.columns[n]] = v[A.columns[n]] + u[A.rows[n]]*A.value[n];
            }

        }
    }
    else{
        for (int n=0;n<A.value.size();n++){

            v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];



        }
    }


}



void Sum( Mat_1_doub &temp1, double a1, Mat_1_doub &temp2, double a2){

    assert(temp1.size()==temp2.size());
    for(int k=0;k<temp1.size();k++){

        temp1[k]=a1*temp1[k] + a2*temp2[k];

    }



}


void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2){

    assert(temp1.size()==temp2.size());
    for(int k=0;k<temp1.size();k++){

        temp1[k]=temp1[k] - x*temp2[k];

    }

}


void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_1_doub & vecG,int lanc_iter){


    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.size(),X.size());
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=sqrt(Y2[i]);
                Ham_(j,i)=sqrt(Y2[i]);
            }
            else{
                Ham_(i,j)=0.0;
            }
        }
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    int Target_state=0;

    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(X.size());
    for(int i=0;i<X.size();i++){
        vecG[i] = Ham_(i,Target_state);
    }

    Ham_.clear();

}


void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_2_doub & vecG,int lanc_iter ,int few_){


    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.size(),X.size());
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=sqrt(Y2[i]);
                Ham_(j,i)=sqrt(Y2[i]);
            }
            else{
                Ham_(i,j)=0.0;
            }
        }
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    int Target_state=0;

    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(few_);

    for(int i=0;i<few_;i++){
        vecG[i].resize(X.size());
    }

    for(int i=0;i<few_;i++){
        for(int j=0;j<X.size();j++){
            vecG[i][j] =  Ham_(j,i);
        }
    }

    Ham_.clear();


}


void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_doub & EG, Mat_2_doub & vecG,int lanc_iter ,int few_, Mat_1_int states_to_look){

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.size(),X.size());
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=sqrt(Y2[i]);
                Ham_(j,i)=sqrt(Y2[i]);
            }
            else{
                Ham_(i,j)=0.0;
            }
        }
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    int Target_state=0;

    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(few_);

    for(int i=0;i<few_;i++){
        vecG[i].resize(X.size());
    }

    for(int i=0;i<few_;i++){
        Target_state=states_to_look[i];
        for(int j=0;j<X.size();j++){
            vecG[i][j] = Ham_(j,Target_state);//mat[j*X.size()+Target_state];
        }
    }

    Ham_.clear();

}


void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG){

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.size(),X.size());
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=sqrt(Y2[i]);
                Ham_(j,i)=sqrt(Y2[i]);
            }
            else{
                Ham_(i,j)=0.0;
            }
        }
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}

    int T_s=0;
    int Target_state=0;

    EG=eigs_[T_s];
    eigs_.clear();
    vecG.clear();
    vecG.resize(X.size());
    for(int i=0;i<X.size();i++){
        vecG[i] = Ham_(i,Target_state); //mat[i*X.size()+Target_state];
    }

    Ham_.clear();

}


void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , double & EG, Mat_1_doub & vecG,
                 Mat_2_doub & Unit_E_vecs, Mat_1_real & Evals_Lanczos){



    Evals_Lanczos.clear();
    Evals_Lanczos.resize(X.size());

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.size(),X.size());
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=sqrt(Y2[i]);
                Ham_(j,i)=sqrt(Y2[i]);
            }
            else{
                Ham_(i,j)=0.0;
            }
        }
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    int T_s=0;
    int Target_state=0;
    /*  if(Target_state>X.size()-1){
        T_s=X.size()-1;
    }
    else{
        T_s=Target_state;
    }*/
    EG=eigs_[T_s];

    for(int i=0;i<X.size();i++){
        Evals_Lanczos[i]=eigs_[i];
    }

    eigs_.clear();
    vecG.clear();
    vecG.resize(X.size());
    Unit_E_vecs.clear();
    Unit_E_vecs.resize(X.size());

    for(int i=0;i<X.size();i++){
        Unit_E_vecs[i].resize(X.size());}

    for(int i=0;i<X.size();i++){

        vecG[i] = Ham_(i,Target_state);//mat[i*X.size()+Target_state];

        for(int j_state=0;j_state<X.size();j_state++){

            Unit_E_vecs[j_state][i] = Ham_(i,j_state);//mat[i*X.size()+j_state];

        }
    }

    Ham_.clear();

}

void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_1_doub & vecG){

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.nrows,X.ncols);

    for(int i=0;i<X.nrows;i++){
        for(int j=i;j<X.ncols;j++){
            Ham_(i,j) = 0.0;
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
       Ham_(r,c) = X.value[i];
       Ham_(c,r) = X.value[i];
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    EVALS.resize(X.nrows);
    for(int i=0;i<EVALS.size();i++){
        EVALS[i]=eigs_[i];
    }
    eigs_.clear();
    vecG.clear();
    vecG.resize(X.nrows);
    for(int i=0;i<X.nrows;i++){
        vecG[i] = Ham_(i,0);//mat[i*X.nrows];
    }

    Ham_.clear();
}

void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_2_doub & vecs){

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.nrows,X.ncols);

    for(int i=0;i<X.nrows;i++){
        for(int j=i;j<X.ncols;j++){
            Ham_(i,j) = 0.0;
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
       Ham_(r,c) = X.value[i];
       Ham_(c,r) = X.value[i];
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    EVALS.resize(X.nrows);
    for(int i=0;i<EVALS.size();i++){
        EVALS[i]=eigs_[i];
    }
    eigs_.clear();

    vecs.clear();
    vecs.resize(X.nrows);
    for(int j=0;j<vecs.size();j++){
        vecs[j].resize(vecs.size());
    for(int i=0;i<X.nrows;i++){
        vecs[j][i] =Ham_(i,j);//mat[i*X.nrows];
    }
    }

    Ham_.clear();
}


void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG){

    Mat_1_real eigs_;
    Matrix<double> Ham_;
    Ham_.resize(X.nrows,X.ncols);

    for(int i=0;i<X.nrows;i++){
        for(int j=i;j<X.ncols;j++){
            Ham_(i,j) = 0.0;
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
       Ham_(r,c) = X.value[i];
       Ham_(c,r) = X.value[i];
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<double> work(3);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dsyev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    EG=eigs_[0];
    eigs_.clear();
    vecG.clear();
    vecG.resize(X.nrows);
    for(int i=0;i<X.nrows;i++){
        vecG[i] =Ham_(i,0);//mat[i*X.nrows];
    }

    Ham_.clear();

}



void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, double value1, double value2){

    int a_i=0;
    int b_j=0;
    int row_a_i,col_a_i, row_b_j,col_b_j;
    Matrix_COO temp;

    if((A.nrows == B.nrows) && (A.ncols == B.ncols)){

        temp.nrows = A.nrows;
        temp.ncols = A.ncols;

        temp.value.clear();
        temp.rows.clear();
        temp.columns.clear();

        while(a_i<A.value.size() || b_j<B.value.size()){

            if(a_i<A.value.size()){
                row_a_i=A.rows[a_i];
                col_a_i=A.columns[a_i];
            }
            else{
                row_a_i=A.nrows+1;
                col_a_i=A.ncols+1;
            }

            if(b_j<B.value.size()){
                row_b_j=B.rows[b_j];
                col_b_j=B.columns[b_j];
            }
            else{
                row_b_j=B.nrows+1;
                col_b_j=B.ncols+1;
            }

            //if element of B comes before element of A
            if( (row_b_j < row_a_i)  ||
                    (row_b_j == row_a_i && col_b_j < col_a_i)
                    )
            {

                temp.value.push_back(value2*B.value[b_j]);

                temp.rows.push_back(row_b_j);
                temp.columns.push_back(col_b_j);
                b_j=b_j+1;

            }
            //if element of A comes before element of B
            if( (row_b_j > row_a_i) ||
                    (row_b_j == row_a_i && col_b_j > col_a_i)
                    )

            {
                temp.value.push_back(value1*A.value[a_i]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;

            }
            //if elements of B, A comes at same place
            if(row_b_j==row_a_i && col_b_j==col_a_i){
                temp.value.push_back(value1*A.value[a_i] + value2*B.value[b_j]);
                temp.rows.push_back(row_a_i);
                temp.columns.push_back(col_a_i);
                a_i=a_i+1;
                b_j=b_j+1;
            }
        }

    }
    else{cout<<"Error in doing Sum"<<endl;}

    C=temp;
    temp.value.clear();
    temp.rows.clear();
    temp.columns.clear();
}

void Calculate_recursive_GF(Mat_1_doub A, Mat_1_doub B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy){

    complex<double> temp_n;
    temp_n.imag(0);temp_n.real(0);

    double val1,val2;
    omega = omega + GS_energy;



    for (int i=A.size()-1;i>=0;i--){
        if(i==(A.size()-1)){
            val1 = omega - A[i];
            val2 = eta;

            temp_n.real(val1/( (val1*val1) + (val2*val2) ) );
            temp_n.imag(-val2/( (val1*val1) + (val2*val2) ) );

        }
        else{
            val1= omega -A[i] - (B2[i+1]*temp_n.real());
            val2= eta - (B2[i+1]*temp_n.imag());

            temp_n.real(val1/( (val1*val1) + (val2*val2) ) );
            temp_n.imag(-val2/( (val1*val1) + (val2*val2) ) );

        }

    }


    Recursive_GF = temp_n;


}
#endif

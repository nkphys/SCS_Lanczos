#ifdef USE_COMPLEX
#ifdef _OPENMP
#include <omp.h>
#endif
#include "functions_complex.h"
#define PI_ 3.14159265

extern "C" void zheev_(char *,char *,int *, complex<double> *, int *, double *,
                       complex<double> *,int *, double *, int *);
//extern "C" void dsyev_(char * , char * , int * , double * , int *, double *, double *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);


int Find_int_in_intarray_using_bisection_(bool &found_, ulli num, Mat_1_ullint &array){

    //assuming increasing order array

    int pos;
    bool not_found=true;
    int ind=0;
    int m_min, m_max;
    m_min=0;m_max=array.size();

    int ind_old=-1000;

    while(not_found){

        ind=(m_min+m_max)/2;
        if(ind==ind_old){
            break;
        }

        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }

        if(num<array[ind]){
            m_max=ind;
        }
        if(num>array[ind]){
            m_min=ind;
        }



        ind_old=ind;
        //ind++;
    }

    found_ = (!not_found);
    return pos;
}


int Find_int_in_intarray(bool &found, ulli num, Mat_1_ullint &array){

    int pos;
    bool not_found=true;
    int ind=0;
    while(not_found){
        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }
        if(ind==(array.size()-1)){
            break;
        }
        ind++;
    }
    found = !not_found;
    return pos;
}

int Find_int_in_intarray_using_bisection(ulli num, Mat_1_ullint &array){

    //assuming increasing order array

    int pos;
    bool not_found=true;
    int ind=0;
    int m_min, m_max;
    m_min=0;m_max=array.size();



    while(not_found){

        ind=(m_min+m_max)/2;
        assert(ind<array.size());
        if(num==array[ind]){
            pos=ind;
            not_found = false;
        }

        if(num<array[ind]){
            m_max=ind;
        }
        if(num>array[ind]){
            m_min=ind;
        }


        //ind++;
    }

    return pos;
}


void DirectProduct_IndicesSet_by_recursion(Mat_1_int &TempState, Mat_1_int &VecSizes, int & Total_no_elements ,Mat_2_int &IndicesSet){


    Mat_1_int Diff;
    Total_no_elements +=1;
    int size = TempState.size();

    if(Total_no_elements==0){
        for(int i=0;i<TempState.size();i++){
            assert(TempState[i]==0);
        }
    }

    /*
for(int i=0;i<TempState.size();i++){
cout<<TempState[i]<<" ";
}
cout<<endl;
*/

    IndicesSet.push_back(TempState);

    //TempState[size-1] +=1;

    Diff.resize(TempState.size());
    for(int i=0;i<TempState.size();i++){
        Diff[i] = (VecSizes[i]-1)-TempState[i];
    }

    int Sum_diff=0;
    int n_to_grow;
    for(int n=(size-1);n>=0;n--){
        Sum_diff +=Diff[n];
        if(Sum_diff!=0){
            n_to_grow=n;
            break;
        }
    }

    bool proceed=(Sum_diff!=0);
    for(int n=0;n<size;n++){
        proceed = proceed && (TempState[n]<VecSizes[n]);
    }

    if(proceed){
        TempState[n_to_grow] +=1;

        for(int n=n_to_grow+1; n<size;n++){
            TempState[n] = 0;
        }

        DirectProduct_IndicesSet_by_recursion(TempState, VecSizes, Total_no_elements, IndicesSet);
    }

}



void Matrix_vector_multiplication(Mat_2_int &A, Mat_1_int &u, Mat_1_int &v){

    v.resize(A.size());

    for(int i=0;i<v.size();i++){
        v[i]=0;
        assert(A[i].size()==u.size());
        for(int j=0;j<u.size();j++){
            v[i] += A[i][j]*u[j];
        }
    }

}



void Perform_SVD(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


    char jobz='A'; //A,S,O,N

    int m=A_.n_row();
    int n=A_.n_col();
    int lda=A_.n_row();
    int ldu=A_.n_row();
    int ldvt=n;

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,m);

    VT_.resize(ldvt,n);


    vector<complex<double>> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));
    int lrwork = max( (5*min(m,n)*min(m,n)) + 5*min(m,n), (2*max(m,n)*min(m,n)) + (2*min(m,n)*min(m,n)) + min(m,n) );
    vector<double> rwork(lrwork);

    // query:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info<0.\n");
        }
    }

    // Ham_.print();



}



int Find_int_in_part_of_intarray(ulli num, Mat_1_ullint &array, int min_i, int max_i){

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

complex<double> conjugate(complex<double> val){

    return (conj(val));
}


double Lorentzian(double eta, double x){
    double val;
    val = (1.0/PI_)*( (eta/1.0) / (  (x*x) + ((eta*eta)/1.0)   ) );
    return val;
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
    {cout<<"Unable to open input file while in the functions_complex.cpp"<<endl;}

    stringstream Dynamics_dos_operators_(Dynamics_dos_operators);
    Dynamics_dos_operators_ >> no_of_fermionic_oprs_string;

    no_of_fermionic_oprs=atoi(no_of_fermionic_oprs_string.c_str());

    TRIO_VEC.resize(no_of_fermionic_oprs);
    values_.resize(no_of_fermionic_oprs);

    for(int i=0;i<no_of_fermionic_oprs;i++){

        //real part
        Dynamics_dos_operators_ >> temp;
        temp_double=atof(temp.c_str());
        values_[i].real(temp_double);

        //imag part
        Dynamics_dos_operators_ >> temp;
        temp_double=atof(temp.c_str());
        values_[i].imag(temp_double);

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

complex<double> reading_pair(string pair_str){

    complex<double> value;
    double part_one, part_two;
    int pair_str_Length=pair_str.length()-1;
    int pos = pair_str.find(",");
    string one_str=pair_str.substr(1,pos-1);
    string two_str=pair_str.substr(pos+1,pair_str_Length-1-pos);

    stringstream onestream(one_str);
    stringstream twostream(two_str);
    onestream>>part_one;
    twostream>>part_two;

    value.real(part_one);
    value.imag(part_two);

    return value;

}

void Read_matrix_from_file(string filepath,
                           Mat_2_doub &Mat, int row, int column){

    ifstream infile(filepath.c_str());

    Mat.resize(row);
    for(int i=0;i<row;i++){
        Mat[i].resize(column);
    }

    string tmp_str;

    for (int i =0;i<row;i++){
        for (int j =0;j<column;j++){

            infile>>tmp_str;
            Mat[i][j]=reading_pair(tmp_str);
        }

    }

}

complex<double> divide(complex<double> z1, complex<double> z2){

    complex<double> z;
    double x1,x2,y1,y2;
    x1=z1.real();
    x2=z2.real();
    y1=z1.imag();
    y2=z2.imag();

    z.real(  (x1*x2 + y1*y2)/(x2*x2  + y2*y2) );
    z.imag(  (y1*x2 - x1*y2)/(x2*x2  + y2*y2) );

    return z;
}

void Print_vector_in_file(Mat_1_doub vec, string filename){

    ofstream outfile(filename.c_str());

    for(int j=0;j<vec.size();j++){
        outfile<<vec[j]<<endl;
    }

}


void Print_Matrix_COO_in_file(Matrix_COO &A, string filename){

    ofstream outfile(filename.c_str());

    outfile<<scientific<<setprecision(16);

    Mat_2_doub B;
    B.resize(A.nrows);
    for(int i=0;i<B.size();i++){
        B[i].resize(A.nrows);
        for(int j=0;j<A.nrows;j++){
            B[i][j]=0.0;
        }
    }

    for(int i=0;i<A.value.size();i++){
        B[A.rows[i]][A.columns[i]] +=A.value[i];
    }

    for(int i=0;i<B.size();i++){
        for(int j=0;j<B.size();j++){

            outfile<<B[i][j]<<"  ";
        }
        outfile<<endl;
    }

}

void Sort_vector_in_decreasing_order_in_file(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old){

    //ofstream outfile(filename.c_str());

    Mat_1_doub Vec_temp;
    Vec_temp=vec;
    Vec_new.clear();
    Index_old.clear();
    Vec_new.resize(Vec_temp.size());
    Index_old.resize(Vec_temp.size());

    int count=0;
    complex<double> max_temp=0.0;
    int index_max;

    while(count<Vec_temp.size()){
        max_temp=0.0;
        for(int j=0;j<Vec_temp.size();j++){
            if(abs(Vec_temp[j])>=abs(max_temp)){
                max_temp=Vec_temp[j];
                index_max = j;
            }

        }
        //outfile<<index_max<<"   "<<Vec_temp[index_max]<<endl;
        Vec_new[count]=Vec_temp[index_max];
        Index_old[count]=index_max;
        count++;
        Vec_temp[index_max]=0.0;
    }


}



void Sort_vector_in_decreasing_order(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old, int n_basis){


    Mat_1_doub Vec_temp;
    Vec_temp=vec;
    Vec_new.clear();
    Index_old.clear();
    Vec_new.resize(n_basis);
    Index_old.resize(n_basis);

    int count=0;
    complex<double> max_temp=0.0;
    int index_max;

    while(count<n_basis){
        max_temp=0.0;

        for(int j=0;j<Vec_temp.size();j++){
            if(abs(Vec_temp[j])>=abs(max_temp)){
                max_temp=Vec_temp[j];
                index_max = j;
            }
        }

        Vec_new[count]=Vec_temp[index_max];
        Index_old[count]=index_max;
        count++;
        Vec_temp[index_max]=0.0;
    }


}


void Sort_vector_in_decreasing_order(Mat_1_doub vec, Mat_1_doub &Vec_new, Mat_1_int &Index_old, Mat_1_int &Index_new, int n_basis){


    Mat_1_doub Vec_temp;
    Vec_temp=vec;
    Vec_new.clear();
    Index_new.clear();
    Vec_new.resize(n_basis);
    Index_new.resize(n_basis);

    int count=0;
    complex<double> max_temp=0.0;
    int index_max;

    while(count<n_basis){
        max_temp=0.0;

        for(int j=0;j<Vec_temp.size();j++){
            if(abs(Vec_temp[j])>=abs(max_temp)){
                max_temp=Vec_temp[j];
                index_max = j;
            }
        }

        Vec_new[count]=Vec_temp[index_max];
        Index_new[count]=Index_old[index_max];
        count++;
        Vec_temp[index_max]=0.0;
    }




}


void Print_file_in_vector(Mat_1_doub &vec, string filename, int rows){

    cout<<"Vector is being read from : '"<<filename<<"'"<<endl;
    vec.clear();
    ifstream infile(filename.c_str());
    //string line;
    double_type temp_double;


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




int Find_intpair_in_intarraypair(int num1, int num2 ,Mat_1_int &array1, Mat_1_int &array2,
                                 int num1_sector, Mat_1_intpair &sectors_offset){


    int pos;
    bool not_found=true;
    int ind, ind_max;

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

int Find_int_in_intarray(ulli num, Mat_1_ullint &array){

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

void Remove_repetitions(Mat_1_int & index_array, Mat_1_doub & val_array, Mat_1_int & index_new_array, Mat_1_doub & val_new_array){


    index_new_array.clear();
    val_new_array.clear();

    int m_min=-100;
    int mi;
    double_type val;
    for(int i=0;i<index_array.size();i++){

        mi = index_array[i];

        if(mi!=m_min){
            index_new_array.push_back(mi);
            val=val_array[i];

            for(int j=i+1;j<index_array.size();j++){

                if(index_array[j]==mi){
                    val +=val_array[j];

                    index_array[j]=-100;
                }

            }

            index_array[i]=-100;

            val_new_array.push_back(val);
        }
    }

}

int partition(Mat_1_int &a, int s, int e)
{
    int piviot = a[e];
    int pind = s;
    int i, t;

    for (i = s; i < e; i++) {
        if (a[i] <= piviot) {
            t = a[i];
            a[i] = a[pind];
            a[pind] = t;
            pind++;
        }
    }

    t = a[e];
    a[e] = a[pind];
    a[pind] = t;

    return pind;
}

void quicksort(Mat_1_int &a, int s, int e)
{
    if (s < e) {
        int pind = partition(a, s, e);
        quicksort(a, s, pind - 1);
        quicksort(a, pind + 1, e);
    }
}


void Remove_repetitions(Mat_1_ullint & index_array, Mat_1_doub & val_array, Mat_1_ullint & index_new_array, Mat_1_doub & val_new_array){


    index_new_array.clear();
    val_new_array.clear();

    Mat_1_int check;
    check.resize(index_array.size());
    for(int i=0;i<check.size();i++){
        check[i]=0;
    }


    int mi;
    double_type val;
    for(int i=0;i<index_array.size();i++){

        mi = index_array[i];

        if(check[i]==0){
            index_new_array.push_back(mi);
            val=val_array[i];

            for(int j=i+1;j<index_array.size();j++){

                if(index_array[j]==mi){
                    val +=val_array[j];

                    check[j]=1;
                }

            }

            check[i]=1;

            val_new_array.push_back(val);
        }
    }

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


void Print_Matrix_COO_format(Matrix_COO &A){


    cout<<"--------------------PRINTING THE MATRIX COO FORMAT:--------------------"<<endl;
    for(int i=0;i<A.rows.size();i++){
            cout<<A.rows[i]<<"  "<<A.columns[i]<<"  "<<A.value[i]<<endl;
    }
    cout<<"-------------------------------------------------------------"<<endl;


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
        B[A.rows[i]][A.columns[i]] +=A.value[i];
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

void Print_Matrix(Mat_2_doub &A){

    for(int i=0;i<A.size();i++){
        for(int j=0;j<A.size();j++){
            cout<<A[i][j]<<"  ";
        }
        cout<<endl;
    }


}


void Print_Matrix(Mat_2_real &A){

    for(int i=0;i<A.size();i++){
        for(int j=0;j<A.size();j++){
            cout<<A[i][j]<<"  ";
        }
        cout<<endl;
    }


}


void Normalize_vec(Mat_1_doub &vec_in){

    complex<double> val_comp;
    double val;

    val_comp = dot_product(vec_in, vec_in);
    val = val_comp.real();
    val = sqrt(val);

    for(int i=0;i<vec_in.size();i++){
        vec_in[i] = vec_in[i]*(1.0/val);
    }

}

void value_multiply_vector(complex<double> value, Mat_1_doub &vec_in){
    //NEVER PARALLELIZE THIS
    for(int i=0;i<vec_in.size();i++){
        vec_in[i] = vec_in[i]*value;
    }
}

complex<double> dot_product(Mat_1_doub &vec1, Mat_1_doub &vec2){
    complex<double> temp;
    complex<double> temp1=zero;
    double temp1_real=0.0;
    double temp1_imag=0.0;
    assert(vec1.size()==vec2.size());
    // bool Parallelize_dot_product;
    // Parallelize_dot_product=true;

    // if(!Parallelize_dot_product){goto skiploop_143;}
    //#pragma omp parallel for default(shared) reduction(+:temp1)
    //skiploop_143:
    //#ifndef _OPENMP
    for(int i=0;i<vec1.size();i++){
        temp1 = temp1 + (vec1[i])*(conj(vec2[i])); //<vec2|vec1>
    }
    //#endif

    /*
#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:temp1_real,temp1_imag)
    for(int i=0;i<vec1.size();i++){
        temp1_real += ((vec1[i])*(conj(vec2[i])) ).real(); //<vec2|vec1>
        temp1_imag += ((vec1[i])*(conj(vec2[i])) ).imag();
    }
    temp1 =complex<double> (temp1_real, temp1_imag);
#endif

 */

    temp = temp1;

    return temp;

}

double Norm(Mat_1_doub &vec1){
    complex<double> temp;
    complex<double> temp1=zero;
    double temp1_real=0.0;
    double temp1_imag=0.0;

#ifndef _OPENMP
    for(int i=0;i<vec1.size();i++){
        temp1 = temp1 + (vec1[i])*(conj(vec1[i])); //<vec1|vec1>
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+:temp1_real,temp1_imag)
    for(int i=0;i<vec1.size();i++){
        temp1_real += ((vec1[i])*(conj(vec1[i]))).real(); //<vec1|vec1>
        temp1_imag += ((vec1[i])*(conj(vec1[i]))).imag();
    }

    temp1=complex<double> (temp1_real, temp1_imag);
#endif


    temp = temp1;

    assert(temp.imag()<1e-6);

    return temp.real();

}

void Matrix_COO_vector_multiplication(string COO_type, Matrix_COO & A,Mat_1_doub &u,Mat_1_doub &v){

    v.clear();
    v.resize(A.nrows);
    assert(A.ncols==u.size());

   // #ifndef _OPENMP
    for (int i=0;i<v.size();i++){
        v[i]=zero;}

    if(COO_type =="U"){

        for (int n=0;n<A.value.size();n++){

            if(A.rows[n]==A.columns[n]){
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
            }
            else{
                v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
                v[A.columns[n]] = v[A.columns[n]] + u[A.rows[n]]*conj(A.value[n]);
            }

        }
    }
    else{

        for (int n=0;n<A.value.size();n++){
            v[A.rows[n]] = v[A.rows[n]] + u[A.columns[n]]*A.value[n];
        }
    }
  //  #endif



#ifdef _OPENMP_NOTUSED
    Mat_2_doub v_temp;
    int thread_id;
    int N_threads;
#pragma omp parallel
    {
        N_threads=omp_get_num_threads();
        //cout<<"N_threads = "<<N_threads<<endl;
    }

    v_temp.resize(N_threads);

#pragma omp parallel for default(shared)
    for(int tid=0;tid<N_threads;tid++){
        v_temp[tid].resize(A.nrows);
        for (int i=0;i<A.nrows;i++){
            v_temp[tid][i]=zero;
        }
    }

    //cout<<"here 1, No_of_threads = "<<N_threads<<endl;

    //assert(false);
    if(COO_type =="U"){

#pragma omp parallel for default(shared) private(thread_id)
        for (int n=0;n<A.value.size();n++){

            thread_id = omp_get_thread_num();

            //            if(thread_id>=N_threads){
            //                cout<<"thread_id = "<<thread_id<<",  N_threads = "<<N_threads<<endl;
            //            }
            //cout<<"n="<<n<<"  thread_id="<<thread_id<<endl;
            if(A.rows[n]==A.columns[n]){
                v_temp[thread_id][A.rows[n]] = v_temp[thread_id][A.rows[n]] + u[A.columns[n]]*A.value[n];
            }
            else{
                v_temp[thread_id][A.rows[n]] = v_temp[thread_id][A.rows[n]] + u[A.columns[n]]*A.value[n];
                v_temp[thread_id][A.columns[n]] = v_temp[thread_id][A.columns[n]] + u[A.rows[n]]*conj(A.value[n]);
            }

        }
    }
    else{
#pragma omp parallel for default(shared) private(thread_id)
        for (int n=0;n<A.value.size();n++){
            thread_id = omp_get_thread_num();
            v_temp[thread_id][A.rows[n]] = v_temp[thread_id][A.rows[n]] + u[A.columns[n]]*A.value[n];
        }
    }


    // cout<<"here 2"<<endl;

    double sum_real, sum_imag;
#pragma omp parallel for default(shared) private (sum_real, sum_imag)
    for(int i=0;i<A.nrows;i++){
        sum_real=0.0;sum_imag=0.0;
        //#pragma omp parallel for default(shared) reduction(+:sum_real, sum_imag)
        for(int tid=0;tid<N_threads;tid++){
            sum_real += v_temp[tid][i].real();
            sum_imag += v_temp[tid][i].imag();
        }
        v[i] = complex<double> (sum_real, sum_imag);
    }

    // cout<<"here 3"<<endl;

    //    sum = 0;
    //    #pragma omp parallel for shared(sum, a) reduction(+: sum)
    //    for (auto i = 0; i < 10; i++)
    //    {
    //        sum += a[i]
    //    }



#pragma omp parallel for default(shared)
    for(int tid=0;tid<N_threads;tid++){
        vector< double_type >().swap( v_temp[tid] );
    }
    v_temp.clear();

    //cout<<"here 4"<<endl;

#endif


    //  assert(false);
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

        B.value[i] = conj(A.value[i]);
        B.rows[i] = A.columns[i];
        B.columns[i] = A.rows[i];

        i=i+1;
    }
    return B;
}


Matrix_COO Identity_COO(int rows_no, int cols_no){
    Matrix_COO B;
    B.nrows=rows_no;
    B.ncols=cols_no;

    B.rows.resize(rows_no);
    B.columns.resize(cols_no);
    B.value.resize(rows_no);


    for(int i=0;i<rows_no;i++){

        B.value[i] = one;
        B.rows[i] = i;
        B.columns[i] = i;

    }

    return B;

}



void Subtract( Mat_1_doub &temp1, double_type x, Mat_1_doub &temp2){

    assert(temp1.size()==temp2.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int k=0;k<temp1.size();k++){
        temp1[k]=temp1[k] - x*temp2[k];
    }

}


void Subtract( Mat_1_doub &temp1, double x, Mat_1_doub &temp2){

    assert(temp1.size()==temp2.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int k=0;k<temp1.size();k++){
        temp1[k]=temp1[k] - x*temp2[k];
    }

}



void Diagonalize(Matrix_COO &X, double & EG, Mat_1_doub & vecG){

    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.nrows,X.ncols);


    for(int i=0;i<X.nrows;i++){
        for(int j=i;j<X.ncols;j++){
            Ham_(i,j) = complex<double>(0.0,0.0);
        }
    }


    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
        Ham_(r,c) += X.value[i];
        if(r!=c){
            assert(c>r);
            Ham_(c,r) += conjugate(X.value[i]);}
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
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

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.nrows;i++){
        vecG[i] =Ham_(i,0);//mat[i*X.nrows];
    }

    Ham_.clear();

}



void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_1_doub & vecG){


    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.nrows,X.ncols);


    for(int i=0;i<X.nrows;i++){
        for(int j=0;j<X.ncols;j++){
            Ham_(i,j) = complex<double>(0.0,0.0);
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
        Ham_(r,c) += X.value[i];
        if(r!=c){
            assert(c>r);
            Ham_(c,r) += conjugate(X.value[i]);}
    }




    //    Mat_1_real eigs_;
    //    Matrix< complex<double> > Ham_;
    //    Ham_.resize(36,36);

    //    for(int i=0;i<36;i++){
    //        for(int j=0;j<36;j++){
    //            Ham_(i,j) = complex<double>(0.0,0.0);
    //        }
    //    }

    //    Ham_(5,10)=complex<double>(0.0,0.5);
    //    Ham_(5,15)=complex<double>(0.0,0.5);
    //    Ham_(5,20)=complex<double>(0.0,0.5);
    //    Ham_(10,15)=complex<double>(0.0,0.5);
    //    Ham_(10,20)=complex<double>(0.0,0.5);
    //    Ham_(15,25)=complex<double>(0.0,-0.5);
    //    Ham_(15,30)=complex<double>(0.0,-0.5);
    //    Ham_(20,25)=complex<double>(0.0,-0.5);
    //    Ham_(20,30)=complex<double>(0.0,-0.5);
    //    Ham_(25,30)=complex<double>(0.0,-0.5);

    //    Ham_(10,5)=complex<double>(0.0,-0.5);
    //    Ham_(15,5)=complex<double>(0.0,-0.5);
    //    Ham_(20,5)=complex<double>(0.0,-0.5);
    //    Ham_(15,10)=complex<double>(0.0,-0.5);
    //    Ham_(20,10)=complex<double>(0.0,-0.5);
    //    Ham_(25,15)=complex<double>(0.0,0.5);
    //    Ham_(30,15)=complex<double>(0.0,0.5);
    //    Ham_(25,20)=complex<double>(0.0,0.5);
    //    Ham_(30,20)=complex<double>(0.0,0.5);
    //    Ham_(30,25)=complex<double>(0.0,0.5);


    //    Ham_(5,10)=complex<double>(-1.0,0.0);
    //    Ham_(5,15)=complex<double>(-1.0,0.0);
    //    Ham_(5,20)=complex<double>(-1.0,0.0);
    //    Ham_(10,15)=complex<double>(-1.0,0.0);
    //    Ham_(10,20)=complex<double>(-1.0,0.0);
    //    Ham_(15,25)=complex<double>(1.0,0.0);
    //    Ham_(15,30)=complex<double>(1.0,0.0);
    //    Ham_(20,25)=complex<double>(1.0,0.0);
    //    Ham_(20,30)=complex<double>(1.0,0.0);
    //    Ham_(25,30)=complex<double>(1.0,0.0);

    //    Ham_(10,5)=complex<double>(-1.0,0.0);
    //    Ham_(15,5)=complex<double>(-1.0,0.0);
    //    Ham_(20,5)=complex<double>(-1.0,0.0);
    //    Ham_(15,10)=complex<double>(-1.0,0.0);
    //    Ham_(20,10)=complex<double>(-1.0,0.0);
    //    Ham_(25,15)=complex<double>(1.0,0.0);
    //    Ham_(30,15)=complex<double>(1.0,0.0);
    //    Ham_(25,20)=complex<double>(1.0,0.0);
    //    Ham_(30,20)=complex<double>(1.0,0.0);
    //    Ham_(30,25)=complex<double>(1.0,0.0);



    //    cout<<scientific<<setprecision(2);
    //    Ham_.print();

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}



    EVALS.resize(eigs_.size());
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<eigs_.size();i++){
        EVALS[i]=eigs_[i];
    }

    eigs_.clear();
    vecG.clear();
    vecG.resize(X.nrows);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.nrows;i++){
        vecG[i] =Ham_(i,0);//mat[i*X.nrows];
    }

    Ham_.clear();

}


void Diagonalize(Matrix_COO &X, Mat_1_real & EVALS, Mat_2_doub & vecs){

    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.nrows,X.ncols);

    for(int i=0;i<X.nrows;i++){
        for(int j=0;j<X.ncols;j++){
            Ham_(i,j) = complex<double>(0.0,0.0);
        }
    }

    for(int i=0;i<X.value.size();i++){
        int r=X.rows[i];
        int c=X.columns[i];
        Ham_(r,c) += X.value[i];
        if(r!=c){
            assert(c>r);
            Ham_(c,r) += conjugate(X.value[i]);}
    }

    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    EVALS.resize(eigs_.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<eigs_.size();i++){
        EVALS[i]=eigs_[i];
    }

    eigs_.clear();
    vecs.clear();
    vecs.resize(X.nrows);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int j=0;j<vecs.size();j++){
        vecs[j].resize(vecs.size());
        for(int i=0;i<X.nrows;i++){
            vecs[j][i] =Ham_(i,j);//mat[i*X.nrows];
        }
    }

    Ham_.clear();

}



void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_1_doub & vecG,int lanc_iter){


    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.size(),X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=complex<double>(sqrt(Y2[i]),0.0);
                Ham_(j,i)=complex<double>(sqrt(Y2[i]),0.0);
            }
            else{
                Ham_(i,j)=complex<double>(0.0,0.0);
            }
        }
    }

    //------------------------------------------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //--------------------------------------------------

    int Target_state=0;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        vecG[i] = Ham_(i,Target_state);
    }

    Ham_.clear();

}



void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_2_doub & vecG,int lanc_iter ,int few_){

    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.size(),X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=complex<double>(sqrt(Y2[i]),0.0);
                Ham_(j,i)=complex<double>(sqrt(Y2[i]),0.0);
            }
            else{
                Ham_(i,j)=complex<double>(0.0,0.0);
            }
        }
    }


    //------------------------------------------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //--------------------------------------------------



#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(few_);

    for(int i=0;i<few_;i++){
        vecG[i].resize(X.size());
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<few_;i++){
        for(int j=0;j<X.size();j++){
            vecG[i][j] =  Ham_(j,i);
        }
    }

    Ham_.clear();


}



void Diagonalize(Mat_1_doub &X ,Mat_1_real &Y2 , Mat_2_real & EG, Mat_2_doub & vecG,int lanc_iter ,int few_, Mat_1_int states_to_look){

    Mat_1_real eigs_;
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.size(),X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=complex<double>(sqrt(Y2[i]),0.0);
                Ham_(j,i)=complex<double>(sqrt(Y2[i]),0.0);
            }
            else{
                Ham_(i,j)=complex<double>(0.0,0.0);
            }
        }
    }

    //------------------------------------------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //--------------------------------------------------

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


    int Target_state=0;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        EG[lanc_iter][i]=eigs_[i];
    }

    eigs_.clear();

    vecG.clear();
    vecG.resize(few_);

    for(int i=0;i<few_;i++){
        vecG[i].resize(X.size());
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
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
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.size(),X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=complex<double>(sqrt(Y2[i]),0.0);
                Ham_(j,i)=complex<double>(sqrt(Y2[i]),0.0);
            }
            else{
                Ham_(i,j)=complex<double>(0.0,0.0);
            }
        }
    }

    //------------------------------------------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //--------------------------------------------------

    int T_s=0;
    int Target_state=0;

    EG=eigs_[T_s];
    eigs_.clear();
    vecG.clear();
    vecG.resize(X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
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
    Matrix< complex<double> > Ham_;
    Ham_.resize(X.size(),X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        for(int j=0;j<=i;j++){
            if(j==i){
                Ham_(i,j)=X[j];
            }
            else if(j==i-1){
                Ham_(i,j)=complex<double>(sqrt(Y2[i]),0.0);
                Ham_(j,i)=complex<double>(sqrt(Y2[i]),0.0);
            }
            else{
                Ham_(i,j)=complex<double>(0.0,0.0);
            }
        }
    }

    //------------------------------------------------
    char jobz='V';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector< complex<double> > work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //--------------------------------------------------


    int T_s=0;
    int Target_state=0;
    /*  if(Target_state>X.size()-1){
        T_s=X.size()-1;
    }
    else{
        T_s=Target_state;
    }*/
    EG=eigs_[T_s];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        Evals_Lanczos[i]=eigs_[i];
    }

    eigs_.clear();
    vecG.clear();
    vecG.resize(X.size());
    Unit_E_vecs.clear();
    Unit_E_vecs.resize(X.size());

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){
        Unit_E_vecs[i].resize(X.size());}

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int i=0;i<X.size();i++){

        vecG[i] = Ham_(i,Target_state);//mat[i*X.size()+Target_state];

        for(int j_state=0;j_state<X.size();j_state++){

            Unit_E_vecs[j_state][i] = Ham_(i,j_state);//mat[i*X.size()+j_state];

        }
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



//C=value1*A + value2*B
void Sum(Matrix_COO A, Matrix_COO B, Matrix_COO & C, complex<double> value1, complex<double> value2){


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

void Calculate_recursive_GF(Mat_1_doub A, Mat_1_real B2, complex<double> &Recursive_GF, double omega,
                            double eta, double GS_energy){


    complex<double> iota(0,1);
    complex<double> temp_n;
    temp_n.imag(0);temp_n.real(0);

    double_type val1;
    omega = omega + GS_energy;



    for (int i=A.size()-1;i>=0;i--){
        if(i==(A.size()-1)){
            val1 = eta*iota + (omega*one) - A[i];
        }
        else{
            val1= eta*iota + omega*one - A[i] - (B2[i+1]*temp_n);
        }
        temp_n = divide(one, val1);

    }

    Recursive_GF = temp_n;


}
#endif

#include <cassert>

#ifndef Matrix_class
#define Matrix_class

template<typename T>
class  Matrix  {
public:
    typedef T value_type;

    //set all elements to zero
    Matrix()
        : nrow_(0), ncol_(0)
    {}

    //allocate number of row col and elements
    Matrix(int nrow,int ncol)
        : nrow_(nrow),ncol_(ncol),data_(nrow*ncol)
    {}

    // copy constructor
    Matrix(const Matrix<T>& m) {
        nrow_=m.nrow_;
        ncol_=m.ncol_;
        data_=m.data_;
    }

    const T& operator()(int i,int j) const;
    T& operator()(int i, int j);
    void print();
    void resize(int newrow, int newcol);
    int n_row();
    int n_col();
    void fill(T val);
    void clear();

private:
    int nrow_,ncol_;
    std::vector<T> data_;
};



/*  
 * ***********
 *  Functions in Class Matrix ------
 *  ***********
*/ 

template<class T>
int Matrix<T>::n_row() {
    return nrow_;
} // ----------

template<class T>
int Matrix<T>::n_col() {
    return ncol_;
} // ----------

template<class T>
void Matrix<T>::fill(T val) {
    std::fill(data_.begin(),data_.end(),val);
} // ----------

template<class T>
void Matrix<T>::resize(int newrow, int newcol) {
    //		assert(i<nrow_ && j<ncol_);
    //		assert(i+j*nrow_<data_.size());
    //~ cout<<"newsize "<< newrow*newcol << endl;
    nrow_=newrow;
    ncol_=newcol;
    data_.clear();
    data_.resize(newrow*newcol);
} // ----------

template<class T>
void Matrix<T>::clear() {
    //		assert(i<nrow_ && j<ncol_);
    //		assert(i+j*nrow_<data_.size());
    //~ cout<<"newsize "<< newrow*newcol << endl;
    nrow_=0;
    ncol_=0;
    data_.clear();
} // ----------


template<class T>
const T& Matrix<T>::operator()(int i, int j) const{
    //		assert(i<nrow_ && j<ncol_);
    //		assert(i+j*nrow_<data_.size());
    return data_[i+j*nrow_];
} // ----------


template<class T>
T& Matrix<T>::operator()(int i,int j){
    //		assert(i<nrow_ && j<ncol_);
    //		assert(i+j*nrow_<data_.size());
    return data_[i+j*nrow_];
} // ----------


template<class T>
void Matrix<T>::print(){
    std::cout<<"shape:= ("<<nrow_<<","<<ncol_<<")"<<std::endl;
    for(int i=0; i<nrow_; i++) {
        for(int j=0; j<ncol_; j++) {
            std::cout << data_[i+j*nrow_] << "\t";
        }
        std::cout << std::endl;
    }
    return;
} // ----------


#endif





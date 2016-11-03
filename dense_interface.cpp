#include <complex>
#include <algorithm>
#include <vector>
#include <iostream>

#include "dense_interface.hpp"
#include "dense_matrix_functions.hpp"
#include "arpack_interface.hpp"

namespace ajaj {

  DenseMatrix::DenseMatrix(Denseint param_rows,Denseint param_cols) : size_(param_rows*param_cols), nrows(param_rows), ncols(param_cols), m_array(size_>0 ? new std::complex<double>[size_] : nullptr){}
  
  DenseMatrix::DenseMatrix(Denseint param_rows, Denseint param_cols,std::complex<double> val ) : size_(param_rows*param_cols), nrows(param_rows), ncols(param_cols),    m_array((size_) ? new std::complex<double>[size_] : nullptr){
    for (Denseint i=0;i<size_;++i){
      m_array[i]=val;
    }
  }

  DenseMatrix::DenseMatrix(Denseint param_rows,Denseint param_cols, std::complex<double>** array_ptr ) : size_(param_rows*param_cols), nrows(param_rows), ncols(param_cols), m_array(nullptr){
    std::swap(m_array,*array_ptr);
  }

  //rare that this should be needed?
  /*DenseMatrix::DenseMatrix(const DenseMatrix& other) : nrows(other.nrows),ncols(other.ncols),m_array(nrows >0 && ncols > 0 ? new std::complex<double>[other.nrows*other.ncols] : nullptr){ //COPY constructor, does a DEEP copy
    std::copy(other.m_array,other.m_array+(nrows*ncols),m_array);

  }*/

  DenseMatrix::DenseMatrix(DenseMatrix&& other) noexcept : m_array(nullptr) {swap(*this,other);}//move

  DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other) {
    swap(*this,other);
    return *this;
  }

  DenseMatrix::DenseMatrix(Denseint param_rows,Denseint param_cols, const std::complex<double>* array ) : size_(param_rows*param_cols),nrows(param_rows), ncols(param_cols), m_array(size_>0 ? new std::complex<double>[nrows*ncols] : nullptr){
    std::copy(array,array+(nrows*ncols),m_array);
    /*for (Denseint d=0;d<nrows*ncols;++d){
      m_array[d]=array[d];
      }*/
  }

  DenseMatrix::~DenseMatrix(){
    delete[] m_array;
  }

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  void DenseMatrix::print() const {
    for (Denseint m=0;m<nrows;++m){
      for (Denseint n=0;n<ncols;++n){
	cout << m_array[m+nrows*n] << " ";
      }
      cout << endl;
    }
  };

  DenseMatrix DenseMatrix::inverse(){
    //currently only square inverse defined
    if (nrows!=ncols){cout << "Not square, can't do inverse" << endl;exit(1);}
    DenseMatrix ans(nrows,nrows);
    densefuncs::inverse_with_lapack(nrows, m_array, ans.m_array);
    return ans;
  }

  DenseMatrix& DenseMatrix::operator*=(const DenseMatrix& rhs)
  {
    //BLAS tends to overwrite things, so need a copy of the right hand side
    std::complex<double>* rhs_copy_array=new std::complex<double>[rhs.nrows*rhs.ncols];
    for (size_t i=0;i<static_cast<size_t>(rhs.nrows*rhs.ncols);++i){
      rhs_copy_array[i]=rhs.m_array[i];
    }
    //new array for the answer
    std::complex<double>* new_lhs=new std::complex<double>[nrows*rhs.ncols];

    densefuncs::gemm(nrows, ncols, rhs.ncols, std::complex<double>(1.0,0.0), std::complex<double>(0.0,0.0), m_array, rhs_copy_array, new_lhs);
    delete[] m_array; //destroy the original array, which now contains junk
    delete[] rhs_copy_array;
    this->m_array=new_lhs;
    this->ncols=rhs.ncols;
    return *this;
  }

  void swap (DenseMatrix& first, DenseMatrix& second)
  {
    using std::swap;
    swap(first.m_array,second.m_array);
    swap(first.size_,second.size_);
    swap(first.nrows,second.nrows);
    swap(first.ncols,second.ncols);
  }

  inline DenseMatrix operator*(DenseMatrix lhs, const DenseMatrix& rhs)
  {
    lhs *= rhs;
    return lhs;
  }

  DenseDecompositionBase<std::complex<double> > DenseMatrix::NHEigenvalues(){
    //check square!
    if (nrows!=ncols){std::cout << "Not square, can't do nh eigenvalue decomposition" << std::endl;exit(1);}
    DenseDecompositionBase<std::complex<double> > ans(nrows);
    densefuncs::diagonalise_with_lapack_nh(nrows, m_array,ans.Values);
    //trashes the original array
    //clear();
    return ans;
  }

  //Hermitian eigen decomposition
  DenseHED DenseMatrix::HED(){
    //check square!
    if (nrows!=ncols){std::cout << "Not square, can't do eigenvalue decomposition " << nrows << " " << ncols << std::endl;exit(1);}
    DenseHED ans(nrows);
    densefuncs::diagonalise_with_lapack(nrows, m_array, ans.EigenVectors.m_array,ans.Values,1,nrows);
    //zheevr trashes the original array
    //clear();
    return ans;
  }

  DenseHED DenseMatrix::HED(Denseint numevals,char which[3]){
    //check square!
    if (nrows!=ncols){std::cout << "Not square, can't do eigenvalue decomposition " << nrows << " " << ncols << std::endl; this->print(); exit(1);}
    static char SMALLESTREAL[]={'S','R','\n'};
    static char LARGESTTREAL[]={'L','R','\n'};
    DenseHED ans(nrows,numevals);
    //std::cout << "nrows " << nrows << ", " << "numevals "  << numevals <<std::endl;
    if (which[0]==SMALLESTREAL[0]){
      densefuncs::diagonalise_with_lapack(nrows, m_array, ans.EigenVectors.m_array,ans.Values,1,numevals);
    }
    else if (which[0]==LARGESTTREAL[0]){
      densefuncs::diagonalise_with_lapack(nrows, m_array, ans.EigenVectors.m_array,ans.Values,nrows-numevals+1,nrows);
    }
    else {
      std::cout << "Incorrect argument to LAPACK " << std::endl; exit(1);
    }
    //clear();
    return ans;
  }

  //Non Hermitian eigen decomposition
  DenseNHED DenseMatrix::NHED(){
    //check square!
    if (nrows!=ncols){std::cout << "Not square, can't do nh eigenvalue decomposition" << std::endl;exit(1);}
    DenseNHED ans(nrows);
    densefuncs::diagonalise_with_lapack_nh(nrows, m_array, ans.RightEigenVectors.m_array,ans.Values);
    //trashes the original array
    //clear();
    return ans;
  }

  //Singular value decomposition
  DenseSVD DenseMatrix::SVD(){
    DenseSVD ans(nrows,ncols);
    densefuncs::svd_with_lapack(nrows, ncols,m_array, ans.U.m_array,ans.Vdagger.m_array,ans.Values);
    return ans;
  }

  std::vector<double> DenseMatrix::SingularValues(){
    DenseDecompositionBase<double> svals(nrows > ncols ? ncols : nrows);
    std::vector<double> ans(nrows > ncols ? ncols : nrows);
    densefuncs::svd_with_lapack(nrows,ncols,m_array,svals.Values);
    //clear(); //dump trash
    std::copy(svals.Values,svals.Values+svals.ValuesSize(),ans.begin());
    return ans;
  }

  DenseHED::DenseHED(Denseint L) : DenseDecompositionBase<double>(L),EigenVectors(L,L) {}     
  DenseHED::DenseHED(Denseint L,Denseint numevals) : DenseDecompositionBase<double>(numevals),EigenVectors(L,numevals) {
    if (numevals>L){std::cout << "Error, too many evals requested" << std::endl; exit(1);}
  }

  DenseNHED::DenseNHED(Denseint L) : DenseDecompositionBase<std::complex<double> >(L),RightEigenVectors(L,L) {}

  DenseSVD::DenseSVD(Denseint M, Denseint N) : DenseDecompositionBase<double>(M>N ? N : M),leftdim(M),rightdim(N),U(M,lineardim),Vdagger(lineardim,N) {}


}

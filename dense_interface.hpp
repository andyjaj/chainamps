#ifndef DMAT_H
#define DMAT_H

#include <complex>
//#include <algorithm>
#include <vector>
#include <iostream>
#include <cs.h>

#include "ajaj_common.hpp" //defines int types
#include "dense_matrix_functions.hpp"

namespace ajaj {

  typedef SuiteSparse_long Denseint;
  //////////////////////////////////////////////////////////////////////////
  template <typename T>
  class DenseDecompositionBase;
  class DenseMatrix;
  class DenseHED;
  class DenseSVD;
  class DenseNHED;

  void swap (DenseMatrix& first, DenseMatrix& second);
  void move_to_dumb_array(DenseMatrix& dm, std::complex<double>* array);
  //////////////////////////////////////////////////////////////////////////
  class DenseMatrix {
  private:
    std::complex<double>* m_array; //default
    Denseint nrows;
    Denseint ncols;
    void clear() {delete[] m_array;m_array=0;nrows=0;ncols=0;}

  public:
    //When first formed matrix is allocated as transpose! When finalised it is transposed again, which orders the rows.
    DenseMatrix(const Denseint& param_rows,const Denseint& param_cols); //create array without initialising entries...
    DenseMatrix(const Denseint& param_rows,const Denseint& param_cols, const std::complex<double>& val ); //initialise entries
    DenseMatrix(const Denseint& param_rows,const Denseint& param_cols, std::complex<double>** array_ptr ); //initialise entries
    DenseMatrix(const DenseMatrix& other); //copy constructor, must DEEP copy
    DenseMatrix(DenseMatrix&& other) noexcept;//move
    DenseMatrix() : m_array(0), nrows(0),ncols(0) {};//useful when creating an instance that needs to be replaced later
    DenseMatrix(const Denseint& param_rows,const Denseint& param_cols, const std::complex<double>* array );
    ~DenseMatrix();

    inline Denseint rows() const {return nrows;}
    inline Denseint cols() const {return ncols;}

    void print() const;
    bool empty() const {return m_array ? 0 : 1;}

    //array operations
    complex<double> Value(Denseint i, Denseint j) const;

    void entry(Denseint i, Denseint j, std::complex<double> value); //insert a value into the array
    void finalise(){};
    DenseMatrix inverse(); //returns inverse, or fails spectacularly if singular
    DenseMatrix& operator=(DenseMatrix other);
    DenseMatrix& operator*=(const DenseMatrix& rhs);
    //SparseMatrix sparse(const double& RELATIVETOL); //make a sparse version by dropping elements that are small relative to the largest element.
    friend void swap(DenseMatrix& first, DenseMatrix& second); //useful for copy-swap
    friend void move_to_dumb_array(DenseMatrix& dm, std::complex<double>* array);
    //decompositions
    //these destroy the contents of the matrix, so be careful
    DenseHED HED(); //
    DenseHED HED(Denseint numevals, char which[3]);
    DenseSVD SVD();
    std::vector<double> SingularValues();
    DenseNHED NHED();
    DenseDecompositionBase<std::complex<double> > NHEigenvalues();
  };

  inline complex<double> DenseMatrix::Value(Denseint i, Denseint j) const {
    return m_array[i+nrows*j];
  }

  inline void DenseMatrix::entry(Denseint i, Denseint j, std::complex<double> value){
    m_array[i+nrows*j]=value;
  }

  inline void move_to_dumb_array(DenseMatrix& dm, std::complex<double>* array){
    std::copy(dm.m_array,dm.m_array+dm.nrows*dm.ncols,array);
    dm.clear();
  }

  /*inline void DenseMatrix::zero(){
    for (Denseint i=0;i<nrows*ncols;++i){
      m_array[i]=0.0;
    }
    }*/

  inline DenseMatrix& DenseMatrix::operator=(DenseMatrix other){
    swap(*this, other);
    return *this;
  }

  inline void DenseMatrixToDumbArray(DenseMatrix* dmvec, std::complex<double>* array){
    //requires that all allowed entries are nonzero!!!
    move_to_dumb_array(*dmvec,array);
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// 
  template <typename T>
  class DenseDecompositionBase {
  private:
    const Denseint lineardim;
  public:
    T* Values;
    Denseint ValuesSize() const {return lineardim;}
    /*DenseDecompositionBase() : lineardim(0){
      Values=0;
      }*/
    DenseDecompositionBase(Denseint L) : lineardim(L){
      Values=new T[L];
    }

    DenseDecompositionBase(const DenseDecompositionBase& other) : lineardim(other.lineardim){
      Values=new T[lineardim];
      for (Denseint i=0;i<lineardim;++i){
	Values[i]=other.Values[i];
      }
    }

    ~DenseDecompositionBase(){
      delete[] Values;
    }
    void printValues() const {for (Denseint i=0;i<lineardim;++i){std::cout << Values[i] << " ";} std::cout << std::endl;}
  };
  //////////////////////////////////////////////////////////////////////////
  class DenseHED : public DenseDecompositionBase<double>{
  public:
    DenseMatrix EigenVectors;
    //DenseHED(){};
    ~DenseHED(){};
    DenseHED(Denseint L) : DenseDecompositionBase<double>(L) {
      EigenVectors=DenseMatrix(L,L);
    };
    DenseHED(Denseint L,Denseint numevals) : DenseDecompositionBase<double>(numevals) {
      if (numevals>L){std::cout << "Error, too many evals requested" << std::endl; exit(1);}
      EigenVectors=DenseMatrix(L,numevals);
    };

  };
  //////////////////////////////////////////////////////////////////////////
  class DenseNHED : public DenseDecompositionBase<std::complex<double> >{
  public:
    DenseMatrix RightEigenVectors;
    //DenseNHED(){};
    ~DenseNHED(){};
    DenseNHED(Denseint L) : DenseDecompositionBase<std::complex<double> >(L) {
      RightEigenVectors=DenseMatrix(L,L);
    };
  };
  //////////////////////////////////////////////////////////////////////////
  class DenseSVD :public DenseDecompositionBase<double> {
  private:
    const Denseint leftdim;
    const Denseint rightdim;
  public:
    DenseMatrix U;
    DenseMatrix Vdagger;
    ~DenseSVD(){};
    DenseSVD(Denseint M, Denseint N) : DenseDecompositionBase<double>(M>N ? N : M),leftdim(M),rightdim(N) {
      U=DenseMatrix(M,M);
      Vdagger=DenseMatrix(N,N);
    };
  };
  //////////////////////////////////////////////////////////////////////////

}
#endif

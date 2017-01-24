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
    Denseint size_;
    Denseint nrows;
    Denseint ncols;
    std::complex<double>* m_array; //default

  public:
    //When first formed matrix is allocated as transpose! When finalised it is transposed again, which orders the rows.
    DenseMatrix(Denseint param_rows, Denseint param_cols); //create array without initialising entries...
    DenseMatrix(Denseint param_rows, Denseint param_cols, std::complex<double> val ); //initialise entries
    DenseMatrix(Denseint param_rows, Denseint param_cols, std::complex<double>** array_ptr ); //initialise entries
    DenseMatrix(Denseint param_rows,Denseint param_cols, const std::complex<double>* array );
    DenseMatrix(const DenseMatrix& other) =delete; //copy constructor, must DEEP copy
    DenseMatrix(DenseMatrix&& other) noexcept;//move
    DenseMatrix() : m_array(nullptr), nrows(0), ncols(0) {};//useful when creating an instance that needs to be replaced later
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
    //DenseMatrix& operator=(DenseMatrix other);
    DenseMatrix& operator=(DenseMatrix&& other);

    DenseMatrix& operator*=(const DenseMatrix& rhs);
    //SparseMatrix sparse(const double& RELATIVETOL); //make a sparse version by dropping elements that are small relative to the largest element.
    friend void swap(DenseMatrix& first, DenseMatrix& second); //useful for copy-swap
    friend void move_to_dumb_array(DenseMatrix& dm, std::complex<double>* array);
    //decompositions
    //these destroy the contents of the matrix, so be careful
    DenseHED HED(); //
    DenseHED HED(Denseint numevals, char which[3]);
    DenseHED HED(Denseint numevals, Denseint il);

    DenseSVD SVD();
    std::vector<double> SingularValues();
    DenseNHED NHED();
    DenseDecompositionBase<std::complex<double> > NHEigenvalues();
  };

  inline std::complex<double> DenseMatrix::Value(Denseint i, Denseint j) const {
    return m_array[i+nrows*j];
  }

  inline void DenseMatrix::entry(Denseint i, Denseint j, std::complex<double> value){
    m_array[i+nrows*j]=value;
  }

  inline void move_to_dumb_array(DenseMatrix& dm, std::complex<double>* array){
    std::copy(dm.m_array,dm.m_array+dm.nrows*dm.ncols,array);
  }

  /*inline void DenseMatrix::zero(){
    for (Denseint i=0;i<nrows*ncols;++i){
      m_array[i]=0.0;
    }
    }*/

  /*inline DenseMatrix& DenseMatrix::operator=(DenseMatrix other){
    swap(*this, other);
    return *this;
    }*/

  /*inline DenseMatrix& DenseMatrix::operator=(DenseMatrix other){ //makes copy?
    swap(*this, other);
    return *this;
  }*/

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
  protected:
    Denseint lineardim;
  public:
    T* Values;
    
    DenseDecompositionBase(Denseint L) : lineardim(L){
      Values=new T[L];
    }
    DenseDecompositionBase(const DenseDecompositionBase& other) : lineardim(other.lineardim),Values(lineardim > 0 ? new T[lineardim] : nullptr){
      /*Values=new T[lineardim];
      for (Denseint i=0;i<lineardim;++i){
	Values[i]=other.Values[i];
      }*/
      std::copy(other.Values,other.Values+lineardim,Values);
    }
    DenseDecompositionBase(DenseDecompositionBase&& other) : lineardim(other.lineardim),Values(other.Values){
      other.lineardim=0;
      other.Values=nullptr;
    }
    ~DenseDecompositionBase(){delete[] Values;}
    void printValues() const {for (Denseint i=0;i<lineardim;++i){std::cout << Values[i] << " ";} std::cout << std::endl;}
    size_t ValuesSize() const {return lineardim;}
  };
  //////////////////////////////////////////////////////////////////////////
  class DenseHED : public DenseDecompositionBase<double>{
  public:
    DenseMatrix EigenVectors;
    ~DenseHED(){};
    DenseHED(Denseint L);
    DenseHED(Denseint L,Denseint numevals);
    DenseHED(DenseHED&&)=default;
  };
  //////////////////////////////////////////////////////////////////////////
  class DenseNHED : public DenseDecompositionBase<std::complex<double> >{
  public:
    DenseMatrix RightEigenVectors;
    ~DenseNHED(){};
    DenseNHED(Denseint L);
    DenseNHED(DenseNHED&&)=default;

  };
  //////////////////////////////////////////////////////////////////////////
  class DenseSVD :public DenseDecompositionBase<double> {
  private:
    Denseint leftdim;
    Denseint rightdim;
  public:
    DenseMatrix U;
    DenseMatrix Vdagger;
    DenseSVD(Denseint M, Denseint N);
    DenseSVD(DenseSVD&&)=default;

    ~DenseSVD(){};

  };
  //////////////////////////////////////////////////////////////////////////

}
#endif

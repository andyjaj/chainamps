#ifndef SMAT_H
#define SMAT_H

#include <complex>
#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include <iostream>
#include <cs.h>
#include "ajaj_common.hpp"
#include "dense_interface.hpp"
#include "arpack_interface.hpp"

#define FILEPREC 16 //used when outputting as text
#define SPARSETOL std::numeric_limits<double>::epsilon()

namespace ajaj { 
  typedef SuiteSparse_long Sparseint;
  //typedef MPXInt Sparseint;
  typedef cs_cl SparseType;
  //////////////////////////////////////////////////////////////////////////
  template <typename T>
  class SparseDecompositionBase;
  //class SparseQR;
  //class SparseLQ;
  class SparseLambdaB;
  class SparseALambda;
  class SparseSVD;
  class SparseHED;
  class SparseED;
  class SparseMatrix;
  struct SingularValues;

  template <typename T>
  class TranslationBlock;
 //////////////////////////////////////////////////////////////////////////
  std::complex<double> FixComplexPrecision(std::complex<double>& C);
  double SquareSum(const std::vector<double>& Values);
  double Sum(const std::vector<double>& Values);
  //std::vector<double>& SquareSumRescale(std::vector<double>& Values, double newsquaresum);
  double SquareSumRescale(std::vector<double>& Values, double newsquaresum);
  double entropy(const std::vector<double>& Values);
  //////////////////////////////////////////////////////////////////////////
  //friend forward declarations
  SparseMatrix NoTransMultiply(const SparseMatrix& A,const SparseMatrix& B);
  SparseMatrix addmult(const SparseMatrix& A,const SparseMatrix& B,const std::complex<double> alpha,const std::complex<double> beta); //will finalise matrices if necessary
  void swap(SparseMatrix& first, SparseMatrix& second) noexcept;
  //SparseMatrix reshape(const SparseMatrix& old,Sparseint newrows,Sparseint newcols);
  SparseMatrix reshape(const SparseMatrix& old,Sparseint newrows);
  SparseMatrix reshape(const SparseMatrix& old,Sparseint old_num_row_idxs,Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order, const bool conjugate);
  SparseMatrix cheap_reshape(const SparseMatrix& old,Sparseint old_num_row_idxs,Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order, const bool conjugate);

  SparseMatrix permute(const SparseMatrix& A, const std::vector<Sparseint>& column_permutations);
  std::vector<Sparseint> GetNonZeroColumns(const SparseMatrix& M);
  SparseMatrix load_SparseMatrix_binary(std::ifstream& infile);
  SparseMatrix copy(const SparseMatrix& other);

  SparseMatrix operator*(const SparseMatrix& lhs, const SparseMatrix& rhs);//why make the copy if it is only about to be destroyed?
  SparseMatrix operator+(const SparseMatrix& lhs, const SparseMatrix& rhs);
  SparseMatrix operator-(const SparseMatrix& lhs, const SparseMatrix& rhs);

  SparseMatrix sparse_multiply_ajaj(const SparseMatrix& lhs, const SparseMatrix& rhs, bool NoSort=0);
  SparseMatrix dense_dense_multiply(const SparseMatrix& lhs, const SparseMatrix& rhs);

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //SparseMatrix class
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  inline void index_decompose(Sparseint idx,Sparseint numdims, const Sparseint* dims_ptr,Sparseint* subidxarray);
  
  class SparseMatrix {
  private:
    SparseType* m_array;
    bool m_finalised;
    SparseMatrix(SparseType* m_array_to_use,const bool is_it_finalised) noexcept;//constructor used by internals only
    SparseMatrix(const SparseMatrix& other); //copy constructor, must DEEP copy
    SparseMatrix(Sparseint param_rows,Sparseint param_cols,Sparseint param_nzmax, bool compressed);
  public:
    //When first formed matrix is allocated as transpose! When finalised it is transposed again, which orders the rows.
    SparseMatrix(Sparseint param_rows,Sparseint param_cols,Sparseint param_nzmax); // will allocate space for the transpose!
    SparseMatrix(Sparseint param_rows,Sparseint param_cols); // will allocate space for the transpose!
    SparseMatrix(SparseMatrix&& other) noexcept : SparseMatrix() {swap(*this,other);} 
    //SparseMatrix(SparseMatrix&& other) noexcept : m_array(other.m_array), m_finalised(other.m_finalised) {other.m_array=nullptr;other.m_finalised=0;std::cout << "Sparse MOVE constructor" << std::endl;} 
    SparseMatrix(const std::vector<complex<double> >& diag,bool inverse=0);  //create diagonal matrix
    SparseMatrix(const std::vector<double>& diag,bool inverse=0);  //create diagonal matrix
    SparseMatrix() noexcept;//useful when creating an instance that needs to be replaced later
    ~SparseMatrix();

    //assignment, unary arithmetic
    SparseMatrix& operator=(SparseMatrix other);//could be and should be better
    SparseMatrix& operator*=(const SparseMatrix& rhs);
    SparseMatrix& operator+=(const SparseMatrix& rhs);
    SparseMatrix& operator-=(const SparseMatrix& rhs);

    SparseMatrix copy_transpose() const; //make a copy and transpose
    SparseMatrix copy_dagger() const; //make a copy and transpose conjugate
    //SparseMatrix copy() const {return SparseMatrix(*this);}
    Sparseint rows() const;
    Sparseint cols() const;
    Sparseint nz() const;
    bool is_finalised() const;
    bool is_dense() const; //*< check to see if the array is actually dense
    Sparseint get_p(Sparseint col) const; //gives p[col]
    Sparseint get_i(Sparseint p) const;//gives row[p]
    std::complex<double> get_x(Sparseint p) const;//gives value[p]
    Sparseint& put_p(Sparseint col);
    Sparseint& put_i(Sparseint p);
    std::complex<double>& put_x(Sparseint p);
    void set_x(const Sparseint p,complex<double> x) {m_array->x[p]=x;}
    void print_emptyness() const;
    void print_sparse_info() const {std::cout << rows() << "*" <<cols() << std::endl << "Non zeros: " << nz() << " "; print_emptyness(); std::cout << "Approx size: " << double(nz()*sizeof(std::complex<double>))/(1024.0*1024.0) << " megabytes" << std::endl; }
    double norm() const;
    double norm(const std::vector<Sparseint>& cols) const;
    double sum_column_square_norms(const std::vector<Sparseint>& cols) const;
    //array operations
    void purge(); //the purpose of this is to empty the array, but to keep the dimensions the same and free it up for new entries
    void entry(Sparseint i, Sparseint j, std::complex<double> value);
    void zero_entry(Sparseint i,Sparseint j);
    void cheap_entry(Sparseint i, Sparseint j, std::complex<double> value);

    SparseMatrix&& finalise(); //compress, then sum up duplicate entries and transpose to order rows
    SparseMatrix&& cheap_finalise(); //compress and transpose only
    SparseMatrix&& cheap_no_transpose_finalise(); //compress only
    SparseMatrix&& resize_finalise(Sparseint x, Sparseint y); //compress and reduce number of rows/cols
    SparseMatrix& transpose();
    SparseMatrix& dagger();
    SparseMatrix& conjugate();
    SparseMatrix& rescale(std::complex<double> factor);
    SparseMatrix& massage();
    SparseMatrix&& order_rows();

    //Sparseint droptol(const double tol);
    Sparseint dropzeros();
    Sparseint droprelative();
    std::complex<double> trace() const; //calculate the trace
    std::complex<double> orthogonality_check();
    std::complex<double> element(const Sparseint i, const Sparseint j) const; //very inefficient for sparse types!
    std::vector<std::complex<double> > diagonal() const;
    DenseMatrix to_dense() const;
    void to_dense(std::complex<double>* array) const;
    //void to_dense_vector(std::complex<double>* array) const;

    //read and write actions
    void print() const;
    bool fprint(std::ofstream& outfile) const;
    bool fprint_binary(std::ofstream& outfile) const;

    void loopy(void (*funcptr)(Sparseint i, Sparseint p, std::complex<double> x));
    SparseMatrix ExtractSubMatrix(const Sparseint old_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<std::pair<Sparseint,Sparseint> >& IndexVal,const bool conjugate) const;
    SparseMatrix ExtractSubMatrix(const std::pair<Sparseint,Sparseint>& RowRange, const std::pair<Sparseint,Sparseint>& ColRange) const;
    SparseMatrix ExtractColumns(const std::vector<Sparseint>& cols) const;
    SparseMatrix ExtractColumnsAndPad(const std::vector<Sparseint>& cols,Sparseint extrarows,Sparseint extracols) const;
    SparseMatrix ZeroLastColumns(Sparseint c) const; //drop all finite values in the last c columns

    SparseSVD SVD(const std::vector<std::vector<Sparseint> >& B,size_t D=0, double min_s_val=std::numeric_limits<double>::epsilon()) const;
    std::vector<double> SVD() const;
    //SparseLambdaB LambdaB(const std::vector<std::vector<Sparseint> >& B) const;
    //SparseALambda ALambda(const std::vector<std::vector<Sparseint> >& B) const;
    SparseHED HED(const std::vector<std::vector<Sparseint> >& B) const;
    SparseHED HED(const std::vector<Sparseint>& B, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseHED HED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseED LeftED(const std::vector<Sparseint>& RowB,const std::vector<Sparseint>& ColB, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseED LeftED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseED RightED(const std::vector<Sparseint>& RowB,const std::vector<Sparseint>& ColB, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseED RightED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;

    SparseED ED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;

    const SparseType* GetArrayPtr() const {return m_array;}

    friend SparseMatrix NoTransMultiply(const SparseMatrix& A,const SparseMatrix& B);
    friend SparseMatrix addmult(const SparseMatrix& A,const SparseMatrix& B,const std::complex<double> alpha,const std::complex<double> beta); //will finalise matrices if necessary
    friend void swap(SparseMatrix& first, SparseMatrix& second) noexcept;
    friend SparseMatrix reshape(const SparseMatrix& old,const Sparseint newrows);
    friend SparseMatrix reshape(const SparseMatrix& old, const Sparseint old_num_row_idxs, const Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order, const bool conjugate);
    friend SparseMatrix cheap_reshape(const SparseMatrix& old, const Sparseint old_num_row_idxs, const Sparseint new_num_row_idxs,const std::vector<Sparseint>& old_idx_dims,const std::vector<Sparseint>& new_idx_order, const bool conjugate);
    friend SparseMatrix permute(const SparseMatrix& A, const std::vector<Sparseint>& column_permutations); //column permutations don't mess up row order.
    //friend SparseMatrix permute(const SparseMatrix& A,  const std::vector<Sparseint>& row_permutations, const std::vector<Sparseint>& column_permutations);
    friend std::vector<Sparseint> GetNonZeroColumns(const SparseMatrix& M);
    friend SparseMatrix load_SparseMatrix_binary(std::ifstream& infile);
    friend SparseMatrix operator*(const SparseMatrix& lhs, const SparseMatrix& rhs);//why make the copy if it is only about to be destroyed?
    friend SparseMatrix operator+(const SparseMatrix& lhs, const SparseMatrix& rhs);
    friend SparseMatrix operator-(const SparseMatrix& lhs, const SparseMatrix& rhs);
    friend SparseMatrix copy(const SparseMatrix& other);
    friend SparseMatrix sparse_multiply_ajaj(const SparseMatrix& lhs, const SparseMatrix& rhs, bool NoSort);
    friend SparseMatrix dense_dense_multiply(const SparseMatrix& lhs, const SparseMatrix& rhs); //Assumes the arrays are actually dense, and row ordered
  };

  inline SparseMatrix NoTransMultiply(const SparseMatrix& A,const SparseMatrix& B){
    //return SparseMatrix(cs_cl_multiply(A.m_array,B.m_array),1);
    return sparse_multiply_ajaj(A,B,1); //1 flag means don't sort rows
  }

  inline Sparseint SparseMatrix::rows() const {return(m_finalised ? m_array->m : m_array->n) ;}
  inline Sparseint SparseMatrix::cols() const {return(m_finalised ? m_array->n : m_array->m) ;}
  inline Sparseint SparseMatrix::nz() const {return(m_array->nz==-1 ? m_array->p[m_array->n] : m_array->nz);}
  inline bool SparseMatrix::is_finalised() const {return m_finalised;}
  inline bool SparseMatrix::is_dense() const {return (nz() == rows() * cols());} 

  inline Sparseint SparseMatrix::get_p(Sparseint col) const {return m_array->p[col];}
  inline Sparseint SparseMatrix::get_i(Sparseint p) const {return m_array->i[p];}
  inline std::complex<double> SparseMatrix::get_x(Sparseint p) const {return m_array->x[p];}
  inline Sparseint& SparseMatrix::put_p(Sparseint col) {return m_array->p[col];}
  inline Sparseint& SparseMatrix::put_i(Sparseint p) {return m_array->i[p];}
  inline std::complex<double>& SparseMatrix::put_x(Sparseint p) {return m_array->x[p];}

  inline Sparseint SparseMatrix::dropzeros() {return cs_cl_dropzeros(m_array);}

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  // Not friends or members
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  // Data containers
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  template <typename T>
  class SparseDecompositionBase {
  protected:
    const Sparseint lineardim;
  public:
    std::vector<T> Values;
    Sparseint ValuesSize() const {return Values.size();}
    SparseDecompositionBase() : lineardim(0){}
    SparseDecompositionBase(const Sparseint& L) : lineardim(L){Values.reserve(L);}
    SparseDecompositionBase(std::vector<T>&& vals) : lineardim(vals.size()),Values(std::move(vals)){}
    void printValues() const {
      for (typename std::vector<T>::const_iterator it=Values.begin();it!=Values.end();++it){std::cout << *it << " ";} 
      std::cout << std::endl;
    }
  };

  class SparseSVD :public SparseDecompositionBase<double> {
  private:
    const Sparseint leftdim;
    const Sparseint rightdim;
    double m_kept_weight;
    double m_discarded_weight;
  public:
    SparseMatrix U;
    SparseMatrix Vdagger;
    SparseSVD(const Sparseint& M, const Sparseint& N) : SparseDecompositionBase<double>(M>N ? N : M),leftdim(M),rightdim(N),m_kept_weight(0.0),m_discarded_weight(0.0) {
    }
    SparseSVD(std::vector<double>&& vals, SparseMatrix&& umatrix, SparseMatrix&& vdmatrix, double keptweight=0.0,double discardedweight=0.0) : SparseDecompositionBase<double>(std::move(vals)),leftdim(umatrix.rows()),rightdim(vdmatrix.cols()),m_kept_weight(keptweight),m_discarded_weight(discardedweight),U(std::move(umatrix)),Vdagger(std::move(vdmatrix)){};
    //~SparseSVD(){};
    void set_weights(double total_weight,double kept_weight){m_kept_weight=kept_weight; m_discarded_weight=total_weight-kept_weight;}
    double kept_weight() const { return m_kept_weight;}
    double discarded_weight() const { return m_discarded_weight;}
  };

  class SparseHED :public SparseDecompositionBase<double> {
  private:
    const Sparseint m_dim;
  public:
    SparseMatrix EigenVectors;
    SparseHED(const Sparseint& dim, const Sparseint& numevals) : SparseDecompositionBase<double>(numevals),m_dim(dim),EigenVectors(SparseMatrix(m_dim,numevals)) {
      //EigenVectors=SparseMatrix(m_dim,numevals);
    }
    SparseHED(std::vector<double>&& values,SparseMatrix&& vecs) : SparseDecompositionBase<double>(std::move(values)),m_dim(values.size()),EigenVectors(std::move(vecs)) {}
    //~SparseHED(){};
  };

  class SparseED :public SparseDecompositionBase<complex<double> > {
  private:
    const Sparseint m_dim;
  public:
    SparseMatrix EigenVectors;
    SparseED(const Sparseint& dim, const Sparseint& numevals) : SparseDecompositionBase<std::complex<double> >(numevals),m_dim(dim),EigenVectors(SparseMatrix(m_dim,numevals)) {
      //EigenVectors=SparseMatrix(m_dim,numevals);
    }
    SparseED(SparseED&& other) noexcept : SparseDecompositionBase<std::complex<double> >(std::move(other.Values)), m_dim(other.m_dim), EigenVectors(std::move(other.EigenVectors)){}
    //~SparseED(){};
  };

  //for translating between a larger sparse matrix and a smaller sub(block) matrix
  //stupid overloads
  inline void setblock(DenseMatrix& b,const Sparseint m, const Sparseint n){b=DenseMatrix(m,n,0.0);}
  inline void setblock(SparseMatrix& b,const Sparseint m, const Sparseint n){b=SparseMatrix(m,n);}
  inline void finishblock(DenseMatrix& b){b.finalise();}
  inline void finishblock(SparseMatrix& b){b.finalise();}

  inline DenseMatrix SparseMatrix::to_dense() const{
    DenseMatrix ans(this->rows(),this->cols(),0.0);
    for (Sparseint col=0;col<this->cols();++col){
      for(Sparseint p=this->get_p(col);p<this->get_p(col+1);++p){
	ans.entry(this->get_i(p),col,this->get_x(p));
      }
    }
    return ans;
  }

  inline void SparseMatrix::to_dense(std::complex<double>* array) const{
    std::fill(array,array+rows()*cols(),0.0);
    for (Sparseint col=0;col<this->cols();++col){
      for(Sparseint p=this->get_p(col);p<this->get_p(col+1);++p){
	array[this->get_i(p)+rows()*col]=this->get_x(p);
      }
    }
  }

  template <typename T>
  class TranslationBlock{
  protected:
    //const SparseMatrix& S;
    const Sparseint orig_rows;
    const Sparseint orig_cols;
    //void FromSparse(SparseMatrix& M);
    void FromSparse(SparseMatrix&& M);

  public:
    T Block;
    std::vector<Sparseint> RowLookups; //for telling us which row of the original sparse matrix it corresponds to
    std::vector<Sparseint> ColLookups; //gets copied in by constructor
    TranslationBlock(const SparseMatrix& sp, const std::vector<Sparseint>& sparsecols) : orig_rows(sp.rows()), orig_cols(sp.cols()),ColLookups(sparsecols) {
      SparseMatrix temp(std::move(sp.ExtractColumns(ColLookups).transpose()));
      //temp.transpose();
      //now find the rows
      for (Sparseint i=0;i<temp.cols();++i){
	if (temp.get_p(i)!=temp.get_p(i+1)){this->RowLookups.push_back(i);}
      }
      this->FromSparse(std::move(temp.ExtractColumns(this->RowLookups).transpose()));
      //Block.print();
    }
    //for ED, with possible zero eigenvalues
    TranslationBlock(const SparseMatrix& sp, const std::vector<Sparseint>& sparserows,const std::vector<Sparseint>& sparsecols) : orig_rows(sp.rows()), orig_cols(sp.cols()),RowLookups(sparserows),ColLookups(sparsecols) {
      this->FromSparse(std::move(sp.ExtractColumns(ColLookups).transpose().ExtractColumns(this->RowLookups).transpose()));
      //Block.print();
    }
    SparseMatrix TranslateRows(const T& Input);
    T ReverseTranslateRows(const T& Input);
  };

  template<>
  inline void TranslationBlock<SparseMatrix>::FromSparse(SparseMatrix&& M){
    //pad with zeros for zero columns and rows
    this->Block=std::move(M);
    //std::cout << "Sparse Method" << std::endl;
  }

  template<>
  inline void TranslationBlock<DenseMatrix>::FromSparse(SparseMatrix&& M){
    //pad with zeros for zero columns and rows
    //M.print();
    this->Block=M.to_dense();
  }

  template<>
  inline SparseMatrix TranslationBlock<SparseMatrix>::TranslateRows(const SparseMatrix& Input){
    SparseMatrix ans(orig_rows,Input.cols(),Input.nz());
    if (Input.rows() > static_cast<Sparseint>(RowLookups.size())){std::cout << "Incorrect matrix dimensions for translating back!" << std::endl; exit(1);}
    for (Sparseint c=0;c<Input.cols();++c){
      //ans.Values.push_back(SpBlockans.Values.at(c));
      for (Sparseint p=Input.get_p(c);p<Input.get_p(c+1);++p){
	ans.entry(RowLookups[Input.get_i(p)],c,Input.get_x(p));
      }
    }
    return ans.finalise();
  }

  template<>
  inline SparseMatrix TranslationBlock<DenseMatrix>::TranslateRows(const DenseMatrix& Input){
    //    DenseMatrix ans(orig_rows,Input.cols());
    SparseMatrix ans(orig_rows,Input.cols());
    if (Input.rows() > static_cast<Sparseint>(RowLookups.size())){std::cout << "Incorrect matrix dimensions for translating back!" << std::endl; exit(1);}
    for (Denseint c=0;c<Input.cols();++c){
      for (Denseint r=0;r<Input.rows();++r){
	if (abs(Input.Value(r,c))>=SPARSETOL)
	  ans.entry(RowLookups[r],c,Input.Value(r,c));
      }
    }
    return ans.finalise();
  }

  template<>
  inline SparseMatrix TranslationBlock<SparseMatrix>::ReverseTranslateRows(const SparseMatrix& Input){
    return SparseMatrix(std::move(Input.copy_transpose().ExtractColumns(RowLookups).transpose()));
  }

  std::pair<SparseMatrix,SparseMatrix> SparseHermitianDecomposition(DenseHED decomp);

  SparseMatrix Exponentiate(const SparseHED& decomp, std::complex<double> factor);

  inline void ConvertSparseVecToDenseVec(SparseMatrix* sparsevec, std::complex<double>* densevec){
    sparsevec->to_dense(densevec);
  }

  struct SparseVectorWithRestriction{
    const SparseMatrix* SPtr;
    const std::vector<Sparseint>* AllowedPtr;
    SparseVectorWithRestriction(const SparseMatrix* s, const std::vector<Sparseint>* a) : SPtr(s), AllowedPtr(a) {if (SPtr) if (SPtr->cols()!=1) {std::cout << "Not a sparse VECTOR!" << std::endl; exit(1);}}
  };

  //Given a sparse vector and a list of allowed entries,
  //produce a reduced vector and pad any allowed empty entries in the sparse vector with zero
  void DumbExtractWithZeros(const SparseMatrix& S,const std::vector<Sparseint>& Allowed, std::complex<double>* out);
  void DumbExtractUpdate(const SparseMatrix& S,const std::vector<Sparseint>& Allowed, std::complex<double>* out, std::complex<double> weight=std::complex<double>(1.0,0.0));

  DenseMatrix ExtractWithZeros(const SparseMatrix& S,const std::vector<Sparseint>& Allowed);

  inline void ConvertSparseVectorWithRestriction(SparseVectorWithRestriction* swr, std::complex<double>* array){
    DumbExtractWithZeros(*(swr->SPtr),*(swr->AllowedPtr),array);
  }

}


#endif

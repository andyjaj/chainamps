/** @file MPX_matrix.hpp
 * MPX_matrix is the base class for MPS and MPO matrices.
 * Essentially a container for a SparseMatrix, with lots of bells and whistles so that they can be blocked by good quantum numbers, and reshaped according
 * to their indices.
 */
#ifndef MPX_MATRIX_H
#define MPX_MATRIX_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp"

namespace ajaj{

  class MPX_matrix;
  class MPXDecomposition;

  MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs);
  template <typename T>
  TranslationBlock<T> contract_conditional(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB,const std::vector<MPXPair>& contractidxs,const std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> >& condition);
  MPX_matrix reorder(const MPX_matrix& A, bool conjA, const std::vector<Sparseint>& newindexorder, const Sparseint numrows); //reorders indices (and reshape sparse array)
  std::pair<MPX_matrix,MPX_matrix> MakeXandXinv(const MPX_matrix& M);
  std::pair<std::vector<double>,MPX_matrix> SqrtDR(const MPX_matrix& M); /**<Uses eigendecomposition, so R^dagger is inverse of R*/
  SparseMatrix reshape_to_vector(const MPX_matrix& A);
  SparseMatrix contract_to_sparse(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs);
  void swap(MPX_matrix& A, MPX_matrix& B);
  MPX_matrix copy(const MPX_matrix& A);
  /** Sometimes after a decomposition we need to figure out the new charges of all the rows, which we can do this using the columns and the sparse structure.*/
  StateArray DeduceFromColumns(const SparseMatrix& array, const MPXIndex& colindices);
  /** Load an MPX_matrix from binary storage format*/
  MPX_matrix load_MPX_matrix_binary(const std::string& filename,const EigenStateArray& spectrum);
  MPX_matrix load_MPX_matrix_binary(std::ifstream& infile,const EigenStateArray& spectrum);
  MPX_matrix load_MPX_matrix(const std::string& filename,const EigenStateArray& spectrum);

  /** BlockStateIndices are used to label a bunch of container indices, (e.g. columns of a matrix, positions in a vector) that correspond to the same charge values (quantum numbers).
      Useful for looping over blocks.*/
  struct BlockStateIndices {
    State BlockState; /**< The actual values of the charges for this group of indices. */
    std::vector<MPXInt> Indices; /**< The indices (e.g. columns) that have this set of charges*/
    BlockStateIndices(const State& is) : BlockState(is) {};
    BlockStateIndices(const State& is,MPXInt c) : BlockState(is) {Indices.push_back(c);}
    BlockStateIndices(const State& is,const std::vector<Sparseint>& idxs) : BlockState(is),Indices(idxs) {}; /**< Constructor using a list of indices.*/
    MPXInt IndicesSize() const {return Indices.size();} /**< How many different indices are there in this instance? */
    void print() const{std::cout << BlockState << ": "; for (std::vector<MPXInt>::const_iterator cit=Indices.begin();cit!=Indices.end();++cit){std::cout << *cit << " ";} std::cout << std::endl;} /**< Print the Quantum numbers and all the indices. */
  };

  //////////////////////////////////
  //////////////////////////////////
  //MPX_matrix class
  //////////////////////////////////
  //////////////////////////////////

  /** MPX_matrix is the base tensor class. It uses a SparseMatrix for storage, knows about the physical spectrum, and keeps record of all of the quantum numbers corresponding to particular values of its MPXIndex indices.*/
  class MPX_matrix {
  protected: //inherited objects need to be able to see these
    const Basis* m_SpectrumPtr; //shouldn't change. Using a ptr to a const means we can use the default equals provided by the compiler
    std::vector<MPXIndex> m_Indices;
    MPXInt m_NumRowIndices;
    SparseMatrix m_Matrix;
  public:
    //1st arg is a ref to the vertex spectrum,
    //2nd arg is a ref to where all the state arrays are stored
    MPX_matrix() noexcept : m_SpectrumPtr(nullptr) {}; 
    MPX_matrix(const Basis& spectrum) noexcept; /**< All tensors need to know about the physical spectrum. This is a design flaw.*/
    MPX_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, const Sparseint numrowindices, const SparseMatrix& matrix); /**< Constructor, indices are listed in order from left to right. numrowindices lets the object know how many of the indices are associated with rows of the SparseMatrix. */
    MPX_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, const Sparseint numrowindices, SparseMatrix&& matrix);
    MPX_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<std::complex<double> >& values,bool inverse=0); /**< Diagonal constructor, with possible inversion of values.*/
    MPX_matrix(const Basis& spectrum, const MPXIndex& Lindex, const MPXIndex& Rindex, const std::vector<std::complex<double> >& values,bool inverse=0); /**< Diagonal constructor, with different left and right indices and possible inversion of values.*/
    MPX_matrix(const Basis& spectrum, const MPXIndex& Lindex, const MPXIndex& Rindex, const std::vector<double>& values,bool inverse=0); /**< Diagonal constructor, with different left and right indices and possible inversion of values.*/

    MPX_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse=0);
    MPX_matrix(MPX_matrix&& rhs) noexcept : m_SpectrumPtr(rhs.m_SpectrumPtr) {swap(*this,rhs);}
    MPX_matrix(const MPX_matrix& rhs) : m_SpectrumPtr(rhs.m_SpectrumPtr), m_Indices(rhs.m_Indices), m_NumRowIndices(rhs.m_NumRowIndices), m_Matrix(copy(rhs.m_Matrix)){}

    bool isEmpty() const {return !(m_SpectrumPtr);} /**< Check to see if null.*/
    bool empty() const {return !(m_SpectrumPtr);} /**< Check to see if null.*/


    const MPXIndex& Index(Sparseint i) const {return m_Indices.at(i);} /**< Lookup the StateArray corresponding to a particular index, interface to storage.*/
    const Basis& GetPhysicalSpectrum() const {return *m_SpectrumPtr;} /**< Return ref to the Physical spectrum, needs renaming to match convention */
    const Basis& getPhysicalSpectrum() const {return *m_SpectrumPtr;} /**< Return ref to the Physical spectrum */
    const Basis& basis() const {return *m_SpectrumPtr;}
    SparseMatrix matrix() {return copy(m_Matrix);}
    bool isConsistent() const; /**< Check dimensions of SparseMatrix match the dimensions of the indices.*/
    bool isHermitian() const;
    std::vector<Sparseint> dimsvector() const; /**< Return a vector containing all the dimensions of the MPXIndex indices, from left to right. */
    void print_matrix() const; /**< Print the SparseMatrix */
    void print_indices() const; /**< Print the dimensions of the indices. Colon indicates how many correspond to rows and how many to columns. */
    void print_indices_values() const;
    void print_sparse_info() const {m_Matrix.print_sparse_info();}
    bool fprint_binary(std::ofstream& outfile) const; /**< Print the MPX_matrix to file in binary format. */
    bool store(const std::string& filename) const; /**< Print the MPX_matrix to file in binary format. */
    //MPX_matrix ExtractSubMPX(const std::vector<MPXPair>& Indexvals) const; /**< Extract a Sub MPX_matrix by giving some indices fixed values according to the pairs in Indexvals . */
    std::vector<BlockStateIndices> GetAllBlockColumns() const; /**< Return a vector of BlockStateIndices corresponding to all the groups of columns with the same quantum numbers.*/
    BlockStateIndices GetBlockColumns(const State& specficstate) const; /**< Return BlockStateIndices for the quantum numbers specified by specificstate. */
    BlockStateIndices GetBlockRows(const State& specficstate) const; /**< Return BlockStateIndices for the quantum numbers specified by specificstate. */
    State CombinationState(std::vector<MPXInt> indices, MPXInt l) const; /**< Combine indices into a single superindex (in order), then find the quantum numbers corresponding to the l-th element. */
    State OutgoingColState(MPXInt j) const; /**< Return the outgoing charges corresponding to column j.*/
    State IngoingRowState(MPXInt i) const; /**< Return the ingoing charges corresponding to row i.*/
    SparseHED Eigs() const; /** Do a full eigenvalue decomposition.*/
    SparseHED Eigs(Sparseint numevals,char which[3],SparseMatrix* initial=NULL) const; /**< Get the numevals eigenvalues and vectors corresponding to which (LM,SM,LR,SR). Can take an optional pointer to an initial guess vector.*/ 
    SparseHED Eigs(const std::vector<Sparseint>& cols, Sparseint numevals,char which[3],SparseMatrix* initial=NULL) const;/**< Considering the block defined by cols, get the numevals eigenvalues and vectors corresponding to which (LM,SM,LR,SR). Can take an optional pointer to an initial guess vector.*/ 
    SparseHED Eigs(const State& blockstate) const; /**< Do a full eigenvalue decomposition for the block corresponding to blockstate. This means all the rows and cols with charges==blockstate, formed into a block. */
    SparseHED Eigs(const State& blockstate, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const; /** Do a decomposition for which numevals eigenvals and vectors, for a specific State. */
    SparseED LeftEigs(const State& blockstate, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const; /**< Left eigen decomposition, for a specific block defined by blockstate. */
    SparseED LeftEigs(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const; /**< Left eigen decomposition. */
    SparseED RightEigs(const State& blockstate, Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;/**< Right eigen decomposition, for a specific block defined by blockstate. */
    SparseED RightEigs(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;/**< Right eigen decomposition. */
    MPXDecomposition SVD(size_t bond_dimension=0,double min_s_val=0.0) const; /**< Singular value decomposition with possible truncation. Row indices become rows of U, col indices become cols of Vdagger. */
    
    std::vector<std::complex<double> > GetDiagonal() {if (!isEmpty()) return m_Matrix.diagonal(); else return std::vector<std::complex<double> >();}
    
    MPX_matrix& Transpose();
    MPX_matrix& RemoveDummyIndices(std::vector<MPXInt> indices_for_removal); /**< Strip off dummy indices. Useful for simplifying open b.c. ends */
    MPX_matrix& MoveDummyIndices(const std::vector<MPXInt>& newindexorder); /**< Move dummy indices only. */
    MPX_matrix& ShiftNumRowIndices(const Sparseint numrows); /**< Move some row MPXIndex indices to columns or vice-versa. */
    MPX_matrix& Rescale(std::complex<double> factor); /**<Scale all values in array by some factor. */
    complex<double> Trace() const {return m_Matrix.trace();} /**< If MPX_matrix is square, find the trace.*/
    MPX_matrix ZeroLastBlock(); /**< Given a particular index and a value of that index, zero the other data*/
    MPX_matrix& CombineSimilarMatrixIndices(bool PhysicalInMiddle=0);
    
    //move assignment operator needed?
    MPX_matrix& operator=(MPX_matrix other){
      swap(*this,other);
      return *this;
    }

    friend bool arrays_equal(const MPX_matrix& A, const MPX_matrix& B,double tol);

    friend MPX_matrix copy(const MPX_matrix& A);
    friend MPX_matrix reorder(const MPX_matrix& A, bool conjA, const std::vector<Sparseint>& newindexorder, const Sparseint numrows); /**< Reorders MPXIndex indices, and reshapes SparseMatrix accordingly. */
    friend MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs); /**< Contracts two MPX_matrix objects using the index pairs given in contractidxs. */
    template <typename T>
    friend TranslationBlock<T> contract_conditional(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB,const std::vector<MPXPair>& contractidxs,const std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> >& condition); /**< Contracts two MPX_matrix objects using the index pairs given in contractidxs, but only include in the result if condition is met. The combined charges of the MPXIndex indices listed in condition must be zero. Useful when forming a superblock Hamiltonian.*/
    friend std::pair<MPX_matrix,MPX_matrix> MakeXandXinv(const MPX_matrix& M); /**< Used only if left orthogonalising, without need for the singular values.*/
    friend std::pair<std::vector<double>,MPX_matrix> SqrtDR(const MPX_matrix& M); /**< Used to decompose into Xdagger sqrt(D) sqrt(D) X form, when orthogonalising. */
    friend SparseMatrix reshape_to_vector(const MPX_matrix& A); /**< Reshape the matrix into a vector (i.e. a single column SparseMatrix). */
    friend SparseMatrix contract_to_sparse(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs); /**< Contract, but only keep the SparseMatrix, not the MPXIndex info. */
    friend void swap(MPX_matrix& A, MPX_matrix& B);
  }; 

  /** Used for storing the results of decompositions. */
  class MPXDecompositionBase {
  public:
    double Truncation;
    MPX_matrix ColumnMatrix;
    std::vector<double> Values;
    MPXDecompositionBase(const Basis& basis) : Truncation(0.0),ColumnMatrix(basis){};
    MPXDecompositionBase(MPX_matrix&& cm, std::vector<double>&& v, double Trunc=0.0) noexcept : Truncation(Trunc),ColumnMatrix(std::move(cm)), Values(std::move(v)){};

    const std::vector<double>& SquareRescale(double sqsum) {
      Truncation/=sqsum;
      SquareSumRescale(Values,sqsum);
      return Values;
    }

    void printValues() const{
      double weight(0.0);
      for (auto v : Values){
	std::cout << v << std::endl;
	weight+=v*v;
      }
      std::cout << "Total squared weight " << weight<<std::endl;
    }

    const Basis& basis() const {
      return ColumnMatrix.basis();
    }

  };

  /** Used for storing the results of decompositions, when both rows and columns are needed. For example in an SVD, ColumnMatrix == U. RowMatrix == Vdagger. */
  class MPXDecomposition : public MPXDecompositionBase {
  public:
    MPX_matrix RowMatrix;
    MPXDecomposition(const Basis& basis) : MPXDecompositionBase(basis), RowMatrix(basis){};
    MPXDecomposition(MPX_matrix&& cm, std::vector<double>&& v, MPX_matrix&& rm, double Trunc=0.0) noexcept : MPXDecompositionBase(std::move(cm),std::move(v),Trunc), RowMatrix(std::move(rm)) {};
  };


  inline bool arrays_equal(const MPX_matrix& A, const MPX_matrix& B,double tol=0.0){
    return check_equal(A.m_Matrix,B.m_Matrix,tol);
  }

  /** Overload for contract, if only one index needs to be contracted. */
  inline MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const MPXPair& contractidx){
    return contract(A,conjA,B,conjB,std::vector<MPXPair>(1,contractidx));
  }

  template <typename T>
  inline TranslationBlock<T> contract_conditional(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB,const std::vector<MPXPair>& contractidxs,const std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> >& condition){
    //first check spectra match
    if (A.m_SpectrumPtr!=B.m_SpectrumPtr){std::cout <<"Physical Indices don't match!" << std::endl;exit(1);}
    //next check contraction indices match sizes match!
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      if (!match(A.m_Indices.at(cit->first),conjA,B.m_Indices.at(cit->second),conjB)){exit(1);}
    }
    //make list of old index order (0 to nA-1 and 0 to nB-1)
    std::vector<Sparseint> neworderA(A.m_Indices.size());
    std::iota(neworderA.begin(),neworderA.end(),0);
    std::vector<Sparseint> neworderB(B.m_Indices.size());
    std::iota(neworderB.begin(),neworderB.end(),0);
    //shift the condition indices to the right
    if (condition.first.begin()==condition.first.end() || condition.second.begin()==condition.second.end()){
      std::cout << "Not enough condition indices?" << std::endl; exit(1);
    }
    std::vector<MPXInt> Adims=A.dimsvector();
    std::vector<MPXInt> Bdims=B.dimsvector();
    if (contractidxs.begin()==contractidxs.end()){
      std::cout << "Not enough contraction indices?" << std::endl; exit(1);
    }
    Sparseint AConditionColumns=1;
    for (std::vector<MPXInt>::const_iterator cit=condition.first.begin();cit!=condition.first.end();++cit){
      *(std::remove(neworderA.begin(),neworderA.end(),abs(*cit)))=abs(*cit);
      AConditionColumns*=Adims.at(*cit);
    }
    Sparseint BConditionColumns=1;
    for (std::vector<MPXInt>::const_iterator cit=condition.second.begin();cit!=condition.second.end();++cit){
      *(std::remove(neworderB.begin(),neworderB.end(),abs(*cit)))=abs(*cit);
      BConditionColumns*=Bdims.at(*cit);
    }
    //shift contraction indices to the right
    Sparseint ContractColumns=1; //there's got to be at least one?
    for (std::vector<MPXPair>::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      *(std::remove(neworderA.begin(),neworderA.end(),cit->first))=cit->first;
      *(std::remove(neworderB.begin(),neworderB.end(),cit->second))=cit->second;
      ContractColumns*=Adims.at(cit->first);
      //container size doesn't change!
    }
    //reshape sparsematrices
    Sparseint Anewrowidxnum=A.m_Indices.size()-contractidxs.size()-condition.first.size();
    Sparseint Bnewrowidxnum=B.m_Indices.size()-contractidxs.size()-condition.second.size();
    SparseMatrix Ar(reshape(A.m_Matrix,A.m_NumRowIndices,Anewrowidxnum,Adims,neworderA,conjA));
    SparseMatrix Br(reshape(B.m_Matrix,B.m_NumRowIndices,Bnewrowidxnum,Bdims,neworderB,conjB));
    SparseMatrix fullsparse(Ar.rows()*Br.rows(),AConditionColumns*BConditionColumns);
    std::vector<MPXInt> extractioncolumns;
    //loop over sub columns of B, decompose into condition indices and remaining contraction indices.
    for (Sparseint Bcdcol=0;Bcdcol<BConditionColumns;++Bcdcol){
      //find charges corresponding to Bcdcol 
      State BCondState(B.CombinationState(condition.second,Bcdcol));
      for (Sparseint Acdcol=0;Acdcol<AConditionColumns;++Acdcol){
	//complete the condition check, and sum over contraction indices if check is true
	if (BCondState==-A.CombinationState(condition.first,Acdcol)){
	  extractioncolumns.push_back(Acdcol+AConditionColumns*Bcdcol);
	  for (Sparseint contractcol=0; contractcol<ContractColumns;++contractcol){
	  //form the actual column of Br
	    Sparseint Bcol=Bcdcol+BConditionColumns*contractcol;
	    if (Br.get_p(Bcol)!=Br.get_p(Bcol+1)){ //not an empty column
	      //form the actual column of Ar
	      Sparseint Acol=Acdcol+AConditionColumns*contractcol;
	      if (Ar.get_p(Acol)!=Ar.get_p(Acol+1)){ //don't mess about with empty columns
		for (Sparseint Bp=Br.get_p(Bcol);Bp<Br.get_p(Bcol+1);++Bp){ //rows of Br
		  for (Sparseint Ap=Ar.get_p(Acol);Ap<Ar.get_p(Acol+1);++Ap){ //rows of Ar
		    //no need to conjugate as the reshapes did that already
		    fullsparse.entry(Ar.get_i(Ap)+Ar.rows()*Br.get_i(Bp),Acdcol+AConditionColumns*Bcdcol,Ar.get_x(Ap)*Br.get_x(Bp));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    fullsparse.finalise();
    //NEEDS to figure out correct columns by charges!!!!!

    return TranslationBlock<T>(fullsparse,extractioncolumns,extractioncolumns);
  }

}
#endif

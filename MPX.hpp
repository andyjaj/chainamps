/** @file MPX.hpp
 * MPX_matrix is the base class for MPS and MPO matrices.
 * Essentially a container for a SparseMatrix, with lots of bells and whistles so that they can be blocked by good quantum numbers, and reshaped according
 * to their indices.
 */
#ifndef MPX_H
#define MPX_H

/** @file MPX.hpp
 * Defines the tensor class and related objects.
 */
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>

#include "ajaj_common.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"

namespace ajaj {
  struct BlockStateIndices;
  class MPXIndex;
  class MPX_matrix;
  class MPO_matrix;
  class MPS_matrix;
  class MPXDecompositionBase;
  class MPXDecomposition;
  class MPSDecomposition;

  //friends for namespace resolution
  bool match(const MPXIndex& A, bool conjA, const MPXIndex& B, bool conjB);
  MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs);
  template <typename T>
  TranslationBlock<T> contract_conditional(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB,const std::vector<MPXPair>& contractidxs,const std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> >& condition);
  MPX_matrix reorder(const MPX_matrix& A, bool conjA, const std::vector<Sparseint>& newindexorder, const Sparseint numrows); //reorders indices (and reshape sparse array)
  std::pair<MPX_matrix,MPX_matrix> MakeXandXinv(const MPX_matrix& M);
  std::pair<std::vector<double>,MPX_matrix> SqrtDR(const MPX_matrix& M);
  SparseMatrix reshape_to_vector(const MPX_matrix& A);
  SparseMatrix contract_to_sparse(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs);
  //void swap(MPXIndex& A, MPXIndex& B);
  void swap(MPX_matrix& A, MPX_matrix& B);
  MPXIndex combine(const MPXIndex& a,const MPXIndex& b);
  void MPS_swap(MPS_matrix& A, MPS_matrix& B);
  MPX_matrix copy(const MPX_matrix& A);


  /** BlockStateIndices are used to label a bunch of container indices, (e.g. columns of a matrix, positions in a vector) that correspond to the same charge values (quantum numbers).*/
  struct BlockStateIndices {
    State BlockState; /**< The actual values of the charges for this group of indices. */
    std::vector<MPXInt> Indices; /**< The indices (e.g. columns) that have this set of charges*/
    BlockStateIndices(const State& is) : BlockState(is) {};
    BlockStateIndices(const State& is,MPXInt c) : BlockState(is) {Indices.push_back(c);}
    BlockStateIndices(const State& is,const std::vector<Sparseint>& idxs) : BlockState(is),Indices(idxs) {}; /**< Constructor using a list of indices.*/
    MPXInt IndicesSize() const {return Indices.size();} /**< How many different indices are there in this instance? */
    void print() const{std::cout << BlockState << ": "; for (std::vector<MPXInt>::const_iterator cit=Indices.begin();cit!=Indices.end();++cit){std::cout << *cit << " ";} std::cout << std::endl;} /**< Print the Quantum numbers and all the indices. */
  };

  /** A directed (ingoing or outgoing) tensor index. Contains a StateArray corresponding to a matrix index, or a reference to an EigenStateArray if it is a physical index. The directional property is necessary to sort out conservation of quantum numbers. */
  class MPXIndex{
  private:
    bool m_isInward;
    bool m_isPhysical;
    StateArray m_IndexStates; //might be empty
    const StateArray* m_IndexStatesPtr; //ptr to the states
  public:
    MPXIndex(bool inward,const StateArray& indexstates) : m_isInward(inward),m_isPhysical(0),m_IndexStates(indexstates),m_IndexStatesPtr(&m_IndexStates){}; /**< Construct a matrix index */
    MPXIndex(bool inward,StateArray&& indexstates) : m_isInward(inward),m_isPhysical(0),m_IndexStates(indexstates),m_IndexStatesPtr(&m_IndexStates){}; /**< Construct a matrix index */

    MPXIndex(bool inward,const EigenStateArray& spectrum) : m_isInward(inward),m_isPhysical(1),m_IndexStatesPtr(&spectrum){}; /**< Construct a physical index. */
    MPXIndex(const MPXIndex& other) : m_isInward(other.m_isInward),m_isPhysical(other.m_isPhysical),m_IndexStates(other.m_IndexStates),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){}; /**< Copy constructor */
    MPXIndex(MPXIndex&& other) noexcept : m_isInward(other.m_isInward),m_isPhysical(other.m_isPhysical),m_IndexStates(std::move(other.m_IndexStates)),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){};/**< Move constructor */
    MPXIndex(bool inward,const MPXIndex& other) : m_isInward(inward),m_isPhysical(other.m_isPhysical),m_IndexStates(other.m_IndexStates),m_IndexStatesPtr(m_isPhysical ? other.m_IndexStatesPtr : &m_IndexStates){}; /**< Specify direction and copy states from another index. */
    //make a trivial MPXIndex with 1 state entry
    MPXIndex(const MPXIndex& other,MPXInt i) : m_isInward(other.m_isInward),m_isPhysical(0),m_IndexStates(StateArray(1,other.at(i))),m_IndexStatesPtr(&m_IndexStates){}; /** Construct a 'dummy' index using just one State from element i of other.*/

    bool Ingoing() const {return m_isInward;} /**< Returns true if index is Ingoing*/
    bool Outgoing() const {return !Ingoing();}
    bool Physical() const {return m_isPhysical;} /**< Returns true if index is Physical */
    State at(Sparseint i) const {return m_IndexStatesPtr->at(i);} /**< Returns a copy of the i-th State in the container. */
    const State& operator[](Sparseint i) const {return (*m_IndexStatesPtr)[i];}

    Sparseint size() const {return m_IndexStatesPtr->size();} /**< How many States are in the container (what is the dimension of the index?) */
    void print() const {for (StateArray::const_iterator cit=m_IndexStatesPtr->begin();cit!=m_IndexStatesPtr->end();++cit){std::cout <<*cit << std::endl;}} /**< Print all the State objects corresponding to the index */
    bool fprint_binary(std::ofstream& outfile) const;

    // MPXIndex& operator=(MPXIndex rhs){swap(*this,rhs);return *this;}
    MPXIndex& operator=(const MPXIndex& rhs){
      m_isInward=rhs.m_isInward;
      m_isPhysical=rhs.m_isPhysical;
      m_IndexStates=rhs.m_IndexStates;//let std::vector default do the work
      m_IndexStatesPtr= m_isPhysical ? rhs.m_IndexStatesPtr : &m_IndexStates;
      return *this;}
    MPXIndex& operator=(MPXIndex&& rhs){
      m_isInward=rhs.m_isInward;
      m_isPhysical=rhs.m_isPhysical;
      m_IndexStates=std::move(rhs.m_IndexStates);//let std::vector default do the work
      m_IndexStatesPtr= m_isPhysical ? rhs.m_IndexStatesPtr : &m_IndexStates;
      return *this;
    }

    MPXIndex flip() const {return MPXIndex((!this->m_isInward),*this);} /**< Switch the direction of the index */

    friend MPXIndex combine(const MPXIndex& a,const MPXIndex& b);

    friend void swap(MPXIndex& A, MPXIndex& B);
    friend bool match(const MPXIndex& A, bool conjA,const MPXIndex& B,bool conjB); /**< Check to see if two indices match up when taking a contraction. One should be Ingoing and one should be Outgoing.*/
  };

  /** MPX_matrix is the base tensor class. It uses a SparseMatrix for storage, knows about the physical spectrum, and keeps record of all of the quantum numbers corresponding to particular values of its MPXIndex indices.*/
  class MPX_matrix {
  protected: //inherited objects need to be able to see these
    const EigenStateArray* m_SpectrumPtr; //shouldn't change. Using a ptr to a const means we can use the default equals provided by the compiler
    std::vector<MPXIndex> m_Indices;
    MPXInt m_NumRowIndices;
    SparseMatrix m_Matrix;
  public:
    //1st arg is a ref to the vertex spectrum,
    //2nd arg is a ref to where all the state arrays are stored
    MPX_matrix() noexcept : m_SpectrumPtr(nullptr) {}; 
    MPX_matrix(const EigenStateArray& spectrum) noexcept; /**< Default constructor. All tensors need to know about the physical spectrum. This might be a design flaw.*/
    MPX_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, const Sparseint numrowindices, const SparseMatrix& matrix); /**< Constructor, indices are listed in order from left to right. numrowindices lets the object know how many of the indices are associated with rows of the SparseMatrix. */
    MPX_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, const Sparseint numrowindices, SparseMatrix&& matrix);
    MPX_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values,bool inverse=0); /**< Diagonal constructor, with possible inversion of values.*/
    MPX_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse=0);
    MPX_matrix(MPX_matrix&& rhs) noexcept : m_SpectrumPtr(rhs.m_SpectrumPtr) {swap(*this,rhs);}
    MPX_matrix(const MPX_matrix& rhs) : m_SpectrumPtr(rhs.m_SpectrumPtr), m_Indices(rhs.m_Indices), m_NumRowIndices(rhs.m_NumRowIndices), m_Matrix(copy(rhs.m_Matrix)){}

    const MPXIndex& Index(Sparseint i) const {return m_Indices.at(i);} /**< Lookup the StateArray corresponding to a particular index, interface to storage.*/
    const EigenStateArray& GetPhysicalSpectrum() const {return *m_SpectrumPtr;} /**< Return ref to the Physical spectrum, needs renaming to match convention */
    const EigenStateArray& getPhysicalSpectrum() const {return *m_SpectrumPtr;} /**< Return ref to the Physical spectrum */
    const Basis& basis() const {return *m_SpectrumPtr;}
    bool isConsistent() const; /**< Check dimensions of SparseMatrix match the dimensions of the indices.*/
    std::vector<Sparseint> dimsvector() const; /**< Return a vector containing all the dimensions of the MPXIndex indices, from left to right. */
    void print_matrix() const; /**< Print the SparseMatrix */
    void print_indices() const; /**< Print the dimensions of the indices. Colon indicates how many correspond to rows and how many to columns. */
    void print_indices_values() const;
    void print_sparse_info() const {m_Matrix.print_sparse_info();}
    bool fprint_binary(std::ofstream& outfile) const; /**< Print the MPX_matrix to file in binary format. */
    bool store(const std::string& filename) const; /**< Print the MPX_matrix to file in binary format. */
    MPX_matrix ExtractSubMPX(const std::vector<MPXPair>& Indexvals) const; /**< Extract a Sub MPX_matrix by giving some indices fixed values according to the pairs in Indexvals . */
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
    MPX_matrix& Transpose();
    MPX_matrix& RemoveDummyIndices(std::vector<MPXInt> indices_for_removal); /**< Strip off dummy indices. Useful for simplifying open b.c. ends */
    MPX_matrix& MoveDummyIndices(const std::vector<MPXInt>& newindexorder); /**< Move dummy indices only. */
    MPX_matrix& ShiftNumRowIndices(const Sparseint numrows); /**< Move some row MPXIndex indices to columns or vice-versa. */
    MPX_matrix& Rescale(std::complex<double> factor); /**<Scale all values in array by some factor. */
    complex<double> Trace() const {return m_Matrix.trace();} /**< If MPX_matrix is square, find the trace.*/
    MPX_matrix RestrictColumnIndex();
    MPX_matrix& CombineSimilarMatrixIndices(bool PhysicalInMiddle=0);



    //void swap(MPX_matrix& other);

    //move assignment operator needed?
    MPX_matrix& operator=(MPX_matrix other){
      swap(*this,other);
      return *this;
    }

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

  /** Inherited class from MPX_matrix. An MPX_matrix with two matrix indices and two physical indices.*/
  class MPO_matrix : public MPX_matrix {
  private:
    void check();
  public:
    MPO_matrix() noexcept : MPX_matrix() {}; 
    MPO_matrix(const EigenStateArray& spectrum);
    MPO_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, const SparseMatrix& matrix);
    MPO_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, SparseMatrix&& matrix);

    //MPO_matrix(const MPX_matrix& MPXref);
    MPO_matrix(MPX_matrix&& MPXref) noexcept;
    MPO_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values,bool inverse=0);
    MPO_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse=0);

    MPO_matrix ExtractMPOBlock(const std::pair<MPXInt,MPXInt>& row_matrix_index_range, const std::pair<MPXInt,MPXInt>& col_matrix_index_range) const;

  };

  struct NamedMPO_matrix {
  public:
    std::string Name;
    ajaj::MPO_matrix Matrix;
    NamedMPO_matrix(const std::string& n,ajaj::MPO_matrix&& m) : Name(n),Matrix(m) {}
  };

  /** Create an identity MPO_matrix. Corresponds to an identity operator.*/
  inline MPO_matrix IdentityMPO_matrix(const EigenStateArray& spectrum){
    MPXIndex dummy(1,ajaj::StateArray(1,spectrum[0]-spectrum[0]));
    return MPO_matrix(spectrum,dummy,std::vector<complex<double> >(spectrum.size(),1.0));  //create diagonal matrix
  }

  MPO_matrix UnitaryTransformMPO_matrix(const Basis&, const std::vector<MPXIndex>&, const SparseMatrix&, size_t, double);


  /** Inherited class from MPX_matrix. An MPX_matrix with two matrix indices and one physical index.*/
  class MPS_matrix : public MPX_matrix {
  private:
    void check();
  public:
    MPS_matrix(const EigenStateArray& spectrum);
    MPS_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, const SparseMatrix& matrix);
    MPS_matrix(const EigenStateArray& spectrum, const std::vector<MPXIndex>& indices, SparseMatrix&& matrix);

    MPS_matrix(MPX_matrix&& MPXref) noexcept;

    const MPXIndex& getInwardMatrixIndex() const {
      const MPXIndex* Indexptr(nullptr);
      for (auto&& i : m_Indices){
	if (i.Ingoing() && !i.Physical()) {Indexptr=&i; break;}
      }
      return *Indexptr;
    }
    const MPXIndex& getOutwardMatrixIndex() const {
      const MPXIndex* Indexptr(nullptr);
      for (auto&& i : m_Indices){
	if (i.Outgoing() && !i.Physical()) {Indexptr=&i; break;}
      }
      return *Indexptr;
    }
    const MPXIndex& getPhysicalIndex() const {
      const MPXIndex* Indexptr(nullptr);
      for (auto&& i : m_Indices){
	if (i.Physical()) {Indexptr=&i; break;}
      }
      return *Indexptr;
    }

    MPXInt InwardMatrixIndexNumber() const {
      for (uMPXInt i=0; i< m_Indices.size() ;++i){
	if (m_Indices[i].Ingoing() && !m_Indices[i].Physical()) {return i;}
      }
      return -1;//should never happen!
    }

    MPXInt OutwardMatrixIndexNumber() const {
      for (uMPXInt i=0; i< m_Indices.size() ;++i){
	if (m_Indices[i].Outgoing() && !m_Indices[i].Physical()) {return i;}
      }
      return -1;
    }
    MPXInt PhysicalIndexNumber() const {
      for (uMPXInt i=0; i< m_Indices.size() ;++i){
	if (m_Indices[i].Physical()) {return i;}
      }
      return -1;
    }

    MPS_matrix&& left_shape() {
      if (m_Indices[0].Physical() && m_Indices[1].Ingoing() && m_Indices[2].Outgoing()){
	return std::move(*this);
      }
      else { //assume it was right shaped
	std::vector<MPXInt> olddims(dimsvector());
	swap(m_Indices[0],m_Indices[1]);
	m_Matrix=std::move(reshape(m_Matrix,m_NumRowIndices,2,olddims,reorder102,0));
	m_NumRowIndices=2;
	return std::move(*this);
      }
    }

    MPS_matrix&& right_shape() {
      if (m_Indices[1].Physical() && m_Indices[0].Ingoing() && m_Indices[2].Outgoing()){
	return std::move(*this);
      }
      else { //assume it was left shaped
	std::vector<MPXInt> olddims(dimsvector());
	swap(m_Indices[0],m_Indices[1]);
	m_Matrix=std::move(reshape(m_Matrix,m_NumRowIndices,1,olddims,reorder102,0));
	m_NumRowIndices=1;
	return std::move(*this);
      }
    }

    friend void MPS_swap(MPS_matrix& A, MPS_matrix& B);
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
  };

  /** Used for storing the results of decompositions, when both rows and columns are needed. For example in an SVD, ColumnMatrix == U. RowMatrix == Vdagger. */
  class MPXDecomposition : public MPXDecompositionBase {
  public:
    MPX_matrix RowMatrix;
    MPXDecomposition(const Basis& basis) : MPXDecompositionBase(basis), RowMatrix(basis){};
    MPXDecomposition(MPX_matrix&& cm, std::vector<double>&& v, MPX_matrix&& rm, double Trunc=0.0) noexcept : MPXDecompositionBase(std::move(cm),std::move(v),Trunc), RowMatrix(std::move(rm)) {};
    
  };

  /** Used for storing the results of decompositions, such as SVD (Schmidt decomposition) that produce two new MPS_matrix objects.*/
  class MPSDecomposition {
  private: 
  public:
    double Truncation;
    std::vector<double> Values;
    MPS_matrix LeftMatrix;
    MPS_matrix RightMatrix;
    MPSDecomposition(MPXDecomposition&& X) noexcept : Truncation(X.Truncation), Values(std::move(X.Values)), LeftMatrix(std::move(X.ColumnMatrix)), RightMatrix(std::move(X.RowMatrix)){};
    MPSDecomposition(const EigenStateArray& spectrum) : Truncation(0.0), Values(), LeftMatrix(spectrum), RightMatrix(spectrum){};
    const std::vector<double>& SquareRescale(double sqsum) {
      Truncation/=sqsum;
      SquareSumRescale(Values,sqsum);
      return Values;
    }
    void printValues() const {
      for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
	std::cout << *cit << std::endl;
      }
    }

    bool store(const std::string& LeftName, const std::string& RightName) const;
    bool store(const std::string& Name, uMPXInt nl,uMPXInt nr) const;
    bool store_left(const std::string& Name, uMPXInt nl) const;
    bool store_right(const std::string& Name,uMPXInt nr) const;

    std::string LeftName(const std::string& s, uMPXInt n) const {
      std::stringstream leftname;
      leftname << s.c_str() << "_Left_" << n << ".MPS_matrix";
      return leftname.str();
    }
    std::string RightName(const std::string& s, uMPXInt n) const {
      std::stringstream rightname;
      rightname << s.c_str() << "_Right_" << n << ".MPS_matrix";
      return rightname.str();
    }

    void OutputPhysicalIndexDensities(std::ofstream& d) const;
  };

  /** A class to hold a unit cell of arbitrary length
   *
   *
   */
  class UnitCell {
  protected:
    const Basis* basis_ptr_;
  public:
    std::vector<MPS_matrix> Matrices;
    std::vector<std::vector<double> > Lambdas;
    //storage format is 'left canonical' based i.e. A_0 A_1 (or lambda_0 Gamma_0 lambda_1 Gamma_1...)
    //even if the matrices are not in fact left canonical
    //so a decomposition is of the form lambda_0 Gamma_0 lambda_1 Gamma_1 lambda_0
    UnitCell()=delete;
    UnitCell(const Basis& spectrum) : basis_ptr_(&spectrum){}
    UnitCell(const Basis& spectrum,const MPS_matrix& m, const std::vector<double>& l,uMPXInt length=1) : basis_ptr_(&spectrum), Matrices(std::vector<MPS_matrix>(length,m)),Lambdas(std::vector<std::vector<double> >(length,l)) {}
    UnitCell(const Basis& spectrum,std::vector<MPS_matrix>&& mvec, std::vector<std::vector<double> >&& lvec) : basis_ptr_(&spectrum), Matrices(mvec),Lambdas(lvec) {}
    UnitCell(const MPSDecomposition& D,const std::vector<double>& PreviousLambda) : basis_ptr_(&D.LeftMatrix.GetPhysicalSpectrum()) {
      Matrices.emplace_back(copy(D.LeftMatrix));
      Matrices.emplace_back(reorder(contract(contract(MPX_matrix(D.LeftMatrix.GetPhysicalSpectrum(),D.RightMatrix.Index(0),D.Values),0,D.RightMatrix,0,contract10),0,MPX_matrix(D.LeftMatrix.GetPhysicalSpectrum(),D.RightMatrix.Index(2),PreviousLambda,1),0,contract20),0,reorder102,2));
      Lambdas.emplace_back(PreviousLambda);
      Lambdas.emplace_back(D.Values);
    }

    void swap(uMPXInt i, uMPXInt j){MPS_swap(Matrices.at(i),Matrices.at(j));std::swap(Lambdas.at(i),Lambdas.at(j));}
    
    /*UnitCell& operator=(UnitCell&& other){
      basis_ptr_=other.basis_ptr_;
      Matrices=std::move(other.Matrices);
      Lambdas=std::move(other.Lambdas);
      return *this;
      }*/

    uMPXInt size() const {return Matrices.size();}
    const MPS_matrix& GetMatrix(size_t i) {return Matrices.at(i);}
    const std::vector<double>& GetLambda(size_t i) {return Lambdas.at(i);}
    void OutputPhysicalIndexDensities(std::ofstream& d) const;
    void OutputOneVertexDensityMatrix(std::ofstream& d) const;
    void OutputOneVertexDensityMatrix(const std::string& name, uMPXInt l) const;
    
    void store(const std::string& filename, uMPXInt l) const; /**< Print UnitCell with to binary file, with an index in filename. */
    void store(const std::string& filename) const; /**< Print the UnitCell to file in binary format. */

    double Entropy() const;
  };

  /** Sometimes after a decomposition we need to figure out the new charges of all the rows, which we can do using the columns and the sparse structure.*/
  StateArray DeduceFromColumns(const SparseMatrix& array, const MPXIndex& colindices);
  /** Load a MPXIndex from binary storage format*/
  MPXIndex load_MPXIndex_binary(std::ifstream& infile,const EigenStateArray& spectrum);
  /** Load a MPX_matrix from binary storage format*/
  MPX_matrix load_MPX_matrix_binary(std::string& filename,const EigenStateArray& spectrum);
  MPX_matrix load_MPX_matrix_binary(std::ifstream& infile,const EigenStateArray& spectrum);
  MPX_matrix load_MPX_matrix(const std::string& filename,const EigenStateArray& spectrum);
  MPO_matrix load_MPO_matrix(const std::string& filename,const EigenStateArray& spectrum);
  MPS_matrix load_MPS_matrix(const std::string& filename,const EigenStateArray& spectrum);
  UnitCell load_UnitCell_binary(std::ifstream& infile, QNVector& charge_rules, Basis& basis);
  UnitCell load_UnitCell_binary(std::ifstream& infile, const QNVector& charge_rules, const Basis& basis);
  MPS_matrix MakeProductState(const EigenStateArray& spectrum, const std::vector<std::pair<uMPXInt,double> >& state_index_vec,State leftstate);
  MPS_matrix MakeProductState(const EigenStateArray& spectrum, uMPXInt state_index,State leftstate);
}

#endif

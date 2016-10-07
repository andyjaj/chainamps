#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <numeric>

#include "ajaj_common.hpp"
#include "sparse_interface.hpp"
#include "arpack_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"

namespace ajaj {
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //friends
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void swap(MPXIndex& A, MPXIndex& B){

    std::swap(A.m_isInward,B.m_isInward);
    std::swap(A.m_isPhysical,B.m_isPhysical);
    std::swap(A.m_IndexStates,B.m_IndexStates);
    std::swap(A.m_IndexStatesPtr,B.m_IndexStatesPtr);
    /*const StateArray* oldAPtr=A.m_IndexStatesPtr;
      const StateArray* oldBPtr=B.m_IndexStatesPtr;*/

    A.m_IndexStatesPtr=A.m_isPhysical ? A.m_IndexStatesPtr : &A.m_IndexStates;
    //oops!!
    B.m_IndexStatesPtr=B.m_isPhysical ? B.m_IndexStatesPtr : &B.m_IndexStates;

  }

  bool match(const MPXIndex& A, bool conjA, const MPXIndex& B, bool conjB){
    bool flag=1;
    if (A.Physical()!=B.Physical()){flag=0;} //both must be of same type
    else if (((A.Outgoing()!=conjA )!=(B.Ingoing()!=conjB))){flag=0;} //A should be outgoing, B should be ingoing
    else if(A.m_IndexStatesPtr->size()!=B.m_IndexStatesPtr->size()){flag=0;} //sizes should match
    if (flag==0){
      std::cout <<"Indices don't match!"<<std::endl;
      std::cout << "Conjugation flags: " << conjA << " " << conjB << std::endl;
      std::cout << "Physical? " << A.Physical() << " " << B.Physical() << std::endl;
      std::cout << "Index Dir? " << A.Outgoing() << " " << B.Ingoing() << std::endl;
      std::cout << "Sizes: " << A.size() << " " << B.size() << std::endl;
    }
#ifndef NDEBUG
    else {
      for (uMPXInt i=0;i<A.m_IndexStatesPtr->size();++i){
	if (A.m_IndexStatesPtr->at(i)!=B.m_IndexStatesPtr->at(i)){
	  std::cout << "Indices don't match because charges are not equal" << std::endl;
	  std::cout << A.m_IndexStatesPtr->at(i) << std::endl;
	  std::cout << B.m_IndexStatesPtr->at(i) << std::endl;
	  exit(1);
	}
      }
    }
#endif
    return flag;
  }

  MPXIndex combine(const MPXIndex& a,const MPXIndex& b){
    StateArray S;

    const StateArray& Ref_to_a_states=*(a.m_IndexStatesPtr);
    const StateArray& Ref_to_b_states=*(b.m_IndexStatesPtr);

    if (a.Ingoing()!=b.Ingoing()){ //different directions...
      for (auto&& b_it : Ref_to_b_states){
	for (auto&& a_it : Ref_to_a_states){
	  S.push_back(a_it-b_it);
	}
      }
    }
    else {
      for (auto&& b_it : Ref_to_b_states){
	for (auto&& a_it : Ref_to_a_states){
	  S.push_back(a_it+b_it);
	}
      }
    }

    return MPXIndex(a.Ingoing(),std::move(S));
  }


  void swap(MPX_matrix& A, MPX_matrix& B){
    std::swap(A.m_SpectrumPtr,B.m_SpectrumPtr);
    swap(A.m_Indices,B.m_Indices);
    std::swap(A.m_NumRowIndices,B.m_NumRowIndices);
    swap(A.m_Matrix,B.m_Matrix);
  }

  void MPS_swap(MPS_matrix& A, MPS_matrix& B){
    std::swap(A.m_SpectrumPtr,B.m_SpectrumPtr);
    swap(A.m_Indices,B.m_Indices);
    std::swap(A.m_NumRowIndices,B.m_NumRowIndices);
    swap(A.m_Matrix,B.m_Matrix);
  }

  MPX_matrix copy(const MPX_matrix& A){
    return MPX_matrix(*(A.m_SpectrumPtr),A.m_Indices,A.m_NumRowIndices,A.m_Matrix);
  }

  //reduces all non physical indices of same direction to one index
  //MPX_matrix indexreduce(const MPX_matrix& A)

  MPX_matrix reorder(const MPX_matrix& A, bool conjA, const std::vector<Sparseint>& newindexorder, const Sparseint numrows){ //reorders indices (and reshape sparse array)
    if (newindexorder.size()!=A.m_Indices.size()){std::cout << "Incorrect MPX reshape params" << std::endl; exit(1);}
    std::vector<MPXIndex> Indices;
    for (std::vector<Sparseint>::const_iterator cit=newindexorder.begin();cit!=newindexorder.end();++cit){ //loop through remaining A indices
      //need ot flip indices if conjugation has occurred
      Indices.push_back(MPXIndex(conjA ? A.m_Indices.at(*cit).flip() : A.m_Indices.at(*cit))); 
    }
    return MPX_matrix(*(A.m_SpectrumPtr),Indices,numrows,reshape(A.m_Matrix,A.m_NumRowIndices,numrows,A.dimsvector(),newindexorder,conjA));
  }

  SparseMatrix contract_to_sparse(const MPX_matrix& A, bool conjA,const MPX_matrix& B,bool conjB, const std::vector<MPXPair>& contractidxs){
#ifndef NDEBUG
    std::cout << "Contract to sparse" << std::endl;
#endif
    if (A.m_SpectrumPtr!=B.m_SpectrumPtr){std::cout <<"Physical Indices don't match!" << std::endl;exit(1);}
    //next check contraction indices match sizes match!
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      if (!match(A.m_Indices.at(cit->first),conjA,B.m_Indices.at(cit->second),conjB)){exit(1);}
    }
    //entries in contract indices tell us which we need to make 'inner'
    //other order should be preserve
    //make list of old index order (0 to nA-1 and 0 to nB-1)
    std::vector<Sparseint> neworderA(A.m_Indices.size());
    std::iota(neworderA.begin(),neworderA.end(),0);
    std::vector<Sparseint> neworderB(B.m_Indices.size());
    std::iota(neworderB.begin(),neworderB.end(),0);
    //we have no guarantee the elements in contractidxs are ordered, even if pair::first is, pair::second might not be
    //therefore only safe thing to use is std::remove which searches for an element and replaces it with the next element
    //end of container is left valid but undefined
    //we then copy in the correct value at the end
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      *(std::remove(neworderA.begin(),neworderA.end(),cit->first))=cit->first;
      *(std::remove(neworderB.begin(),neworderB.end(),cit->second))=cit->second;
      //container size doesn't change!
    }
    //confusingly we move the inner indices to the left of A and the right of B
    //this is because we then only need one transpose to order the rows of the final sparsematrix
    //Last contractidxs.size() elements of A need to become the first, without messing up their individual order!
#if defined(USETBB)
    std::rotate(neworderB.begin(),neworderB.end()-contractidxs.size(),neworderB.end());
    return std::move(reshape(A.m_Matrix,A.m_NumRowIndices,neworderA.size()-contractidxs.size(),A.dimsvector(),neworderA,conjA)*reshape(B.m_Matrix,B.m_NumRowIndices,contractidxs.size(),B.dimsvector(),neworderB,conjB));
#else
    std::rotate(neworderA.begin(),neworderA.end()-contractidxs.size(),neworderA.end());
    return std::move(NoTransMultiply(reshape(B.m_Matrix,B.m_NumRowIndices,neworderB.size()-contractidxs.size(),B.dimsvector(),neworderB,conjB),reshape(A.m_Matrix,A.m_NumRowIndices,contractidxs.size(),A.dimsvector(),neworderA,conjA)).transpose());
#endif
  }

  MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B, bool conjB, const std::vector<MPXPair >&  contractidxs){
#ifndef NDEBUG
    std::cout << "Contract" << std::endl;
#endif
    //first check spectra match
    if (A.m_SpectrumPtr!=B.m_SpectrumPtr){std::cout <<"Physical Indices don't match!" << std::endl;exit(1);}
    //next check contraction indices match sizes match!
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      if (!match(A.m_Indices.at(cit->first),conjA,B.m_Indices.at(cit->second),conjB)){
	exit(1);}
    }
    //entries in contract indices tell us which we need to make 'inner'
    //other order should be preserve
    //make list of old index order (0 to nA-1 and 0 to nB-1)
    std::vector<Sparseint> neworderA(A.m_Indices.size(),0);
    std::iota(neworderA.begin(),neworderA.end(),0);
    std::vector<Sparseint> neworderB(B.m_Indices.size(),0);
    std::iota(neworderB.begin(),neworderB.end(),0);
    
    //we have no guarantee the elements in contractidxs are ordered, even if pair::first is, pair::second might not be
    //therefore only safe thing to use is std::remove which searches for an element and replaces it with the next element
    //end of container is left valid but undefined
    //we then copy in the correct value at the end
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      *(std::remove(neworderA.begin(),neworderA.end(),cit->first))=cit->first;
      *(std::remove(neworderB.begin(),neworderB.end(),cit->second))=cit->second;
      //container size doesn't change!
    }
#if defined(USETBB)
    std::rotate(neworderB.begin(),neworderB.end()-contractidxs.size(),neworderB.end());
    std::vector<MPXIndex> Indices;
    for (std::vector<Sparseint>::const_iterator cit=neworderA.begin();cit!=neworderA.end()-contractidxs.size();++cit){ //loop through remaining A indices
      Indices.emplace_back(conjA ? A.m_Indices.at(*cit).flip() : A.m_Indices.at(*cit));
    }
    for (std::vector<Sparseint>::const_iterator cit=neworderB.begin()+contractidxs.size();cit!=neworderB.end();++cit){ //loop through remaining B indices
      Indices.emplace_back(conjB ? B.m_Indices.at(*cit).flip() : B.m_Indices.at(*cit));
    }
    return MPX_matrix(*(A.m_SpectrumPtr),Indices,neworderA.size()-contractidxs.size(),std::move(reshape(A.m_Matrix,A.m_NumRowIndices,neworderA.size()-contractidxs.size(),A.dimsvector(),neworderA,conjA)*reshape(B.m_Matrix,B.m_NumRowIndices,contractidxs.size(),B.dimsvector(),neworderB,conjB)));

    //return MPX_matrix(*(A.m_SpectrumPtr),Indices,neworderA.size()-contractidxs.size(),std::move(NoTransMultiply(reshape(A.m_Matrix,A.m_NumRowIndices,neworderA.size()-contractidxs.size(),A.dimsvector(),neworderA,conjA),reshape(B.m_Matrix,B.m_NumRowIndices,contractidxs.size(),B.dimsvector(),neworderB,conjB)).order_rows()));

#else
    //confusingly we move the inner indices to the left of A and the right of B
    //this is because we then only need one transpose to order the rows of the final sparsematrix
    //Last contractidxs.size() elements of A need to become the first, without messing up their individual order!
    std::rotate(neworderA.begin(),neworderA.end()-contractidxs.size(),neworderA.end());
    //form new indices
    //preserves order of non contracted indices
    std::vector<MPXIndex> Indices;
    for (std::vector<Sparseint>::const_iterator cit=neworderA.begin()+contractidxs.size();cit!=neworderA.end();++cit){ //loop through remaining A indices
      Indices.emplace_back(conjA ? A.m_Indices.at(*cit).flip() : A.m_Indices.at(*cit));
    }
    for (std::vector<Sparseint>::const_iterator cit=neworderB.begin();cit!=neworderB.end()-contractidxs.size();++cit){ //loop through remaining B indices
      Indices.emplace_back(conjB ? B.m_Indices.at(*cit).flip() : B.m_Indices.at(*cit));
    }
    return MPX_matrix(*(A.m_SpectrumPtr),Indices,neworderA.size()-contractidxs.size(),std::move(NoTransMultiply(reshape(B.m_Matrix,B.m_NumRowIndices,neworderB.size()-contractidxs.size(),B.dimsvector(),neworderB,conjB),reshape(A.m_Matrix,A.m_NumRowIndices,contractidxs.size(),A.dimsvector(),neworderA,conjA)).transpose()));
#endif

  }

  SparseMatrix reshape_to_vector(const MPX_matrix& A){
    return reshape(A.m_Matrix,A.m_Matrix.rows()*A.m_Matrix.cols());
  }

  std::pair<MPX_matrix,MPX_matrix> MakeXandXinv(const MPX_matrix& M){
    //check square, check indices match
    if (M.m_Matrix.rows()!=M.m_Matrix.cols()){std::cout << "MPX_matrix not square, can't decompose!" << std::endl; exit(1);}
    SparseHED EigDecomp(M.Eigs());
    //form SparseMatrices
    std::vector<double> roots;
    std::vector<double> inverseroots;
    for (std::vector<double>::const_iterator cit=EigDecomp.Values.begin();cit!=EigDecomp.Values.end();++cit){
      if (*cit<=0.0){
	std::cout << "Negative eigenvalue: " << *cit << " Matrix isn't positive semi definite, can't decompose into X and X inverse!" << std::endl;
	for (std::vector<double>::const_iterator it=EigDecomp.Values.begin();it!=EigDecomp.Values.end();++it){
	  std::cout << "All eigenvalues" << std::endl;
	  std::cout << *cit << std::endl;
	}
	exit(1);
      }
      roots.push_back(sqrt(*cit));
      inverseroots.push_back(1.0/sqrt(*cit));
    }
    //there is no reason at this stage to expect the new row MPXIndices to match the cols
    //we have to deduce them
    std::vector<MPXIndex> XIndices;
    XIndices.emplace_back(1,DeduceFromColumns(EigDecomp.EigenVectors.dagger(),M.Index(0)));
    XIndices.emplace_back(0,M.Index(0));
    MPX_matrix X(M.GetPhysicalSpectrum(),XIndices,1,SparseMatrix(roots)*EigDecomp.EigenVectors);
    //and here we reuse the deduced index
    std::vector<MPXIndex> XinvIndices;
    XinvIndices.emplace_back(1,M.Index(0));
    XinvIndices.emplace_back(0,X.Index(0));
    return std::pair<MPX_matrix,MPX_matrix>(std::move(X),MPX_matrix(M.GetPhysicalSpectrum(),XinvIndices,1,std::move(EigDecomp.EigenVectors.dagger()*SparseMatrix(inverseroots))));
  }

  std::pair<std::vector<double>,MPX_matrix> SqrtDR(const MPX_matrix& M){
    //Decompose M into Rdagger D R, where D is diagonal
    //return vector for Sqrt(D) and return an MPX_matrix for R with correct indices
    if (M.m_Matrix.rows()!=M.m_Matrix.cols()){std::cout << "MPX_matrix not square, can't decompose!" << std::endl; exit(1);}
    SparseHED EigDecomp(M.Eigs());

    //form SparseMatrices
    std::vector<double> roots;
    //first find largest eval



    //first check the absolute magnitude of the summed evals
    double abssum=0.0;
    for (std::vector<double>::const_iterator it=EigDecomp.Values.begin();it!=EigDecomp.Values.end();++it){
      abssum+=abs(*it);
    }
    //get average magnitude
    abssum*=((1.0e-10)/double(EigDecomp.Values.size()));
    uMPXInt num_finite=0;
    uMPXInt num_negative=0;

    for (std::vector<double>::const_iterator cit=EigDecomp.Values.begin();cit!=EigDecomp.Values.end();++cit){
      if (abs(*cit)<=abssum) {
	roots.push_back(0.0);
      }
      else {
	roots.push_back(sqrt(abs(*cit)));
	++num_finite;
	if(*cit<0.0) ++num_negative; 
      }
    }
    if (num_negative!=0 && num_negative!=num_finite){
      std::cout << "Negative eigenvalue, can't decompose into Rdagger D R!" << std::endl;
      for (std::vector<double>::const_iterator it=EigDecomp.Values.begin();it!=EigDecomp.Values.end();++it){
	std::cout << "All eigenvalues:" << std::endl;
	std::cout << *it << std::endl;
      }
      return std::pair<std::vector<double>,MPX_matrix>(std::vector<double>(),MPX_matrix());
    }
    
    std::vector<MPXIndex> XIndices;
    XIndices.emplace_back(1,DeduceFromColumns(EigDecomp.EigenVectors.dagger(),M.Index(0)));
    XIndices.emplace_back(0,M.Index(0));
    return std::pair<std::vector<double>,MPX_matrix>(roots,MPX_matrix(M.GetPhysicalSpectrum(),XIndices,1,EigDecomp.EigenVectors));
  }

  bool MPXIndex::fprint_binary(std::ofstream& outfile) const{
    if (outfile.is_open()){
      outfile.write(reinterpret_cast<const char*>(&m_isInward),sizeof(bool));
      outfile.write(reinterpret_cast<const char*>(&m_isPhysical),sizeof(bool));
      size_t indexsize(m_IndexStatesPtr->size());
      outfile.write(reinterpret_cast<const char*>(&indexsize),sizeof(size_t));
      if (!m_isPhysical && m_IndexStatesPtr->size()>0){	
	for (StateArray::const_iterator cit=m_IndexStates.begin();cit!=m_IndexStates.end();++cit){
	  cit->fprint_binary(outfile);
	}
      }
      return 0;
    }
    else {
      return 1;
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPX members
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPX_matrix::MPX_matrix(const EigenStateArray& spectrum) noexcept : m_SpectrumPtr(&spectrum){};
  MPX_matrix::MPX_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, const Sparseint numrowindices, const SparseMatrix& matrix) : m_SpectrumPtr(&spectrum),m_Indices(indices),m_NumRowIndices(numrowindices),m_Matrix(copy(matrix)){
    if (!m_Matrix.is_finalised()){std::cout << "Initialising MPX_matrix from unfinalised SparseMatrix is not allowed!" << std::endl; exit(1);}
  };

  MPX_matrix::MPX_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, const Sparseint numrowindices, SparseMatrix&& matrix) : m_SpectrumPtr(&spectrum),m_Indices(indices),m_NumRowIndices(numrowindices),m_Matrix(std::move(matrix)){
    if (!m_Matrix.is_finalised()){std::cout << "Initialising MPX_matrix from unfinalised SparseMatrix is not allowed!" << std::endl; exit(1);}
  };

  MPX_matrix::MPX_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values, bool inverse): m_SpectrumPtr(&spectrum) {
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=1;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPX_matrix::MPX_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse): m_SpectrumPtr(&spectrum) {
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=1;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPX_matrix& MPX_matrix::ShiftNumRowIndices(const Sparseint numrows){ //reorders indices (and reshape sparse array)
#ifndef NDEBUG
    std::cout << "Shift Num Row Indices" << std::endl;
#endif
    if (static_cast<size_t>(numrows)>this->m_Indices.size() || numrows <1){std::cout << "Incorrect MPX reshape params for change rows" << std::endl; exit(1);}
    MPXInt newrowdim=1;
    for (std::vector<MPXIndex>::const_iterator cit=this->m_Indices.begin();cit!=this->m_Indices.begin()+numrows;++cit){
      newrowdim*=cit->size();
    }
    //std::cout << newrowdim << std::endl;
    this->m_Matrix=reshape(this->m_Matrix,newrowdim);
    this->m_NumRowIndices=numrows;
    return *this;
  }

  MPX_matrix& MPX_matrix::MoveDummyIndices(const std::vector<MPXInt>& newindexorder){
    if (newindexorder.size()!=m_Indices.size()){std::cout << "Incorrect MPX reshape params" << std::endl; exit(1);}
    //record the non dummy indices and form the new index vector
    std::vector<MPXInt> nondummyorder;
    std::vector<MPXIndex> Indices;
    for (std::vector<Sparseint>::const_iterator nit=newindexorder.begin();nit!=newindexorder.end();++nit){
      if(m_Indices.at(*nit).size()>1){nondummyorder.push_back(*nit);}
      if (*nit==m_NumRowIndices){
	m_NumRowIndices=std::distance(static_cast<std::vector<Sparseint>::const_iterator>(newindexorder.begin()),nit);
      }
      Indices.push_back(m_Indices.at(*nit));
    }
    //ensure the non dummy indices are not out of order
    for (std::vector<Sparseint>::const_iterator ndit=nondummyorder.begin()+1;ndit!=nondummyorder.end();++ndit){ 
      if(!(*ndit>*(ndit-1))){std::cout << "Non dummy order shouldn't change!" << std::endl; exit(1);}//
    }
    std::swap(m_Indices,Indices);
    return *this;
  }
  
  SparseHED MPX_matrix::Eigs() const {
    //only if we want full decomposition
    //but still uses blocks to reduce the size of the problem
    const std::vector<BlockStateIndices>& blocks=GetAllBlockColumns();
    std::vector<std::vector<MPXInt> > columns;
    for (std::vector<BlockStateIndices>::const_iterator cit=blocks.begin();cit!=blocks.end();++cit){
      columns.push_back(cit->Indices);
    }
    return m_Matrix.HED(columns);
  }

  SparseHED MPX_matrix::Eigs(Sparseint numevals,char which[3],SparseMatrix* initial) const {
    //only if we want full decomposition
    return m_Matrix.HED(numevals,which,initial);
  }

  SparseHED MPX_matrix::Eigs(const std::vector<Sparseint>& cols, Sparseint numevals,char which[3],SparseMatrix* initial) const {
    return m_Matrix.HED(cols,numevals,which,initial);
  }

  SparseHED MPX_matrix::Eigs(const State& blockstate) const { //all evals and vecs in a sector
    const BlockStateIndices& block=GetBlockColumns(blockstate);
    return m_Matrix.HED(std::vector<std::vector<MPXInt> >(1,block.Indices));
  }

  SparseHED MPX_matrix::Eigs(const State& blockstate, Sparseint numevals,char which[3],SparseMatrix* initial) const { //get specific number of which evals in a particular sector
    const BlockStateIndices& block=GetBlockColumns(blockstate);
    return m_Matrix.HED(block.Indices,numevals,which,initial);
  }

  SparseED MPX_matrix::LeftEigs(const State& blockstate, Sparseint numevals,char which[3],SparseMatrix* initial) const { //get specific number of which evals in a particular sector
    const BlockStateIndices& blockcols=GetBlockColumns(blockstate);
    const BlockStateIndices& blockrows=GetBlockRows(blockstate);
    return m_Matrix.LeftED(blockrows.Indices,blockcols.Indices,numevals,which,initial);
  }

  SparseED MPX_matrix::LeftEigs(Sparseint numevals,char which[3],SparseMatrix* initial) const { //get specific number of which evals
    return m_Matrix.LeftED(numevals,which,initial);
  }

  SparseED MPX_matrix::RightEigs(const State& blockstate, Sparseint numevals,char which[3],SparseMatrix* initial) const { //get specific number of which evals in a particular sector
    const BlockStateIndices& blockcols=GetBlockColumns(blockstate);
    const BlockStateIndices& blockrows=GetBlockRows(blockstate);
    return m_Matrix.RightED(blockrows.Indices,blockcols.Indices,numevals,which,initial);
  }

  SparseED MPX_matrix::RightEigs(Sparseint numevals,char which[3],SparseMatrix* initial) const { //get specific number of which evals in a particular sector
    return m_Matrix.RightED(numevals,which,initial);
  }

  MPXDecomposition MPX_matrix::SVD(size_t bond_dimension,double min_s_val) const{
#ifndef NDEBUG
    std::cout << "MPX_matrix.SVD()" << std::endl;
#endif
    std::vector<BlockStateIndices> Blocks(this->GetAllBlockColumns());
    std::vector<std::vector<Sparseint> > JustColumns;
    for (std::vector<BlockStateIndices>::const_iterator cit=Blocks.begin();cit!=Blocks.end();++cit){
      JustColumns.emplace_back(cit->Indices);
    }
#ifndef NDEBUG
    std::cout << "SparseMatrix.SVD()" << std::endl;
#endif
    SparseSVD decomposition(m_Matrix.SVD(JustColumns,bond_dimension,min_s_val));
    //still need to produce new MPX_matrices
    //what is the newly formed MPX_index?
    StateArray newstates;
#ifndef NDEBUG
    std::cout << "Figure out new charges" << std::endl;
#endif
    for (Sparseint c=0;c<decomposition.U.cols();++c){
      //columns better not be empty!
      if (decomposition.U.get_p(c)==decomposition.U.get_p(c+1)) {std::cout << "Empty column! Malformed SVD" << std::endl; exit(1);}
      newstates.emplace_back(IngoingRowState(decomposition.U.get_i(decomposition.U.get_p(c)))); //inward row charges = outward col charges
    }
    //push first num_row_indices then add new index on the end
    std::vector<MPXIndex> leftindices(m_Indices.begin(),m_Indices.begin()+m_NumRowIndices);
    leftindices.emplace_back(0,newstates);
    std::vector<MPXIndex> rightindices;
    rightindices.reserve(m_Indices.size()-m_NumRowIndices+1);
    rightindices.emplace_back(1,newstates);
    rightindices.insert(rightindices.end(),m_Indices.begin()+m_NumRowIndices,m_Indices.end());
#ifndef NDEBUG
    std::cout << "return SVD" << std::endl;
#endif				     
    return MPXDecomposition(MPX_matrix(this->GetPhysicalSpectrum(),leftindices,m_NumRowIndices,std::move(decomposition.U)),std::move(decomposition.Values),MPX_matrix(this->GetPhysicalSpectrum(),rightindices,1,std::move(decomposition.Vdagger)),decomposition.discarded_weight());
  }

  MPX_matrix& MPX_matrix::Transpose(){
    m_Matrix.transpose();
    uMPXInt new_num_row_indices=m_Indices.size()-m_NumRowIndices;
    std::rotate(m_Indices.begin(),m_Indices.begin()+m_NumRowIndices,m_Indices.end());
    for (std::vector<MPXIndex>::iterator it=m_Indices.begin();it!=m_Indices.end();++it){
      if (!it->Physical()){/* flip*/
	it->flip();
      }
    }
    m_NumRowIndices=new_num_row_indices;
    return *this;
  }

  std::vector<Sparseint> MPX_matrix::dimsvector() const {
    std::vector<Sparseint> ans;
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.end();++cit){
      ans.push_back(cit->size());
    }
    return ans;
  }

  void MPX_matrix::print_matrix() const {
    m_Matrix.print();
  }

  bool MPX_matrix::fprint_binary(std::ofstream& outfile) const {
    if (outfile.is_open()){
      size_t numindices(m_Indices.size());
      //print total number of indices and number of row indices
      outfile.write(reinterpret_cast<const char*>(&numindices),sizeof(size_t));
      outfile.write(reinterpret_cast<const char*>(&m_NumRowIndices),sizeof(Sparseint));
      //print the MPXIndex objects
      for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.end();++cit){
	cit->fprint_binary(outfile);
      }
      //print the sparse matrix  
       m_Matrix.fprint_binary(outfile);
      return 0;
    }
    else {
      return 1;
    }
  }

  bool MPX_matrix::store(const std::string& filename) const {
    std::ofstream outfile;
    outfile.open(filename.c_str(),ios::out | ios::trunc | ios::binary);  
    if (outfile.is_open()){
      return fprint_binary(outfile);
    }
    else {
      return 1; //error
    }
  }

  void MPX_matrix::print_indices() const {
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.begin()+m_NumRowIndices;++cit){
      std::cout << (cit->Ingoing() ? "+" : "-");
      std::cout << cit->size() << " ";
    }
    std::cout << ": ";
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin()+m_NumRowIndices;cit!=m_Indices.end();++cit){
      std::cout << (cit->Ingoing() ? "+" : "-");
      std::cout << cit->size()<< " ";
    }
    std::cout << std::endl;
  }

  void MPX_matrix::print_indices_values() const {
    uMPXInt count(0);
    for (auto&& I : m_Indices){
      std::cout << "Index " << count << ", Direction " << (I.Ingoing() ? "+" : "-") <<std::endl;
      I.print();
    }
  }

  bool MPX_matrix::isConsistent() const {
    MPXInt rows=1;
    MPXInt cols=1;
    for (size_t i=0;i<static_cast<size_t>(m_NumRowIndices);++i){
      rows*=m_Indices.at(i).size();
    }
    for (size_t i=m_NumRowIndices;i<m_Indices.size();++i){
      cols*=m_Indices.at(i).size();
      
    }
    if (rows!=m_Matrix.rows()){std::cout << "Incorrect row dimensions: " << rows << " " << m_Matrix.rows() << std::endl; return 0;}
    if (cols!=m_Matrix.cols()){std::cout << "Incorrect col dimensions: " << cols << " " << m_Matrix.cols() << std::endl; return 0;}
    //now check charges are all ok
    for (Sparseint c=0;c<m_Matrix.cols();++c){
      State colstate=OutgoingColState(c);
      for (Sparseint p=m_Matrix.get_p(c);p<m_Matrix.get_p(c+1);++p){
	if (IngoingRowState(m_Matrix.get_i(p))!=colstate){
	  std::cout << "Row and col states don't match " << m_Matrix.get_i(p) << " " << c << " " << m_Matrix.get_x(p) << std::endl;
	  std::cout << IngoingRowState(m_Matrix.get_i(p)) << "!= " << colstate << std::endl;
	  print_matrix();
	  print_indices();
	  
	  getPhysicalSpectrum().print();
	  return 0;
	}
      }
    }
    return 1;
  }

  MPX_matrix MPX_matrix::ExtractSubMPX(const std::vector<MPXPair >& IndexVal) const {
    MPX_matrix ans(*(m_SpectrumPtr),m_Indices,m_NumRowIndices,m_Matrix.ExtractSubMatrix(m_NumRowIndices,dimsvector(),IndexVal,0)); 
    for (std::vector<MPXPair>::const_iterator cit=IndexVal.begin();cit!=IndexVal.end();++cit){
      ans.m_Indices.at(cit->first)=MPXIndex(m_Indices.at(cit->first),cit->second);
    }
    return ans;
  }

  MPX_matrix MPX_matrix::RestrictColumnIndex(){
    return MPX_matrix(*(m_SpectrumPtr),m_Indices,m_NumRowIndices,m_Matrix.ZeroLastColumns(m_SpectrumPtr->size()));
  }

  State MPX_matrix::CombinationState(std::vector<MPXInt> indices, MPXInt l) const {
    State accumulatorstate(m_SpectrumPtr->at(0).getChargeRules());
    for (std::vector<MPXInt>::const_iterator cit=indices.begin();cit!=indices.end();++cit){
      const MPXIndex& indexref= m_Indices.at(*cit);
      if (indexref.Outgoing()){
	  accumulatorstate+=indexref.at(l % indexref.size());
      }
      else { //if not Outgoing then must be incoming
	accumulatorstate-=indexref.at(l % indexref.size());
      }
      l/=indexref.size();
    }
    return accumulatorstate;
  }

  State MPX_matrix::OutgoingColState(MPXInt j) const{
    //calculate net Outgoing charges for a particular col
    State accumulatorstate(m_SpectrumPtr->at(0).getChargeRules());
    //take col index and break apart
    for (size_t cp=0;cp<m_Indices.size()-m_NumRowIndices;++cp){
      const MPXIndex& indexref= m_Indices.at(m_NumRowIndices+cp);
	if (indexref.Outgoing()){
	  accumulatorstate+=indexref.at(j % indexref.size());
	}
	else { //if not Outgoing then must be incoming
	  accumulatorstate-=indexref.at(j % indexref.size());
	}
	j/=indexref.size();
    }
    return accumulatorstate;
  }

  State MPX_matrix::IngoingRowState(MPXInt i) const {
    //calculate net Incoming charges for a particular row
    State accumulatorstate(m_SpectrumPtr->at(0).getChargeRules());
    //take row index and break apart
    for (Sparseint rp=0;rp<m_NumRowIndices;++rp){
      const MPXIndex& indexref= m_Indices.at(rp);
	if (indexref.Ingoing()){
	  accumulatorstate+=indexref.at(i % indexref.size());
	}
	else { //if not incoming must be outgoing
	  accumulatorstate-=indexref.at(i % indexref.size());
	}
	i/=indexref.size();
    }
    return accumulatorstate;
  }

  std::vector<BlockStateIndices> MPX_matrix::GetAllBlockColumns() const{
    std::vector<BlockStateIndices> ans;
    for (Sparseint col=0;col<m_Matrix.cols();++col){
      State teststate=OutgoingColState(col);
      //now test accumulator state
      bool flag=0;
      for (std::vector<BlockStateIndices>::iterator it=ans.begin();it!=ans.end();++it){
	//check to see if we already have this sector
	if (teststate==it->BlockState){
	  //if yes, then update flag and add state to sector
	  flag=1;
	  it->Indices.push_back(col);
	  break;
	}
      }
      if (flag==0){
	ans.push_back(BlockStateIndices(teststate,col)); //start new group
      }
    }
    return ans;
  }

  BlockStateIndices MPX_matrix::GetBlockColumns(const State& specificstate) const{
    BlockStateIndices ans(specificstate);
    for (Sparseint col=0;col<m_Matrix.cols();++col){
      State teststate=OutgoingColState(col);
      //now test accumulator state
      //check to see if we already have this sector
      if (teststate==ans.BlockState){
	//if yes, then add to sector
	ans.Indices.push_back(col);
      }
    }
    return ans;
  }

  BlockStateIndices MPX_matrix::GetBlockRows(const State& specificstate) const{
    BlockStateIndices ans(specificstate);
    for (Sparseint row=0;row<m_Matrix.rows();++row){
      State teststate=IngoingRowState(row);
      //now test accumulator state
      //check to see if we already have this sector
      if (teststate==ans.BlockState){
	//if yes, then add to sector
	ans.Indices.push_back(row);
      }
    }
    return ans;
  }

  MPX_matrix& MPX_matrix::RemoveDummyIndices(std::vector<MPXInt> indices_for_removal){
    //sort indices so they are ascending
    std::sort(indices_for_removal.begin(),indices_for_removal.end());
    //MPXInt rows_to_remove(0);
    //using reverse iterator, move through the vector
    for (std::vector<MPXInt>::const_reverse_iterator rit=indices_for_removal.rbegin();rit!=indices_for_removal.rend();++rit){
      if (*rit<this->m_NumRowIndices){--m_NumRowIndices;}
      this->m_Indices.erase(m_Indices.begin()+*rit,m_Indices.begin()+*rit+1); //would be better to rotate bad values to the end then erase
    }
    //this->m_NumRowIndices-=rows_to_remove;
    return *this;
  }

  MPX_matrix& MPX_matrix::Rescale(std::complex<double> factor){
    m_Matrix.rescale(factor);
    return *this;
  }

  MPX_matrix& MPX_matrix::CombineSimilarMatrixIndices(bool PhysicalInMiddle){
    //find all indices that are ingoing and outgoing
    std::vector<Sparseint> IngoingMatrixIndices;
    std::vector<Sparseint> OutgoingMatrixIndices;
    std::vector<Sparseint> PhysicalIndices;
    for (uMPXInt m=0;m<m_Indices.size();++m){
      if (!m_Indices[m].Physical()){
	m_Indices[m].Ingoing() ? IngoingMatrixIndices.push_back(m) : OutgoingMatrixIndices.push_back(m);
      }
      else {PhysicalIndices.push_back(m);}
    }
    //for (auto l : IngoingMatrixIndices) std::cout << "Ingoing " << l << std::endl;
    //for (auto l : OutgoingMatrixIndices) std::cout << "Outgoing " << l << std::endl;
    //lists in hand we need to reshape
    //for which we need one long list of the new order
    std::vector<Sparseint> new_order;
    new_order.reserve(m_Indices.size());
    if (!PhysicalInMiddle){
      new_order.insert(new_order.end(),PhysicalIndices.begin(),PhysicalIndices.end());
      new_order.insert(new_order.end(),IngoingMatrixIndices.begin(),IngoingMatrixIndices.end());
    }
    else {
      new_order.insert(new_order.end(),IngoingMatrixIndices.begin(),IngoingMatrixIndices.end());
      new_order.insert(new_order.end(),PhysicalIndices.begin(),PhysicalIndices.end());
    }
    new_order.insert(new_order.end(),OutgoingMatrixIndices.begin(),OutgoingMatrixIndices.end());

    Sparseint reshaperowindices = PhysicalInMiddle ? IngoingMatrixIndices.size() : IngoingMatrixIndices.size()+PhysicalIndices.size();
    if (!reshaperowindices) {std::cout << "Error: incorrect reshape and combine parameters, no row indices" << std::endl; exit(1);}

    //reshape matrix
    m_Matrix=std::move(reshape(m_Matrix,m_NumRowIndices,reshaperowindices,dimsvector(),new_order,0));

    //set new num row indices
    m_NumRowIndices=0;
    if (!PhysicalInMiddle){
      if (PhysicalIndices.size()) ++m_NumRowIndices;
    }
    if (IngoingMatrixIndices.size()) ++m_NumRowIndices;
    if (!m_NumRowIndices) {std::cout << "Error: incorrect reshape and combine parameters, no row indices" << std::endl; exit(1);}

    //now combine MPXIndex's 

    std::vector<MPXIndex> NewIndexVector;
    if (!PhysicalInMiddle){
      for (auto i : PhysicalIndices){ 
	NewIndexVector.emplace_back(m_Indices[i]);
      }
      //now do ingoing matrix indices
      if (IngoingMatrixIndices.size()>0){
	MPXIndex In(m_Indices[IngoingMatrixIndices[0]]);
	for (uMPXInt i=1; i<IngoingMatrixIndices.size();++i){
	  In=combine(In,m_Indices[IngoingMatrixIndices[i]]);   
	}
	NewIndexVector.emplace_back(In);
      }
    }

    else {
      if (IngoingMatrixIndices.size()>0){
	MPXIndex In(m_Indices[IngoingMatrixIndices[0]]);
	for (uMPXInt i=1; i<IngoingMatrixIndices.size();++i){
	  In=combine(In,m_Indices[IngoingMatrixIndices[i]]);   
	}
	NewIndexVector.emplace_back(In);
      }
      for (auto i : PhysicalIndices){ 
	NewIndexVector.emplace_back(m_Indices[i]);
      }
    }

    if (OutgoingMatrixIndices.size()>0){
      MPXIndex Out(m_Indices[OutgoingMatrixIndices[0]]);
      for (uMPXInt i=1; i<OutgoingMatrixIndices.size();++i){
	Out=combine(Out,m_Indices[OutgoingMatrixIndices[i]]);   
      }
      NewIndexVector.emplace_back(Out);
    }

    m_Indices=std::move(NewIndexVector);
    return *this;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPO members
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void MPO_matrix::check(){
    if (m_Indices.size()!=4/* || m_NumRowIndices!=2*/){std::cout << "Incorrect number of total indices for MPO_matrix" << std::endl; exit(1);}
    MPXInt numphysical=0;
    MPXInt nummatrix=0;
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.end();++cit){
      if (cit->Physical()){++numphysical;}
      else {++nummatrix;}
    }
    if (numphysical!=2){std::cout << "Incorrect number of physical indices for MPO_matrix: " << numphysical << std::endl; exit(1);}
    if (nummatrix!=2){std::cout << "Incorrect number of matrix indices for MPO_matrix: " << numphysical << std::endl; exit(1);}
  }

  MPO_matrix::MPO_matrix(const EigenStateArray& spectrum) : MPX_matrix(spectrum){};

  MPO_matrix::MPO_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, const SparseMatrix& matrix) : MPX_matrix(spectrum,indices,2,matrix)
  {
    check();
  }

  MPO_matrix::MPO_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, SparseMatrix&& matrix) : MPX_matrix(spectrum,indices,2,std::move(matrix))
  {
    check();
  }

  MPO_matrix::MPO_matrix(MPX_matrix&& MPXref) noexcept : MPX_matrix(std::move(MPXref))
  {
    check();
  }

  MPO_matrix::MPO_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values, bool inverse): MPX_matrix(spectrum) {
    m_Indices.reserve(4);
    m_Indices.emplace_back(1,spectrum);
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,spectrum);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=2;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPO_matrix::MPO_matrix(const EigenStateArray& spectrum, const MPXIndex& index, const std::vector<double>& values, bool inverse): MPX_matrix(spectrum) {
    m_Indices.reserve(4);
    m_Indices.emplace_back(1,spectrum);
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,spectrum);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=2;
    m_Matrix=SparseMatrix(values,inverse);
  }

  MPO_matrix MPO_matrix::ExtractMPOBlock(const std::pair<MPXInt,MPXInt>& matrix_index_row_range, const std::pair<MPXInt,MPXInt>& matrix_index_col_range) const {
    //get row range and col range
    std::pair<MPXInt,MPXInt> array_row_range(matrix_index_row_range.first*basis().size(),(matrix_index_row_range.second+1)*basis().size()-1);
    std::pair<MPXInt,MPXInt> array_col_range(matrix_index_col_range.first*basis().size(),(matrix_index_col_range.second+1)*basis().size()-1);

   StateArray row_matrix_index;
   for (MPXInt r=matrix_index_row_range.first;r<=matrix_index_row_range.second;++r){
      row_matrix_index.emplace_back(Index(1)[r]);
    }

   StateArray col_matrix_index;
   for (MPXInt c=matrix_index_col_range.first;c<=matrix_index_col_range.second;++c){
      col_matrix_index.emplace_back(Index(3)[c]);
    }
    
    return MPO_matrix(basis(),std::vector<MPXIndex>({{1,basis()},{1,std::move(row_matrix_index)},{0,basis()},{0,std::move(col_matrix_index)}}),m_Matrix.ExtractSubMatrix(array_row_range,array_col_range));

  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPS members
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void MPS_matrix::check(){
    //carry out check if constructing from already existing parts
    if (m_Indices.size()!=3){std::cout << "Incorrect number of indices for MPS_matrix: " << m_Indices.size() << std::endl; exit(1);}
    MPXInt numphysical=0;
    MPXInt nummatrix=0;
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.end();++cit){
      if (cit->Physical()){++numphysical;}
      else {++nummatrix;}
    }
    if (numphysical!=1){std::cout << "Incorrect number of physical indices for MPS_matrix: " << numphysical << std::endl; exit(1);}
    if (nummatrix!=2){std::cout << "Incorrect number of matrix indices for MPS_matrix: " << numphysical << std::endl; exit(1);}
  };

  MPS_matrix::MPS_matrix(const EigenStateArray& spectrum) : MPX_matrix(spectrum){};

  //assume two row indices so left type storage
  MPS_matrix::MPS_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, const SparseMatrix& matrix) : MPX_matrix(spectrum,indices,2,matrix)
  {
    check();
  }
  MPS_matrix::MPS_matrix(const EigenStateArray& spectrum,const std::vector<MPXIndex>& indices, SparseMatrix&& matrix) : MPX_matrix(spectrum,indices,2,std::move(matrix))
  {
    check();
  }
  /*MPS_matrix::MPS_matrix(const MPX_matrix& MPXref) : MPX_matrix(MPXref)
  {
    check();
    }*/

  MPS_matrix::MPS_matrix(MPX_matrix&& MPXref) noexcept : MPX_matrix(std::move(MPXref))
  {
    check();
  }

  /*MPS_matrix::MPS_matrix(MPS_matrix&& rhs) noexcept : MPS_matrix()
  {
    MPS_swap(*this,rhs);
    }*/

  bool MPSDecomposition::store(const std::string& LeftName, const std::string& RightName) const {
    std::ofstream LeftOutfile;
    LeftOutfile.open(LeftName.c_str(),ios::out | ios::trunc | ios::binary);    
    if (LeftMatrix.fprint_binary(LeftOutfile)) return 1;
    else {
      std::ofstream RightOutfile;
      RightOutfile.open(RightName.c_str(),ios::out | ios::trunc | ios::binary);
      if (RightMatrix.fprint_binary(RightOutfile)) return 1;
      else return 0;
    }
  }

  bool MPSDecomposition::store(const std::string& Name, uMPXInt nl,uMPXInt nr) const {
    std::ofstream LeftOutfile;
    LeftOutfile.open(LeftName(Name,nl).c_str(),ios::out | ios::trunc | ios::binary);    
    if (LeftMatrix.fprint_binary(LeftOutfile)) return 1;
    else {
      std::ofstream RightOutfile;
      RightOutfile.open(RightName(Name,nr).c_str(),ios::out | ios::trunc | ios::binary);
      if (RightMatrix.fprint_binary(RightOutfile)) return 1;
      else {
	std::stringstream LambdaName;
	LambdaName << Name.c_str() << "_Lambda_" << nl << "_" << nr << ".MPX_matrix";
	std::ofstream LambdaOutfile;
	LambdaOutfile.open(LambdaName.str().c_str(),ios::out | ios::trunc | ios::binary);
	if (MPX_matrix(LeftMatrix.GetPhysicalSpectrum(),LeftMatrix.Index(2),Values).fprint_binary(LambdaOutfile)){return 1;}
	else return 0;
      }
    }
  }

  bool MPSDecomposition::store_left(const std::string& Name, uMPXInt nl) const{
    return LeftMatrix.store(LeftName(Name,nl));
  }

  bool MPSDecomposition::store_right(const std::string& Name, uMPXInt nr) const{
      return RightMatrix.store(RightName(Name,nr));
  }

  void MPSDecomposition::OutputPhysicalIndexDensities(std::ofstream& d) const {
    if (d.is_open()){
      std::vector<std::complex<double> > densities;
      {
	MPX_matrix temp(contract(LeftMatrix,0,MPX_matrix(LeftMatrix.GetPhysicalSpectrum(),LeftMatrix.Index(2),Values),0,contract20));
	contract_to_sparse(temp,1,temp,0,contract1122).diagonal();
      }      
      for (auto &&i : densities){
	d << real(i) << " ";
      }
      d << std::endl;
    }
  }

  void UnitCell::OutputPhysicalIndexDensities(std::ofstream& d) const {
    if (d.is_open() && size()){
      std::vector<std::complex<double> > densities;
      MPX_matrix accumulator(contract(Matrices.at(0),1,Matrices.at(0),0,contract0011));
      for (uMPXInt i=1;i<Matrices.size()-1;++i){
	accumulator=std::move(contract(contract(accumulator,0,Matrices.at(i),1,contract01),0,Matrices.at(i),0,contract0110));
      }
      {
	const MPS_matrix& A=Matrices.back();
	const std::vector<double>& L=Lambdas.at(0);
	MPX_matrix temp(contract(A,0,MPX_matrix(A.GetPhysicalSpectrum(),A.Index(2),L),0,contract20));
	densities=contract_to_sparse(contract(temp,1,accumulator,0,contract10),0,temp,0,contract1221).diagonal();
      }
      for (auto &&i : densities){
	d << real(i) << " ";
      }
      d << std::endl;
    }
  }
  
  void UnitCell::OutputOneVertexDensityMatrix(std::ofstream& d) const {
    if (d.is_open() && size()){
      MPX_matrix accumulator(contract(Matrices.at(0),1,Matrices.at(0),0,contract0011));
      for (uMPXInt i=1;i<Matrices.size()-1;++i){
	accumulator=std::move(contract(contract(accumulator,0,Matrices.at(i),1,contract01),0,Matrices.at(i),0,contract0110));
      }
      {
	const MPS_matrix& A=Matrices.back();
	const std::vector<double>& L=Lambdas.at(0);
	MPX_matrix temp(contract(A,0,MPX_matrix(A.GetPhysicalSpectrum(),A.Index(2),L),0,contract20));
	SparseMatrix S(contract_to_sparse(contract(temp,1,accumulator,0,contract10),0,temp,0,contract1221));
	S.fprint(d);
      }
    }
  }

  void UnitCell::OutputOneVertexDensityMatrix(const std::string& name, uMPXInt l) const {
    std::stringstream DensityMatrixNameStream;
    DensityMatrixNameStream << name << "_" << l << ".dat";
    std::ofstream DensityMatrixFileStream;
    DensityMatrixFileStream.open(DensityMatrixNameStream.str().c_str(),ios::out | ios::trunc);
    OutputOneVertexDensityMatrix(DensityMatrixFileStream);
    DensityMatrixFileStream.close();
  }

  void UnitCell::store(const std::string& filename, uMPXInt l) const {
    std::stringstream FileNameStream;
    FileNameStream << filename << "_" << l;
    store(FileNameStream.str());
  }

  void UnitCell::store(const std::string& filename) const {
    //first do a consistency check
    if (Matrices.size()!=Lambdas.size()){
      std::cout << "Invalid Unit Cell, skipping store." << std::endl; 
    }
    else {
      std::stringstream FileNameStream;
      FileNameStream << filename << ".UNITCELL"; 
      std::ofstream outfile;
      outfile.open(FileNameStream.str().c_str(),ios::out | ios::trunc | ios::binary);
      /*//first output basis
      size_t basis_size(basis_ptr_->size());
      for (auto&& s : (*basis_ptr_)){
	s.fprint_binary(outfile);
	}*/      
      size_t num_in_unit_cell(Matrices.size());
      outfile.write(reinterpret_cast<const char*>(&num_in_unit_cell),sizeof(size_t));
      for (auto&& m : Matrices){
	m.fprint_binary(outfile);
      }
      for (auto&& l : Lambdas){
	size_t lambda_size(l.size());
	outfile.write(reinterpret_cast<const char*>(&lambda_size),sizeof(size_t));	
	for (auto v : l){
	  outfile.write(reinterpret_cast<const char*>(&v),sizeof(double));
	}
      }
    }
  }


  double UnitCell::Entropy() const {
    if (Lambdas.size()){
      return entropy(Lambdas[0]);
    }
    else {
      std::cout << "UnitCell illformed. Can't compute entanglement entropy, returning -1" <<std::endl;
      return -1.0; //indicates failure
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  StateArray DeduceFromColumns(const SparseMatrix& array, const MPXIndex& colindices){
    StateArray rowstates;
    SparseMatrix trans(array.copy_transpose());
    for (Sparseint col=0;col<trans.cols();++col){
      if (trans.get_p(col)==trans.get_p(col+1)){ //empty
	rowstates.push_back(colindices.at(0)-colindices.at(0));
      }
      else {
	double testvalue=0.0;
	Sparseint i=-1;//dummy value
	for (Sparseint p=trans.get_p(col);p<trans.get_p(col+1);++p){
	  //first pass find largest abs
	  if (abs(trans.get_x(p))>testvalue){
	    testvalue=abs(trans.get_x(p));
	    i=trans.get_i(p);
	  }
	}
	rowstates.push_back(colindices.at(i));
	for (Sparseint p=trans.get_p(col);p<trans.get_p(col+1);++p){
	  if (rowstates[col]!=colindices.at(trans.get_i(p))){
	    //shouldn't be zero
	    /*if (abs(trans.get_x(p))<testvalue*SPARSETOL){
	      trans.set_x(p,0.0);
	    }
	    else {*/
	      std::cout << "Large state mismatch value" << trans.get_x(p) << std::endl; exit(1);
	      //}
	  }
	}
      }
    }
    trans.dropzeros();
    return rowstates;
    }

  MPXIndex load_MPXIndex_binary(std::ifstream& infile,const EigenStateArray& spectrum) {
    if (infile.is_open()){
      bool Inward;
      bool Physical;
      size_t indexsize;
      infile.read(reinterpret_cast<char*>(&Inward),sizeof(bool));
      infile.read(reinterpret_cast<char*>(&Physical),sizeof(bool));
      infile.read(reinterpret_cast<char*>(&indexsize),sizeof(size_t));
      if (Physical){ return MPXIndex(Inward,spectrum);}
      else {
	//we get charge rules from the physical spectrum
	const QNVector& cr(spectrum[0].getChargeRules());
	QuantumNumberInt ncr(cr.size());
	QuantumNumberInt* buffer=new QuantumNumberInt[ncr];
	StateArray loadedstates;
	loadedstates.reserve(indexsize);
	for (size_t s=0;s<indexsize;++s){
	  infile.read(reinterpret_cast<char*>(buffer),sizeof(QuantumNumberInt)*ncr);
	  loadedstates.push_back(State(cr,QNVector(buffer,buffer+ncr)));
	}
	delete[] buffer;
	return MPXIndex(Inward,loadedstates);
      }
    }
    else {
      return MPXIndex(0,spectrum);
    }
  }

  MPX_matrix load_MPX_matrix_binary(const std::string& filename,const EigenStateArray& spectrum){
    std::cout << "LOADING " << filename << std::endl; 
    std::ifstream infile;
    infile.open(filename.c_str(),ios::in | ios::binary);
    return load_MPX_matrix_binary(infile,spectrum);
  }

  MPX_matrix load_MPX_matrix_binary(std::ifstream& infile,const EigenStateArray& spectrum){
    if (infile.is_open()){
      size_t numindices(0);
      Sparseint numrowindices(0);
      infile.read(reinterpret_cast<char*>(&numindices),sizeof(size_t));
      infile.read(reinterpret_cast<char*>(&numrowindices),sizeof(Sparseint));
      std::vector<MPXIndex> indices;
      indices.reserve(numindices);
      for (size_t i=0;i<numindices;++i){
	indices.push_back(load_MPXIndex_binary(infile,spectrum));
      }
      MPX_matrix ans(spectrum,indices,numrowindices,load_SparseMatrix_binary(infile));
      //ans.print_indices();
      return ans;
    }
    else return MPX_matrix(spectrum);
  }

  MPX_matrix load_MPX_matrix(const std::string& filename,const EigenStateArray& spectrum){
    return load_MPX_matrix_binary(filename,spectrum);
  }
  

  MPO_matrix load_MPO_matrix(const std::string& filename,const EigenStateArray& spectrum){
    return MPO_matrix(std::move(load_MPX_matrix_binary(filename,spectrum)));
  }

  MPS_matrix load_MPS_matrix(const std::string& filename,const EigenStateArray& spectrum){
    return MPS_matrix(std::move(load_MPX_matrix_binary(filename,spectrum)));
  }

  MPS_matrix MakeProductState(const EigenStateArray& spectrum, uMPXInt state_index){
    std::vector<MPXIndex> indices;
    indices.emplace_back(1,spectrum);
    indices.emplace_back(1,StateArray(1,State(spectrum.at(0).getChargeRules())));
    indices.emplace_back(0,StateArray(1,State(spectrum.at(state_index))));
    SparseMatrix array(spectrum.size(),1,1);
    array.entry(state_index,0,1.0);
    array.cheap_finalise();
    return MPS_matrix(spectrum,indices,array);
  }

  UnitCell load_UnitCell_binary(std::ifstream& infile, const Basis& basis){
    UnitCell ans(basis);
    if (infile.is_open()){
      //read in size of unitcell
      size_t num_in_cell(0);
      infile.read(reinterpret_cast<char*>(&num_in_cell),sizeof(size_t));
      for (size_t n=0;n<num_in_cell;++n){
	ans.Matrices.emplace_back(MPS_matrix(load_MPX_matrix_binary(infile,basis)));
      }
      for (size_t n=0;n<num_in_cell;++n){
	std::vector<double> lambda;
	size_t num_in_lambda(0);
	infile.read(reinterpret_cast<char*>(&num_in_lambda),sizeof(size_t));
	for (size_t l=0;l<num_in_lambda;++l){
	  double val(0.0);
	  infile.read(reinterpret_cast<char*>(&val),sizeof(double));
	  lambda.push_back(val);
	}
	ans.Lambdas.emplace_back(std::move(lambda));
      }
      return ans;
    }
    else {
      return ans;
    }
  }
}

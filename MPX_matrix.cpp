#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <numeric>

#include "sparse_interface.hpp"
//#include "arpack_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp"
#include "MPX_matrix.hpp"

namespace ajaj{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPX members
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPX_matrix::MPX_matrix(const Basis& spectrum) noexcept : m_SpectrumPtr(&spectrum){};
  MPX_matrix::MPX_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, const Sparseint numrowindices, const SparseMatrix& matrix) : m_SpectrumPtr(&spectrum),m_Indices(indices),m_NumRowIndices(numrowindices),m_Matrix(copy(matrix)){
    if (!m_Matrix.is_finalised()){std::cout << "Initialising MPX_matrix from unfinalised SparseMatrix is not allowed!" << std::endl; exit(1);}
  };

  MPX_matrix::MPX_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, const Sparseint numrowindices, SparseMatrix&& matrix) : m_SpectrumPtr(&spectrum),m_Indices(indices),m_NumRowIndices(numrowindices),m_Matrix(std::move(matrix)){
    if (!m_Matrix.is_finalised()){std::cout << "Initialising MPX_matrix from unfinalised SparseMatrix is not allowed!" << std::endl; exit(1);}
  };

  MPX_matrix::MPX_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<std::complex<double> >& values, bool inverse): m_SpectrumPtr(&spectrum) {
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=1;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPX_matrix::MPX_matrix(const Basis& spectrum, const MPXIndex& Lindex, const MPXIndex& Rindex, const std::vector<std::complex<double> >& values, bool inverse): m_SpectrumPtr(&spectrum) {
    if (Lindex.size()!=Rindex.size()){
      std::cout << "Malformed diagonal MPX_matrix, left and right indices are of different lengths!" <<std::endl; exit(1);
    }
    m_Indices.emplace_back(1,Lindex);
    m_Indices.emplace_back(0,Rindex);
    m_NumRowIndices=1;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPX_matrix::MPX_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse): m_SpectrumPtr(&spectrum) {
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=1;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPX_matrix::MPX_matrix(const Basis& spectrum, const MPXIndex& Lindex, const MPXIndex& Rindex, const std::vector<double>& values, bool inverse): m_SpectrumPtr(&spectrum) {
    if (Lindex.size()!=Rindex.size()){
      std::cout << "Malformed diagonal MPX_matrix, left and right indices are of different lengths!" <<std::endl; exit(1);
    }
    m_Indices.emplace_back(1,Lindex);
    m_Indices.emplace_back(0,Rindex);
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

  bool MPX_matrix::isHermitian() const {
    if (m_Matrix.rows()!= m_Matrix.cols()) return 0;
    else return check_equal(m_Matrix,m_Matrix.copy_dagger(),1.0e-14);
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
	  
	  basis().print();
	  return 0;
	}
      }
    }
    return 1;
  }

  /*MPX_matrix MPX_matrix::ExtractSubMPX(const std::vector<MPXPair >& IndexVal) const {
    MPX_matrix ans(*(m_SpectrumPtr),m_Indices,m_NumRowIndices,m_Matrix.ExtractSubMatrix(m_NumRowIndices,dimsvector(),IndexVal,0)); 
    for (std::vector<MPXPair>::const_iterator cit=IndexVal.begin();cit!=IndexVal.end();++cit){
      ans.m_Indices.at(cit->first)=MPXIndex(m_Indices.at(cit->first),cit->second);
    }
    return ans;
    }*/

  MPX_matrix MPX_matrix::ZeroLastBlock(){
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
    if (indices_for_removal.size()>=m_Indices.size()-1){
      std::cout << "Too many indices for removal!" << std::endl;
	*this=MPX_matrix();
	return *this;
    }

    //check dummies
    for (auto i : indices_for_removal){
      if (m_Indices[i].size()!=1 || i >=m_Indices.size()) {
	std::cout << "Not a dummy index!" << std::endl;
	*this=MPX_matrix();
	return *this;
      }
    }

    MPXInt num_removed_row_indices(0); //how many of the removed indices are row indices?
    for (auto i : indices_for_removal){
      if (i<m_NumRowIndices) ++num_removed_row_indices;
    }
    if (num_removed_row_indices==m_NumRowIndices){
      std::cout <<"All row indices would be removed, requiring a reshape!"<<std::endl;
      exit(1);
    }

    std::vector<MPXInt> newindices(m_Indices.size());
    std::iota(newindices.begin(),newindices.end(),0); //0,1,2...,m_Indices.size()-1
    auto nend(newindices.end());
    for (auto i : indices_for_removal){ //go through and move the ones we want to get rid of
      nend=std::remove(newindices.begin(),nend,i);
    }
    newindices.erase(nend,newindices.end()); //erase the ones we moved

    std::vector<MPXIndex> new_m_Indices;
    for (auto n : newindices){
      new_m_Indices.emplace_back(m_Indices[n]);
    }

    swap(new_m_Indices,m_Indices);
    m_NumRowIndices-=num_removed_row_indices;

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


  ///
  //Friends
  ///

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

  MPX_matrix load_MPX_matrix_binary(const std::string& filename,const Basis& spectrum){
    std::cout << "LOADING " << filename << std::endl; 
    std::ifstream infile;
    infile.open(filename.c_str(),ios::in | ios::binary);
    return load_MPX_matrix_binary(infile,spectrum);
  }

  MPX_matrix load_MPX_matrix_binary(std::ifstream& infile,const Basis& spectrum){
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
      return MPX_matrix(spectrum,indices,numrowindices,load_SparseMatrix_binary(infile));
    }
    return MPX_matrix(); //return null if it didn't work
  }

  MPX_matrix load_MPX_matrix(const std::string& filename,const Basis& spectrum){
    return load_MPX_matrix_binary(filename,spectrum);
  }
  
  //Ancilliary helpers

  void swap(MPX_matrix& A, MPX_matrix& B){
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
    std::cout << "##################" << std::endl;
    std::cout << "Contract to sparse" << std::endl;
    std::cout << "Left info" << std::endl;
    A.print_sparse_info();
    std::cout << "Right info" << std::endl;
    B.print_sparse_info();
    std::cout << "##################" << std::endl;
#endif
    if (A.m_SpectrumPtr!=B.m_SpectrumPtr){std::cout <<"Physical Indices don't match!" << std::endl;exit(1);}
    //next check contraction indices match sizes match!
    for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
      if (!match(A.m_Indices.at(cit->first),conjA,B.m_Indices.at(cit->second),conjB)){exit(1);}
    }

    //what if the arrays are already suitably shaped?
    bool readyA(1);
    bool readyB(1);

    if (A.m_Indices.size()-A.m_NumRowIndices!=contractidxs.size()) readyA=0;
    if (B.m_NumRowIndices!=contractidxs.size()) readyB=0;
    
    if (readyA){
      MPXInt checkA=-1;
      for (auto&& ci : contractidxs){
	if (ci.first<A.m_NumRowIndices) {//not a column index
	  readyA=0;
	  break;
	}
	if (ci.first<checkA){ //not ordered!
	  readyA=0;
	  break;
	}
	else if (ci.first==checkA) {std::cout << "Illegal contraction, repeated indices!" <<std::endl; exit(1);}
	checkA=ci.first;
      }
    }

    if (readyB){
      MPXInt checkB=-1;
      for (auto&& ci : contractidxs){
	if (ci.second>=B.m_NumRowIndices) {//not a column index
	  readyB=0;
	  break;
	}
	if (ci.second<checkB){ //not ordered!
	  readyB=0;
	  break;
	}
	else if (ci.second==checkB) {std::cout << "Illegal contraction, repeated indices!" <<std::endl; exit(1);}
	checkB=ci.second;
      }
    }

    //std::cout << "A" << readyA << "\t B" << readyB <<std::endl; 
    if (readyA && readyB){
      return conj_multiply(A.m_Matrix,conjA,B.m_Matrix,conjB);
    }
    else if (readyA){
      std::vector<Sparseint> neworderB(B.m_Indices.size());
      std::iota(neworderB.begin(),neworderB.end(),0);
      for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
	*(std::remove(neworderB.begin(),neworderB.end(),cit->second))=cit->second;
      }
      std::rotate(neworderB.begin(),neworderB.end()-contractidxs.size(),neworderB.end());
      return conj_multiply(A.m_Matrix,conjA,reshape(B.m_Matrix,B.m_NumRowIndices,contractidxs.size(),B.dimsvector(),neworderB,conjB),0);
    }
    else if (readyB){
      std::vector<Sparseint> neworderA(A.m_Indices.size());
      std::iota(neworderA.begin(),neworderA.end(),0);
      for (std::vector<MPXPair >::const_iterator cit=contractidxs.begin();cit!=contractidxs.end();++cit){
	*(std::remove(neworderA.begin(),neworderA.end(),cit->first))=cit->first;
      }
      return conj_multiply(reshape(A.m_Matrix,A.m_NumRowIndices,A.m_Indices.size()-contractidxs.size(),A.dimsvector(),neworderA,conjA),0,B.m_Matrix,conjB);
    }
    else { //neither are currently shaped correctly, and we can avoid a transpose by including it in the reshape...
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
#if defined(USETBB) //no avoidance of extra transpose...
    std::rotate(neworderB.begin(),neworderB.end()-contractidxs.size(),neworderB.end());
    return std::move(reshape(A.m_Matrix,A.m_NumRowIndices,neworderA.size()-contractidxs.size(),A.dimsvector(),neworderA,conjA)*reshape(B.m_Matrix,B.m_NumRowIndices,contractidxs.size(),B.dimsvector(),neworderB,conjB));
#else
    //avoid extra transpose through reshape
    std::rotate(neworderA.begin(),neworderA.end()-contractidxs.size(),neworderA.end());
    return std::move(NoTransMultiply(reshape(B.m_Matrix,B.m_NumRowIndices,neworderB.size()-contractidxs.size(),B.dimsvector(),neworderB,conjB),reshape(A.m_Matrix,A.m_NumRowIndices,contractidxs.size(),A.dimsvector(),neworderA,conjA)).transpose());
#endif
    }
  }

  MPX_matrix contract(const MPX_matrix& A, bool conjA,const MPX_matrix& B, bool conjB, const std::vector<MPXPair>& contractidxs){
#ifndef NDEBUG
    std::cout << "##################" << std::endl;
    std::cout << "Contract" << std::endl;
    std::cout << "Left info" << std::endl;
    A.print_sparse_info();
    std::cout << "Right info" << std::endl;
    B.print_sparse_info();
    std::cout << "##################" << std::endl;
#endif
    //figure out new indices
    std::vector<MPXIndex> Indices;

    for (size_t a=0;a<A.m_Indices.size();++a){
      bool drop(0);
      for (auto&& ci : contractidxs){
	if (ci.first==a) {drop=1; break;}
      }
      if (!drop){ //if we are keeping this index
	Indices.emplace_back(conjA ? A.m_Indices.at(a).flip() : A.m_Indices.at(a));
      }
    }

    MPXInt new_num_row_indices(Indices.size());

    for (size_t b=0;b<B.m_Indices.size();++b){
      bool drop(0);
      for (auto&& ci : contractidxs){
	if (ci.second==b) {drop=1; break;}
      }
      if (!drop){
	Indices.emplace_back(conjB ? B.m_Indices.at(b).flip() : B.m_Indices.at(b));
      }
    }

    if (Indices.size()!=A.m_Indices.size()+B.m_Indices.size()-2*contractidxs.size()){
      for (auto&& i : contractidxs){
	std::cout << "(" << i.first << "," << i.second << ") "; 
      }
      std::cout << std::endl << new_num_row_indices << " " << Indices.size() <<std::endl;

      A.print_indices();
      B.print_indices();
      std::cout <<"ERROR assembling post contraction indices..." <<std::endl; exit(1);
    }

    return MPX_matrix(A.basis(),Indices,new_num_row_indices,contract_to_sparse(A,conjA,B,conjB,contractidxs));
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

    std::vector<double> roots;

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

}

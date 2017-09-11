#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp" 

namespace ajaj{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //friends
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void MPS_swap(MPS_matrix& A, MPS_matrix& B){
    std::swap(A.m_SpectrumPtr,B.m_SpectrumPtr);
    swap(A.m_Indices,B.m_Indices);
    std::swap(A.m_NumRowIndices,B.m_NumRowIndices);
    swap(A.m_Matrix,B.m_Matrix);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPS members
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

  MPS_matrix::MPS_matrix(const Basis& spectrum) : MPX_matrix(spectrum){};

  //assume two row indices so left type storage
  MPS_matrix::MPS_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, const SparseMatrix& matrix) : MPX_matrix(spectrum,indices,2,matrix)
  {
    check();
  }
  MPS_matrix::MPS_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, SparseMatrix&& matrix) : MPX_matrix(spectrum,indices,2,std::move(matrix))
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

  const MPXIndex& MPS_matrix::getInwardMatrixIndex() const {
    const MPXIndex* Indexptr(nullptr);
    for (auto&& i : m_Indices){
      if (i.Ingoing() && !i.Physical()) {Indexptr=&i; break;}
    }
    return *Indexptr;
  }
  const MPXIndex& MPS_matrix::getOutwardMatrixIndex() const {
    const MPXIndex* Indexptr(nullptr);
    for (auto&& i : m_Indices){
      if (i.Outgoing() && !i.Physical()) {Indexptr=&i; break;}
    }
    return *Indexptr;
  }
  const MPXIndex& MPS_matrix::getPhysicalIndex() const {
    const MPXIndex* Indexptr(nullptr);
    for (auto&& i : m_Indices){
      if (i.Physical()) {Indexptr=&i; break;}
    }
    return *Indexptr;
  }

  MPXInt MPS_matrix::InwardMatrixIndexNumber() const {
    for (uMPXInt i=0; i< m_Indices.size() ;++i){
      if (m_Indices[i].Ingoing() && !m_Indices[i].Physical()) {return i;}
    }
    return -1;//should never happen!
  }

  MPXInt MPS_matrix::OutwardMatrixIndexNumber() const {
    for (uMPXInt i=0; i< m_Indices.size() ;++i){
      if (m_Indices[i].Outgoing() && !m_Indices[i].Physical()) {return i;}
    }
    return -1;//should never happen!
  }
  MPXInt MPS_matrix::PhysicalIndexNumber() const {
    for (uMPXInt i=0; i< m_Indices.size() ;++i){
      if (m_Indices[i].Physical()) {return i;}
    }
    return -1;//should never happen!
  }

  MPS_matrix&& MPS_matrix::left_shape() {
    if (is_left_shape()){
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

  MPS_matrix&& MPS_matrix::right_shape() {
    if (is_right_shape()){
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

  bool MPSDecomposition::store(const std::string& LeftName, const std::string& RightName) const {
    std::ofstream LeftOutfile;
    LeftOutfile.open(LeftName.c_str(),ios::out | ios::trunc | ios::binary);    
    if (LeftMatrix.fprint_binary(LeftOutfile)) return 1;
    else {
      std::ofstream RightOutfile;
      RightOutfile.open(RightName.c_str(),ios::out | ios::trunc | ios::binary);
      if (RightMatrix.fprint_binary(RightOutfile)) return 1;
      else {
	return 0;
      }
    }
  }

  bool MPSDecomposition::store(const std::string& Name, uMPXInt nl,uMPXInt nr, bool StoreLambda) const {
    std::ofstream LeftOutfile;
    LeftOutfile.open(LeftName(Name,nl).c_str(),ios::out | ios::trunc | ios::binary);    
    if (LeftMatrix.fprint_binary(LeftOutfile)) return 1;
    else {
      std::ofstream RightOutfile;
      RightOutfile.open(RightName(Name,nr).c_str(),ios::out | ios::trunc | ios::binary);
      if (RightMatrix.fprint_binary(RightOutfile)) return 1;
      else {
	if (StoreLambda){
	  std::stringstream LambdaName;
	  LambdaName << Name.c_str() << "_Lambda_" << nl << "_" << nr << ".MPX_matrix";
	  std::ofstream LambdaOutfile;
	  LambdaOutfile.open(LambdaName.str().c_str(),ios::out | ios::trunc | ios::binary);
	  if (MPX_matrix(LeftMatrix.GetPhysicalSpectrum(),LeftMatrix.Index(2),Values).fprint_binary(LambdaOutfile)){return 1;}
	  else {
	    //now do one vertex density matrix
	    std::ofstream rhofile;
	    std::stringstream rhoName;
	    rhoName << "fDMRGRho_"<<nl <<".dat";
	    rhofile.open(rhoName.str().c_str(),ios::out | ios::trunc);
	    if (rhofile.is_open()){
	      OutputOneVertexDensityMatrix(rhofile);
	      return 0;
	    }
	    else return 1;
	  }
	}
	return 0;
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

  void MPSDecomposition::OutputOneVertexDensityMatrix(std::ofstream& d) const {
    if (d.is_open()){
      MPX_matrix L(LeftMatrix.basis(),LeftMatrix.Index(2),Values);
      MPX_matrix temp(contract(LeftMatrix,0,L,0,contract20));
      contract_to_sparse(temp,1,temp,0,contract1122).fprint(d);
    }
  }

  MPS_matrix load_MPS_matrix(const std::string& filename,const Basis& spectrum){
    return MPS_matrix(std::move(load_MPX_matrix_binary(filename,spectrum)));
  }

  MPS_matrix MakeProductState(const Basis& spectrum, const c_specifier_vector& state_index_vec,State leftstate){
    if (state_index_vec.size()){
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,spectrum);
      State matrix_index(spectrum.at(state_index_vec.back().first));
      indices.emplace_back(1,StateArray(1,leftstate));
      indices.emplace_back(0,StateArray(1,matrix_index+leftstate));
      SparseMatrix array(spectrum.size(),1,state_index_vec.size());
      for (auto&& s :state_index_vec){
	if (spectrum[s.first]==matrix_index)
	  array.entry(s.first,0,s.second);
	else {
	  std::cout << "Charges not consistent for initial state!" << std::endl;
	  std::cout << matrix_index << " and " << s.first << "," << spectrum[s.first] <<std::endl;
	  return MPS_matrix(spectrum);
	}
      }
      array.finalise();
      return MPS_matrix(spectrum,indices,array);
    }
    else return MPS_matrix(spectrum);
  }

  MPS_matrix MakeProductState(const Basis& spectrum, uMPXInt state_index, State leftstate){
    return MakeProductState(spectrum,c_specifier_vector(1,c_specifier(state_index,1.0)),leftstate);
  }

  c_specifier_array LoadCNumbers(const std::string& filename){
    std::ifstream infile;
    infile.open(filename.c_str(),ios::in);
    if (infile.is_open()){
      c_specifier_array c;
      std::string s;
      while (getline(infile,s)){
	if (s.empty()) continue;
	std::stringstream ss(s);
	if (ss.peek()=='#') continue;
	c.push_back(c_specifier_vector());
	uMPXInt idx;
	std::complex<double> value;
	while (ss >> idx >> value){ //take in pairs
	  c.back().push_back(c_specifier(idx,value));
	}
      }
      return c;
    }
    return c_specifier_array();
  }

}

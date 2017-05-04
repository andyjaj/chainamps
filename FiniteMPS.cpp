#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp"
#include "FiniteMPS.hpp"

namespace ajaj{

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt mix_idx): Basis_(model_basis),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(canon),MixPoint_(mix_idx){}

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs) : Basis_(model_basis),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(0),MixPoint_(num+1) {
    //if coeffs is empty, then nothing is changed
    if (coeffs.size()){
      for (auto&& c_vect : coeffs){
	if (c_vect.size()==0){
	  std::cout << "Error: Empty c number specification for vertex!" <<std::endl; exit(1);
	}
	State test_state(Basis_[c_vect.begin()->first]);
	for (auto&& c_pair : c_vect){
	  if (Basis_[c_pair.first]!=test_state){ //check for charge consistency (i.e. preservation of good quantum numbers)
	    std::cout << "Error: Inconsistent c numbers for vertex state conservation laws: " << c_vect.begin()->first << " " << c_pair.first  << std::endl; exit(1);
	  }
	  else {std::cout << "(" << c_pair.first << " : " << c_pair.second << ") ";}
	}
	std::cout << std::endl;
      }
      
      State LeftState(Basis_.getChargeRules()); //dummy state (all zeros)
      size_t counter=0;
      for (uMPXInt n=1;n<=NumVertices_;++n){
	auto& c_vec=coeffs[counter++];
	Current_=std::pair<uMPXInt,MPS_matrix>(n,MakeProductState(Basis_,c_vec,LeftState));
	LeftState+=Basis_[c_vec.begin()->first];
	if(counter==coeffs.size()){
	  counter=0;//back to beginning of c number defns
	}
	Current_.second.store(filename(position(),1/*Left*/));
      }
      //can't guarantee that we are canonical
    }
  }
			      
  void FiniteMPS::fetch_matrix(uMPXInt i, bool Left){
    //check buffer
    if (i <= NumVertices_ && i >= 1){ //inside range
      if (i!=position())
	Current_=std::pair<uMPXInt,MPS_matrix>(i,load_MPS_matrix(filename(i,Left),Basis_));
      else if (Left)
	Current_.second.left_shape();//does nothing if left shape already
      else 
	Current_.second.right_shape();

      if (matrix().empty()) {//failed to fetch
	std::cout << "Failed to fetch " << filename(i,Left) <<std::endl;
	exit(1);
      }
    }
	
    else {
      std::cout <<"Outside range of finite MPS!"<<std::endl;
      exit(1);
    }
  }

  std::string FiniteMPS::filename(uMPXInt i, bool Left, const std::string& name) const{
    std::stringstream namestream;
    if (Left)
      namestream << (name.empty() ? MPSName_ : name) << "_Left_" << i << ".MPS_matrix";
    else
      namestream << (name.empty() ? MPSName_ : name) << "_Right_" << NumVertices_-i+1 << ".MPS_matrix";
    
    return namestream.str();

  }

  void FiniteMPS::store_current(){

    if (!matrix().is_left_shape() && !matrix().is_right_shape()){
      std::cout << "Misshapen MPS_matrix, left shaping and storing..." <<std::endl;
      Current_.second.left_shape();
    }
    
    Current_.second.store(filename(position(),matrix().is_left_shape()));

  }																		 

  std::complex<double> FiniteMPS::makeLC(const std::string& name){
    std::cout << "Left canonising initial state..." <<std::endl;
    //do a check
    if (CheckFilesExist()==CanonicalType::Error){
      return 0.0;
    }

    MPX_matrix Vd;
    if (Canonical_ && MixPoint_ < NumVertices_ && MixPoint_ > 0){ //lambda should exist and be used
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << MixPoint_ << "_" << NumVertices_-MixPoint_ << ".MPX_matrix";
      Vd=load_MPX_matrix(Lambdanamestream.str(),Basis_);
    }

    const std::vector<MPXPair>& contractids=Canonical_ ? contract10 : contract11;

    for (uMPXInt p=1;p<=NumVertices_;++p){
      fetch_matrix(p,!Canonical_ || p<=MixPoint_); //if not canonical, load left, otherwise load left if left before mix point
      
      if (!Canonical_ || p>MixPoint_){ //if not canonical, or on right side of mix point
	//contract Vd on and decompose, if initial step for right canonical or non canonical, Vd will be empty
	MPXDecomposition decomp((Vd.empty() ? Current_.second : MPS_matrix(contract(Vd,0,matrix(),0,contractids))).left_shape().SVD());
	//record row vectors (Vdagger part) as new Vd
	Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
	Current_.second=std::move(decomp.ColumnMatrix);
	store_current(); //store new left canonical matrix
      }

      //at end of step store copy if requested
      if (!name.empty() && name!=MPSName_)
	matrix().store(filename(position(),matrix().is_left_shape(),name));
    }

    Canonical_=1;
    MixPoint_=NumVertices_;
    if (Vd.empty()){
      return 1.0;
    }
    else {
      return Vd.Trace();
    }
  }

  CanonicalType FiniteMPS::CheckFilesExist(){ //checks files exist for left, right or mixed, doesn't establish canon
    //assume nothing about defined Canonical_ or MixPoint_

    //check for a lambda file
    bool LAMBDA(0);
    uMPXInt MP(NumVertices_+1);

    for (uMPXInt L=1;L<=NumVertices_;++L){
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << L << "_" << NumVertices_-L << ".MPX_matrix";
      std::ifstream lambdafile;
      lambdafile.open(Lambdanamestream.str().c_str(),ios::in);
      if (lambdafile.is_open()){
	LAMBDA=1;
	MP=L;
      }
    }

    //then check for left and right parts
    uMPXInt LMAX(0);
    uMPXInt RMAX(0);

    for (uMPXInt L=1;L<=NumVertices_;++L){
      std::ifstream Lfile;
      Lfile.open(filename(L,1));
      if (Lfile.is_open()){
	LMAX=L;
      }
      else {
	break;
      }
    }

    for (uMPXInt R=1;R<=NumVertices_;++R){
      std::ifstream Rfile;
      Rfile.open(filename(R,0));
      if (Rfile.is_open()){
	RMAX=R;
      }
      else {
	break;
      }
    }

    //if lambda exists, then we assume mixed over all other types

    if (LAMBDA && LMAX>=MP && RMAX >=NumVertices_-MP){
      //by assumption...
      Canonical_=1;
      MixPoint_=MP;
      return CanonicalType::Mixed;
    }
    else if (LMAX==NumVertices_){
      //assume nothing
      MixPoint_=NumVertices_+1;
      return CanonicalType::Left;
    }
    else if (RMAX==NumVertices_){ //unlikely case, as drivers don't output right canonical files alone
      //by assumption...
      Canonical_=1;
      MixPoint_=0;
      return CanonicalType::Right;
    }
    else {
      //missing files
      std::cout <<"Error: Missing files for specified Finite MPS state!" <<std::endl;
      return CanonicalType::Error;
    }
    
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

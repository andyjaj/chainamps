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
#include "MPO_matrix.hpp"
#include "FiniteMPS.hpp"

namespace ajaj{

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& oldname, const std::string& newname, uMPXInt num, bool canon, uMPXInt mix_idx):Basis_(model_basis),MPSName_(oldname),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(canon),MixPoint_(mix_idx),Weight_(1.0){
    
    Canonization_=CheckFilesExist(newname);
    if (!(Canonization_==MPSCanonicalType::Error)){
      std::cout << "Found files for MPS " << MPSName_ << std::endl;
    }
    if (!newname.empty()) {//switch storage name if we have specified a new name
      MPSName_=newname;
      std::cout << "Switching storage to copy, " << MPSName_ <<std::endl;
    }
    std::cout << std::endl;
  }
  
  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt mix_idx): FiniteMPS(model_basis,name,std::string(),num,canon,mix_idx){}

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs) : Basis_(model_basis),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(0),MixPoint_(num+1),Canonization_(MPSCanonicalType::Non),Weight_(1.0) {
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
	MatrixCanonizations_.emplace_back(MPS_matrixCanonicalType::Non);
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
    //do a check, which will establish the (claimed by the individual matrices) MPSCanonicalType
    /*if (CheckFilesExist()==MPSCanonicalType::Error){
      return 0.0;
      }*/

    MPX_matrix Vd;
    //if (Canonical_ && MixPoint_ < NumVertices_ && MixPoint_ > 0){ //lambda should exist and be used
    if (Canonization_==MPSCanonicalType::Mixed){ //lambda should exist and be used
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << MixPoint_ << "_" << NumVertices_-MixPoint_ << ".MPX_matrix";
      Vd=load_MPX_matrix(Lambdanamestream.str(),Basis_);
    }

    //const std::vector<MPXPair>& contractids=Canonical_ ? contract10 : contract11;

    for (uMPXInt p=1;p<=NumVertices_;++p){
      //only load if not left canonical, or if we have requested a copy.
      if ((!name.empty() && name!=MPSName_) || MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Left){
	fetch_matrix(p,MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right);
	//if not left canonical, then we need to canonize
	if (MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Left){
	  //contract Vd on and decompose, if initial step for right canonical or non canonical, Vd will be empty
	  MPXDecomposition decomp((Vd.empty() ? Current_.second : MPS_matrix(contract(Vd,0,matrix(),0,MatrixCanonizations_[p-1]==MPS_matrixCanonicalType::Non ? contract11 : contract10))).left_shape().SVD());
	  //record row vectors (Vdagger part) as new Vd
	  Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
	  Current_.second=std::move(decomp.ColumnMatrix);
	  store_current(); //store new left canonical matrix
	  set_matrix_canonization(p,MPS_matrixCanonicalType::Left);
	  //if there is a non empty Vd and the next matrix was left canonical, we need to change the flag to trick it into being loaded.
	  if (p<NumVertices_ && !Vd.empty() && MatrixCanonizations_[p]==MPS_matrixCanonicalType::Left)
	    set_matrix_canonization(p+1,MPS_matrixCanonicalType::Non);

	}
	//if we want to store a copy...
	if (!name.empty() && name!=MPSName_)
	  matrix().store(filename(position(),matrix().is_left_shape(),name));
      }
    }
      /*fetch_matrix(p,!Canonical_ || p<=MixPoint_); //if not canonical, load left, otherwise load left if left before mix point
      
      if (!Canonical_ || p>MixPoint_){ //if not canonical, or on right side of mix point
	//contract Vd on and decompose, if initial step for right canonical or non canonical, Vd will be empty
	MPXDecomposition decomp((Vd.empty() ? Current_.second : MPS_matrix(contract(Vd,0,matrix(),0,contractids))).left_shape().SVD());
	//record row vectors (Vdagger part) as new Vd
	Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
	Current_.second=std::move(decomp.ColumnMatrix);
	store_current(); //store new left canonical matrix
	MatrixCanonizations_[p-1]=MPS_matrixCanonicalType::Left;
      }

      //at end of step store copy if requested
      if (!name.empty() && name!=MPSName_)
      matrix().store(filename(position(),matrix().is_left_shape(),name));
      }*/
    

    Canonical_=1;
    MixPoint_=NumVertices_;
    Canonization_=MPSCanonicalType::Left;
    if (Vd.empty()){
      return Weight_;
    }
    else {
      Weight_*=Vd.Trace();
      return Weight_;
    }
  }

  std::complex<double> FiniteMPS::makeRC(const std::string& name){
    std::cout << "Right canonising initial state..." <<std::endl;
    //do a check
    /*if (CheckFilesExist()==MPSCanonicalType::Error){
      return 0.0;
    }*/

    MPX_matrix U; //temp storage
    //if (Canonical_ && MixPoint_ < NumVertices_ && MixPoint_ > 0){ //lambda should exist and be used
    if (Canonization_==MPSCanonicalType::Mixed){ //lambda should exist and be used
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << MixPoint_ << "_" << NumVertices_-MixPoint_ << ".MPX_matrix";
      U=load_MPX_matrix(Lambdanamestream.str(),Basis_);
    }


    for (uMPXInt p=NumVertices_;p>0;--p){
      //only load if not right canonical, or if we have requested a copy.
      if ((!name.empty() && name!=MPSName_) || MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right){
	fetch_matrix(p,MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right);
	//if not right canonical, then we need to canonize
	if (MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right){
	  //contract U on and decompose, if initial step for left canonical or non canonical, U will be empty
	  MPXDecomposition decomp((U.empty() ? Current_.second : MPS_matrix(contract(matrix(),0,U,0,contract20))).right_shape().SVD());
	  U=std::move(contract(decomp.ColumnMatrix,0,MPX_matrix(Basis_,decomp.ColumnMatrix.Index(1),decomp.Values),0,contract10));
	  Current_.second=std::move(decomp.RowMatrix);
	  store_current(); //store new right canonical matrix
	  set_matrix_canonization(p,MPS_matrixCanonicalType::Right);
	  //if there is a non empty Vd and the next matrix was left canonical, we need to change the flag to trick it into being loaded.
	  if (p>1 && !U.empty() && MatrixCanonizations_[p-2]==MPS_matrixCanonicalType::Right)
	    set_matrix_canonization(p-1,MPS_matrixCanonicalType::Non);	  
	}
	//if we want to store a copy...
	if (!name.empty() && name!=MPSName_)
	  matrix().store(filename(position(),matrix().is_left_shape(),name));
      }
    }
    
    /*for (uMPXInt p=NumVertices_;p>0;--p){
      fetch_matrix(p,!Canonical_ || p<=MixPoint_); //if not canonical, load left, otherwise load left if left before mix point
      
      if (!Canonical_ || p<=MixPoint_){ //if not canonical, or on left side of mix point
	//contract U on and decompose, if initial step for left canonical or non canonical, U will be empty
	MPXDecomposition decomp((U.empty() ? Current_.second : MPS_matrix(contract(matrix(),0,U,0,contract20))).right_shape().SVD());
	U=std::move(contract(decomp.ColumnMatrix,0,MPX_matrix(Basis_,decomp.ColumnMatrix.Index(1),decomp.Values),0,contract10));
	Current_.second=std::move(decomp.RowMatrix);
	store_current(); //store new right canonical matrix
	MatrixCanonizations_[p-1]=MPS_matrixCanonicalType::Right;
      }

      //at end of step store copy if requested
      if (!name.empty() && name!=MPSName_)
	matrix().store(filename(position(),matrix().is_left_shape(),name));
	}*/

    Canonical_=1;
    MixPoint_=0;
    Canonization_=MPSCanonicalType::Right;
    if (U.empty()){
      return Weight_;
    }
    else {
      Weight_*=U.Trace();
      return Weight_;
    }
  }
  
  MPSCanonicalType FiniteMPS::CheckFilesExist(const std::string& newname){ //checks files exist for left, right or mixed, doesn't explicitly establish canonical status of individual matrices
    //assume nothing about defined Canonical_ or MixPoint_
    //check for lambda files, and establish the middle.
    bool LAMBDA(0);
    uMPXInt MP(NumVertices_+1);//initial guess is state is stored left shaped

    //assume even lengths only, and output in mixed form is centred.
    {
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
      std::ifstream lambdafile;
      lambdafile.open(Lambdanamestream.str().c_str(),ios::in);
      if (lambdafile.is_open()){
	LAMBDA=1;
	MP=NumVertices_/2;
	if (!newname.empty()){
	  std::stringstream newLambdanamestream;
	  newLambdanamestream << newname << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
	  //load and re-store lambda
	  load_MPX_matrix(Lambdanamestream.str(),Basis_).store(newLambdanamestream.str());
	}
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
	if (!newname.empty()){
	  load_MPS_matrix(filename(L,1),Basis_).store(filename(L,1,newname));
	}
      }
      else {
	break;
      }
    }

    for (uMPXInt R=NumVertices_;R>=1;--R){
      std::ifstream Rfile;
      Rfile.open(filename(R,0));
      if (Rfile.is_open()){
	RMAX=NumVertices_-R+1;
	if (!newname.empty()){
	  load_MPS_matrix(filename(R,0),Basis_).store(filename(R,0,newname));
	}
      }
      else {
	break;
      }
    }

    //if lambda exists, then we assume mixed over all other types

    if (LAMBDA && LMAX>=MP && RMAX >=NumVertices_-MP){
      //by assumption...
      Canonical_=1; //override previous value
      MixPoint_=MP;

      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(MP,MPS_matrixCanonicalType::Left);
      std::vector<MPS_matrixCanonicalType> Rpart(NumVertices_-MP,MPS_matrixCanonicalType::Right);
      MatrixCanonizations_.insert(MatrixCanonizations_.end(),Rpart.begin(),Rpart.end());
      return MPSCanonicalType::Mixed;
    }
    else if (LMAX==NumVertices_){
      //assume nothing
      MixPoint_=NumVertices_; //Check this?
      MPS_matrixCanonicalType ThisType= Canonical_ ? MPS_matrixCanonicalType::Left : MPS_matrixCanonicalType::Non;
      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(NumVertices_,ThisType);
      return MPSCanonicalType::Left;
    }
    else if (RMAX==NumVertices_){ //unlikely case, as drivers don't output right canonical files alone
      //by assumption...
      Canonical_=1;
      MixPoint_=0;
      MPS_matrixCanonicalType ThisType= Canonical_ ? MPS_matrixCanonicalType::Right : MPS_matrixCanonicalType::Non;
      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(NumVertices_,ThisType);
      return MPSCanonicalType::Right;
    }
    else {
      //missing files
      std::cout <<"Error: Missing files for specified Finite MPS state '" << MPSName_<< "'" << std::endl;
      std::cout << "LAMBDA "<< MP << ", " << LMAX << " " << RMAX << std::endl;

      return MPSCanonicalType::Error;
    }
    
  }

  void FiniteMPS::set_matrix_canonization(uMPXInt i,const MPS_matrixCanonicalType& c){
    MatrixCanonizations_[i-1]=c;
    update_MPS_canonization_status();
  }

  void FiniteMPS::update_MPS_canonization_status(){
    //go through list of matrix canonizations and fix the overall MPS status
    
    MPS_matrixCanonicalType CurrentStatus=MatrixCanonizations_[0];

    Canonization_=MPSCanonicalType::Error;
    
    for (uMPXInt i=1; i<=NumVertices_;++i){
      if (MatrixCanonizations_[i-1]==MPS_matrixCanonicalType::Non){
	Canonization_=MPSCanonicalType::Non;
	break;
      }
      else if (i>1) {
	if ((MatrixCanonizations_[i-1]==MPS_matrixCanonicalType::Left) && (MatrixCanonizations_[i-2]==MPS_matrixCanonicalType::Right)){
	  Canonization_=MPSCanonicalType::Non;
	  break;
	}
	else if ((MatrixCanonizations_[i-1]==MPS_matrixCanonicalType::Right) && (MatrixCanonizations_[i-2]==MPS_matrixCanonicalType::Left)){
	  Canonization_= (MixPoint_== i-1 ? MPSCanonicalType::Mixed : MPSCanonicalType::Error);
	}
	else if (i==NumVertices_ && Canonization_!=MPSCanonicalType::Mixed){
	  if (MatrixCanonizations_[i-1]==MPS_matrixCanonicalType::Left){
	    Canonization_=MPSCanonicalType::Left;
	  }
	  else {
	    Canonization_=MPSCanonicalType::Right;
	  }
	}
      }
    }
  }
  
  std::complex<double> ApplySingleVertexOperatorToMPS(const MPO_matrix& Op, FiniteMPS& F, uMPXInt vertex, const MPSCanonicalType& RequestedCanonization){

    if (F.Canonization_==MPSCanonicalType::Mixed || RequestedCanonization==MPSCanonicalType::Error ){
      std::cout << "Not supported yet." <<std::endl;
      exit(1);
    }
    
    //apply operator at vertex
    F.fetch_matrix(vertex,F.MatrixCanonizations_[vertex-1]!=MPS_matrixCanonicalType::Right);

    F.Current_.second=std::move(contract(Op,0,F.Current_.second,0,F.MatrixCanonizations_[vertex-1]!=MPS_matrixCanonicalType::Right ? contract20 : contract21).RemoveDummyIndices({{1,2}}));
    //now vertex is no longer canonical
    F.Current_.second.ShiftNumRowIndices(2); //make it standard format for a non canonical matrix, will always be stored as 'Left'!  
    F.set_matrix_canonization(vertex,MPS_matrixCanonicalType::Non);
    F.store_current();
    
    if (RequestedCanonization==MPSCanonicalType::Left){
      F.makeLC();
    }
    else if (RequestedCanonization==MPSCanonicalType::Right){
      F.makeRC();
    }
    
    return F.weight();
  }
  
}

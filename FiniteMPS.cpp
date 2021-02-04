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

  //at the moment will only find consistent files if lambda is at centre of system
  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& oldname, const std::string& newname, uMPXInt num, bool canon, uMPXInt lp):Basis_(model_basis),TotalCharges_(Basis_.Identity()),MPSName_(oldname),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),HasLambda_(0),Weight_(1.0){

    //check files exist will set HasLambda_ and will set LambdaPosition_
    if (check_files_exist(canon,newname)==MPSCanonicalType::Error || (HasLambda_ && (LambdaPosition_!=lp))){
      std::cout << "Couldn't find consistent files for this finite MPS" <<std::endl;
      exit(1);
    }
    else {
      std::cout << "Found files for MPS " << MPSName_ << std::endl;

      if (!newname.empty()) {//switch storage name if we have specified a new name
	MPSName_=newname;
	std::cout << "Switching storage to copy, " << MPSName_ <<std::endl;
      }
    }
  }

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt lp): FiniteMPS(model_basis,name,std::string(),num,canon,lp){}

  FiniteMPS::FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs) : Basis_(model_basis),TotalCharges_(Basis_.Identity()),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonization_(MPSCanonicalType::Non),HasLambda_(0),Weight_(1.0) {
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
	MatrixCanonizations_.emplace_back(MPS_matrixCanonicalType::Non);       //can't guarantee that we are canonical
      }
      TotalCharges_=Current_.second.Index(2)[0];       //get TotalCharges_
    }
    else {
      check_files_exist(0,name); //look to load a state. Assume non-canonical i.e. unnormalised to be safe
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

  MPX_matrix FiniteMPS::fetch_lambda() const {
    //if the state has a defined lambda matrix, and the position isn't ill defined
    if (HasLambda_  && LambdaPosition_<=NumVertices_){
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << LambdaPosition_ << "_" << NumVertices_-LambdaPosition_ << ".MPX_matrix";
      return load_MPX_matrix(Lambdanamestream.str(),Basis_);
    }
    else
      return MPX_matrix();
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

  std::pair<MPX_matrix,MPX_matrix> FiniteMPS::left_canonize_to(uMPXInt i,const std::string& name){
    //this internal function will sweep through from the left, left canonizing up to and including i, or the end of the system, whichever is smaller

    uMPXInt start_position=NumVertices_+1; //
    for (uMPXInt p=1;p<=std::min(NumVertices_,i);++p){
      //search through til we hit the first matrix that isn't left canonical
      //once reached all matrices after will need treating
      //only load if not left canonical, or if we have requested a copy.
      if (MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Left){
	start_position=p; //we now know the position of the first matrix that needs to be treated
	break;
      }
      if (!name.empty() && name!=MPSName_){
	fetch_matrix(p,MatrixCanonizations_.at(p-1)!=MPS_matrixCanonicalType::Right);
	matrix().store(filename(position(),matrix().is_left_shape(),name)); //save copy if requested
      }
    }

    MPX_matrix Vd;
    MPX_matrix DecompLambda;

    for (uMPXInt p=start_position;p<=std::min(NumVertices_,i);++p){
      fetch_matrix(p,MatrixCanonizations_.at(p-1)!=MPS_matrixCanonicalType::Right); //get the matrix
      if (!Vd.empty()) { //have a Vd from previous iteration, must therefore also have a DecompLambda. Combine them
	Vd=contract(DecompLambda,0,Vd,0,contract10);
      }
      if (HasLambda_ && (p-1==LambdaPosition_) ) {
	MPX_matrix Lambda=fetch_lambda(); //must also incorporate any mixed state lambda immediately to the left of this matrix
	if (!Vd.empty()) {Vd=contract(Vd,0,Lambda,0,contract10);}
	else {Vd=std::move(Lambda);}
	HasLambda_=0; //have now incorporated lambda, so the stored one is no longer valid
      }

      //now we can do the decomposition      
      MPXDecomposition decomp((Vd.empty() ? Current_.second : MPS_matrix(contract(Vd,0,matrix(),0,MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right ? contract11 : contract10))).left_shape().SVD());
      Vd=std::move(decomp.RowMatrix); //decomp row vectors
      DecompLambda=MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values); //decomp lambda
      Current_.second=std::move(decomp.ColumnMatrix); //update matrix for this position
      store_current(); //store new left canonical matrix 
      set_matrix_canonization(p,MPS_matrixCanonicalType::Left); //record that this matrix is now left canonical

      if (!name.empty() && name!=MPSName_){ //if copy is requested store this matrix again
	matrix().store(filename(position(),matrix().is_left_shape(),name));
      }
      
    }

    if (!Vd.empty()) { //if we have had to do a decomposition, then the state overall must now be non canonical
      Canonization_=MPSCanonicalType::Non;
    }
    return std::pair<MPX_matrix,MPX_matrix>(DecompLambda,Vd);
  }

  std::complex<double> FiniteMPS::makeLC(const std::string& name){
    std::pair<MPX_matrix,MPX_matrix> LambdaVd=left_canonize_to(NumVertices_);
    //Canonical_=1;
    //LambdaPosition_=NumVertices_;
    HasLambda_=0;
    Canonization_=MPSCanonicalType::Left;
    if (LambdaVd.second.empty()){
      return Weight_;
    }
    else {
      Weight_*=contract(LambdaVd.first,0,LambdaVd.second,0,contract10).Trace();
      return Weight_;
    }
  }

  std::pair<MPX_matrix,MPX_matrix> FiniteMPS::right_canonize_to(uMPXInt i,const std::string& name){
    //this internal function will sweep through from the right, right canonizing up to and including i, or the end of the system, whichever is smaller

    uMPXInt start_position=0; //
    for (uMPXInt p=NumVertices_;p>=std::max(i,uMPXInt(1));--p){
      //search through til we hit the first matrix that isn't right canonical
      //once reached all matrices after will need treating
      //only load if not right canonical, or if we have requested a copy.
      if (MatrixCanonizations_[p-1]!=MPS_matrixCanonicalType::Right){
	start_position=p; //we now know the position of the first matrix that needs to be treated
	break;
      }
      if (!name.empty() && name!=MPSName_){
	fetch_matrix(p,MatrixCanonizations_.at(p-1)!=MPS_matrixCanonicalType::Right);
	matrix().store(filename(position(),matrix().is_left_shape(),name)); //save copy if requested
      }
    }

    MPX_matrix U;
    MPX_matrix DecompLambda;
    
    for (uMPXInt p=start_position;p>=std::max(i,uMPXInt(1));--p){
      fetch_matrix(p,MatrixCanonizations_.at(p-1)!=MPS_matrixCanonicalType::Right); //get the matrix
      if (!U.empty()) { //have a Vd from previous iteration, must therefore also have a DecompLambda. Combine them
	U=contract(U,0,DecompLambda,0,contract10);
      }
      if (HasLambda_ && (p-1==LambdaPosition_) ) {
	MPX_matrix Lambda=fetch_lambda(); //must also incorporate any mixed state lambda immediately to the left of this matrix
	if (!U.empty()) {U=contract(Lambda,0,U,0,contract10);}
	else {U=std::move(Lambda);}
	HasLambda_=0; //have now incorporated lambda, so the stored one is no longer valid
      }

      //now we can do the decomposition
      MPXDecomposition decomp((U.empty() ? Current_.second : MPS_matrix(contract(matrix(),0,U,0,contract20))).right_shape().SVD());

      U=std::move(decomp.ColumnMatrix); //decomp row vectors
      DecompLambda=MPX_matrix(Basis_,decomp.ColumnMatrix.Index(1),decomp.Values); //decomp lambda
      Current_.second=std::move(decomp.RowMatrix); //update matrix for this position
      store_current(); //store new left canonical matrix 
      set_matrix_canonization(p,MPS_matrixCanonicalType::Right); //record that this matrix is now left canonical

      if (!name.empty() && name!=MPSName_){ //if copy is requested store this matrix again
	matrix().store(filename(position(),matrix().is_left_shape(),name));
      }
      
    }

    if (!U.empty()) { //if we have had to do a decomposition, then the state overall we must now be non canonical
      Canonization_=MPSCanonicalType::Non;
    }
    return std::pair<MPX_matrix,MPX_matrix>(U,DecompLambda);
  }

  std::complex<double> FiniteMPS::makeRC(const std::string& name){
    std::pair<MPX_matrix,MPX_matrix> ULambda=right_canonize_to(1);
    Canonization_=MPSCanonicalType::Right;
    HasLambda_=0;
    if (ULambda.first.empty()){
      return Weight_;
    }
    else {
      Weight_*=contract(ULambda.first,0,ULambda.second,0,contract10).Trace();
      return Weight_;
    }
  }

  double FiniteMPS::mixed(uMPXInt mixposition,const std::string& name){

    if(mixposition>NumVertices_){
      std::cout << "Invalid position " << mixposition << " > " << NumVertices_ <<std::endl;
      exit(1);
    }

    //needs a check to decide which end to start from if the current state is already mixed canonical
    //this will help sweep up lambda before issues occur

    std::complex<double> weight=makeLC(name); //this should sweep up any existing lambdas, but currently non optimal

    std::pair<MPX_matrix,MPX_matrix> ULambda=right_canonize_to(mixposition+1,name);
    //shift U to left
    fetch_matrix(mixposition); //fetch, should be left canonical default
    if (!ULambda.first.empty()){
      Current_.second=contract(Current_.second,0,ULambda.first,0,contract20);
      store_current(); //store new left canonical matrix    
      //store Lambda
      std::stringstream Lambdanamestream;
      Lambdanamestream << MPSName_ << "_Lambda_" << mixposition << "_" << NumVertices_-mixposition << ".MPX_matrix";
      //load and re-store lambda
      ULambda.second.store(Lambdanamestream.str());
    
      if (!name.empty() && name!=MPSName_){ //store copy if requested
	matrix().store(filename(position(),matrix().is_left_shape(),name));
	std::stringstream newLambdanamestream;
	newLambdanamestream << name << "_Lambda_" << mixposition << "_" << NumVertices_-mixposition << ".MPX_matrix";
	//load and re-store lambda
	ULambda.second.store(newLambdanamestream.str());
      }

      //finish in a well defined mixed state
      HasLambda_=1;
      LambdaPosition_=mixposition;
      Canonization_=MPSCanonicalType::Mixed;
      return ULambda.second.norm()*(weight * conj(weight)).real();
    }
    else {
      return (weight * conj(weight)).real();
    }
  }

  std::string FiniteMPS::canon_type_string() const {
    switch(Canonization_){
    case MPSCanonicalType::Non :
      return(std::string("Non Canonical"));
    case MPSCanonicalType::Mixed :
      return(std::string("Mixed"));
    case MPSCanonicalType::Right :
      return(std::string("Right"));
    case MPSCanonicalType::Left :
      return(std::string("Left"));
    default :
      return(std::string("Error"));
    }
  }
  
  MPSCanonicalType FiniteMPS::check_files_exist(bool canon, const std::string& newname){ //checks files exist for left, right or mixed, doesn't explicitly establish canonical status of individual matrices
    //assume nothing about LambdaPosition_
    //use canon value to decide on status
    //make copy if newname not empty

    //Possible cases: left canonical, Lambda at right end, or none
    //                right canonical, Lambda at left end, or none
    //                mixed canonical, must have a Lambda, favour the one nearest the centre if several files exist
    //                non canonical, may or may not have a Lambda.
    // For a consistent state, any lambda must be sandwiched between L and R type files.

    //Find the lambda files
    //Find the largest L file and the largest R file

    
    //check for lambda files, and establish the middle.
    //bool LAMBDA(0);
    //uMPXInt MP(NumVertices_+1);

    //Find L,R shape files first
    //then check for left and right parts
    uMPXInt LMAX(0);
    uMPXInt RMAX(0);

    for (uMPXInt L=1;L<=NumVertices_;++L){
      std::ifstream Lfile;
      Lfile.open(filename(L,1));
      if (!Lfile.is_open()){
	break;
      }
      else {
	MPS_matrix M=load_MPS_matrix(filename(L,1),Basis_);
	if (M.isConsistent()) {
	  LMAX=L;
	  if (!newname.empty()){
	    M.store(filename(L,1,newname)); 
	  }
	  if (L==NumVertices_){
	    TotalCharges_=M.Index(2).at(0);
	  }
	}
	else {
	  break;
	}
      }
    }

    for (uMPXInt R=NumVertices_;R>=1;--R){
      std::ifstream Rfile;
      Rfile.open(filename(R,0));
      if (!Rfile.is_open()){
	break;
      }
      else {
	MPS_matrix M=load_MPS_matrix(filename(R,0),Basis_);
	if (M.isConsistent()) {
	  LMAX=NumVertices_-R+1;
	  if (!newname.empty()){
	    M.store(filename(R,0,newname)); 
	  }
	  if (R==NumVertices_){
	    if (LMAX==NumVertices_ && (TotalCharges_!=M.Index(2).at(0)) ){
	      std::cout << "Left and right shaped matrices for position " << NumVertices_ << " give different total charges!" << std::endl;
	      TotalCharges_.print();
	      M.Index(2).at(0).print();
	      exit(1);
	    }
	    else {
	      TotalCharges_=M.Index(2).at(0);
	    }
	  }
	}
	else {
	  break;
	}
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

    //This version only looks for lambda files in the middle of the system
    
    std::stringstream Lambdanamestream;
    Lambdanamestream << MPSName_ << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
    std::ifstream lambdafile;
    lambdafile.open(Lambdanamestream.str().c_str(),ios::in);
    if (lambdafile.is_open()){
      HasLambda_=1;
      LambdaPosition_=NumVertices_;
      if (!newname.empty()){
	std::stringstream newLambdanamestream;
	newLambdanamestream << newname << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
	//load and re-store lambda
	load_MPX_matrix(Lambdanamestream.str(),Basis_).store(newLambdanamestream.str());
      }
    }

    //if lambda exists, then we assume mixed over all other types

    if (HasLambda_ && LMAX>=LambdaPosition_ && RMAX >=NumVertices_-LambdaPosition_){

      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(LambdaPosition_,MPS_matrixCanonicalType::Left);
      std::vector<MPS_matrixCanonicalType> Rpart(NumVertices_-LambdaPosition_,MPS_matrixCanonicalType::Right);
      MatrixCanonizations_.insert(MatrixCanonizations_.end(),Rpart.begin(),Rpart.end());
      //update_MPS_canonization_status();
      Canonization_= canon ? MPSCanonicalType::Mixed : MPSCanonicalType::Non; // Just in case we somehow defined a non canonical state with a lambda
    }
    else if (LMAX==NumVertices_){
      //assume nothing
      //LambdaPosition_=NumVertices_; //Check this?
      MPS_matrixCanonicalType ThisType= canon ? MPS_matrixCanonicalType::Left : MPS_matrixCanonicalType::Non;
      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(NumVertices_,ThisType);
      Canonization_= canon ? MPSCanonicalType::Left : MPSCanonicalType::Non;    }
    else if (RMAX==NumVertices_){ //unlikely case, as drivers don't output right canonical files alone
      //by assumption...
      //Canonical_=1;
      //LambdaPosition_=0;
      MPS_matrixCanonicalType ThisType= canon ? MPS_matrixCanonicalType::Right : MPS_matrixCanonicalType::Non;
      MatrixCanonizations_=std::vector<MPS_matrixCanonicalType>(NumVertices_,ThisType);
      Canonization_= canon ? MPSCanonicalType::Right : MPSCanonicalType::Non;    }
    else {
      //missing files
      std::cout <<"Error: Missing files for specified Finite MPS state '" << MPSName_<< "'" << std::endl;
      std::cout << "LAMBDA "<< LambdaPosition_ << ", " << LMAX << " " << RMAX << std::endl;
      //update_MPS_canonization_status();
      Canonization_=MPSCanonicalType::Error;
    }
    return Canonization_;
    
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
	  Canonization_= (LambdaPosition_== i-1 ? MPSCanonicalType::Mixed : MPSCanonicalType::Error);
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

  MPS_matrix ConstFiniteMPS::matrix(uMPXInt i) const {
    //check buffer
    if (i <= NumVertices_ && i >= 1){ //inside range
      if (Canonization_==MPSCanonicalType::Mixed){
	//be careful with lambda matrix?
	//for the moment just fail spectacularly
	std::cout <<"Not supported yet!" <<std::endl;
	exit(1);
	//return MPS_matrix();
      }
      else
	return load_MPS_matrix(filename(i,MatrixCanonizations_.at(i-1)==MPS_matrixCanonicalType::Right ? 0 : 1),Basis_);
    }
	
    else {
      std::cout <<"Outside range of finite MPS! " << i << " : " << NumVertices_ << std::endl;
      exit(1);
    }
  }

  std::string ConstFiniteMPS::filename(uMPXInt i, bool Left, const std::string& name) const{
    std::stringstream namestream;
    if (Left)
      namestream << (name.empty() ? MPSName_ : name) << "_Left_" << i << ".MPS_matrix";
    else
      namestream << (name.empty() ? MPSName_ : name) << "_Right_" << NumVertices_-i+1 << ".MPS_matrix";
    
    return namestream.str();
  }
  
}

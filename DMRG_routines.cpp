#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <iomanip>

#include "arpack_interface.hpp" //arpack
#include "sparse_interface.hpp" //SparseHED
#include "states.hpp" //EigenStateArray etc.
#include "MPX.hpp" //MPX_matrix etc.
#include "DMRG_routines.hpp" //MPX_matrix etc.
#include "data.hpp"

namespace ajaj {
#ifdef TIMING
  static double formation_time_millisecs (0.0);
  static size_t formations (0);

  static double decomp_time_millisecs (0.0);
  static size_t decomps (0);

  static double vec_time_millisecs (0.0);
  static size_t vec_formations (0);

  static double c1_time_ms (0.0);
  static size_t c1_reps(0);

  static double c2_time_ms (0.0);
  static size_t c2_reps(0);

  static double c3_time_ms (0.0);
  static size_t c3_reps(0);

  static double c4_time_ms (0.0);
  static size_t c4_reps(0);

  static size_t num_op_x(0);

  void dmrg_print_time_info() {
    std::cout << std::setprecision(8) << std::endl << "TIME INFO" <<std::endl;
    std::cout << "AVERAGE LEFT AND RIGHT FORMATION TIME=" << formation_time_millisecs/formations << " ms" << std::endl;
    std::cout << "AVERAGE VECTOR FORMATION TIME=" << vec_time_millisecs/vec_formations/1000.0 << " ms" << std::endl;
    std::cout << "AVERAGE Contraction 1 TIME=" << c1_time_ms/c1_reps << " ms" << std::endl;
    std::cout << "AVERAGE Contraction 2 TIME=" << c2_time_ms/c2_reps << " ms" << std::endl;
    std::cout << "AVERAGE Reshape TIME=" << c3_time_ms/c3_reps/1000.0 << " ms" << std::endl;
    std::cout << "AVERAGE Extract TIME=" << c4_time_ms/c4_reps/1000.0 << " ms" << std::endl;
    std::cout << "NUM OP* x for this step=" << num_op_x <<std::endl;
    std::cout << "AVERAGE DECOMP TIME=" << decomp_time_millisecs/decomps/1000.0 << " s" << std::endl;

    std::cout <<std::endl;
  }
#endif

  BlocksStructure::BlocksStructure(const std::string&  Name, const EigenStateArray& Spectrum, uMPXInt num_vertices, uMPXInt numLeft, uMPXInt numMiddle) : Name_(Name), SpectrumPtr_(&Spectrum),LeftBlock(Spectrum),RightBlock(Spectrum),num_vertices_(num_vertices),numLeft_(numLeft),numMiddle_(numMiddle) {
    if (num_vertices % 2 !=0){
      std::cout << "Currently only even length systems are supported!" <<std::endl;
      exit(1);
    }
    if (num_vertices!=0){
      //needs new initialiser that sets things up correctly if sizes aren't zero
      if(load_left_block() || load_right_block()){
	std::cout << "Invalid block files!" << std::endl;
	exit(1);
      }
      ;
    }
  }
  
  bool BlocksStructure::save_left_block(uMPXInt l){
    std::stringstream leftname;
    leftname << Name_ << "_Left_" << l << ".BLOCK";
    return LeftBlock.store(leftname.str());
  }
  bool BlocksStructure::save_right_block(uMPXInt r){
    std::stringstream rightname;
    rightname << Name_ << "_Right_" << r << ".BLOCK";
    return RightBlock.store(rightname.str());
  }
  bool BlocksStructure::save_left_block(){
    std::stringstream leftname;
    leftname << Name_ << "_Left_" << left_size() << ".BLOCK";
    return LeftBlock.store(leftname.str());
  }
  bool BlocksStructure::save_right_block(){
    std::stringstream rightname;
    rightname << Name_ << "_Right_" << right_size() << ".BLOCK";
    return RightBlock.store(rightname.str());
  }
  bool BlocksStructure::save_blocks(){
    return (save_left_block() && save_right_block());
  }
  bool BlocksStructure::load_left_block(){
    std::stringstream leftname;
    leftname << Name_ << "_Left_" << left_size() << ".BLOCK";
    LeftBlock=std::move(load_MPX_matrix(leftname.str(),*SpectrumPtr_));
    return LeftBlock.isEmpty();
  }
  bool BlocksStructure::load_right_block(){
    std::stringstream rightname;
    rightname << Name_ << "_Right_" << right_size() << ".BLOCK";
    RightBlock=std::move(load_MPX_matrix(rightname.str(),*SpectrumPtr_));
    return RightBlock.isEmpty();

  }
  void BlocksStructure::ends(MPX_matrix&& leftend, MPX_matrix&& rightend) {
    std::stringstream leftname;
    leftname << Name_ << "_Left_0" << ".BLOCK";
    std::stringstream rightname;
    rightname << Name_ << "_Right_0" << ".BLOCK";
    leftend.store(leftname.str());
    rightend.store(rightname.str());
    //so store some dummy blocks,
    //don't knwo anything about the middle
  }
  void BlocksStructure::initial_2(MPX_matrix&& newleft, MPX_matrix&& newright) {
    ///special
    //for pre storing the first left and right blocks
    num_vertices_=2;
    numMiddle_=0;
    numLeft_=1;
    LeftBlock=std::move(newleft);
    RightBlock=std::move(newright);
    save_left_block();
    save_right_block();
    //so store some dummy blocks,
    //don't know anything about the middle, but prep
  }
  void BlocksStructure::set_2() {
    if (num_vertices_==2) {
      numMiddle_=2;
      numLeft_=0;
      load_left_block();
      load_right_block();
    }
  }
  void BlocksStructure::insert_2(MPX_matrix&& newleft, MPX_matrix&& newright) {
    numMiddle_=2;
    num_vertices_+=2;
    ++numLeft_;
    LeftBlock=std::move(newleft);
    RightBlock=std::move(newright);
    save_left_block();save_right_block();
  }
  void BlocksStructure::insert_2() {
    numMiddle_=2;
    num_vertices_+=2;
  }
  void BlocksStructure::shift_left(MPX_matrix&& newright){
    if (left_size()>0) {
      --numLeft_;
      RightBlock=std::move(newright);
      save_right_block();
      load_left_block();
    }
  }
  void BlocksStructure::shift_right(MPX_matrix&& newleft){
    if (right_size()>0) {
      ++numLeft_;
      LeftBlock=std::move(newleft);
      save_left_block();
      load_right_block();
    }
  }
  void BlocksStructure::shift_left(const std::string& OtherLeftName, MPX_matrix&& newright){
    if (left_size()>0) {
      --numLeft_;
      RightBlock=std::move(newright);
      save_right_block();
      std::stringstream leftname;
      leftname << OtherLeftName << "_Left_" << left_size() << ".BLOCK";
      LeftBlock=std::move(load_MPX_matrix(leftname.str(),*SpectrumPtr_));
    }
  }
  void BlocksStructure::shift_right(MPX_matrix&& newleft, const std::string& OtherRightName){
    if (right_size()>0) {
      ++numLeft_;
      LeftBlock=std::move(newleft);
      save_left_block();
      std::stringstream rightname;
      rightname << OtherRightName << "_Right_" << right_size() << ".BLOCK";
      RightBlock=std::move(load_MPX_matrix(rightname.str(),*SpectrumPtr_));
    }
  }
  void BlocksStructure::set_blocks(uMPXInt ls, MPX_matrix&& newleft, MPX_matrix&& newright){ //if blocks are generated by an algorithm other than grow
    numLeft_=ls;
    LeftBlock=std::move(newleft);
    RightBlock=std::move(newright);
    save_left_block();
    save_right_block();
  }

  SuperBlock::SuperBlock(const std::string& Name, const MPO_matrix& H, const State& TargetState, uMPXInt num_vertices, uMPXInt numLeft, uMPXInt numMiddle) : BlocksStructure(Name,H.GetPhysicalSpectrum(),num_vertices,numLeft,numMiddle),H_ptr_(&H),TargetState_(TargetState),CentralDecomposition(H.getPhysicalSpectrum()) {
      std::stringstream dnamestream;
      dnamestream << getName();// << "_Density_Matrix.dat";
      DensityFileName_=dnamestream.str();

      if (num_vertices!=0){//we are resuming a run

	//collect the parts from file storage and do a trivial SVD on the lambda matrix to extract the singular values cleanly.
	//what we will lose here, relative to a run that hasn't been stopped and resumed from file, is any knowledge of previous truncations.
	std::stringstream Lambdanamestream;
	Lambdanamestream << Name << "_Lambda_" << num_vertices/2 << "_" << num_vertices/2 << ".MPX_matrix";
	
	CentralDecomposition.Values=std::move(vector_of_reals(load_MPX_matrix(Lambdanamestream.str(),basis()).GetDiagonal()));

	if (CentralDecomposition.Values.size()==0){
	  std::cout << "Invalid file " << Lambdanamestream.str() << std::endl;
	  exit(1);
	}
	
	std::stringstream Lnamestream;
	Lnamestream << Name << "_Left_" << num_vertices/2 << ".MPS_matrix";
	std::stringstream Rnamestream;
	Rnamestream << Name << "_Right_" << num_vertices/2 << ".MPS_matrix";
	
	CentralDecomposition.LeftMatrix=std::move(load_MPS_matrix(Lnamestream.str(),basis()).left_shape());
	CentralDecomposition.RightMatrix=std::move(load_MPS_matrix(Rnamestream.str(),basis()).right_shape());

	if (CentralDecomposition.LeftMatrix.isEmpty()){
	  std::cout << "Invalid file " << Lnamestream.str() << std::endl;
	  exit(1);
	}
	if (CentralDecomposition.RightMatrix.isEmpty()){
	  std::cout << "Invalid file " << Rnamestream.str() << std::endl;
	  exit(1);
	}
	
	//we still require the previous lambda
       
	if (num_vertices==2){ //resuming from the very first step, need to generate trivial previous lambda
	  previous_lambda_=std::vector<double>(1,1.0);
	}
	else {
	  std::stringstream PrevLambdanamestream;
	  PrevLambdanamestream << Name << "_Lambda_" << num_vertices/2 -1 << "_" << num_vertices/2 -1 << ".MPX_matrix";
	  
	  previous_lambda_=std::move(vector_of_reals(load_MPX_matrix(PrevLambdanamestream.str(),basis()).GetDiagonal()));
	  
	  if (previous_lambda_.size()==0){
	    std::cout << "Invalid file " << PrevLambdanamestream.str() << std::endl;
	    exit(1);
	  }
	}

	std::cout << "Info on loaded data:" <<std::endl;
	std::cout << "Lambda size " << CentralDecomposition.Values.size() <<std::endl;
	std::cout << "Previous Lambda size " << previous_lambda_.size() <<std::endl;
	std::cout << "Left MPS" <<std::endl;
	CentralDecomposition.LeftMatrix.print_indices();
	std::cout << "Right MPS" <<std::endl;
	CentralDecomposition.RightMatrix.print_indices();
	
	//build prediction vector and calculate fidelity
	pred_=MakePrediction(CentralDecomposition,previous_lambda_);//allows checking overlap of current state with the previous one.
	fidelity_=CheckConvergence(pred_,previous_lambda_);
	
      }   
  }
  
  Data SuperBlock::initialise(uMPXInt chi, double truncation){
    //Assumes any previous contents of the blocks are junk.
    if(!H_ptr_->isConsistent()){std::cout << "Hamiltonian MPO is malformed. Aborting." << std::endl; H_ptr_->print_indices_values(); exit(1);}
    //first we store some dummy blocks at the left and right ends (position 0), which are useful
    ends(MakeDummyLeftBlock(getH(),getTargetState()),MakeDummyRightBlock(getH(),getTargetState()));
    //we need the left and right MPO forms of H for open boundary conditions
    MPO_matrix LeftH(LeftOpenBCHamiltonian(getH())); //left Hamiltonian MPO
    MPO_matrix RightH(RightOpenBCHamiltonian(getH())); //right Hamiltonian MPO
    
    //why not save these to file as well?
    LeftH.store("LeftH.MPO_matrix");
    RightH.store("RightH.MPO_matrix");
    std::cout << "Initializing DMRG with two vertices..." << std::endl;
    //chi and truncation zero? take as an indication that we should do a full diag...
    uMPXInt bonddim = (chi<1 && truncation <=0.0) ? getH().GetPhysicalSpectrum().size() : chi;
    Data two_vertex_energy;
    //solve
    CentralDecomposition=TwoVertexSVD(TwoVertexInitialWavefunction(LeftH,RightH,TargetState_,two_vertex_energy),bonddim,truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,1);

    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    //at this stage we can use the generated left and right Hamiltonians to save the left and right blocks
    //update superblock values
    initial_2(MakeInitialLeftBlock(LeftH,CentralDecomposition.LeftMatrix),MakeInitialRightBlock(RightH,CentralDecomposition.RightMatrix));
    //at this point there are 2 vertices, and we have pre stored the blocks they contribute to
    previous_lambda_=std::vector<double>(1,1.0);
    pred_=MakePrediction(CentralDecomposition,previous_lambda_);//allows checking overlap of current state with the previous one.
    fidelity_=CheckConvergence(pred_,previous_lambda_);
    double S_E(entropy(CentralDecomposition.Values));
    std::cout << "Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    two_vertex_energy.Real_measurements.push_back(S_E);
    two_vertex_energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    two_vertex_energy.Real_measurements.push_back(fidelity_);
    return two_vertex_energy;
  }

  Data SuperBlock::grow_two_vertex(uMPXInt chi, double truncation){
    if (size()<2){ //haven't initialised yet?
      return initialise(chi,truncation); //store dummies, and first left and right blocks
    }
    else if (size()==2){ //special case, pre stored blocks, but still need to add middle blocks
      insert_2();
    }
    else { /*size()>4, make new left and right blocks, insert dummy vertices ready for solving*/
      insert_2(MakeLeftBlock(getLeftBlock(),getH(),CentralDecomposition.LeftMatrix),MakeRightBlock(getRightBlock(),getH(),CentralDecomposition.RightMatrix));
    }
    //Prediction Next(MakePrediction(CentralDecomposition,previous_lambda_));
    //std::cout << "Checking overlap" << std::endl; //actually we are checking the convergence of the previous iteration...
    //fidelity_=CheckConvergence(Next,previous_lambda_);
    //we don't need the central decomp anymore, so we could steal/swap from it
    previous_lambda_=CentralDecomposition.Values;
    //update numbers
    //solve
    Data energy;
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy,&(pred_.Guess)),chi,truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,1);
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    pred_=MakePrediction(CentralDecomposition,previous_lambda_);
    fidelity_=CheckConvergence(pred_,previous_lambda_);
    double S_E(entropy(CentralDecomposition.Values));
    energy.Real_measurements.push_back(S_E);
    energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    energy.Real_measurements.push_back(fidelity_);
    std::cout <<"1-fidelity: " << fidelity_ << ", Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    //after calculation, insert the two new sites
    return energy;
  }

  Data SuperBlock::move_right_two_vertex(uMPXInt chi, double truncation){
    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();
    std::cout << "Move right" << std::endl;
    //check we can move right
    if (right_size()<1){std::cout << "At right end of superblock already..." << std::endl; exit(1);};
    //are we at the far left end?
    if (left_size()<1){
      shift_right(MakeInitialLeftBlock(load_MPO_matrix("LeftH.MPO_matrix",spectrum),CentralDecomposition.LeftMatrix));
    }
    else {
      shift_right(MakeLeftBlock(getLeftBlock(),getH(),CentralDecomposition.LeftMatrix));
    }
    Data energy;
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;

    previous_lambda_=CentralDecomposition.Values;
    std::stringstream namestream;
    namestream << getName() << "_Right_" << right_size()+1 << ".MPS_matrix";
    pred_=MakeRFinitePrediction(CentralDecomposition,load_MPS_matrix(namestream.str(),basis()));

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy,&(pred_.Guess)),chi,right_size()<=1 ? -0.0 : truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    double S_E(entropy(CentralDecomposition.Values));
    energy.Real_measurements.push_back(S_E);
    energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    energy.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition));

    //energy.Real_measurements.push_back(left_size());
    //energy.Real_measurements.push_back(right_size());
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return energy;
  }

  Data SuperBlock::move_left_two_vertex(uMPXInt chi, double truncation){
    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();

    std::cout << "Move left" << std::endl;
    //check we can move left
    if (left_size()<1){std::cout << "At left end of superblock already..." << std::endl; exit(1);};
    //are we at the far right end?
    if (right_size()<1){
      shift_left(MakeInitialRightBlock(load_MPO_matrix("RightH.MPO_matrix",spectrum),CentralDecomposition.RightMatrix));
    }
    else {
      shift_left(MakeRightBlock(getRightBlock(),getH(),CentralDecomposition.RightMatrix));
    }
    Data energy;
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;

    previous_lambda_=CentralDecomposition.Values;
    std::stringstream namestream;
    namestream << getName() << "_Left_" << left_size()+1 << ".MPS_matrix";
    pred_=MakeLFinitePrediction(CentralDecomposition,load_MPS_matrix(namestream.str(),basis()));

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy,&(pred_.Guess)),chi,left_size()<=1 ? -0.0 : truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    double S_E(entropy(CentralDecomposition.Values));
    energy.Real_measurements.push_back(S_E);
    energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    energy.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition));
    //energy.Real_measurements.push_back(left_size());
    //energy.Real_measurements.push_back(right_size());
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return energy;
  }

  std::pair<MPS_matrix,MPS_matrix> ProjectorBlocks::FetchProjectorStatePair(uMPXInt ls){
    std::stringstream lnstream;
    std::stringstream rnstream;
    uMPXInt rs= size()-ls-middle_size();
    
    if (ls+1>size()/2){ //to right of middle
      lnstream << ProjectorStateName_ << "_Right_"<< rs+middle_size() << ".MPS_matrix"; //B type matrix
    }
    else {
      lnstream << ProjectorStateName_ << "_Left_"<< ls+1 << ".MPS_matrix"; //A type
    }
    if (rs+1>size()/2){
      rnstream << ProjectorStateName_ << "_Left_"<< ls+middle_size() << ".MPS_matrix"; //A type
    }
    else {
      rnstream << ProjectorStateName_ << "_Right_"<< rs+1 << ".MPS_matrix"; //B type
    }
    
    if (ls==rs){
      //on a rightward sweep this is fine, lambda gets included in the left block as we shift right
      std::stringstream cnstream;
      cnstream << ProjectorStateName_ << "_Lambda_"<< ls+1 << "_" << rs+1 << ".MPX_matrix";
      return std::pair<MPS_matrix,MPS_matrix>(std::move(MPS_matrix(contract(load_MPS_matrix(lnstream.str(),getSpectrum()),0,load_MPX_matrix(cnstream.str(),getSpectrum()),0,contract20)).left_shape()),std::move(load_MPS_matrix(rnstream.str(),getSpectrum()).right_shape()));
    }
    else if (ls==size()/2-2){      //on a leftward sweep lambda can be missed when we refetch...
      std::stringstream cnstream;
      cnstream << ProjectorStateName_ << "_Lambda_"<< size()/2 << "_" << size()/2 << ".MPX_matrix";
      return std::pair<MPS_matrix,MPS_matrix>(std::move(load_MPS_matrix(lnstream.str(),getSpectrum()).left_shape()),std::move(MPS_matrix(contract(load_MPS_matrix(rnstream.str(),getSpectrum()),0,load_MPX_matrix(cnstream.str(),getSpectrum()),0,contract20)).right_shape()));
    }
    else {
      return std::pair<MPS_matrix,MPS_matrix>(std::move(load_MPS_matrix(lnstream.str(),getSpectrum()).left_shape()),std::move(load_MPS_matrix(rnstream.str(),getSpectrum()).right_shape()));
    }
  }

  void ProjectorBlocks::ShiftProjectorStatePair(bool dir, uMPXInt new_pos){
    //first swap contents
    std::swap(ProjectorStatePair_.first,ProjectorStatePair_.second);

    std::stringstream nstream;
    bool need_lambda=0;
    
    if (dir==Rightwards_){//moving right
      if (new_pos<=size()/2){ //to right of middle
	nstream << ProjectorStateName_ << "_Right_"<< new_pos << ".MPS_matrix"; //B type matrix
      }
      else {
	nstream << ProjectorStateName_ << "_Left_"<< size()-new_pos+1 << ".MPS_matrix"; //A type
	if (size()-new_pos+1==size()/2) need_lambda=1;
      }
    }
    else {//moving left

      if (new_pos<=size()/2){ //to left of middle
	nstream << ProjectorStateName_ << "_Left_"<< new_pos << ".MPS_matrix"; //A type matrix
	if (new_pos==size()/2) need_lambda=1;
      }
      else {
	nstream << ProjectorStateName_ << "_Right_"<< size()-new_pos+1 << ".MPS_matrix"; //B type
      }
    }

    //when to get lambda too?
    //if we are loading _Left_N/2
    if (need_lambda){
      std::stringstream cnstream;
      cnstream << ProjectorStateName_ << "_Lambda_"<< size()/2 << "_" << size()/2 << ".MPX_matrix";
      if (dir)
	ProjectorStatePair_.second=std::move(MPS_matrix(contract(load_MPS_matrix(nstream.str(),getSpectrum()),0,load_MPX_matrix(cnstream.str(),getSpectrum()),0,contract20)).right_shape());
      else
	ProjectorStatePair_.first=std::move(MPS_matrix(contract(load_MPS_matrix(nstream.str(),getSpectrum()),0,load_MPX_matrix(cnstream.str(),getSpectrum()),0,contract20)).left_shape());
    }
    else { //don't need to contract with lambda
      if (dir)
	ProjectorStatePair_.second=std::move(load_MPS_matrix(nstream.str(),getSpectrum()));
      else
	ProjectorStatePair_.first=std::move(load_MPS_matrix(nstream.str(),getSpectrum()));
    }

    ProjectorStatePair_.first.left_shape();
    ProjectorStatePair_.second.right_shape();
    
  }

  
  TensorWeightPair ProjectorBlocks::makePTensor() const{
    const MPX_matrix& PL=getLeftBlock();
    const MPX_matrix& PR=getRightBlock();
    const MPS_matrix& psi_f=ProjectorStatePair_.first;
    const MPS_matrix& psi_s=ProjectorStatePair_.second;
    //return TensorWeightPair(std::move(reshape_to_vector(contract(contract(contract(psi_f,1,PL,0,contract10),0,psi_s,1,contract10),0,PR,0,contract30)).conjugate()),getWeight());
    return TensorWeightPair(std::move(contract(contract(contract(psi_f,1,PL,0,contract10),0,psi_s,1,contract10),0,PR,0,contract30)),getWeight());
  }

  void ProjectorBlocks::move_right_two_vertex(const MPS_matrix& LeftMatrix) {
    //make new LeftBlock
    MPX_matrix newLeftBlock(contract(ProjectorStatePair_.first,1,contract(getLeftBlock(),0,LeftMatrix,0,contract11),0,contract0110));
    //load in a new projector state pair
    //ProjectorStatePair_=FetchProjectorStatePair(left_size()+1);
    ShiftProjectorStatePair(Rightwards_,right_size());
    
    if (right_size()-1<right_mark_){
      // dummy case
      const MPXIndex& rightindex=ProjectorStatePair_.second.getOutwardMatrixIndex();
      set_blocks(left_size()+1,std::move(newLeftBlock),std::move(MPX_matrix(getSpectrum(),rightindex,std::vector<double>(rightindex.size(),1.0)).Transpose()));
      --right_mark_;
    }
    else{
      //recall from storage
      shift_right(std::move(newLeftBlock));
    }

    //ShiftProjectorStatePair(Rightwards_);
    
  }

  void ProjectorBlocks::move_left_two_vertex(const MPS_matrix& RightMatrix) {
    //make new RightBlock
    MPX_matrix newRightBlock(contract(ProjectorStatePair_.second,1,contract(getRightBlock(),0,RightMatrix,0,contract12),0,contract1220));
    //load in a new projector state pair
    //ProjectorStatePair_=FetchProjectorStatePair(left_size()-1);
    ShiftProjectorStatePair(Leftwards_,left_size());

    //ProjectorStatePair_.first.print_indices();
    //ProjectorStatePair_.second.print_indices();
    if (left_size()-1<left_mark_){
      // dummy case
      const MPXIndex& leftindex=ProjectorStatePair_.first.getInwardMatrixIndex();
      set_blocks(left_size()-1,MPX_matrix(getSpectrum(),leftindex,std::vector<double>(leftindex.size(),1.0)),std::move(newRightBlock));
      --left_mark_;
    }
    else{
      //recall from storage
      shift_left(std::move(newRightBlock));
    }
    //ShiftProjectorStatePair(Leftwards_);
  }

  void iDMRG::run(uMPXInt number_of_steps, double convergence_criterion, uMPXInt chi, double truncation){
    if (number_of_steps==0 && convergence_criterion<=0.0){std::cout << "Need a finite number of steps OR convergence criterion" << std::endl;}
    if (convergence_criterion>=1.0){std::cout << "Need a convergence criterion <1.0" << std::endl;}
    //currently two site only
    double convergence(1.0);
    uMPXInt step_number=0;
    while (step_number<number_of_steps || (convergence>convergence_criterion && convergence_criterion>0.0)){
      Data this_step(grow_two_vertex(chi,truncation));
      //CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
      convergence=this_step.Real_measurements[2];
      ++step_number;
      output_ref_.push(this_step);
    }
    set_2();//if num vertices is only 2, need to adjust so that result is ok as input for other methods.
  }

  void FiniteDMRG::run(uMPXInt num_sweeps, uMPXInt chi, double truncation){
    chi_=chi;
    truncation_=truncation;
    std::cout << "Performing finite sweeps" << std::endl;
    //we start at the middle of the system
    //in case we need excited states later, we store MPS matrices as well as the L and R blocks.
    if (size()==2) {
      //necessary to leave a well formed state to hand to excited states or TEBD
      //would be better to just store lambda, as left and right will have been stored by growth stage
      CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());//store left and right
      std::cout << "Only two vertices. Skipping sweeps and storing lambda." << std::endl;
    }
    else {
      for (uMPXInt n=num_sweeps;n>0;--n){
	double cumulative_truncation=0.0;
	std::cout << std::endl << "Starting sweep: " << num_sweeps-n+1<< std::endl;
	std::cout << left_size() << " " << middle_size() << " " << right_size() << std::endl;
	for (uMPXInt r=right_size();r>0;--r){
	  Data this_step(move_right_two_vertex(chi_,truncation_));
	  //output_ref_.push(this_step);//at midpoint, push output

	  cumulative_truncation+=CentralDecomposition.Truncation;
	}

	for (uMPXInt l=left_size();l>0;--l){
	  Data this_step(move_left_two_vertex(chi_,truncation_));
	  //output_ref_.push(this_step);//at midpoint, push output

	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);//at midpoint, push output
	  }
	  cumulative_truncation+=CentralDecomposition.Truncation;
	}
	//CentralDecomposition.store_left(getName(),left_size()+1);
	for (uMPXInt r=right_size();r>left_size();--r){
	  Data this_step(move_right_two_vertex(chi_,truncation_));
	  //output_ref_.push(this_step);//at midpoint, push output

	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);//at midpoint, push output
	  }
	  cumulative_truncation+=CentralDecomposition.Truncation;
	}
	cumulative_truncation=cumulative_truncation/(2.0*size());
	std::cout << "Average (per vertex) truncation for sweep = " << cumulative_truncation <<std::endl;
      }
    }
  }

  void ExcitedStateFiniteDMRG::init(double chi, double truncation, bool converge){
    //first step
    //LeftBlock and RightBlock are correct? but we need to solve centre and store.
    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();
    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Calculating Excited State " << PBlocks_.size() << std::endl;

    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;

    //not the best guess, as we want something orthogonal to this...
    SparseMatrix Guess(reshape_to_vector(contract(CentralDecomposition.LeftMatrix,0,contract(MPX_matrix(CentralDecomposition.basis(),CentralDecomposition.RightMatrix.Index(0),CentralDecomposition.Values),0,CentralDecomposition.RightMatrix,0,contract10),0,contract20)));

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results,nullptr/*&Guess*/,converge),chi,truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());

    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(1.0); //overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    output_ref_.push(results);
  }

  Data ExcitedStateFiniteDMRG::move_right_two_vertex(uMPXInt chi, double truncation, bool converge){
    //if this is the initial run through, or we need a dummy block only, we get Hamiltonian blocks from GS storage
    std::string StorageName;
    if (right_size()==1){//just need a dummy block
      StorageName=GSName_;
    }
    else if (init_flag_){
      StorageName=HBlocksName_;
    }
    else {
      StorageName=getName();
    }

    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();
    //check we can move right
    if (right_size()<1){std::cout << "At right end of superblock already..." << std::endl; exit(1);};
    //are we at the far left end?
    if (left_size()<1){
      shift_right(MakeInitialLeftBlock(load_MPO_matrix("LeftH.MPO_matrix",spectrum),CentralDecomposition.LeftMatrix),StorageName);
    }
    else {
      shift_right(MakeLeftBlock(getLeftBlock(),getH(),CentralDecomposition.LeftMatrix),StorageName);
    }
    //shift projectors
    for (auto&& pb : PBlocks_){
      pb.move_right_two_vertex(CentralDecomposition.LeftMatrix);
    }

    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;
    previous_lambda_=CentralDecomposition.Values;
    
    std::stringstream namestream;
    namestream << StorageName << "_Right_" << right_size()+1 << ".MPS_matrix";
    pred_=MakeRFinitePrediction(CentralDecomposition,load_MPS_matrix(namestream.str(),basis()));

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results,&(pred_.Guess),converge),chi,right_size()<=1 ? -0.0 : truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());

    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition));
    //results.Real_measurements.push_back(left_size());
    //results.Real_measurements.push_back(right_size());

    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return results;
  }

  Data ExcitedStateFiniteDMRG::move_left_two_vertex(uMPXInt chi, double truncation, bool converge){
    //if this is the initial run through, or we need a dummy block only, we get blocks from GS storage
    //first pass to left, we will have the left blocks for the right hand side only
    std::string StorageName;
    if (left_size()==1){//just need a dummy block, and ground state provides one.
      StorageName=GSName_;
    }
    else if ((left_size()<=right_size()+middle_size()) && init_flag_){//if in init phase we have not formed all the new blocks yet
      StorageName=HBlocksName_;
    }
    else {
      StorageName=getName();
    }

    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();
    //check we can move left
    if (left_size()<1){std::cout << "At left end of superblock already..." << std::endl; exit(1);};
    //are we at the far right end?
    if (right_size()<1){
      shift_left(StorageName,MakeInitialRightBlock(load_MPO_matrix("RightH.MPO_matrix",spectrum),CentralDecomposition.RightMatrix));
    }
    else {
      shift_left(StorageName,MakeRightBlock(getRightBlock(),getH(),CentralDecomposition.RightMatrix));
    }

    for (auto&& pb : PBlocks_){
      pb.move_left_two_vertex(CentralDecomposition.RightMatrix);
    }

    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;
    previous_lambda_=CentralDecomposition.Values;
    //mismatch between storage and prediction names at midpoint, while in init phase,
    //reaching the mid point from the right requires block with HBlocksName, but we have an MPS_matrix with the new name that we must use!
    std::stringstream namestream;
    std::string PredictionName = ((left_size()==right_size()) && init_flag_) ? getName() : StorageName;
    namestream << PredictionName << "_Left_" << left_size()+1 << ".MPS_matrix";
    pred_=MakeLFinitePrediction(CentralDecomposition,load_MPS_matrix(namestream.str(),basis()));

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results,&(pred_.Guess),converge),chi,left_size()<=1 ? -0.0 : truncation);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1,left_size()==right_size());

    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition));
    //results.Real_measurements.push_back(left_size());
    //results.Real_measurements.push_back(right_size());
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return results;
  }

  void ExcitedStateFiniteDMRG::run(uMPXInt num_sweeps, uMPXInt chi, double truncation){
    
    if (init_flag_) init(chi, truncation, 1);

    //    if (init_flag_) init(chi,(num_sweeps>1 && size()!=2) ? 0.0 : truncation, (num_sweeps>1 && size()!=2) ? 0 : 1);
    std::cout << "Performing finite sweeps" << std::endl;
    //we start at the middle of the system
    //in case we need excited states later, we store MPS matrices as well as the L and R blocks.
    if (size()==2) {
      //CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
      std::cout << "Skipping finite sweeps, only two vertices..." << std::endl;
    }
    else {
      bool converge=num_sweeps>1 ? 0 : 1;

      for (uMPXInt n=num_sweeps;n>0;--n){//sweep towards right
	uMPXInt chi_local = chi;
	
	double cumulative_truncation=0.0;
	std::cout << std::endl << "Starting sweep: " << num_sweeps-n+1<< std::endl;
	for (uMPXInt r=right_size();r>0;--r){
	  Data this_step(move_right_two_vertex(chi_local,truncation,converge));

	  cumulative_truncation+=CentralDecomposition.Truncation;
	}
	for (uMPXInt l=left_size();l>0;--l){
	  Data this_step(move_left_two_vertex(chi_local,truncation,converge));
	    
	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);
	  }
	  cumulative_truncation+=CentralDecomposition.Truncation;
	}
	
	if (n==num_sweeps) init_flag_=0; //first time through, turn off init here as we have formed all blocks once...
	
	for (uMPXInt r=right_size();r>left_size();--r){
	  Data this_step(move_right_two_vertex(chi_local,truncation,converge));

	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);
	  }
	  cumulative_truncation+=CentralDecomposition.Truncation;
	}
	cumulative_truncation=cumulative_truncation/(2.0*size()); //factor of two for sweep back and forth
	std::cout << "Average (per vertex) truncation for sweep = " << cumulative_truncation <<std::endl;
      }
    }
  }

  MPX_matrix MakeInitialLeftBlock(const MPO_matrix& LeftH, const MPS_matrix& A){
    std::cout << "Forming initial Left Block" << std::endl;

    return contract(A,1,contract(LeftH,0,A,0,std::vector<MPXPair>(1,MPXPair(2,0))),0,std::vector<MPXPair>(1,MPXPair(0,0)));
  }
  MPX_matrix MakeInitialRightBlock(const MPO_matrix& RightH, const MPS_matrix& B){
    std::cout << "Forming initial Right Block" << std::endl;

    return contract(B,1,contract(RightH,0,B,0,std::vector<MPXPair>(1,MPXPair(2,1))),0,std::vector<MPXPair>(1,MPXPair(1,0)));
  }

  MPX_matrix MakeLeftBlock(const MPX_matrix& LB, const MPX_matrix& H, const MPS_matrix& A){
    /*static MPXInt dints[]={2,0,3,1,4,5};
      std::vector<MPXInt> reorder203145(dints,dints+6);*/
    std::vector<MPXPair> secondcontraction(1,MPXPair(1,3));
    secondcontraction.push_back(MPXPair(2,5));
    std::vector<MPXPair> thirdcontraction(1,MPXPair(0,0));
    thirdcontraction.push_back(MPXPair(1,3));
    return std::move(contract(A,1,contract(H,0,contract(LB,0,A,0,MPXPair(5,1)),0,secondcontraction),0,thirdcontraction).MoveDummyIndices(reorder203145)); //203145
  }
  MPX_matrix MakeRightBlock(const MPX_matrix& RB, const MPX_matrix& H, const MPS_matrix& B){
    /*static MPXInt dints[]={0,2,1,3,5,4};
      std::vector<MPXInt> reorder021354(dints,dints+6);*/
    std::vector<MPXPair> secondcontraction(1,MPXPair(2,6));
    secondcontraction.push_back(MPXPair(3,2));
    std::vector<MPXPair> thirdcontraction(1,MPXPair(1,0));
    thirdcontraction.push_back(MPXPair(2,2));
    return std::move(contract(B,1,contract(H,0,contract(RB,0,B,0,MPXPair(4,2)),0,secondcontraction),0,thirdcontraction).MoveDummyIndices(reorder021354)); //021354
  }

  MPX_matrix MakeDummyLeftBlock(const MPO_matrix& H, const State& TargetState){
    //need six mpxindices, all of size 1
    //outgoing dummy
    //ingoing dummy
    //ingoing dummy
    //outgoing , selects last matrix row of H MPO
    //ingoing dummy
    //outgoing dummy
    std::vector<MPXIndex> indices;
    indices.emplace_back(0,StateArray(1,TargetState.Identity()));
    indices.emplace_back(1,StateArray(1,TargetState.Identity()));
    indices.emplace_back(1,StateArray(1,TargetState.Identity()));
    indices.emplace_back(0,H.Index(1));
    indices.emplace_back(1,StateArray(1,TargetState.Identity()));
    indices.emplace_back(0,StateArray(1,TargetState.Identity()));
    SparseMatrix array(H.Index(1).size(),1,1);
    array.cheap_entry(H.Index(1).size()-1,0,1.0);
    array.cheap_finalise();
    return MPX_matrix(H.GetPhysicalSpectrum(),indices,4,array);
  }

  MPX_matrix MakeDummyRightBlock(const MPO_matrix& H, const State& TargetState){
    //need six mpxindices, all of size 1
    //outgoing targetstate
    //ingoing dummy
    //ingoing , selects first matrix column of H MPO
    //outgoing dummy
    //ingoing target state
    //outgoing dummy
    std::vector<MPXIndex> indices;
    indices.emplace_back(0,StateArray(1,TargetState));
    indices.emplace_back(1,StateArray(1,TargetState.Identity()));
    indices.emplace_back(1,H.Index(3));
    indices.emplace_back(0,StateArray(1,TargetState.Identity()));
    indices.emplace_back(1,StateArray(1,TargetState));
    indices.emplace_back(0,StateArray(1,TargetState.Identity()));
    SparseMatrix array(H.Index(3).size(),1,1);
    array.cheap_entry(0,0,1.0);
    array.cheap_finalise();
    return MPX_matrix(H.GetPhysicalSpectrum(),indices,3,array);
  }

  MPX_matrix TwoVertexInitialHamiltonian(const MPO_matrix& LeftH, const MPO_matrix& RightH){
    MPX_matrix ans(reorder(contract(LeftH,0,RightH,0,std::vector<MPXPair>(1,MPXPair(3,1))),0,reorder032415,2));
    return ans;
  }

  MPX_matrix TwoVertexInitialWavefunction(const MPO_matrix& LeftH, const MPO_matrix& RightH, const State& TargetSector, Data& result){
    const MPX_matrix H2(TwoVertexInitialHamiltonian(LeftH,RightH));

#ifndef DNDEBUG
    if (H2.basis().size()<=5){
      std::cout <<"*****************************" <<std::endl;
      std::cout <<"TWO VERTEX HAMILTONIAN OUTPUT" <<std::endl;
      std::cout <<"*****************************" <<std::endl;
      H2.print_matrix();
    }
#endif

    static char SMALLESTREAL[]={'S','R','\n'}; //lowest real part for energies

    std::cout << "Eigensolver starting for two vertex wavefunction..." << std::endl;    

    SparseHED decomp(H2.Eigs(TargetSector,1,SMALLESTREAL));
    std::cout << "Lowest energy/Number of Vertices: " << decomp.Values[0]/2.0 <<std::endl;

    std::cout << "Raw eigenvalues" <<std::endl;
    decomp.printValues();

    result.Real_measurements.push_back(decomp.Values[0]);
    result.Real_measurements.push_back(decomp.Values[0]/2.0);
    std::vector<MPXIndex> wfindices;
    wfindices.push_back(MPXIndex(1,H2.GetPhysicalSpectrum())); //ingoing
    wfindices.push_back(MPXIndex(1,StateArray(1,TargetSector-TargetSector))); //ingoing
    wfindices.push_back(MPXIndex(1,H2.GetPhysicalSpectrum())); //ingoing
    wfindices.push_back(MPXIndex(0,StateArray(1,TargetSector))); //outgoing

    return MPX_matrix(H2.GetPhysicalSpectrum(),wfindices,2,reshape(decomp.EigenVectors.ExtractColumns(std::vector<MPXInt>(1,0)),H2.GetPhysicalSpectrum().size()));
  }

  Prediction MakePrediction(const MPSDecomposition& Decomp, const std::vector<double>& PreviousLambda){
    std::cout <<"Making prediction vector" <<std::endl;
    const Basis& spectrum(Decomp.LeftMatrix.GetPhysicalSpectrum());
    Prediction ans;
    MPX_matrix svals(spectrum,Decomp.LeftMatrix.Index(2),Decomp.Values);
    MPX_matrix LambdaLB(reorder(contract(Decomp.LeftMatrix,0,svals,0,contract20),0,reorder102,1));
    //MPXDecomposition RotDecomp(LambdaLB.SVD());
    //ans.LambdaL=contract_to_sparse(RotDecomp.ColumnMatrix,0,MPX_matrix(spectrum,RotDecomp.ColumnMatrix.Index(1),RotDecomp.Values),0,contract10);
    MPX_matrix ALambdaR(reorder(contract(svals,0,Decomp.RightMatrix,0,contract10),0,reorder102,2));
    MPXDecomposition RotDecomp=ALambdaR.SVD();
    ans.Auxiliary=contract_to_sparse(MPX_matrix(spectrum,RotDecomp.RowMatrix.Index(0),RotDecomp.Values),0,RotDecomp.RowMatrix,0,contract10);
    MPX_matrix InversePreviousLambda(spectrum,Decomp.RightMatrix.Index(2),Decomp.LeftMatrix.Index(1),PreviousLambda,1); //1 means take inverse values    
    ans.Guess=reshape_to_vector(contract(contract(ALambdaR,0,InversePreviousLambda,0,contract20),0,LambdaLB,0,contract20));
    return ans;
  }

  Prediction MakeLFinitePrediction(const MPSDecomposition& Decomp, const MPS_matrix& A){
    std::cout <<"Making prediction vector" <<std::endl;
    return Prediction(reshape_to_vector(contract(A,0,contract(Decomp.LeftMatrix,0,MPX_matrix(Decomp.LeftMatrix.basis(),Decomp.LeftMatrix.Index(2),Decomp.Values),0,contract20),0,contract21)));
  }
  Prediction MakeRFinitePrediction(const MPSDecomposition& Decomp, const MPS_matrix& B){
    std::cout <<"Making prediction vector" <<std::endl;
    return Prediction(reshape_to_vector(contract(reorder(contract(MPX_matrix(Decomp.RightMatrix.basis(),Decomp.LeftMatrix.Index(2),Decomp.Values),0,Decomp.RightMatrix,0,contract10),0,reorder102,2),0,B,0,contract20)));
  }

  MPX_matrix TwoVertexWavefunction(const MPX_matrix& LB, const MPO_matrix& H, const MPX_matrix& RB, const std::vector<ProjectorBlocks>* ProjectorBlocksPtr, MPXInt NumVertices, Data& result, SparseMatrix* guessptr, bool converge){
    //Need to figure out 'length' then choose strategy appropriately.
    //form LeftPart, RightPart
#ifdef TIMING
    num_op_x=0; //reset
    auto t1 = std::chrono::high_resolution_clock::now();
#endif

    TwoVertexComponents stuff(LB,H,RB,ProjectorBlocksPtr);

#ifdef TIMING
    auto t2 = std::chrono::high_resolution_clock::now();
    formation_time_millisecs+=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    formations++;
#endif

    MPXInt num_to_find = 1;

#ifdef TIMING
    auto tD1 = std::chrono::high_resolution_clock::now();
#endif

    SparseHED decomp(stuff.HED(num_to_find,SMALLESTREAL,guessptr,converge));
#ifdef TIMING
    auto tD2 = std::chrono::high_resolution_clock::now();
    decomp_time_millisecs+=std::chrono::duration_cast<std::chrono::milliseconds>(tD2-tD1).count();
    decomps++;
#endif
    if (!ProjectorBlocksPtr){ //are we looking for excited states?
      std::cout << "Lowest energy/Number of Vertices: " << decomp.Values[0]/double(NumVertices) <<std::endl;
    }
    std::cout << "Lowest eigenvalue: " << decomp.Values[0] <<std::endl;
    if (decomp.ValuesSize()>1) {
      std::cout << "Next lowest eigenvalues: ";
      for (MPXInt n=1;n<num_to_find;++n) {
	std::cout << decomp.Values[n] << " ";
      }
      std::cout << std::endl;
    }
    if (ProjectorBlocksPtr){ //looking for excited states, then don't divide
      result.Real_measurements.push_back(decomp.Values[0]);
    }
    else {
      result.Real_measurements.push_back(decomp.Values[0]);
      result.Real_measurements.push_back(decomp.Values[0]/double(NumVertices));
    }
    //need to reshape from a vector
    //(sigma a), (sigma a) matrix form
    return MPX_matrix(LB.GetPhysicalSpectrum(),stuff.indices,2,reshape(decomp.EigenVectors.ExtractColumns(std::vector<MPXInt>(1,0)),LB.GetPhysicalSpectrum().size()*LB.Index(5).size()));
  }

  double CheckConvergence(const Prediction& guess,const std::vector<double>& PreviousLambda){
    return 1.0-Sum((guess.Auxiliary*SparseMatrix(PreviousLambda)).SVD());
  }

  double CheckConvergence(const Prediction& g,const MPSDecomposition& D){
    //g.Guess is just a sparse matrix at this point...
    return 1.0-abs((reshape(g.Guess,D.basis().size()*D.LeftMatrix.Index(1).size()).dagger()*contract_to_sparse(contract(D.LeftMatrix,0,MPX_matrix(D.basis(),D.RightMatrix.Index(0),D.Values),0,contract20),0,D.RightMatrix,0,contract20)).trace());
  }


  TwoVertexComponents::TwoVertexComponents(const MPX_matrix& L, const MPO_matrix& HMPO, const MPX_matrix& R, const std::vector<ProjectorBlocks>* P, const State* StatePtr) : LeftBlock(L),H(HMPO),RightBlock(R),ProjectorsPtr(P),TargetStatePtr(StatePtr),m_length(L.Index(5).size()*H.Index(2).size()*R.Index(4).size()*H.Index(2).size()),LeftPart(reorder(contract(H,0,LeftBlock,0,contract13).RemoveDummyIndices(std::vector<MPXInt>({{3,5,6}})),0,reorder03214,3)),RightPart(reorder(contract(H,0,RightBlock,0,contract32).RemoveDummyIndices(std::vector<MPXInt>({{4,5,7}})),0,reorder12403,3)) {
    
#ifndef NDEBUG
    std::cout << "H MPO" << std::endl;
    HMPO.print_indices();
    std::cout << "Left and right blocks" << std::endl;
    LeftBlock.print_indices();
    RightBlock.print_indices();

    std::cout << "Forming left and right parts" << std::endl;
    LeftPart.print_indices();
    RightPart.print_indices();
    LeftPart.print_sparse_info();
    RightPart.print_sparse_info();
#endif

    indices.emplace_back(1,HMPO.Index(2));
    indices.emplace_back(1,LeftBlock.Index(5));
    indices.emplace_back(1,HMPO.Index(2));
    indices.emplace_back(0,RightBlock.Index(4));

    const MPXIndex& sigma1(indices[0]);
    const MPXIndex& a0(indices[1]);
    const MPXIndex& sigma2(indices[2]); 
    const MPXIndex& a2(indices[3]);

    vrows=sigma1.size()*a0.size();
    vcols=sigma2.size()*a2.size();

    for (uMPXInt c_a2=0;c_a2<a2.size();++c_a2){

      for (MPXInt c_s2=0;c_s2<sigma2.size();++c_s2){
	State loop2acc(a2[c_a2]-sigma2[c_s2]);
	MPXInt colpart(c_s2+sigma2.size()*c_a2);
	for (MPXInt r_a0=0;r_a0<a0.size();++r_a0){
	  State loop3acc(loop2acc-a0[r_a0]);
	  for (MPXInt r_s1=0;r_s1<sigma1.size();++r_s1){
	    if (sigma1[r_s1]==loop3acc){ //check for consistency of charges
	      
	      allowed_indices.push_back(r_s1+sigma1.size()*r_a0+sigma1.size()*a0.size()*colpart);
	      rows_and_cols.push_back(std::array<MPXInt,2> {{r_s1+sigma1.size()*r_a0,colpart}});
	    }
	  }
	}
      }
    }

    if (ProjectorsPtr){
      const std::vector<ProjectorBlocks>& PBvec(*ProjectorsPtr);
      for (auto&& pb : PBvec){
	std::cout << "Making PTensor" <<std::endl;
	ProjectorTensors.emplace_back(pb.makePTensor());	
      }
    }
  };

  SparseHED TwoVertexComponents::HED(MPXInt numevals, char which[3],const SparseMatrix* initial,bool converge) const {
    //set up workspace (Evals, Evecs) SparseHED
    uMPXInt fulldim(this->length());
    std::cout <<"Eigensolver for matrix of length " << fulldim << std::endl;
    std::cout <<"Using reduced subspace of length " << allowed_indices.size() <<std::endl;
    //need to guard against numevals being 0 or larger than matrix size
    numevals=(numevals/*+ProjectorTensors.size()*/ > allowed_indices.size() ? allowed_indices.size()/*ProjectorTensors.size()*/ : numevals);
    if (numevals==0){
      //std::cout << "All excited states in sector have been found!" <<std::endl;
      std::cout << "Number of requested eigenvalues=0!" <<std::endl;
      return SparseHED(fulldim,0);
    }
    SparseVectorWithRestriction guess_struct(initial,&allowed_indices);

    //internal contractions may be large.
    //however if allowed_indices is 'small' we should use a dense method

    if (allowed_indices.size()<=SPARSE_THRESHOLD && /*allowed_indices.size()>*/!ProjectorTensors.size()){
      //just make the (dense) matrix and use lapack
      static std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> > condition={{1,7},{2,6}};
      TranslationBlock<DenseMatrix> TB(contract_conditional<DenseMatrix>(contract(H,0,LeftBlock,0,contract13),0,contract(H,0,RightBlock,0,contract32),0,contract21,condition));
      //use lapack
      DenseHED dense_ans(TB.Block.HED(1,ProjectorTensors.size()));
      //DenseHED dense_ans(TB.Block.HED(1,0)); //find smallest real eval
      return SparseHED(std::vector<double>(dense_ans.Values,dense_ans.Values+dense_ans.ValuesSize()),TB.TranslateRows(dense_ans.EigenVectors));
    }
    else {
      std::cout <<"Allocate storage for evals and evecs" << std::endl;
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      SparseHED ans(fulldim,numevals);
      std::cout <<"Calling arpack_eigs()..." << std::endl;
      arpack::arpack_eigs<TwoVertexComponents,SparseVectorWithRestriction> eigensystem(this,&TwoVertexMPOMPSMultiply,allowed_indices.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs,converge);
      if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}

      for (size_t v=0;v<static_cast<size_t>(numevals);++v){
	std::cout << "Eval: " << Evals[v] <<std::endl;
	ans.Values.push_back(Evals[v].real()); //arpack doesn't assume a Hermitian form
	double maxabsvalue=0.0;
	std::complex<double> maxvalue=0.0;
	//first pass find largest abs value
	for (MPXInt i=0;i<allowed_indices.size();++i){
	  if (abs(Evecs[i+v*allowed_indices.size()])>maxabsvalue){
	    maxabsvalue = abs(maxvalue=Evecs[i+v*allowed_indices.size()]);
	  }
	}
	maxvalue=maxabsvalue/maxvalue;
	//now get rid of possible phase
	for(MPXInt i=0;i<allowed_indices.size();++i){
	  std::complex<double> value=Evecs[i+v*allowed_indices.size()]*maxvalue;
	  if (abs(value)>SPARSETOL*maxabsvalue){
	    //convert back to real rows
	    ans.EigenVectors.entry(allowed_indices[i],v,FixComplexPrecision(value));
	  }
	}
      }
      delete[] Evecs;
      delete[] Evals;
      ans.EigenVectors.finalise();
      return ans;
    }
  }

  void TwoVertexMPOMPSMultiply(const TwoVertexComponents* arraystuff, std::complex<double> *in, std::complex<double> *out){
    SparseMatrix V(arraystuff->vcols,arraystuff->vrows,arraystuff->vrows);//undocumented behaviour of SparseMatrix, using transposed rows and cols to avoid unnecessary row ordering...
#ifdef TIMING
    auto vt1 = std::chrono::high_resolution_clock::now();
#endif

    for (struct {std::vector<std::array<MPXInt,2> >::const_iterator cit; MPXInt idx;} itstruct ={arraystuff->rows_and_cols.begin(), 0};itstruct.cit!=arraystuff->rows_and_cols.end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    MPX_matrix Vector(arraystuff->LeftPart.GetPhysicalSpectrum(),arraystuff->indices,2,V.cheap_no_transpose_finalise());

#ifdef TIMING
    auto vt2 = std::chrono::high_resolution_clock::now();
    vec_time_millisecs+=std::chrono::duration_cast<std::chrono::microseconds>(vt2-vt1).count();
    vec_formations++;

    auto c1t1 = std::chrono::high_resolution_clock::now();
    //MPX_matrix M(contract(arraystuff->LeftPart,0,Vector,0,contract1071));
    MPX_matrix M(std::move(contract(arraystuff->LeftPart,0,Vector,0,contract3041).ShiftNumRowIndices(2)));
    auto c1t2 = std::chrono::high_resolution_clock::now();
    c1_time_ms+=std::chrono::duration_cast<std::chrono::milliseconds>(c1t2-c1t1).count();
    c1_reps++;

    auto c2t1 = std::chrono::high_resolution_clock::now();
    //SparseMatrix SM(contract_to_sparse(M,0,arraystuff->RightPart,0,contract116276));
    SparseMatrix SM(contract_to_sparse(M,0,arraystuff->RightPart,0,contract203142));
    auto c2t2 = std::chrono::high_resolution_clock::now();
    c2_time_ms+=std::chrono::duration_cast<std::chrono::milliseconds>(c2t2-c2t1).count();
    c2_reps++;

    num_op_x++;

    auto c3t1 = std::chrono::high_resolution_clock::now();
    SM=std::move(reshape(SM,arraystuff->length()));
    auto c3t2 = std::chrono::high_resolution_clock::now();
    c3_time_ms+=std::chrono::duration_cast<std::chrono::microseconds>(c3t2-c3t1).count();
    c3_reps++;

    auto c4t1 = std::chrono::high_resolution_clock::now();
    DumbExtractWithZeros(SM,arraystuff->allowed_indices,out);
    auto c4t2 = std::chrono::high_resolution_clock::now();
    c4_time_ms+=std::chrono::duration_cast<std::chrono::microseconds>(c4t2-c4t1).count();
    c4_reps++;
#else
 
    DumbExtractWithZeros(reshape(contract_to_sparse(std::move(contract(arraystuff->LeftPart,0,Vector,0,contract3041).ShiftNumRowIndices(2)),0,arraystuff->RightPart,0,contract203142),arraystuff->length()),arraystuff->allowed_indices,out);

#endif

    //projector bits
    const std::vector<TensorWeightPair>& PTensors(arraystuff->ProjectorTensors);
    //loop through list
    for (auto&& PT : PTensors){
      //contract_to_sparse with in vector and multiply by weight
      //std::cout << "Projector vector contraction" << std::endl;
      //PT.first.print_indices();
      //Vector.print_indices();

      SparseMatrix s(contract_to_sparse(PT.first,0,Vector,0,contract00112233));//really should be a scalar
      //std:: cout << "s: " << s.get_x(0) <<std::endl;
      //std::cout << std::endl << std::endl;
      //s.print();
      //std::cout << std::endl << std::endl;

      if (s.rows()!=1 || s.cols()!=1){
	std::cout << "Error contracting projector down to scalar: " << s.rows() << " " << s.cols() << std::endl;
      }
      else if (s.nz()!=0){ //otherwise memory errors...
	//std::complex<double> w=PT.second*(contract_to_sparse(PT.first,0,Vector,0,contract00112233).get_x(0));
	std::complex<double> w=PT.second*s.get_x(0); 
	//std::cout << "Done Projector vector contraction" <<std::endl;
	//update out vector
	DumbExtractUpdate(std::move(reshape_to_vector(PT.first).conjugate().rescale(w)),arraystuff->allowed_indices,out);
	//std::cout << "Done conjugation and rescale contraction" <<std::endl;
      }
    }
  }

}

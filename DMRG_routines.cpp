#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include "ajaj_common.hpp" //MPXInt MPXPair
#include "sparse_interface.hpp" //SparseHED
#include "states.hpp" //EigenStateArray etc.
#include "MPX.hpp" //MPX_matrix etc.
#include "DMRG_routines.hpp" //MPX_matrix etc.
#include "data.hpp"

namespace ajaj {

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
  void BlocksStructure::load_left_block(){
    std::stringstream leftname;
    leftname << Name_ << "_Left_" << left_size() << ".BLOCK";
    LeftBlock=std::move(load_MPX_matrix(leftname.str(),*SpectrumPtr_));
  }
  void BlocksStructure::load_right_block(){
    std::stringstream rightname;
    rightname << Name_ << "_Right_" << right_size() << ".BLOCK";
    RightBlock=std::move(load_MPX_matrix(rightname.str(),*SpectrumPtr_));
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

  void SuperBlock::push_density() const {
    std::ofstream DensityFileStream;
    DensityFileStream.open(DensityFileName_.c_str(),ios::out);
    CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream);
    DensityFileStream.close();
  }

  Data SuperBlock::initialise(uMPXInt chi, double smin){
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
    //chi and smin zero? take as an indication that we should do a full diag...
    uMPXInt bonddim = (chi<1 && smin <=0.0) ? getH().GetPhysicalSpectrum().size() : chi;
    Data two_vertex_energy;
    //solve
    CentralDecomposition=TwoVertexSVD(TwoVertexInitialWavefunction(LeftH,RightH,TargetState_,two_vertex_energy),bonddim,smin);
    CentralDecomposition.SquareRescale(1.0);
    push_density();
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    //at this stage we can use the generated left and right Hamiltonians to save the left and right blocks
    //update superblock values
    initial_2(MakeInitialLeftBlock(LeftH,CentralDecomposition.LeftMatrix),MakeInitialRightBlock(RightH,CentralDecomposition.RightMatrix));
    //at this point there are 2 vertices, and we have pre stored the blocks they contribute to
    previous_lambda_=std::vector<double>(1,1.0);
    pred_=MakePrediction(CentralDecomposition,previous_lambda_);//give use one way of checking overlap of current state with the previous one.
    fidelity_=CheckConvergence(pred_,previous_lambda_);
    double S_E(entropy(CentralDecomposition.Values));
    std::cout << "Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    two_vertex_energy.Real_measurements.push_back(S_E);
    two_vertex_energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    two_vertex_energy.Real_measurements.push_back(fidelity_);
    return two_vertex_energy;
  }

  Data SuperBlock::grow_two_vertex(uMPXInt chi, double smin){
    if (size()<2){ //haven't initialised yet?
      return initialise(chi,smin); //store dummies, and first left and right blocks
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
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy,&(pred_.Guess)),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    push_density();
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

  Data SuperBlock::move_right_two_vertex(uMPXInt chi, double smin){
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
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    push_density();
    //pred_=MakePrediction(CentralDecomposition,previous_lambda_);
    //fidelity_=CheckConvergence(pred_,previous_lambda_);
    double S_E(entropy(CentralDecomposition.Values));
    energy.Real_measurements.push_back(S_E);
    energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    energy.Real_measurements.push_back(1.0);//overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return energy;
  }

  Data SuperBlock::move_left_two_vertex(uMPXInt chi, double smin){
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
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),nullptr,size(),energy),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    //CentralDecomposition.OutputPhysicalIndexDensities(DensityFileStream_);
    push_density();
    //pred_=MakePrediction(CentralDecomposition,previous_lambda_);
    //fidelity_=CheckConvergence(pred_,previous_lambda_);
    double S_E(entropy(CentralDecomposition.Values));
    energy.Real_measurements.push_back(S_E);
    energy.Real_measurements.push_back(CentralDecomposition.Truncation);
    energy.Real_measurements.push_back(1.0);//overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return energy;
  }

  std::pair<MPS_matrix,MPS_matrix> ProjectorBlocks::FetchProjectorStatePair(uMPXInt ls){
    std::stringstream lnstream;
    std::stringstream rnstream;
    uMPXInt rs= size()-ls-middle_size();
    if (ls+1>size()/2){
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
      std::stringstream cnstream;
      cnstream << ProjectorStateName_ << "_Lambda_"<< ls+1 << "_" << rs+1 << ".MPX_matrix";
      return std::pair<MPS_matrix,MPS_matrix>(std::move(MPS_matrix(contract(load_MPS_matrix(lnstream.str(),getSpectrum()),0,load_MPX_matrix(cnstream.str(),getSpectrum()),0,contract20)).left_shape()),std::move(load_MPS_matrix(rnstream.str(),getSpectrum()).right_shape()));
    }
    else {
      return std::pair<MPS_matrix,MPS_matrix>(std::move(load_MPS_matrix(lnstream.str(),getSpectrum()).left_shape()),std::move(load_MPS_matrix(rnstream.str(),getSpectrum()).right_shape()));
    }
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
    ProjectorStatePair_=FetchProjectorStatePair(left_size()+1);
    //ProjectorStatePair_.first.print_indices();
    //ProjectorStatePair_.second.print_indices();
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
  }

  void ProjectorBlocks::move_left_two_vertex(const MPS_matrix& RightMatrix) {
    //make new RightBlock
    MPX_matrix newRightBlock(contract(ProjectorStatePair_.second,1,contract(getRightBlock(),0,RightMatrix,0,contract12),0,contract1220));
    //load in a new projector state pair
    ProjectorStatePair_=FetchProjectorStatePair(left_size()-1);
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
  }

  void iDMRG::run(uMPXInt number_of_steps, double convergence_criterion, uMPXInt chi, double smin){
    if (number_of_steps==0 && convergence_criterion<=0.0){std::cout << "Need a finite number of steps OR convergence criterion" << std::endl;}
    if (convergence_criterion>=1.0){std::cout << "Need a convergence criterion <1.0" << std::endl;}
    //currently two site only
    double convergence(1.0);
    uMPXInt step_number=0;
    while (step_number<number_of_steps || (convergence>convergence_criterion && convergence_criterion>0.0)){
      Data this_step(grow_two_vertex(chi,smin));
      convergence=this_step.Real_measurements[2];
      ++step_number;
      output_ref_.push(this_step);
    }
    set_2();//if num vertices is only 2, need to adjust so that result is ok as input for other methods.
  }

  void FiniteDMRG::run(uMPXInt num_sweeps, uMPXInt chi, double smin){
    chi_=chi;
    smin_=smin;
    std::cout << "Performing finite sweeps" << std::endl;
    //we start at the middle of the system
    //in case we need excited states later, we store MPS matrices as well as the L and R blocks.
    if (size()==2) {
      CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
      std::cout << "Skipping finite sweeps, only two vertices..." << std::endl;
    }
    else {
      for (uMPXInt n=num_sweeps;n>0;--n){
	std::cout << std::endl << "Starting sweep: " << num_sweeps-n+1<< std::endl;
	std::cout << left_size() << " " << middle_size() << " " << right_size() << std::endl;
	for (uMPXInt r=right_size();r>0;--r){
	  move_right_two_vertex(chi_,smin_);
	}
	CentralDecomposition.store_right(getName(),right_size()+1);
	for (uMPXInt l=left_size();l>0;--l){
	  Data this_step(move_left_two_vertex(chi_,smin_));
	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);//at midpoint, push output
	    CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
	  }
	  else if (left_size()>right_size()) {CentralDecomposition.store_right(getName(),right_size()+1);} //just store right
	}
	CentralDecomposition.store_left(getName(),left_size()+1);
	for (uMPXInt r=right_size();r>left_size();--r){
	  Data this_step(move_right_two_vertex(chi_,smin_));
	  if (left_size()==right_size()) {
	    output_ref_.push(this_step);
	    CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
	  }
	  else if (left_size()<right_size()) {CentralDecomposition.store_left(getName(),left_size()+1);} //just store left
	}
      }
    }
  }

  void ExcitedStateFiniteDMRG::init(double chi, double smin){
    //first step
    //LeftBlock and RightBlock are correct? but we need to solve centre and store.
    const EigenStateArray& spectrum=getH().GetPhysicalSpectrum();
    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Calculating Excited State" << std::endl;

    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    CentralDecomposition.store(getName(),left_size()+1,right_size()+1);

    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(1.0); //overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    output_ref_.push(results);
  }

  Data ExcitedStateFiniteDMRG::move_right_two_vertex(uMPXInt chi, double smin){
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
    // PROBLEM IS HERE
    //  INCORRECT STORAGE NAME FOR HIGHER EXCITED STATES?

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
      /*std::cout << "Shifting" <<std::endl;
      std::cout << pb.getName() <<std::endl;
      pb.getLeftBlock().print_indices();
      pb.getRightBlock().print_indices();*/
    }

    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;
    previous_lambda_=CentralDecomposition.Values;
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(1.0);//overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return results;
  }

  Data ExcitedStateFiniteDMRG::move_left_two_vertex(uMPXInt chi, double smin){
    //if this is the initial run through, or we need a dummy block only, we get blocks from GS storage
    //first pass to left, we will have the left blocks for the right hand side only
    std::string StorageName;
    if (left_size()==1){//just need a dummy block
      StorageName=GSName_;
    }
    else if ((left_size()<=right_size()+middle_size()) && init_flag_){
      StorageName=HBlocksName_;
    }
    else {
      StorageName=getName();
    }
    //std::string StorageName = (((left_size()<=right_size()+middle_size()) && init_flag_) || left_size()==1) ? HBlocksName_ : getName();
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
      //std::cout << "Shifting" <<std::endl;
      //pb.getLeftBlock().print_indices();
      //pb.getRightBlock().print_indices();
    }

    Data results;
    results.Real_measurements.push_back(double(PBlocks_.size()));
    std::cout << "Position: " << left_size() << " " << right_size() << std::endl;
    previous_lambda_=CentralDecomposition.Values;
    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results),chi,smin);
    CentralDecomposition.SquareRescale(1.0);
    double S_E(entropy(CentralDecomposition.Values));
    results.Real_measurements.push_back(S_E);
    results.Real_measurements.push_back(CentralDecomposition.Truncation);
    results.Real_measurements.push_back(1.0);//overlap not calculated at moment
    std::cout <<"Current Bond Dimension: " << CentralDecomposition.Values.size() << ", Entropy: " << S_E << std::endl;
    return results;
  }

  void ExcitedStateFiniteDMRG::run(uMPXInt num_sweeps, uMPXInt chi, double smin){
    if (init_flag_) init(chi,smin);
    std::cout << "Performing finite sweeps" << std::endl;
    //we start at the middle of the system
    //in case we need excited states later, we store MPS matrices as well as the L and R blocks.
    for (uMPXInt n=num_sweeps;n>0;--n){//sweep towards right
      std::cout << std::endl << "Starting sweep: " << num_sweeps-n+1<< std::endl;
      for (uMPXInt r=right_size();r>0;--r){
	move_right_two_vertex(chi,smin);
      }
      CentralDecomposition.store_right(getName(),right_size()+1);
      for (uMPXInt l=left_size();l>0;--l){
	Data this_step(move_left_two_vertex(chi,smin));
	if (left_size()==right_size()) {
	  output_ref_.push(this_step);//at midpoint, push output
	  CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
	}
	else if (left_size()>right_size()) {CentralDecomposition.store_right(getName(),right_size()+1);} //just store right
      }
      CentralDecomposition.store_left(getName(),left_size()+1);

      if (n==num_sweeps) init_flag_=0; //first time through, turn off init here as we have formed all blocks once...

      for (uMPXInt r=right_size();r>left_size();--r){
	Data this_step(move_right_two_vertex(chi,smin));
	if (left_size()==right_size()) {
	  output_ref_.push(this_step);
	  CentralDecomposition.store(getName(),left_size()+1,right_size()+1);//store left and right
	}
	else if (left_size()<right_size()) {CentralDecomposition.store_left(getName(),left_size()+1);} //just store left
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
    static char SMALLESTREAL[]={'S','R','\n'}; //lowest real part for energies
    //SparseMatrix InitialGuess(LeftH.GetPhysicalSpectrum().size()*RightH.GetPhysicalSpectrum().size(),1,1);
    //a good initial guess would be equal weighting to all combinations in target sector?
    //InitialGuess.entry(0,0,1.0);
    //InitialGuess.finalise();
    std::cout << "Eigensolver starting for two vertex wavefunction..." << std::endl;
    SparseHED decomp(H2.Eigs(TargetSector,2,SMALLESTREAL));//uses arpack, finds two eigenvals for fun
    decomp.printValues();
    std::cout << "Lowest energy/Number of Vertices: " << decomp.Values[0]/2.0 <<std::endl;
    std::cout << "Next lowest energy/Number of Vertices: " << decomp.Values[1]/2.0 <<std::endl;
    result.Real_measurements.push_back(decomp.Values[0]/2.0);
    std::vector<MPXIndex> wfindices;
    wfindices.push_back(MPXIndex(1,H2.GetPhysicalSpectrum())); //ingoing
    wfindices.push_back(MPXIndex(1,StateArray(1,TargetSector-TargetSector))); //ingoing
    wfindices.push_back(MPXIndex(1,H2.GetPhysicalSpectrum())); //ingoing
    wfindices.push_back(MPXIndex(0,StateArray(1,TargetSector))); //outgoing
    return MPX_matrix(H2.GetPhysicalSpectrum(),wfindices,2,reshape(decomp.EigenVectors.ExtractColumns(std::vector<MPXInt>(1,0)),H2.GetPhysicalSpectrum().size()));
  }

  Prediction MakePrediction(const MPSDecomposition& Decomp, const std::vector<double>& PreviousLambda){
    const EigenStateArray& spectrum(Decomp.LeftMatrix.GetPhysicalSpectrum());
    Prediction ans;
    std::cout << "Forming State Prediction Vector" << std::endl;
    MPX_matrix svals(spectrum,Decomp.LeftMatrix.Index(2),Decomp.Values);
    MPX_matrix LambdaLB(reorder(contract(Decomp.LeftMatrix,0,svals,0,contract20),0,reorder102,1));
    MPXDecomposition RotDecomp(LambdaLB.SVD());
    ans.LambdaL=contract_to_sparse(RotDecomp.ColumnMatrix,0,MPX_matrix(spectrum,RotDecomp.ColumnMatrix.Index(1),RotDecomp.Values),0,contract10);
    MPX_matrix ALambdaR(reorder(contract(svals,0,Decomp.RightMatrix,0,contract10),0,reorder102,2));
    RotDecomp=ALambdaR.SVD();
    ans.LambdaR=contract_to_sparse(MPX_matrix(spectrum,RotDecomp.RowMatrix.Index(0),RotDecomp.Values),0,RotDecomp.RowMatrix,0,contract10);
    MPX_matrix InversePreviousLambda(spectrum,Decomp.LeftMatrix.Index(1),PreviousLambda,1); //1 means take inverse values    
    ans.Guess=reshape_to_vector(contract(contract(ALambdaR,0,InversePreviousLambda,0,contract20),0,LambdaLB,0,contract20));
    //ans.Guess.rescale(sqrt(1.0/ans.Guess.square_norm()));

    return ans;
  }

  MPX_matrix TwoVertexWavefunction(const MPX_matrix& LB, const MPO_matrix& H, const MPX_matrix& RB, const std::vector<ProjectorBlocks>* ProjectorBlocksPtr, MPXInt NumVertices, Data& result, SparseMatrix* guessptr){
    //Need to figure out 'length' then choose strategy appropriately.
    //form LeftPart, RightPart
    TwoVertexComponents stuff(LB,H,RB,ProjectorBlocksPtr);
    MPXInt num_to_find = stuff.length() >= 5 ? 4 : 1;
    SparseHED decomp(stuff.HED(num_to_find,SMALLESTREAL,guessptr));
    std::cout << "Lowest energy/Number of Vertices: " << decomp.Values[0]/double(NumVertices) <<std::endl;
    if (num_to_find>1) {
      std::cout << "Next lowest eigenvalues/Number of Vertices: ";
      for (MPXInt n=1;n<num_to_find;++n) {
	std::cout << decomp.Values[n]/double(NumVertices)  << " ";
      }
      std::cout << std::endl;
    }
    result.Real_measurements.push_back(decomp.Values[0]/double(NumVertices));
    //need to reshape from a vector
    //(sigma a), (sigma a) matrix form
    return MPX_matrix(LB.GetPhysicalSpectrum(),stuff.indices,2,reshape(decomp.EigenVectors.ExtractColumns(std::vector<MPXInt>(1,0)),LB.GetPhysicalSpectrum().size()*LB.Index(5).size()));
  }

  double CheckConvergence(const Prediction& guess,const std::vector<double>& PreviousLambda){
    return 1.0-Sum((guess.LambdaR*SparseMatrix(PreviousLambda)).SVD());
  }

  TwoVertexComponents::TwoVertexComponents(const MPX_matrix& L, const MPO_matrix& HMPO, const MPX_matrix& R, const std::vector<ProjectorBlocks>* P, const State* StatePtr) : LeftBlock(L),H(HMPO),RightBlock(R),ProjectorsPtr(P),TargetStatePtr(StatePtr),m_length(L.Index(5).size()*H.Index(2).size()*R.Index(4).size()*H.Index(2).size()),LeftPart(contract(H,0,LeftBlock,0,contract13)),RightPart(contract(H,0,RightBlock,0,contract32)) {
    
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

  SparseHED TwoVertexComponents::HED(MPXInt numevals, char which[3],const SparseMatrix* initial) const {
    //set up workspace (Evals, Evecs) SparseHED
    uMPXInt fulldim(this->length());
    std::cout <<"Eigensolver for matrix of length " << fulldim << std::endl;
    std::cout <<"Using reduced subspace of length " << allowed_indices.size() <<std::endl;

    SparseVectorWithRestriction guess_struct(initial,&allowed_indices);

    //internal contractions may be large.
    //however if allowed_indices is small we should use a dense method

    if (allowed_indices.size()<400){
      //just make the (dense) matrix and use lapack
      static std::pair<const std::vector<MPXInt>,const std::vector<MPXInt> > condition={{1,7},{2,6}};
      //TranslationBlock<DenseMatrix> TB(contract_conditional<DenseMatrix>(LeftPart,0,RightPart,0,contract21,condition));
      TranslationBlock<DenseMatrix> TB(contract_conditional<DenseMatrix>(contract(H,0,LeftBlock,0,contract13),0,contract(H,0,RightBlock,0,contract32),0,contract21,condition));
      //use lapack
      DenseHED dense_ans(TB.Block.HED(numevals,which));
      return SparseHED(std::vector<double>(dense_ans.Values,dense_ans.Values+numevals),TB.TranslateRows(dense_ans.EigenVectors));
    }
    else {
      std::cout <<"Allocate storage for evals and evecs" << std::endl;
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      SparseHED ans(fulldim,numevals);
      std::cout <<"Calling arpack_eigs()..." << std::endl;
      arpack::arpack_eigs<TwoVertexComponents,SparseVectorWithRestriction> eigensystem(this,&TwoVertexMPOMPSMultiply,allowed_indices.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
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
    for (struct {std::vector<std::array<MPXInt,2> >::const_iterator cit; MPXInt idx;} itstruct ={arraystuff->rows_and_cols.begin(), 0};itstruct.cit!=arraystuff->rows_and_cols.end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    MPX_matrix Vector(arraystuff->LeftPart.GetPhysicalSpectrum(),arraystuff->indices,2,V.cheap_no_transpose_finalise());
    DumbExtractWithZeros(reshape(contract_to_sparse(contract(arraystuff->LeftPart,0,Vector,0,contract1071),0,arraystuff->RightPart,0,contract116276),arraystuff->length()),arraystuff->allowed_indices,out);
    //projector bits
    const std::vector<TensorWeightPair>& PTensors(arraystuff->ProjectorTensors);
    //loop through list
    for (auto&& PT : PTensors){
      //contract_to_sparse with in vector and multiply by weight
      //std::cout << "Projector vector contraction" << std::endl;
      //PT.first.print_indices();
      //Vector.print_indices();

      SparseMatrix s(contract_to_sparse(PT.first,0,Vector,0,contract00112233));//really should be a scalar
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

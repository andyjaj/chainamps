#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <iomanip>

//#include "arpack_interface.hpp" //arpack
#include "sparse_interface.hpp" //SparseHED
#include "states.hpp" //EigenStateArray etc.
#include "MPX.hpp" //MPX_matrix etc.
//#include "DMRG_routines.hpp"
#include "fDMRG_routines.hpp"
//#include "data.hpp"

namespace ajaj {

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

    CentralDecomposition=TwoVertexSVD(TwoVertexWavefunction(getLeftBlock(),getH(),getRightBlock(),&PBlocks_,size(),results,&Guess,converge),chi,truncation);
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
    //results.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition)); 
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
    //results.Real_measurements.push_back(CheckConvergence(pred_,CentralDecomposition));
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

}

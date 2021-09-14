/**
 *@file TEBD_DYN_MEASURE.cpp Driver file for dynamical measurements made on TEBD files.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>
#include <vector>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "TEBD_routines.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){

  ajaj::TEBD_DYN_Args RuntimeArgs(argc,argv); //Process the command line input
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs)); //Create the model, based on input
    ajaj::uMPXInt CHI(RuntimeArgs.chi()); //bond dimension
    double trunc(RuntimeArgs.trunc()); //allowed truncation error
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices()); //number of vertices/nodes
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order()); //set the trotter decomposition order

    //TEBD_DYN_MEASURE can only measure one pair of operators per run, because different operators lead to different
    //time evolution for unequal time correlators.
    //The required operators need their MPOs built, the below is probably overkill, as it is borrowed from other drivers
    //with more measurements.
    std::vector<ajaj::NamedMPO_matrix> generated_MPOs; //storage for MPOs
    std::vector<std::string> ops({RuntimeArgs.Op1_name(), RuntimeArgs.Op2_name()});
    std::vector<ajaj::uMPXInt> operator_storage_position;//indices for the operator's position in storage
    
    for (auto&& op : ops){
      bool found(0);
      size_t index;
      for (size_t m=0;m< generated_MPOs.size();++m){//search through generated MPOs (initial size is zero, so skipped on first run)
	if (op==generated_MPOs.at(m).Name){ //if found in the list, then skip it
	  found=1;
	  index=m;
	  break;//found
	}
      }

      if (!found){ //if not already generated
	if (myModel.vertex.operator_exists(op)){ //do matrix elements already exist in vertex def?
	  std::cout << "Operator " << op << " found in vertex definition." << std::endl << std::endl;
	  found=1;
	  index=generated_MPOs.size();
	  generated_MPOs.emplace_back(op,myModel.vertex.make_one_site_operator(op));
	}
	else {
	  std::cout <<"Assuming name includes position info" <<std::endl;
	  ajaj::ShiftedOperatorInfo trial(op);
	  if (myModel.vertex.operator_exists(trial.Name)){
	    std::cout << std::endl;
	    std::cout << "Operator " << trial.Name << " found in vertex definition." << std::endl;
	    std::cout << "Requested unitary transformation using basis quantum number " << trial.WhichCharge << " and factor " << trial.Factor <<std::endl << std::endl;
	    if (trial.WhichCharge<myModel.basis().getChargeRules().size()){
	      found=1;
	      index=generated_MPOs.size();
	      generated_MPOs.emplace_back(op,myModel.vertex.make_one_site_operator(trial));
	    }
	  }
	}
      }
      //at this point if we still can't find the operator definition, we give up
      if(!found) {
	std::cout << "Operator " << op << " couldn't be found or created from predefined matrix elements." <<std::endl;
	std::cout << "Note that the name should be a name defined in your operators file, or by a built-in model." <<std::endl;
	std::cout << "It should NOT be a .SPARSEMATRIX file name!" <<std::endl;
	return 1;
      } //not found

      operator_storage_position.emplace_back(index);
      
    }

    //The model and operators are now loaded
    //Next load the info from the time file.
    //The name should be INPUTFILE_TEBD_NUMVERTICES_INITIALSTATENAME_Evolution.dat
    std::stringstream tfnamestream;
    tfnamestream << RuntimeArgs.filename() << "_TEBD_" << RuntimeArgs.num_vertices() << "_" << RuntimeArgs.initial_state_name() << "_Evolution.dat";
    //From this need the indices and times
    //Note that there is no way to discern what Trotter order was used!
    typedef std::pair<size_t,double> idx_dble;//typedef for convenience
    std::vector<idx_dble> idx_times; //storage for timeslice index and time
    
    idx_times.emplace_back(0,0.0);
    
    if (tfnamestream.str().empty()){
      std::cout << "Time filename empty? Something is wrong." << std::endl;
      return 1;
    }
    else {
      std::ifstream time_file;
      time_file.open(tfnamestream.str().c_str(),ios::in);
      
      if (time_file.is_open()){
	std::cout << "Opened " << tfnamestream.str() << "." << std::endl;
	size_t idx;
	double time;
	std::string sdump;
	time_file >> std::ws;
	if (time_file.peek()=='#'){getline(time_file,sdump);}//dump the comment line
	while (time_file >> idx >> time){
	  getline(time_file,sdump); //dump the rest of the line
	  idx_times.emplace_back(idx,time);
	}

	std::cout << "There are " << idx_times.size() << " (index,time) pairs." <<std::endl;
	std::cout << "Note, index 0 at time 0.0 is created by default. " <<std::endl << std::endl;
	
      }
      else {
	std::cout <<"Could not open specified time data file '" << tfnamestream.str() << "'." <<std::endl;
	std::cout <<"File may be missing or you needed to specify an initial state name other than '" << RuntimeArgs.initial_state_name() << "'" << std::endl ;

	return 1;
      }
    }


    //make a list of the timeslices at which the Hamiltonian changes (includes 0)
    std::vector<ajaj::uMPXInt> H_MPO_indices; 
    for (auto& t : idx_times){
      if (ajaj::check_for_H_MPO_file(ajaj::SAVEALLNAME,t.first)){
	H_MPO_indices.push_back(t.first);
      }
      else H_MPO_indices.push_back(H_MPO_indices.back());
    }
    std::cout << "H_MPO index at each time slice " <<std::endl;
    ajaj::uMPXInt hi=0;
    for (auto &h : H_MPO_indices){
      std::cout << idx_times[hi++].second << " " << h <<std::endl;
    }
    
    //Before continuing need to check for and load H_MPO for time slice 0.
    /*std::stringstream FirstHMPOnamestream;
    FirstHMPOnamestream << ajaj::SAVEALLNAME << "_0_H.MPO_matrix";
    ajaj::MPO_matrix H_MPO(load_MPO_matrix(FirstHMPOnamestream.str(),myModel.basis()));
    if (!H_MPO.isConsistent()){
      std::cout << "Problem with " << FirstHMPOnamestream.str() << std::endl;
      return 1;
    }
    else {
      std::cout << "Successfully loaded " << FirstHMPOnamestream.str() << std::endl;
      H_MPO.print_indices();
      }*/
    

    //set up results file
    std::ostringstream commentlinestream;
    commentlinestream << "Index, t1, t2";
    
    for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1) {
      commentlinestream << ", Re(<" << RuntimeArgs.Op1_name() <<"(" << y1 << ",t1)" << RuntimeArgs.Op2_name() <<"(" << RuntimeArgs.y2() << ",t2)>)";
      commentlinestream << ", Im(<" << RuntimeArgs.Op1_name() <<"(" << y1 << ",t1)" << RuntimeArgs.Op2_name() <<"(" << RuntimeArgs.y2() << ",t2)>)";
    }
    
    std::stringstream outfilenamestream;
    outfilenamestream<<RuntimeArgs.filename()<<"_DYN_"<<number_of_vertices<<"_"<<RuntimeArgs.initial_state_name()<<"_" << RuntimeArgs.Op1_name() << "_y1_" << RuntimeArgs.Op2_name() << "_"<< RuntimeArgs.y2() <<".dat";
    ajaj::DataOutput results(outfilenamestream.str(),commentlinestream.str());
    ajaj::DataOutput dummyresults;

    std::string WorkingName("TempMPS");
    
    //loop over t2

    ajaj::uMPXInt H_MPO_idx=0;
    ajaj::uMPXInt results_index=0;

    //Don't want to be recreating the bond operators needlessly
    
    
    for (const auto& t2p : idx_times){
      std::cout << t2p.first << " " << t2p.second <<std::endl;

      //We will need to load the MPS at each timeslice
      //Its files should be named SAVEALLNAME_timesliceidx_INITIALSTATENAME_Left_i.MPS_matrix
      //SAVEALLNAME is defined in TEBD_routines.hpp
      
      std::stringstream mpsrootnamestream;
      mpsrootnamestream << ajaj::SAVEALLNAME << "_" << t2p.first << "_" << RuntimeArgs.initial_state_name();
      
      ajaj::FiniteMPS F(myModel.basis(),mpsrootnamestream.str(),WorkingName,number_of_vertices,1/*should be canonical*/,number_of_vertices); //load finite MPS from storage, and save a working copy 'F'.

      //equal time bit
      //removed for testing
      {
	ajaj::ConstFiniteMPS CF(F);
	std::vector<std::complex<double> > equal_time_results;
	for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
	  equal_time_results.emplace_back(ajaj::GeneralisedOverlap(CF,std::vector<ajaj::meas_pair>{{y1,&(generated_MPOs[operator_storage_position[0]].Matrix)},{RuntimeArgs.y2(),&(generated_MPOs[operator_storage_position[1]].Matrix)}}));
	}
	results.push(++results_index,ajaj::Data(std::vector<double>{t2p.second,t2p.second},equal_time_results));
      }
      
      //apply operator_2 to F and canonize new F
      std::complex<double> op2weight=ApplySingleVertexOperatorToMPS(generated_MPOs[operator_storage_position[1]].Matrix,F,RuntimeArgs.y2(),ajaj::MPSCanonicalType::Left);

      //need to setup a TEBD object with an initial MPO and step size

      //create TEBD object for inner loop
      std::vector<ajaj::MultiVertexMeasurement> dummy_measurements; //TEBD object requires a measurement object, even if empty.
      //time ordered part
      {
	ajaj::TEBD finrun(myModel.H_MPO,F,dummyresults); //doesn't mess about creating bond operators.
	for (ajaj::uMPXInt tidx=t2p.first; tidx<idx_times.back().first; ++tidx){
	  //check H for this step and update if necessary
	  if (tidx==t2p.first || H_MPO_idx!=H_MPO_indices[tidx]){
	    double stepsize= idx_times[tidx+1].second-idx_times[tidx].second;
	    H_MPO_idx=H_MPO_indices[tidx];
	    std::cout << "Time slice " << tidx << ". Updating H_MPO idx to " << H_MPO_indices[H_MPO_idx] <<std::endl;
	    std::stringstream tSliceHMPOnamestream;
	    tSliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << H_MPO_idx << "_H.MPO_matrix";
	    finrun.change_bond_operator(myModel.update_H_MPO_from_file(tSliceHMPOnamestream.str()),stepsize,trotter_order);
	  }
	  //evolve one step
	  //std::cout << "Evolve one step " <<std::endl;
	  finrun.evolve(1,dummy_measurements,CHI,trunc);
	  ajaj::ConstFiniteMPS TEBDCF(finrun.GetEvolvingState());
	  std::stringstream mpst1namestream;
	  mpst1namestream << ajaj::SAVEALLNAME << "_" << tidx+1 << "_" << RuntimeArgs.initial_state_name();
	  //std::cout << "Loading bra state" <<std::endl;
	  ajaj::ConstFiniteMPS BraState(myModel.basis(),mpst1namestream.str(),number_of_vertices);
	  std::vector<std::complex<double> > unequal_time_results;
	  for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
	    //std::cout << "measurement" <<std::endl;
	    
	    unequal_time_results.emplace_back(op2weight*ajaj::GeneralisedOverlap(BraState,TEBDCF,generated_MPOs[operator_storage_position[0]].Matrix,y1));
	  }
	  results.push(++results_index,ajaj::Data(std::vector<double>{idx_times[tidx+1].second,t2p.second},unequal_time_results));
	}
      }

      //anti-time ordered
      if (RuntimeArgs.include_reverse()){
	ajaj::TEBD finrun_rev(myModel.H_MPO,F,dummyresults); //doesn't mess about creating bond operators.

	for (ajaj::uMPXInt tidx=t2p.first; tidx>0; --tidx){
	  //check H for this step and update if necessary
	  if (tidx==t2p.first || H_MPO_idx!=H_MPO_indices[tidx-1]){
	    double stepsize= idx_times[tidx-1].second-idx_times[tidx].second;
	    H_MPO_idx=H_MPO_indices[tidx-1];
	    std::cout << "Time slice " << tidx << ". Updating H_MPO idx to " << H_MPO_indices[H_MPO_idx] <<std::endl;	   
	    std::stringstream tSliceHMPOnamestream;
	    tSliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << H_MPO_idx << "_H.MPO_matrix";
	    finrun_rev.change_bond_operator(myModel.update_H_MPO_from_file(tSliceHMPOnamestream.str()),stepsize,trotter_order);
	  }
	  //evolve one step
	  finrun_rev.evolve(1,dummy_measurements,CHI,trunc);
	  ajaj::ConstFiniteMPS TEBDCF(finrun_rev.GetEvolvingState());
	  std::stringstream mpst1namestream;
	  mpst1namestream << ajaj::SAVEALLNAME << "_" << tidx-1 << "_" << RuntimeArgs.initial_state_name();
	  //std::cout << "Loading bra state" <<std::endl;
	  ajaj::ConstFiniteMPS BraState(myModel.basis(),mpst1namestream.str(),number_of_vertices);
	  std::vector<std::complex<double> > unequal_time_results;
	  for (ajaj::uMPXInt y1=1;y1<=number_of_vertices;++y1){
	    //std::cout << "measurement" <<std::endl;
	    
	    unequal_time_results.emplace_back(op2weight*ajaj::GeneralisedOverlap(BraState,TEBDCF,generated_MPOs[operator_storage_position[0]].Matrix,y1));
	  }
	  results.push(++results_index,ajaj::Data(std::vector<double>{idx_times[tidx-1].second,t2p.second},unequal_time_results));
	}
      }
      
    }
    
    return 0;
    
  }

  return 1;
}

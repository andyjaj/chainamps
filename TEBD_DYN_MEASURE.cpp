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

  ajaj::TEBD_DYN_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double trunc(RuntimeArgs.trunc());
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());

    std::vector<ajaj::NamedMPO_matrix> generated_MPOs; //actual storage for MPOs
    typedef std::vector<std::vector<std::pair<size_t,ajaj::uMPXInt> > > MPOIndexVertexPairs;
    MPOIndexVertexPairs measurement_lookups;
    
    std::ostringstream measurednames;
    //build all required measurement MPOs (no repeats) and index them

    std::vector<std::string> ops({RuntimeArgs.Op1_name(), RuntimeArgs.Op2_name()});

    std::vector<ajaj::uMPXInt> operator_storage_position;
    
    for (auto&& op : ops){
      bool found(0);
      size_t index;
      for (size_t m=0;m< generated_MPOs.size();++m){//search through generated MPOs
	if (op==generated_MPOs.at(m).Name){
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

    //We have loaded the model and operators
    //We now need to load the info from the time file.
    //The name will be INPUTFILE_TEBD_NUMVERTICES_INITIALSTATENAME_Evolution.dat
    std::stringstream tfnamestream;
    tfnamestream << RuntimeArgs.filename() << "_TEBD_" << RuntimeArgs.num_vertices() << "_" << RuntimeArgs.initial_state_name() << "_Evolution.dat";
    
    //From this we will need the indices and times
    typedef std::pair<size_t,double> idx_dble;
    std::vector<idx_dble> idx_times;
    
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
    std::vector<ajaj::uMPXInt> H_MPO_change_indices; 
    for (auto& t : idx_times){
      if (ajaj::check_for_H_MPO_file(ajaj::SAVEALLNAME,t.first)){
	H_MPO_change_indices.push_back(t.first);
      }
    }

    std::cout << "There are " << H_MPO_change_indices.size() << " time slices at which the Hamiltonian changes (including t=0)." <<std::endl;

    
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
    
    for (auto y1=1;y1<=number_of_vertices;++y1)
      commentlinestream << ", <" << RuntimeArgs.Op1_name() <<"(" << y1 << ",t1)" << RuntimeArgs.Op2_name() <<"(" << RuntimeArgs.y2() << ",t2)>";
    
    std::stringstream outfilenamestream;
    outfilenamestream<<RuntimeArgs.filename()<<"_DYN_"<<number_of_vertices<<"_"<<RuntimeArgs.initial_state_name()<<"_" << RuntimeArgs.Op1_name() << "_y1_" << RuntimeArgs.Op2_name() << "_"<< RuntimeArgs.y2() <<".dat";
    ajaj::DataOutput results(outfilenamestream.str(),commentlinestream.str());
    ajaj::DataOutput dummyresults;

    std::string WorkingName("TempMPS");
    
    //loop over t2

    ajaj::uMPXInt current_H_index=0;
    ajaj::uMPXInt results_index=0;
    
    for (const auto& t2p : idx_times){
      std::cout << t2p.first << " " << t2p.second <<std::endl;

     
      //We will need to load the MPS at each timeslice (i.e. t=0)
      //Its files should be named should be SAVEALLNAME_timesliceidx_INITIALSTATENAME_Left_i.MPS_matrix
      //SAVEALLNAME is defined in TEBD_routines.hpp
      
      std::stringstream mpsrootnamestream;
      mpsrootnamestream << ajaj::SAVEALLNAME << "_" << t2p.first << "_" << RuntimeArgs.initial_state_name();
      ajaj::FiniteMPS F(myModel.basis(),mpsrootnamestream.str(),WorkingName,number_of_vertices,1/*should be canonical*/,number_of_vertices);
      //we have loaded and checked files, and stored a working copy which is what we will use	
      //apply operator and canonize
      std::complex<double> op2weight=ApplySingleVertexOperatorToMPS(generated_MPOs[operator_storage_position[1]].Matrix,F,RuntimeArgs.y2(),ajaj::MPSCanonicalType::Left);
      std::cout << op2weight <<std::endl;

      std::vector<ajaj::MultiVertexMeasurement> dummy_measurements; //TEBD object requires a measurement object, even if empty.
      
      //time ordered
      {

	//establish the correct start slice index for our H_MPO
	if (t2p.first!=current_H_index){
	  for (ajaj::uMPXInt h=H_MPO_change_indices.size()-1;h<H_MPO_change_indices.size();--h){
	    if (H_MPO_change_indices[h]<=t2p.first){
	      if (H_MPO_change_indices[h]!=current_H_index){
		current_H_index=H_MPO_change_indices[h];
		std::cout << "Updating H_MPO to that defined at timeslice " << current_H_index <<std::endl <<std::endl;
		std::stringstream t2SliceHMPOnamestream;
		t2SliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << current_H_index << "_H.MPO_matrix";
		myModel.update_H_MPO_from_file(t2SliceHMPOnamestream.str());
	      }
	      break;
	    }
	  }
	}

	//set the initial step size, necessary to set up the TEBD object.
	results.push(++results_index,ajaj::Data(std::vector<double>{t2p.second,t2p.second}));

	if (t2p.first<idx_times.size()-1){
	  double stepsize= idx_times[t2p.first+1].second-t2p.second;//idx_times[t2p.first].second;

	  //calculate the equal time piece first
	  {
	    ajaj::ConstFiniteMPS CF(myModel.basis(),WorkingName,number_of_vertices);
	  }
	  
	  //evolve one step
	  ajaj::TEBD finrun(myModel.H_MPO,F,stepsize,dummyresults,trotter_order,nullptr,0);
	  finrun.evolve(1,dummy_measurements,CHI,trunc,1);

	  //inner loop over timeslices, from t2p.first to the max time doing TEBD evolution steps
	  for (ajaj::uMPXInt t1sliceindex=t2p.first+1; t1sliceindex<=idx_times.back().first; ++t1sliceindex){
	    //now measure...
	    ajaj::ConstFiniteMPS TEBDCF(finrun.GetEvolvingState());
	    //std::cout << TEBDCF.name() <<std::endl;
	    results.push(++results_index,ajaj::Data(std::vector<double>{idx_times[t1sliceindex].second,t2p.second}));

	    
	    if (t1sliceindex!=idx_times.back().first){
	      //do updates of step size and H_MPO
	      std::cout << idx_times[t1sliceindex].second << " " << t2p.second <<std::endl<<std::endl;
	      
	      stepsize=idx_times[t1sliceindex].second-t2p.second;
	      
	      if (ajaj::check_for_H_MPO_file(ajaj::SAVEALLNAME, t1sliceindex) && (t1sliceindex != current_H_index)){ //if this is an update step, and we don't already have the right H
		std::cout << "Updating H_MPO from timeslice " << current_H_index << " to " << t1sliceindex <<std::endl <<std::endl;
		current_H_index=t1sliceindex;
		std::stringstream t1SliceHMPOnamestream;
		t1SliceHMPOnamestream << ajaj::SAVEALLNAME << "_" << current_H_index<< "_H.MPO_matrix";
		finrun.change_bond_operator(myModel.update_H_MPO_from_file(t1SliceHMPOnamestream.str()),stepsize);
	      }

	      //finrun.evolve();
	      
	    }
	  }
	}
	
      }
      //anti timeordered
      if (RuntimeArgs.include_reverse()){
	std::stringstream mpsrootnamestream;
	mpsrootnamestream << ajaj::SAVEALLNAME << "_t2p.first_" << RuntimeArgs.initial_state_name();
	ajaj::FiniteMPS F(myModel.basis(),mpsrootnamestream.str(),number_of_vertices,1/*is canonical*/,number_of_vertices);
      }

      
      
    }

    
    return 0;
  }

  return 1;
}

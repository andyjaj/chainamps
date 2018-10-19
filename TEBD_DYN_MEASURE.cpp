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
    commentlinestream << "Index1, t1, Index2, t2, y1, y2";
    commentlinestream << ", <" << RuntimeArgs.Op1_name() <<"(y1,t1)" << RuntimeArgs.Op2_name() <<"(y2,t2)>";
    std::stringstream outfilenamestream;
    outfilenamestream<<RuntimeArgs.filename()<<"_DYN_"<<number_of_vertices<<"_"<<RuntimeArgs.initial_state_name()<<"_" << RuntimeArgs.Op1_name() << "_y1_" << RuntimeArgs.Op2_name() << "_"<< RuntimeArgs.y2() <<".dat";
    ajaj::DataOutput results(outfilenamestream.str(),commentlinestream.str());

    std::string WorkingName("TempMPS");
    
    //loop over t2
    
    for (const auto& t2p : idx_times){
      std::cout << t2p.first << " " << t2p.second <<std::endl;

      //time ordered
      {
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

	//loop over timeslices, from t2p.first to the max time doing TEBD evolution steps
	//// evaluate Op1 on evolving state at timeslice.
	
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

/**
 *@file TEBD_DRV.cpp Driver file for TEBD.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>

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

  ajaj::TEBD_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(1.0e-14);
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());
    ajaj::uMPXInt number_of_time_steps(RuntimeArgs.number_of_steps());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());
    ajaj::uMPXInt measurement_interval(RuntimeArgs.measurement_interval());
    double time_step(RuntimeArgs.step_size());

    std::vector<ajaj::NamedMPO_matrix> generated_MPOs; //actual storage for MPOs
    typedef std::vector<std::vector<std::pair<size_t,ajaj::uMPXInt> > > MPOIndexVertexPairs;
    MPOIndexVertexPairs measurement_lookups;
    std::ostringstream measurednames;
    //build all required measurement MPOs (no repeats) and index them
    for (auto&& fm : RuntimeArgs.finite_measurements()){ //loop over all measurements
      std::ostringstream currentname;
      std::vector<std::pair<size_t,ajaj::uMPXInt> > temp_measurement_vec;
      for (auto&& op : fm){//for all operators in an individual measurement
	//first see if we already generated this
	bool found(0);
	size_t index;
	for (size_t m=0;m< generated_MPOs.size();++m){//search through generated MPOs
	  if (op.first==generated_MPOs.at(m).Name){
	    found=1;
	    index=m;
	    break;//found
	  }
	}
	if (!found){ //if not already generated
	  if (myModel.vertex.operator_exists(op.first)){ //do matrix elements already exist in vertex def?
	    std::cout << "Operator " << op.first << " found in vertex definition." << std::endl;
	    found=1;
	    index=generated_MPOs.size();
	    generated_MPOs.emplace_back(op.first,myModel.vertex.make_one_site_operator(op.first));
	  }
	  else {
	    std::cout <<"Assuming name includes position info" <<std::endl;
	    ajaj::ShiftedOperatorInfo trial(op.first);
	    if (myModel.vertex.operator_exists(trial.Name)){
	      std::cout << "Operator " << trial.Name << " found in vertex definition." << std::endl;
	      std::cout << "Requested unitary transformation using basis quantum number " << trial.WhichCharge << " and factor " << trial.Factor <<std::endl;
	      if (trial.WhichCharge<myModel.basis().getChargeRules().size()){
		found=1;
		index=generated_MPOs.size();
		generated_MPOs.emplace_back(op.first,myModel.vertex.make_one_site_operator(trial));
	      }
	    }
	  }
	}
	if (found){
	  if (temp_measurement_vec.size()>0){
	    currentname << ":";
	  }
	  currentname << op.first <<"," <<op.second;
	  temp_measurement_vec.emplace_back(index,op.second);
	}
	else {
	  std::cout << "Operator " << op.first << "couldn't be found or created from predefined matrix elements." <<std::endl;
	  std::cout << "Note that the name should be a name defined in your operators file, or by a built-in model," <<std::endl;
	  std::cout << "NOT a .SPARSEMATRIX file name!" <<std::endl;
	  return 0;
	} //not found
      }
      if (temp_measurement_vec.size()==fm.size()){ //if we correctly generated enough mpos for measurement...
	measurement_lookups.emplace_back(temp_measurement_vec);
	measurednames << ", Re(" << currentname.str() << "), Im(" << currentname.str() << ") ";
      }
    }

    //go through and make pointers

    std::vector<ajaj::MPO_matrix> measured_operators(1,myModel.vertex.make_one_site_operator(1)); //for Ising this is the fermion occupation number on a chain
    //ajaj::DataOutput results(ajaj::OutputName(RuntimeArgs.filename(),"Evolution.dat"),"Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap), Re(Op1), Im(Op1), ...");
    std::ostringstream commentline;
    commentline << "Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap)";
    if (measurement_lookups.size()) commentline << measurednames.str();  
    ajaj::DataOutput results(ajaj::OutputName(RuntimeArgs.filename(),"Evolution.dat"),commentline.str());
    //measurements.emplace_back(ajaj::MultiVertexMeasurement(number_of_vertices/2/*first measurement position*/,&measured_operators[0],number_of_vertices/2+1/* second measurement position*/,&measured_operators[0]));

    std::vector<ajaj::MultiVertexMeasurement> measurements;
    for (auto&& m : measurement_lookups){
      if (m.size()==2){
	measurements.emplace_back(ajaj::MultiVertexMeasurement(m[0].second/*first measurement position*/,&(generated_MPOs[m[0].first].Matrix),m[1].second/* second measurement position*/,&(generated_MPOs[m[1].first].Matrix)));
      }
      else if (m.size()==1){
	measurements.emplace_back(ajaj::MultiVertexMeasurement(m[0].second/*first measurement position*/,&(generated_MPOs[m[0].first].Matrix)));
      }
      else {
	std::cout <<"Unsupported measurement!" <<std::endl;
	return 0;
      }
    }

    //TEBD, select between files, product state definition by c numbers, or default

    std::string StateName;
    ajaj::c_specifier_array CSpec;

    if (RuntimeArgs.initial_state_name()!=""){
      StateName=RuntimeArgs.initial_state_name();
    }
    else if (RuntimeArgs.c_number_filename()!=""){
      StateName=RuntimeArgs.c_number_filename();
      CSpec=ajaj::LoadCNumbers(RuntimeArgs.c_number_filename());
    }
    else {
      //default case, product state of vertex states in
      std::cout << "No initial state specified." << std::endl <<"Defaulting to product state of local basis [0] states." <<std::endl; 
      StateName="DefaultState";
      CSpec=ajaj::c_specifier_array(1,ajaj::c_specifier_vector(1,ajaj::c_specifier(0,1.0)));
    }

    std::cout << "Using initial state '" << StateName << "'." <<std::endl; 

    ajaj::FiniteMPS F(myModel.basis(),StateName,number_of_vertices,CSpec); //if CSpec is empty, nothing is changed.
    ajaj::TEBD finrun(myModel.H_MPO,F,time_step,results,trotter_order);
    finrun.evolve(number_of_time_steps,measurements,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
    if (finrun.good()){
      return 0;
    }
  }
  return 1;
}

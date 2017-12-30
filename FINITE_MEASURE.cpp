/**
 *@file FINITE_MEASURE.cpp Driver file for measurements on finite states.
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
  
  ajaj::fMEAS_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());


    //We need to generate the operators and store them
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
    commentline << "Index, Entropy"; //refers to order states are processed
    if (measurement_lookups.size()) commentline << measurednames.str(); 

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

    ajaj::DataOutput results("FINITE_RESULTS.dat",commentline.str());
    
    size_t Index(0);
    if (!RuntimeArgs.fdmrg_mode()){
      for (auto&& StateName : RuntimeArgs.state_names()){
	std::cout << StateName <<std::endl;
	ajaj::FiniteMPS F(myModel.basis(),StateName,number_of_vertices);
	ajaj::TEBD finrun(myModel.H_MPO,F,results);
	finrun.evolve(Index,measurements);
	if (!finrun.good()){
	  return 1;
	}
	++Index;
      }
      return 0;
    }
    else {
      bool gs=1;
      bool bad_name=0;
      while (!bad_name){
	std::string StateName;
	if (gs) {StateName="GroundState"; gs=0;}
	else {
	  std::stringstream fnss;
	  fnss << "Excited_" << Index;
	  StateName=fnss.str();
	}
	ajaj::FiniteMPS F(myModel.basis(),StateName,number_of_vertices);
	if (!F.valid_files()) {bad_name=1; continue;}
	
	std::cout << StateName <<std::endl;

	ajaj::TEBD finrun(myModel.H_MPO,F,results);
	finrun.evolve(Index,measurements);
	if (!finrun.good()){
	  return 1;
	}
	++Index;  
      }
      return 0;
    }
  }
  return 1;
}

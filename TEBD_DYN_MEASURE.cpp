/**
 *@file TEBD_DYN_MEASURE.cpp Driver file for dynamical measurements made on TEBD files.
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
    for (auto&& fm : RuntimeArgs.finite_measurements()){ //loop over all measurements
      std::ostringstream currentname;
      std::vector<std::pair<size_t,ajaj::uMPXInt> > temp_measurement_vec;
      //std::vector<ajaj::meas_pair> temp_mp_vec;
      
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
	    //generated_MPOs.back().Matrix.print_matrix();
	    //generated_MPOs.back().Matrix.print_indices(); /**< Print the dimensions of the indices. Colon indicates how many correspond to rows and how many to columns. */
	    //generated_MPOs.back().Matrix.print_indices_values();
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
	  std::cout << "Operator " << op.first << " couldn't be found or created from predefined matrix elements." <<std::endl;
	  std::cout << "Note that the name should be a name defined in your operators file, or by a built-in model." <<std::endl;
	  std::cout << "It should NOT be a .SPARSEMATRIX file name!" <<std::endl;
	  return 0;
	} //not found
      }
      if (temp_measurement_vec.size()==fm.size()){ //if we correctly generated enough mpos for measurement...
	measurement_lookups.emplace_back(std::move(temp_measurement_vec));
	measurednames << ", Re(" << currentname.str() << "), Im(" << currentname.str() << ") ";
      }
    }

    std::ostringstream commentline;
    commentline << "Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap)";
    if (measurement_lookups.size()) commentline << measurednames.str(); 

    std::vector<ajaj::MultiVertexMeasurement> measurements;
    
    for (auto&& m : measurement_lookups){//loop over all measurements
      std::vector<ajaj::meas_pair> mvmd;
      for (auto&& idx_v : m){//loop through all operators in measurement
	if (m.size()){
	  mvmd.emplace_back(idx_v.second,&(generated_MPOs[idx_v.first].Matrix));
	}
	else {
	  std::cout <<"Unsupported measurement!" <<std::endl;
	  return 0;
	}
      }

      measurements.emplace_back(ajaj::MultiVertexMeasurement(std::move(mvmd)));

    }

    return 0;
  }

  return 1;
}

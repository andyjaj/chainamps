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
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double trunc(RuntimeArgs.trunc());
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());
    ajaj::uMPXInt number_of_time_steps(RuntimeArgs.number_of_steps());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());
    ajaj::uMPXInt measurement_interval(RuntimeArgs.measurement_interval());
    double time_step_size(RuntimeArgs.step_size());
    bool save_all_flag(RuntimeArgs.save_all_flag());

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
	  std::cout << std::endl;
	  std::cout << "Operator " << op.first << " couldn't be found or created from predefined matrix elements." <<std::endl;
	  std::cout << "Note that the name should be a name defined in your operators file, or by a built-in model." <<std::endl;
	  std::cout << "It should NOT be a .SPARSEMATRIX file name!" <<std::endl;
	  return 1;
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

    //TEBD, select between files, product state definition by c numbers, or default
    std::string StateName;
    ajaj::c_specifier_array CSpec;

    if (!RuntimeArgs.initial_state_name().empty()){
      StateName=RuntimeArgs.initial_state_name();
    }
    else if (!RuntimeArgs.c_number_filename().empty()){
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

    std::stringstream Rss;
    Rss<<RuntimeArgs.filename()<<"_TEBD_"<<number_of_vertices<<"_"<<StateName;
    ajaj::DataOutput results(ajaj::OutputName(Rss.str(),"Evolution.dat"),commentline.str());
    ajaj::FiniteMPS F(myModel.basis(),StateName,number_of_vertices,CSpec); //if CSpec is empty, nothing is changed.
    
    //if we don't have time dependent couplings...
    if (!myModel.times().size()){ //need this for builtin models, with old style coupling params
      std::cout <<"Evolution hamiltonian is static." <<std::endl;
      ajaj::TEBD finrun(myModel.H_MPO,F,time_step_size,results,trotter_order,nullptr,save_all_flag);
      
      finrun.evolve(number_of_time_steps,measurements,CHI/*bond dimension*/,trunc,measurement_interval);
      if (finrun.good()){
	return 0;
      }
    }
    //if we do have time dep couplings...
    else {
      std::cout <<"Evolution Hamiltonian is time dependent." <<std::endl;
      
      //if ramp step size is smaller than step size, then use that until ramp over (check each time)
      //if step size smaller than ramp step size then use size that is commensurate with ramp step, but smaller than step size.
      ajaj::uMPXInt ramp_step=1;
      //do the first explicitly in order to create the TEBD object
      double ramp_step_size_1=myModel.times()[ramp_step]-myModel.times()[ramp_step-1];
      double current_step_size_1=ramp_step_size_1;
      ajaj::uMPXInt num_1=1;
      while (current_step_size_1>time_step_size){current_step_size_1=ramp_step_size_1/(++num_1);}
      if (num_1>number_of_time_steps) num_1=number_of_time_steps;

      ajaj::TEBD finrun(myModel.H_MPO,F,current_step_size_1,results,trotter_order,nullptr,save_all_flag);     
      finrun.evolve(num_1,measurements,CHI/*bond dimension*/,trunc/*min s val*/,measurement_interval);
      number_of_time_steps-=num_1;
      ++ramp_step;

      //now continue with time dep params
      while (number_of_time_steps>0 && ramp_step < myModel.times().size() && finrun.good()){
	double ramp_step_size=myModel.times()[ramp_step]-myModel.times()[ramp_step-1];
	double current_step_size=ramp_step_size;
	ajaj::uMPXInt num=1;
	while (current_step_size>time_step_size){current_step_size=ramp_step_size/(++num);}
	if (num>number_of_time_steps) num=number_of_time_steps;

	finrun.change_bond_operator(myModel.change_H_MPO(myModel.coupling_arrays()[ramp_step-1]),current_step_size);
	
	finrun.evolve(num,measurements,CHI/*bond dimension*/,trunc/*min s val*/,measurement_interval);
	number_of_time_steps-=num;
	++ramp_step;
      }
      
      if (number_of_time_steps>0 && finrun.good()){
	std::cout <<"End of time dependent hamiltonian stage." <<std::endl;
	finrun.change_bond_operator(myModel.change_H_MPO(myModel.coupling_arrays()[ramp_step-1]),time_step_size);
	finrun.evolve(number_of_time_steps,measurements,CHI/*bond dimension*/,trunc/*min s val*/,measurement_interval);
      }

      if (finrun.good()){
	return 0;
      }
    }
  }
  return 1;
}

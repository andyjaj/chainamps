/**
 *@file iTEBD_DRV.cpp Driver file for iTEBD.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>

#include "common_defs.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "TEBD_routines.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){

  ajaj::iTEBD_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    myModel.basis().print();

#ifndef DNDEBUG
    myModel.H_MPO.print_matrix();
#endif

    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(1.0e-14);
    ajaj::uMPXInt number_of_time_steps(RuntimeArgs.number_of_steps());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());
    ajaj::uMPXInt measurement_interval(RuntimeArgs.measurement_interval());
    double time_step_size(RuntimeArgs.step_size());
    std::vector<ajaj::MPO_matrix> measured_operators;// can insert runtime measurements here

    ajaj::UnitCell Initial(myModel.basis());
    std::string InitialStateName;

    if (!RuntimeArgs.initial_unit_cell().empty()){
      InitialStateName=RuntimeArgs.initial_unit_cell();
      std::ifstream infile;
      infile.open(RuntimeArgs.initial_unit_cell().c_str(),ios::in | ios::binary);
      if (infile.is_open()){
	Initial=ajaj::load_UnitCell_binary(infile,myModel.basis().getChargeRules(),myModel.basis());//populates basis
      }
      else {
	std::cout << "Failed to open " << RuntimeArgs.initial_unit_cell() << std::endl;
	return 1;
      }
    }
    else if (!RuntimeArgs.c_number_filename().empty()) {
      InitialStateName=RuntimeArgs.c_number_filename();
      Initial=MakeProductStateUnitCell(myModel.basis(),ajaj::LoadCNumbers(RuntimeArgs.c_number_filename()),ajaj::State(myModel.basis().getChargeRules()));
    }
    else {
      std::cout << "No initial state specified, generating default state." << std::endl;
      InitialStateName="DefaultState";
      Initial=MakeProductStateUnitCell(myModel.basis(),0,ajaj::State(myModel.basis().getChargeRules()));
    }

    std::istringstream fnss(InitialStateName);
    getline(fnss,InitialStateName,'.');
    InitialStateName=ajaj::StripName(InitialStateName);
    if (!RuntimeArgs.initial_unit_cell().empty()){
      Initial.store(InitialStateName);
    }

    std::stringstream Rss;
    Rss<<RuntimeArgs.filename()<<"_iTEBD_"<<InitialStateName;

    std::cout << "Using initial state '" << InitialStateName << "'" <<std::endl;

    ajaj::DataOutput results(ajaj::OutputName(Rss.str(),"Evolution.dat"),"Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap)");

    //do we have time dep couplings?
    if (!myModel.times().size()){ //need this for builtin models, with old style coupling params
      std::cout <<"Evolution hamiltonian is static." <<std::endl;
      ajaj::iTEBD infrun(myModel.H_MPO,Initial,time_step_size,results,InitialStateName,trotter_order);
      infrun.evolve(number_of_time_steps,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
      return 0;
    }
    else {

      std::cout <<"Evolution hamiltonian is time dependent." <<std::endl;
      //if ramp step size is smaller than step size, then use that until ramp over (check each time)
      //if step size smaller than ramp step size then use size that is commensurate with ramp step, but smaller than step size.
      ajaj::uMPXInt ramp_step=1;
      //do the first explicitly in order to creat the TEBD object
      double ramp_step_size_1=myModel.times()[ramp_step]-myModel.times()[ramp_step-1];
      double current_step_size_1=ramp_step_size_1;
      ajaj::uMPXInt num_1=1;
      while (current_step_size_1>time_step_size){current_step_size_1=ramp_step_size_1/(++num_1);}
      if (num_1>number_of_time_steps) num_1=number_of_time_steps;
      ajaj::iTEBD infrun(myModel.H_MPO,Initial,current_step_size_1,results,InitialStateName,trotter_order);
      infrun.evolve(num_1,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
      number_of_time_steps-=num_1;
      ++ramp_step;

      //now continue with time dep params
      while (number_of_time_steps>0 && ramp_step < myModel.times().size()){
	double ramp_step_size=myModel.times()[ramp_step]-myModel.times()[ramp_step-1];
	double current_step_size=ramp_step_size;
	ajaj::uMPXInt num=1;
	while (current_step_size>time_step_size){current_step_size=ramp_step_size/(++num);}
	if (num>number_of_time_steps) num=number_of_time_steps;

	infrun.change_bond_operator(myModel.change_H_MPO(myModel.coupling_arrays()[ramp_step-1]),current_step_size);
	infrun.evolve(num,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
	number_of_time_steps-=num;
	++ramp_step;
      }
  
      if (number_of_time_steps>0){
	std::cout <<"End of time dependent hamiltonian stage." <<std::endl;
	infrun.change_bond_operator(myModel.change_H_MPO(myModel.coupling_arrays()[ramp_step-1]),time_step_size);
	infrun.evolve(number_of_time_steps,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
      }

      return 0;

    }
  }
  else {
    return 1;
  }
}


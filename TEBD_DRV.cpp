/**
 *@file TEBD_DRV.cpp Driver file for TEBD.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>

#include "ajaj_common.hpp"
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
    std::vector<ajaj::MPO_matrix> measured_operators(1,myModel.vertex.make_one_site_operator(1)); //for Ising this is the fermion occupation number on a chain
    ajaj::DataOutput results("Evolution.dat");
    std::vector<ajaj::MultiVertexMeasurement> measurements;
    measurements.emplace_back(ajaj::MultiVertexMeasurement(number_of_vertices/2/*first measurement position*/ ,&measured_operators[0],number_of_vertices/2+1/* second measurement position*/,&measured_operators[0]));
    //TEBD stuff
    if (RuntimeArgs.initial_state_name()!=""){
      if (ajaj::CheckMPSFilesExist(RuntimeArgs.initial_state_name(),number_of_vertices)){
	ajaj::TEBD finrun(myModel.H_MPO,RuntimeArgs.initial_state_name(),number_of_vertices,time_step,results,trotter_order);
	finrun.evolve(number_of_time_steps,measurements,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
      }
      else {
	return 1;
      }
    }
    else {
      ajaj::SparseMatrix SpM(myModel.basis().size(),1,1);
      SpM.entry(0,0,1.0);
      SpM.cheap_finalise();
      ajaj::MPS_matrix A(myModel.basis(),std::vector<ajaj::MPXIndex>({{ajaj::MPXIndex(1,myModel.basis()),ajaj::MPXIndex(1,ajaj::StateArray(1,myModel.basis()[0])),ajaj::MPXIndex(0,ajaj::StateArray(1,myModel.basis()[0]))}}),SpM);
      ajaj::TEBD finrun(myModel.H_MPO,std::string("TEBDState"),A,number_of_vertices,time_step,results,trotter_order); //use the provided MPS_matrix A, give state the name "TEBDState"
      finrun.evolve(number_of_time_steps,measurements,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
    }
    return 0;
  }
  else {
    return 1;
  }

}

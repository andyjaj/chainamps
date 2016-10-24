/**
 *@file iTEBD_DRV.cpp Driver file for iTEBD.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>

#include "ajaj_common.hpp"
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
    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    myModel.basis().print();
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(1.0e-14);
    ajaj::uMPXInt number_of_time_steps(RuntimeArgs.number_of_steps());
    ajaj::uMPXInt trotter_order(RuntimeArgs.trotter_order());
    ajaj::uMPXInt measurement_interval(RuntimeArgs.measurement_interval());
    double time_step(RuntimeArgs.step_size());
    std::vector<ajaj::MPO_matrix> measured_operators;// can insert runtime measurements here
    ajaj::DataOutput results(ajaj::OutputName(RuntimeArgs.filename(),"Evolution.dat"),"Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap)");

    //symmetric superposition starting state
    /*ajaj::StateArray sa(1,myvertex.Spectrum[0]);
      ajaj::SparseMatrix SpM(myvertex.Spectrum.size(),1,2);
      SpM.entry(0,0,sqrt(0.5));
      SpM.entry(1,0,sqrt(0.5));
      SpM.cheap_finalise();
      ajaj::MPS_matrix A(myvertex.Spectrum,std::vector<ajaj::MPXIndex>({{ajaj::MPXIndex(1,myvertex.Spectrum),ajaj::MPXIndex(1,sa),ajaj::MPXIndex(0,sa)}}),SpM);
      ajaj::UnitCell Superposition(myvertex.Spectrum);
      Superposition.Matrices.push_back(A);
      Superposition.Matrices.push_back(A);
      Superposition.Lambdas.push_back(std::vector<double>({1.0}));
      Superposition.Lambdas.push_back(std::vector<double>({1.0}));
      ajaj::iTEBD infrun(H,Superposition,time_step,results,trotter_order);*/

    //starting state is product state of ground states
    const ajaj::UnitCell Initial(MakeProductStateUnitCell(myModel.basis(),0/*Currently this should be set to the index of a state that has zero for all charges*/));
    ajaj::iTEBD infrun(myModel.H_MPO,Initial,time_step,results,trotter_order);
    infrun.evolve(number_of_time_steps,measured_operators,CHI/*bond dimension*/,minS/*min s val*/,measurement_interval);
    return 0;
  }
  else {
    return 1;
  }
}


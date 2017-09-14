/**
 *@file fDMRG_DRV.cpp Driver file for finite size DMRG.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
#include <limits>
#include <chrono>

#include "common_defs.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "DMRG_routines.hpp"
#include "measurement.hpp"
#include "data.hpp"
#include "command_line_input.hpp"
#include "model.hpp"
#include "make_model.hpp"

int main(int argc, char** argv){
  ajaj::fDMRG_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));
    myModel.basis().print();

    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(SPARSETOL);
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());
    ajaj::uMPXInt number_of_excited_states(RuntimeArgs.num_excited());
    ajaj::uMPXInt number_of_finite_vol_sweeps(RuntimeArgs.num_sweeps());
    double weight_factor(RuntimeArgs.weight_factor());

    ajaj::State TargetState(myModel.make_target(RuntimeArgs.target()));
    ajaj::DataOutput results(ajaj::OutputName(RuntimeArgs.filename(),"Energies.dat"),"Index, Energy, Energy/vertices, Entropy, Truncation, 1-Fidelity"); //open file for output
    ajaj::iDMRG infvol(std::string("GroundState"),myModel.H_MPO,TargetState,results); //create an infinite sweep simulation object
    auto t1 = std::chrono::high_resolution_clock::now();
    infvol.run(number_of_vertices/2,-1.0,CHI,minS); //convergence criterion: not used if negative, number of vertices used instead
    auto t2 = std::chrono::high_resolution_clock::now();

    ajaj::FiniteDMRG fvol(infvol,results); //create a finite sweep object from the output of infinite volume algorithm
    auto t3 = std::chrono::high_resolution_clock::now();
    fvol.run(number_of_finite_vol_sweeps,CHI,minS); //run for a set number of sweeps
    auto t4 = std::chrono::high_resolution_clock::now();
    if (number_of_excited_states>0){
      ajaj::DataOutput ex_results(ajaj::OutputName(RuntimeArgs.filename(),"Excited_Energies.dat"),"Sweep Index, State Index, Energy, Entropy, Truncation, 1-abs(guess overlap)"); //open file for output
      
      if (number_of_finite_vol_sweeps==0){
	fvol.run(1,CHI,minS); //do a finite vol sweep to finish off.
      }
      ajaj::ExcitedStateFiniteDMRG Exfvol(std::string("Excited"),fvol,weight_factor,ex_results);
      Exfvol.run(number_of_finite_vol_sweeps,CHI,minS);
      for (ajaj::uMPXInt l=1;l<number_of_excited_states;++l){
	Exfvol.next_state(weight_factor);
	Exfvol.run(number_of_finite_vol_sweeps,CHI,minS);
      }
    }
    

    std::cout << "Infinite Volume part took "<< std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<< " milliseconds" << std::endl;
    std::cout << "Finite Volume part took "<< std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count()<< " milliseconds" << std::endl;

    return 0;
  }
  else {
    return 1;
  }
}

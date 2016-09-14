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

#include "ajaj_common.hpp"
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
    double minS(1.0e-14);
    ajaj::uMPXInt number_of_vertices(RuntimeArgs.num_vertices());
    ajaj::uMPXInt number_of_excited_states(RuntimeArgs.num_excited());
    ajaj::uMPXInt number_of_finite_vol_sweeps(RuntimeArgs.num_sweeps());
    double weight_factor(RuntimeArgs.weight_factor());

    //set precision used in printing output
    //std::cout << std::setprecision(16);
    //ajaj::VertexParameterArray params;
    //ajaj::VertexParameterArray couplings;
    //#include "params.hpp"
    //now generate the vertex using definitions in ./vertexdefs/vertex_generator.hpp file
    //static const ajaj::Vertex myvertex(ajaj::VertexGenerator(params));
    //myvertex.Spectrum.print(); //print out spectrum
    //generate a Hamiltonian MPO using ./vertexdefs/vertex_generator.hpp file
    //const ajaj::MPO_matrix H(ajaj::MakeHamiltonian(myvertex,couplings));
    ajaj::State TargetState(myModel.basis().getChargeRules()); //state with all charges=0
    //ajaj::State TargetState(myvertex.Spectrum[0]); //target quantum numbers of the chain ground state
    ajaj::DataOutput results("Energies.dat"); //open file for output
    ajaj::iDMRG infvol(std::string("GroundState"),myModel.H_MPO,TargetState,results); //create an infinite sweep simulation object
    auto t1 = std::chrono::high_resolution_clock::now();
    infvol.run(number_of_vertices/2,-1.0,CHI,minS); //convergence criterion: not used if negative, number of vertices used instead
    auto t2 = std::chrono::high_resolution_clock::now();
    ajaj::FiniteDMRG fvol(infvol,results); //create a finite sweep object from the output of infinite volume algorithm
    auto t3 = std::chrono::high_resolution_clock::now();
    fvol.run(number_of_finite_vol_sweeps,CHI,minS); //run for a set number of sweeps
    auto t4 = std::chrono::high_resolution_clock::now();
    ajaj::DataOutput ex_results("Excited_Energies.dat"); //open file for output
    if (number_of_excited_states>0){
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

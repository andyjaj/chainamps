/**
 *@file iDMRG_DRV.cpp Driver file for iDMRG.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>

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
  ajaj::iDMRG_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){

    static double convergence_test=-0.0;
    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));

    ajaj::State TargetState(myModel.make_target(RuntimeArgs.target()));
    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double trunc(RuntimeArgs.trunc());
    ajaj::uMPXInt steps(RuntimeArgs.number_of_steps());

    std::string PreName("");
    
    if (RuntimeArgs.resume()){
      std::cout << "Resuming a previous iDMRG run at step " << RuntimeArgs.resume() <<std::endl;
      std::stringstream ns;
      ns << "Resumed_Step_" <<  RuntimeArgs.resume()<< "_";
      PreName=ns.str();
    }

    ajaj::DataOutput results(ajaj::OutputName(PreName+RuntimeArgs.filename(),"Energies.dat"),"Index, Energy, Energy/vertex, Entropy, Truncation, Fidelity",RuntimeArgs.resume());

    
    ajaj::iDMRG infvol(std::string("GroundState"),myModel.H_MPO,TargetState,results,2*RuntimeArgs.resume());

    //iDMRG requires at least two previous growth steps. If less than two then we just do an infinite volume step and skip the orthogonalisation bits later.
    ajaj::uMPXInt warmup_steps=steps+RuntimeArgs.resume() <= 2 ? steps : (RuntimeArgs.resume() <2 ? 2-RuntimeArgs.resume() : 0);

    std::cout << "Need to do " << warmup_steps << " warmup steps." <<std::endl; 
    
    infvol.run(warmup_steps,convergence_test,CHI,trunc);
    //infvol.run(steps+RuntimeArgs.resume() > 2 ? 2 : steps ,convergence_test,CHI,trunc);
    
    const ajaj::MPO_matrix H1(myModel.vertex.make_one_site_operator("Vertex_Hamiltonian")); //form the on-vertex part of a Hamiltonian
    const ajaj::MPO_matrix ColX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(1).size()-2),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0)));
    const ajaj::MPO_matrix RowX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(myModel.H_MPO.Index(1).size()-1,myModel.H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(3).size()-2)));

    std::vector<double> iDMRGEnergies;

    ajaj::DataOutput infvolresults(PreName+ajaj::OutputName(RuntimeArgs.filename(),"iDMRGEnergies.dat"),"Index, Energy/vertex, Entropy",RuntimeArgs.resume()>2 ? RuntimeArgs.resume()-2 : 0);

    ajaj::uMPXInt VarCHI(CHI);
    
    //for (ajaj::uMPXInt r=0;r< (steps>2 ? steps-2 : 0) ;++r){
    for (ajaj::uMPXInt r=0;r< steps-warmup_steps ;++r){
      infvol.run(1,-0.0,VarCHI,trunc);

#ifdef TIMING
      ajaj::dmrg_print_time_info();
#endif

      ajaj::UnitCell Ortho(Orthogonalise(infvol.getCentralDecomposition(),infvol.getPreviousLambda()));
      if(Ortho.size()) Ortho.store("Ortho",r); //don't store if the unitcell couldn't be formed
      else {
	std::cout << "Couldn't form unit cell, aborting." <<std::endl;
	exit(1);
      }
      iDMRGEnergies.push_back(real(ajaj::iTwoVertexEnergy(ColX,RowX,H1,Ortho)));

      std::cout << "iDMRG energy per vertex " << iDMRGEnergies.back() << std::endl;
      ajaj::Data inf_data(std::vector<double>({{iDMRGEnergies.back(),Ortho.Entropy()}}));
      infvolresults.push(inf_data);
      Ortho.OutputOneVertexDensityMatrix("OneVertexRho",r);
    }

    std::cout << "SUMMARY of iDMRG thermodynamic energy per vertex" << std::endl;
    std::cout << std::setprecision(16);
    for (auto iE : iDMRGEnergies){
      std::cout << iE << std::endl;
    }    
    return 0;
  }
  return 1;
}


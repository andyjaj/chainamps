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
    ajaj::DataOutput results(ajaj::OutputName(RuntimeArgs.filename(),"Energies.dat"),"Index, Energy, Energy/vertex, Entropy, Truncation, Fidelity");
    ajaj::iDMRG infvol(std::string("GroundState"),myModel.H_MPO,TargetState,results);

    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(0.0);
    ajaj::uMPXInt steps(RuntimeArgs.number_of_steps());

    infvol.run(steps > 2 ? 2 : steps ,convergence_test,CHI,minS);
    //done for at least "steps" number of steps
    //now do orthogs and averaging
    //left orthogonalise the unit cell and produce the new singular values for the infinite state
    //const ajaj::MPXIndex dummy(1,ajaj::StateArray(1,myModel.basis().getChargeRules())); 
    const ajaj::MPO_matrix H1(myModel.vertex.make_one_site_operator("Vertex_Hamiltonian")); //form the on-vertex part of a Hamiltonian
    const ajaj::MPO_matrix ColX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(1).size()-2),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0)));
    const ajaj::MPO_matrix RowX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(myModel.H_MPO.Index(1).size()-1,myModel.H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(3).size()-2)));

    std::vector<double> iDMRGEnergies;

    ajaj::DataOutput infvolresults(ajaj::OutputName(RuntimeArgs.filename(),"iDMRGEnergies.dat"),"Index, Energy/vertex, Entropy");

    ajaj::uMPXInt VarCHI(CHI);

    for (ajaj::uMPXInt r=0;r< (steps>2 ? steps-2 : 0) ;++r){
      infvol.run(1,-0.0,VarCHI,minS);

#ifdef TIMING
      ajaj::dmrg_print_time_info();
#endif

      ajaj::UnitCell Ortho(Orthogonalise(infvol.getCentralDecomposition(),infvol.getPreviousLambda()));
      if(Ortho.size()) Ortho.store("Ortho",r); //don't store if the unitcell couldn't be formed
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


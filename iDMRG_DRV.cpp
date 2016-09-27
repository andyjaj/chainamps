/**
 *@file iDMRG_DRV.cpp Driver file for iDMRG.
 */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <numeric>
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
  ajaj::iDMRG_Args RuntimeArgs(argc,argv);
  if (RuntimeArgs.is_valid()){
    static double convergence_test=-0.0;
    const ajaj::Model myModel(ajaj::MakeModelFromArgs(RuntimeArgs));

    ajaj::State TargetState(myModel.basis().getChargeRules());
    ajaj::DataOutput results("Energies.dat");
    ajaj::iDMRG infvol(std::string("GroundState"),myModel.H_MPO,TargetState,results);

    ajaj::uMPXInt CHI(RuntimeArgs.chi());
    double minS(1.0e-14);
    ajaj::uMPXInt steps(RuntimeArgs.number_of_steps());

    auto t1 = std::chrono::high_resolution_clock::now();
    infvol.run(steps > 2 ? 2 : steps ,convergence_test,CHI,minS);
    //done for at least "steps" number of steps
    //now do orthogs and averaging
    //left orthogonalise the unit cell and produce the new singular values for the infinite state
    const ajaj::MPXIndex dummy(1,ajaj::StateArray(1,myModel.basis().getChargeRules())); 
    const ajaj::MPO_matrix H1(myModel.vertex.make_one_site_operator("Vertex_Hamiltonian")); //form the on-vertex part of a Hamiltonian
    //const ajaj::MPO_matrix I(ajaj::IdentityMPO_matrix(myModel.basis())); //form an identity MPO. Useful for some measurements
    //const ajaj::MPO_matrix LeftH(LeftOpenBCHamiltonian(myModel.H_MPO));
    //const ajaj::MPO_matrix RightH(RightOpenBCHamiltonian(myModel.H_MPO));

    const ajaj::MPO_matrix ColX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(1).size()-2),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0)));
    const ajaj::MPO_matrix RowX(myModel.H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(myModel.H_MPO.Index(1).size()-1,myModel.H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,myModel.H_MPO.Index(3).size()-2)));

    //std::vector<double> iDMRGEnergies_simple;
    std::vector<double> iDMRGEnergies_accurate;

    ajaj::DataOutput infvolresults("iDMRGEnergies.dat");
    for (ajaj::uMPXInt r=0;r< (steps>2 ? steps-2 : 0) ;++r){
      infvol.run(1,-0.0,CHI,minS);
      ajaj::UnitCell Ortho(Orthogonalise(infvol.getCentralDecomposition(),infvol.getPreviousLambda()));
      //at the moment measuring the energy is done in an inelegant way
      //iDMRGEnergies_simple.push_back(real(ajaj::SimpleEnergy(LeftH,RightH,H1,I,Ortho)));
      iDMRGEnergies_accurate.push_back(real(ajaj::SophisticatedEnergy(ColX,RowX,H1,Ortho)));

      std::cout << "iDMRG energy per vertex " << iDMRGEnergies_accurate.back() << std::endl;
      ajaj::Data inf_data(std::vector<double>({{iDMRGEnergies_accurate.back(),ajaj::entropy(*(Ortho.Lambdas.begin()))}}));
      infvolresults.push(inf_data);
      {
	std::stringstream DensityMatrixNameStream;
	DensityMatrixNameStream << "DensityMatrix_" << r << ".dat";
	std::ofstream DensityMatrixFileStream;
	DensityMatrixFileStream.open(DensityMatrixNameStream.str().c_str(),ios::out | ios::trunc);
	Ortho.OutputOneVertexDensityMatrix(DensityMatrixFileStream);
	DensityMatrixFileStream.close();
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "SUMMARY of iDMRG thermodynamic energy per vertex" << std::endl;
    std::cout << std::setprecision(16);
    for (auto iE : iDMRGEnergies_accurate){
      std::cout << iE << std::endl;
    }
    std::cout << "Run took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " milliseconds" << std::endl;
    
    return 0;
  }
  return 1;
}


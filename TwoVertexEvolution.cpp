#include <sstream>

#include "states.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"
#include "TwoVertexEvolution.hpp"

namespace ajaj{

  SparseHED BondHamiltonianDecomp(const MPO_matrix& HMPO,FiniteMPS& F){
     if (F.size()==2) {
      MPXIndex RightIndex=F.matrix(2,1).Index(2);
      MPXIndex LeftIndex=F.matrix(1,1).Index(1);
      std::cout << LeftIndex.size() << " " << RightIndex.size() <<std::endl;
      if (RightIndex.size()==1 && LeftIndex.size()==1){
	std::cout <<"Two vertex only special case" <<std::endl;
	return MakeBondHamiltonian(HMPO).Eigs(RightIndex.at(0)-LeftIndex.at(0));
      }
    }
    return MakeBondHamiltonian(HMPO).Eigs();
  }
  
  TwoVE::TwoVE(const MPO_matrix& HMPO, FiniteMPS& F, double time_step_size, DataOutput& results) : HMPO_(HMPO),F_(F),InitialMPSName_(F_.name()),CurrentMPSName_(make_evolving_name()),InitialWeight_(F_.makeLC(CurrentMPSName_)),CurrentStep_(0),StepSize_(time_step_size),Results_(results),HDecomp_(BondHamiltonianDecomp(HMPO_,F_)) {

    std::cout << "Initial state weight was " << InitialWeight_ <<std::endl;

    if (InitialWeight_!=0.0)
      GoodInitial_=1;
    else
      GoodInitial_=0;
    //get blockstate

    std::cout << "Number of eigenvectors/values found: " << HDecomp_.ValuesSize() <<std::endl;
    
  }
  
  void TwoVE::evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements){
    
  }

  std::string TwoVE::make_evolving_name(){
    std::stringstream Evolvingnamestream;
    Evolvingnamestream << "Evolving_" << InitialMPSName_;
    return Evolvingnamestream.str();
  }

  
}

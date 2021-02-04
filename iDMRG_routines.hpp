#ifndef IDMRGROUTINES_H
#define IDMRGROUTINES_H

#include "DMRG_routines.hpp" //MPX_matrix etc.

namespace ajaj {

  class iDMRG : public SuperBlock {
  private:
    DataOutput& output_ref_;
  public:
    //iDMRG(const std::string& Name, const MPO_matrix& H, const State& TargetState, DataOutput& resultsref) : SuperBlock(Name,H,TargetState),output_ref_(resultsref) {};
    iDMRG(const std::string& Name, const MPO_matrix& H, const State& TargetState, DataOutput& resultsref, uMPXInt num_vertices=0) : SuperBlock(Name,H,TargetState,num_vertices,num_vertices >2 ? num_vertices/2-1 : (num_vertices==2 ? 1 : 0), num_vertices > 2 ? 2 : 0),output_ref_(resultsref) {}
    iDMRG(SuperBlock& SB, DataOutput& resultsref) : SuperBlock(std::move(SB)),output_ref_(resultsref) {}
    //note above the special case of 2 vertices, where the blocks for the next step have been precalculated, and need to be fetched (as they are expected to be in storage).
    
    void run(uMPXInt number_of_steps=0, double convergence_criterion=0.0,  uMPXInt chi=0, double truncation=0.0); /**< Perform infinite algorithm growth steps*/

    double fidelity() const {return fidelity_;}
  };


}
#endif

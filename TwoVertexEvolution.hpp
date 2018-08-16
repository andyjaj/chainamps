/** @file TwoVertexOnly.hpp
 * Routines relating to the special case of systems with only two vertices
 * Here the MPS machinery is largely redundant, and we want to reduce to basic TSA
 */

#ifndef TwoVertex_H
#define TwoVertex_H

#include "common_defs.hpp"
#include "states.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"

namespace ajaj {

  class BondHamiltonianDecomp{
  private:
    const MPO_matrix* Hptr_;
  public:
    BondHamiltonianDecomp(const MPO_matrix& HMPO) : Hptr_(&HMPO) {}
  };
  
  class TwoVE{
  private:
    uMPXInt CurrentStep_;
    double StepSize_;
    DataOutput& Results_;
    const std::string MPSName_;
    const Basis& Basis_;

    BondHamiltonianDecomp HDecomp_;
    
    bool GoodInitial_;
   
    
  public:
    TwoVE(const MPO_matrix& H, FiniteMPS& F, double time_step_size, DataOutput& results) : CurrentStep_(0),StepSize_(time_step_size),Results_(results),MPSName_(F.name()),Basis_(H.basis()),HDecomp_(H),GoodInitial_(0) {}

    double time_step_size() const {return StepSize_;}
    
    void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements){}// at the end of evolve, the MPS is actually updated

    
  };
  
}

#endif

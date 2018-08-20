/** @file TwoVertexOnly.hpp
 * Routines relating to the special case of systems with only two vertices
 * Here the MPS machinery is largely redundant, and we want to reduce to basic TSA
 */

#ifndef TwoVertex_H
#define TwoVertex_H

#include <utility>
#include <string>


#include "common_defs.hpp"
#include "states.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"
#include "TEBD_routines.hpp"

namespace ajaj {

  SparseHED BondHamiltonianDecomp(const MPO_matrix& HMPO,FiniteMPS& F);
  
  class TwoVE{
  private:
    const MPO_matrix& HMPO_;
    FiniteMPS& F_;
    const std::string InitialMPSName_;
    std::string CurrentMPSName_;
    const std::complex<double> InitialWeight_;

    uMPXInt CurrentStep_;  
    double StepSize_;
    DataOutput& Results_;
    SparseHED HDecomp_;   
    bool GoodInitial_;


    std::string make_evolving_name();

    
    
  public:
    TwoVE(const MPO_matrix& HMPO, FiniteMPS& F, double time_step_size, DataOutput& results);
    void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements);// at the end of evolve, the MPS is actually updated
    double time_step_size() const {return StepSize_;} 

    
  };
  
}

#endif

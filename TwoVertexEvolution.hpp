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
  
  MPX_matrix Make2VBondHamiltonian(const MPO_matrix& HMPO); //Make the two vertex hamiltonian
  MPX_matrix Make2VEvolutionOperator(const MPX_matrix& BondH, double timestep,const State& blockstate);

  class TwoVE : public TimeBase {
  private:
    const std::string MPSName_;
    const EigenStateArray& Basis_;
    MPX_matrix EvolutionOperator_;
    bool GoodInitial_;
    State BlockState_;
    
    void apply();
    void left_canonise_measure(std::vector<MultiVertexMeasurement>& measurement);
    
  public:
    TwoVE(const MPO_matrix& HMPO, FiniteMPS& F, double time_step_size, DataOutput& results, const State& blockstate); //use FiniteMPS class   
    void change_bond_operator(const MPO_matrix& HMPO, double time_step_size);   
    void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements);
    bool good() const {return GoodInitial_;}
  };
  
}

#endif

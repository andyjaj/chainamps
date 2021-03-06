/** @file TEBD_routines.hpp
 * Driver functions for TEBD
 *
 */
#ifndef TEBD_H
#define TEBD_H

#include <vector>
#include <utility>


#include "common_defs.hpp"
#include "states.hpp"
#include "vertex.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"


namespace ajaj {

  class TrotterDecomposition{
  private:
    const MPO_matrix* m_H_ptr;
    double m_time_step_size;
    uMPXInt m_order;
    std::vector<MPX_matrix> BondOperators;
    //std::vector<MPO_matrix*> Measurement_Operator_ptrs;
  public:
    std::vector<const MPX_matrix*> OrderedOperatorPtrs;
    TrotterDecomposition(const MPO_matrix& H,double time_step_size,uMPXInt order=1);
    uMPXInt order() const {return m_order;}
    double time_step_size() const {return m_time_step_size;}
  };

  class SeparatedTrotterDecomposition{
  private:
    static constexpr double BONDDECOMPTOL=1.0e-15;//needs setting externally
    const MPO_matrix& H_;
    double TimeStepSize_;
    uMPXInt Order_;
    std::vector<MPX_matrix> UOperators_;
    std::vector<MPX_matrix> UbarOperators_;
    std::vector<uMPXInt> OperatorOrder_;
  public:
    const MPX_matrix& getU(uMPXInt i) const {return UOperators_.at(OperatorOrder_.at(i));}
    const MPX_matrix& getUbar(uMPXInt i) const {return UbarOperators_.at(OperatorOrder_.at(i));}

    SeparatedTrotterDecomposition(const MPO_matrix& H,double time_step_size,uMPXInt order=1);
    uMPXInt order() const {return Order_;}
    double time_step_size() const {return TimeStepSize_;}
  };

  class TimeBase {
  protected:
    uMPXInt m_current_time_step;
    double m_time_step_size;
    DataOutput& m_results;
    double m_truncation;
    double m_current_time;
    double update_time() {++m_current_time_step; return m_current_time+=m_time_step_size;}
  public:
    TimeBase(double time_step_size, DataOutput& results) : m_current_time_step(0),m_time_step_size(time_step_size),m_results(results),m_truncation(0.0),m_current_time(0.0){}

    double time_step_size() const {return m_time_step_size;}
    double current_time() const {return m_current_time;}
  };

  class iTEBD : public TimeBase {
  private:
    TrotterDecomposition m_EvolutionOperators;
    const UnitCell& m_initial_unit;
    UnitCell m_unit;
    std::string Name_;
    const UnitCell& apply_and_decompose(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
    void do_measurements(const UnitCell& ortho, const std::vector<MPO_matrix>& measuredMPOs);
  public:
    iTEBD(const MPO_matrix& H,const UnitCell& C, double time_step_size, DataOutput& results, const std::string& Name, uMPXInt order=1);
    const UnitCell& evolve(uMPXInt num_steps, const std::vector<MPO_matrix>& measuredMPOs, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);
    uMPXInt order() const {return m_EvolutionOperators.order();}
    void change_bond_operator(const MPO_matrix& H, double time_step_size);

  };

class TEBD : public TimeBase {
private:
  const std::string MPSName_;
  const EigenStateArray& Basis_;
  const uMPXInt NumVertices_;
  MPO_matrix SingleVertexOp_; //for open boundary conditions
  TrotterDecomposition m_EvolutionOperators;
  bool GoodInitial_;

  void apply_to_odd_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  void apply_to_even_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  void left_canonise(uMPXInt chi=0,double minS=0);
  void left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements,uMPXInt chi=0,double minS=0,bool overlap_requested=1);
  void left_canonise_measure_special(std::vector<MultiVertexMeasurement>& measurements, uMPXInt Index=0);

  
  double max_truncation_=0.0;

public:
  TEBD(const MPO_matrix& H, FiniteMPS& F, DataOutput& results);
  TEBD(const MPO_matrix& H, FiniteMPS& F, double time_step_size, DataOutput& results, uMPXInt order=1); //use FiniteMPS class

  void change_bond_operator(const MPO_matrix& H, double time_step_size);

  void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);
  void left_info();
  void right_info();
  bool good() const {return GoodInitial_;}
};

  MPX_matrix MakeBondHamiltonian(const MPO_matrix& H); //Returns part of hamiltonian that refers to a single bond (no double counting of vertex part).
  MPX_matrix MakeBondEvolutionOperator(const MPX_matrix& BondH, double timestep);//Forms exponential bond update operator
  MPO_matrix MakeSingleSiteEvolutionOperatorFromLowTriMPO(const MPO_matrix& H_MPO, double timestep);

}

#endif

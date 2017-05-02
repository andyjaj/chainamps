/** @file TEBD_routines.hpp
 * Driver functions for TEBD
 *
 */
#ifndef TEBD_H
#define TEBD_H

#include <vector>
#include <utility>


#include "ajaj_common.hpp"
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
  public:
    TimeBase(double time_step_size, DataOutput& results) : m_current_time_step(0),m_time_step_size(time_step_size),m_results(results),m_truncation(0.0){}

    //uMPXInt order() const {return m_EvolutionOperators.order();}
    double time_step_size() const {return m_time_step_size;}
    double current_time() const {return m_current_time_step*m_time_step_size;}
  };

  class iTEBD : public TimeBase {
  private:
    const TrotterDecomposition m_EvolutionOperators;
    const UnitCell& m_initial_unit;
    UnitCell m_unit;
    std::string Name_;
    const UnitCell& apply_and_decompose(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
    void do_measurements(const UnitCell& ortho, const std::vector<MPO_matrix>& measuredMPOs);
  public:
    iTEBD(const MPO_matrix& H,const UnitCell& C, double time_step_size, DataOutput& results, const std::string& Name, uMPXInt order=1);
    const UnitCell& evolve(uMPXInt num_steps, const std::vector<MPO_matrix>& measuredMPOs, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);
    uMPXInt order() const {return m_EvolutionOperators.order();}
  };

class TEBD : public TimeBase {
private:
  const std::string MPSName_;
  const EigenStateArray& Basis_;
  const uMPXInt NumVertices_;
  const MPX_matrix SingleVertexOp_; //for open boundary conditions
  const TrotterDecomposition m_EvolutionOperators;
  void apply_to_odd_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  void apply_to_even_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  void left_canonise(uMPXInt chi=0,double minS=0);
  void left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements,uMPXInt chi=0,double minS=0);

  double max_truncation_=0.0;

public:
  TEBD(const MPO_matrix& H, const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //load from file
  TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //use repeating simple MPS_matrix defined by InitialMPS_matrix
  void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);
  void left_info();
  void right_info();
};

  UnitCell MakeProductStateUnitCell(const Basis& basis, const std::vector<std::pair<uMPXInt,double> >& state_index_vec, State leftstate, uMPXInt length=2);
  UnitCell MakeProductStateUnitCell(const Basis& basis, uMPXInt state_index, State leftstate, uMPXInt length=2);

  MPX_matrix MakeBondHamiltonian(const MPO_matrix& H); //Returns part of hamiltonian that refers to a single bond (no double counting of vertex part).
  MPX_matrix MakeBondEvolutionOperator(const MPX_matrix& BondH, double timestep);//Forms exponential bond update operator
  MPX_matrix MakeSingleSiteEvolutionOperator(const MPO_matrix& H, double timestep);

  bool CheckMPSFilesExist(const std::string& MPSName,uMPXInt NumVertices);

  //#include "temp_TEBD.hpp"

}

#endif

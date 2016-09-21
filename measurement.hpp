/** @file Measurement.hpp
 * Routines for measuring expectation values, and orthogonalising infinite states.
 *
 */

#ifndef MEASURE_H
#define MEASURE_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <complex>


#include "ajaj_common.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"

namespace ajaj {

  class meas_pair{
  private:
    uMPXInt VertexNumber_;
    const MPO_matrix* OperatorPtr_;
  public:
    meas_pair(uMPXInt vn, const MPO_matrix* op) : VertexNumber_(vn), OperatorPtr_(op) {}
    uMPXInt position() const {return VertexNumber_;}
    const MPO_matrix* MPO_ptr() const {return OperatorPtr_;} 
  };
  
  //class designed to be used for measurement on a single site sweep
  class MultiVertexMeasurement{
  private:
    std::vector<meas_pair> VertexOperatorPtrs_;
    MPX_matrix T_;
    std::complex<double> Result_;

    void link_(const MPO_matrix* Op, const MPX_matrix& A);
    void start_chain_(const MPO_matrix* Op, const MPX_matrix& A);
    void finish_chain_(const MPX_matrix& Lambda);
  public:
    MultiVertexMeasurement(uMPXInt start, const MPO_matrix* Op1Ptr, uMPXInt finish, const MPO_matrix* Op2Ptr) : VertexOperatorPtrs_({{meas_pair(start,Op1Ptr), meas_pair(finish,Op2Ptr)}}), T_(Op1Ptr->GetPhysicalSpectrum()),Result_(0.0) {
      if (start>=finish || start<1){
	std::cout << "Incorrectly defined MultiVertexMeasurement positions, start, finish: " << start << "," << finish <<std::endl;
	std::cout << "start must be >=1 and < finish. finish must be <= Number of Vertices"<<std::endl;
	exit(1);
      }
    }
    //MultiVertexMeasurement(uMPXInt start, uMPXInt finish) : VertexOperatorPtrs_({{meas_pair(start,nullptr), meas_pair(finish,nullptr)}}), Result_(0.0) {}
    //MultiVertexMeasurement() : VertexOperatorPtrs_({{meas_pair(0,nullptr), meas_pair(0,nullptr)}}), Result_(0.0) {}
    const meas_pair& operator[](uMPXInt i) const {return VertexOperatorPtrs_.at(i);}
    uMPXInt start() const {return VertexOperatorPtrs_.begin()->position();}
    uMPXInt finish() const {return VertexOperatorPtrs_.back().position();}
    void update(uMPXInt v, const MPXDecomposition& D) {
      if (v<finish() && v>start()){
	const MPO_matrix* op_p=nullptr;
	for (auto&& vop : VertexOperatorPtrs_){
	  if (v==vop.position()) {op_p=vop.MPO_ptr(); break;}	
	}
	link_(op_p,D.ColumnMatrix);
      }
      else {
	//depending on v, start measurement chain, or continue it
	if (v==start()){
	 start_chain_(VertexOperatorPtrs_.begin()->MPO_ptr(),D.ColumnMatrix); //contract first measurement etc.
	}
	if (v==finish()){
	  //apply and finish off with lambda^2
	  link_(VertexOperatorPtrs_.back().MPO_ptr(),D.ColumnMatrix);
	  finish_chain_(MPX_matrix(D.ColumnMatrix.GetPhysicalSpectrum(),D.ColumnMatrix.Index(2),D.Values));
	}
	else { //v<start() or v>finish()
	  //do nothing
	}
      }
    }
    std::complex<double> result() const {
      return Result_;
    }
  };

  class TransferMatrixParts{
  public:
    const UnitCell& BraCell;
    const UnitCell& KetCell;
    const State* TargetStatePtr;
    std::vector<MPXIndex> left_indices;
    std::vector<MPXIndex> right_indices;

    std::vector<Sparseint> allowed_indices;
    std::vector<Sparseint> allowed_indices_dagger;

    std::vector<std::array<Sparseint,2> > rows_and_cols;
    uMPXInt vrows;
    uMPXInt vcols;

    TransferMatrixParts(const UnitCell& B,const UnitCell& K,const State* T=NULL);
    TransferMatrixParts(const UnitCell& C,const State* T=NULL);
    uMPXInt length(){return m_length;}
    SparseED LeftED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL);
    SparseED RightED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL);

  private:
    uMPXInt m_length;
    void m_init();
  };

  void LeftTransferMatrixMultiply(TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out);
  void RightTransferMatrixMultiply(TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out);

  /** Orthogonalise a new decomposition, A L_n B inverse L_(n-1). Return a left orthogonal unit cell, and new lambda*/
  UnitCell Orthogonalise(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda);
  UnitCell Orthogonalise(const UnitCell& C);
  UnitCell OrthogonaliseInversionSymmetric(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda);
  UnitCell OrthogonaliseInversionSymmetric(const UnitCell& C);

  /** Sandwich two MPO_matrix objects, to make a transfer matrix. Can be used to find thermodynamic limit eigenvalues, if MPO_matrices are simple (NOT LOWER TRIANGULAR) */
  MPX_matrix MakeMeasurementTransferMatrix(const MPO_matrix& W1, const MPO_matrix& W2, const MPS_matrix& A1, const MPS_matrix& A2);
  /** Form a transfer matrix for two possibly different wavefunctions.*/
  MPX_matrix MakeTransferMatrix(const MPS_matrix& A1bra, const MPS_matrix& A2bra, const MPS_matrix& A1ket, const MPS_matrix& A2ket);
  MPX_matrix MakeLTransferMatrix(const UnitCell& bra, const UnitCell& ket);
  void ShiftLTransferMatrixToR(MPX_matrix& LTransferMatrix, const std::vector<double>& Lambda0);
  void ShiftLTransferMatrixToR(MPX_matrix& LTransferMatrix);

  std::complex<double> Overlap(const UnitCell& bra, const UnitCell& ket);

  /** Measure an operator that spans only two vertices. */
  std::complex<double> TwoVertexMeasurement(const MPO_matrix& W1, const MPO_matrix& W2, const MPS_matrix& A1ket, const MPS_matrix& A2ket, const MPX_matrix& Lambda);
  std::complex<double> OneVertexMeasurement(const MPO_matrix& W, const UnitCell& U);

  /** Special routine to measure the energy per vertex, because Hamiltonian is a lower triangular matrix product operator */
  std::complex<double> SimpleEnergy(const MPO_matrix& LeftH,const MPO_matrix& RightH,const MPO_matrix& H1,const MPO_matrix& I,const UnitCell& Ortho);
  std::complex<double> SophisticatedEnergy(const MPO_matrix& ColX,const MPO_matrix& RowX,const MPO_matrix& H1,const UnitCell& Ortho);

}

#endif

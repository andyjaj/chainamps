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
#include <array>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"

namespace ajaj {

  class meas_pair{
  private:
    uMPXInt VertexNumber_;
    const MPO_matrix* OperatorPtr_;
  public:
    meas_pair() : VertexNumber_(0), OperatorPtr_(nullptr) {}
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
    uMPXInt start_; //should be set by find_start_and_finish_()
    uMPXInt finish_; //should be set as above

    void link_(const std::vector<const MPO_matrix*>& ops, const MPX_matrix& A);
    void start_chain_(const MPX_matrix& A);
    void finish_chain_(const MPX_matrix& Lambda);
    std::vector<const MPO_matrix*> get_ops(uMPXInt v);
    void find_start_and_finish_();
    
  public:
    MultiVertexMeasurement(uMPXInt start, const MPO_matrix* OpPtr) : VertexOperatorPtrs_({{meas_pair(start,OpPtr)}}), T_(OpPtr->basis()),Result_(0.0) {find_start_and_finish_();} //special single measurement case
    MultiVertexMeasurement(uMPXInt start, const MPO_matrix* Op1Ptr, uMPXInt finish, const MPO_matrix* Op2Ptr);//two measurement case
    MultiVertexMeasurement(const std::vector<meas_pair>& VOPs) : VertexOperatorPtrs_(VOPs) {find_start_and_finish_();}
    
    const meas_pair& operator[](uMPXInt i) const {return VertexOperatorPtrs_.at(i);}
    uMPXInt start() const {return start_;}//argh!
    uMPXInt finish() const {return finish_;}
    void update(uMPXInt v, const MPXDecomposition& D); //update measurement with next vertex
    std::complex<double> result() const {
      return Result_;
    }
  };

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
  std::vector<std::complex<double> > Overlap(const UnitCell& bra, const UnitCell& ket, uMPXInt nev);
  std::vector<std::complex<double> > TransferMatrixEigs(const UnitCell& ket, uMPXInt nev,const State& S);

  /** Measure an operator that spans only two vertices. */
  std::complex<double> TwoVertexMeasurement(const MPO_matrix& W1, const MPO_matrix& W2, const MPS_matrix& A1ket, const MPS_matrix& A2ket, const MPX_matrix& Lambda);
  std::complex<double> OneVertexMeasurement(const MPO_matrix& W, const UnitCell& U);
  std::complex<double> TwoVertexMeasurement(const MPO_matrix& W1,const MPO_matrix& W2,const UnitCell& U,uMPXInt separation=0);

  /** Special routine to measure the energy per vertex, because Hamiltonian is a lower triangular matrix product operator */
  std::complex<double> SimpleEnergy(const MPO_matrix& LeftH,const MPO_matrix& RightH,const MPO_matrix& H1,const MPO_matrix& I,const UnitCell& Ortho);
  std::complex<double> iTwoVertexEnergy(const MPO_matrix& H_MPO,const UnitCell& Ortho);
  std::complex<double> iTwoVertexEnergy(const MPO_matrix& ColX,const MPO_matrix& RowX,const MPO_matrix& H1,const UnitCell& Ortho);

  class TransferMatrixComponents {
  private:
    const State TargetState_;
    const uMPXInt CellSize_;
    const bool Hermitian_answer_; //setting this true insists that the final eigenvector can be reshaped into a (obviously square) Hermitian matrix, and assumes that dominant e val is real...

    const UnitCell& KetCell_;
    const UnitCell& BraCell_;

    std::vector<std::pair<const MPS_matrix*,const MPS_matrix*> > BraKetMatrixPtrs_;

    uMPXInt length_;
    std::vector<MPXIndex> left_indices_;
    std::vector<MPXIndex> right_indices_;
    std::vector<Sparseint> allowed_indices_;
    std::vector<Sparseint> allowed_indices_dagger_;//useful in order to ensure eigenvector can be made reshaped into Hermitian matrix to high precision
    std::vector<std::array<Sparseint,2> > rows_and_cols_;
    uMPXInt vrows_;
    uMPXInt vcols_;

    MPX_matrix accumulate_() const;

  public:

    //TransferMatrixComponents(const std::vector<const MPS_matrix*>& BraPtrs, const std::vector<const MPS_matrix*>& KetPtrs, bool HV, const State S);
    //TransferMatrixComponents(const std::vector<const MPS_matrix*>& KetPtrs, bool HV, const State S);

    TransferMatrixComponents(const UnitCell& BraCell, const UnitCell& KetCell, bool HV, const State S);
    TransferMatrixComponents(const UnitCell& KetCell, bool HV, const State S);

    const MPS_matrix& BraMatrix(uMPXInt i) const {return *(BraKetMatrixPtrs_.at(i).first);} //public interface
    const MPS_matrix& KetMatrix(uMPXInt i) const {return *(BraKetMatrixPtrs_.at(i).second);}
    uMPXInt last() const {return CellSize_-1;}
    uMPXInt length() const {return length_;}

    SparseED LeftED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;
    SparseED RightED(Sparseint numevals, char which[3],SparseMatrix* initial=NULL) const;

    const std::vector<Sparseint>& allowed_indices() const {return allowed_indices_;}
    const std::vector<std::array<Sparseint,2> >& rows_and_cols() const {return rows_and_cols_;}
    uMPXInt vrows() const {return vrows_;}
    uMPXInt vcols() const {return vcols_;}
    const std::vector<MPXIndex>& left_indices() const {return left_indices_;}
    const std::vector<MPXIndex>& right_indices() const {return right_indices_;}

  };

  void LeftComponentsMultiply(const TransferMatrixComponents* stuff, std::complex<double> *in, std::complex<double> *out);
  void RightComponentsMultiply(const TransferMatrixComponents* stuff, std::complex<double> *in, std::complex<double> *out);

  double MultiVertexEntropy(const UnitCell& U,uMPXInt v=1);

  bool check_for_H_MPO_file(const std::string& name, uMPXInt sliceindex);

  std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Bra, const ConstFiniteMPS& Ket, const std::vector<meas_pair>& ops);
  std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Ket, const std::vector<meas_pair>& ops);

  std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Bra, const ConstFiniteMPS& Ket, const MPO_matrix& Op=MPO_matrix(), uMPXInt v=0);
  
  //inline std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Ket, const MPO_matrix& Op=MPO_matrix(), uMPXInt v=0) {return GeneralisedOverlap(Ket,Ket,Op,v);}

}

#endif

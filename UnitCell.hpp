/** @file UnitCell.hpp
 * Container for a repeating unit of multiple MPS_matrix objects 
 */
#ifndef UNITCELL_H
#define UNITCELL_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <cstdlib>
#include <numeric>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp" //specialisations, MPS

namespace ajaj {

  /** A class to hold a unit cell of arbitrary length
   *
   *
   */
  class UnitCell {
  protected:
    const Basis* basis_ptr_;
  public:
    std::vector<MPS_matrix> Matrices;
    std::vector<std::vector<double> > Lambdas;
    //storage format is 'left canonical' based i.e. A_0 A_1 (or lambda_0 Gamma_0 lambda_1 Gamma_1...)
    //even if the matrices are not in fact left canonical
    //so a decomposition is of the form lambda_0 Gamma_0 lambda_1 Gamma_1 lambda_0
    UnitCell()=delete;
    UnitCell(const Basis& spectrum) : basis_ptr_(&spectrum){}
    UnitCell(const Basis& spectrum,const MPS_matrix& m, const std::vector<double>& l,uMPXInt length=1) : basis_ptr_(&spectrum), Matrices(std::vector<MPS_matrix>(length,m)),Lambdas(std::vector<std::vector<double> >(length,l)) {}
    UnitCell(const Basis& spectrum,std::vector<MPS_matrix>&& mvec, std::vector<std::vector<double> >&& lvec) : basis_ptr_(&spectrum), Matrices(mvec),Lambdas(lvec) {}
    UnitCell(const MPSDecomposition& D,const std::vector<double>& PreviousLambda) : basis_ptr_(&D.LeftMatrix.GetPhysicalSpectrum()) {
      Matrices.emplace_back(copy(D.LeftMatrix));
      Matrices.emplace_back(reorder(contract(contract(MPX_matrix(D.LeftMatrix.GetPhysicalSpectrum(),D.RightMatrix.Index(0),D.Values),0,D.RightMatrix,0,contract10),0,MPX_matrix(D.LeftMatrix.GetPhysicalSpectrum(),D.RightMatrix.Index(2),PreviousLambda,1),0,contract20),0,reorder102,2));
      Lambdas.emplace_back(PreviousLambda);
      Lambdas.emplace_back(D.Values);
    }

    void swap(uMPXInt i, uMPXInt j){MPS_swap(Matrices.at(i),Matrices.at(j));std::swap(Lambdas.at(i),Lambdas.at(j));}
    
    /*UnitCell& operator=(UnitCell&& other){
      basis_ptr_=other.basis_ptr_;
      Matrices=std::move(other.Matrices);
      Lambdas=std::move(other.Lambdas);
      return *this;
      }*/

    uMPXInt size() const {return Matrices.size();}
    const MPS_matrix& GetMatrix(size_t i) {return Matrices.at(i);}
    const std::vector<double>& GetLambda(size_t i) {return Lambdas.at(i);}
    void OutputPhysicalIndexDensities(std::ofstream& d) const;
    void OutputOneVertexDensityMatrix(std::ofstream& d) const;
    void OutputOneVertexDensityMatrix(const std::string& name, uMPXInt l) const;
    
    void store(const std::string& filename, uMPXInt l) const; /**< Print UnitCell with to binary file, with an index in filename. */
    void store(const std::string& filename) const; /**< Print the UnitCell to file in binary format. */

    double Entropy() const;
    const Basis& basis() const {
      return *basis_ptr_;
    }
  };

  UnitCell load_UnitCell_binary(std::ifstream& infile, QNVector& charge_rules, Basis& basis);
  UnitCell load_UnitCell_binary(std::ifstream& infile, const QNVector& charge_rules, const Basis& basis);

  //UnitCell MakeProductStateUnitCell(const Basis& basis, const std::vector<std::pair<uMPXInt,std::complex<double> > >& state_index_vec, State leftstate, uMPXInt length=2);
  UnitCell MakeProductStateUnitCell(const Basis& basis, const c_specifier_array& c_spec, State leftstate, uMPXInt length=2);
  UnitCell MakeProductStateUnitCell(const Basis& basis, uMPXInt state_index, State leftstate, uMPXInt length=2);


}

#endif

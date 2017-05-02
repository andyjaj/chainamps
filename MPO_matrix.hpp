/** @file MPO_matrix.hpp
 * MPO_matrix class, derived from MPX_matrix.
 * Essentially a container for a SparseMatrix, with lots of bells and whistles so that they can be blocked by good quantum numbers, and reshaped according
 * to their indices.
 */
#ifndef MPO_MATRIX_H
#define MPO_MATRIX_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class

namespace ajaj{

  //friends for namespace resolution
  //void MPS_swap(MPS_matrix& A, MPS_matrix& B);

  //classes

  /** Inherited class from MPX_matrix. An MPX_matrix with two matrix indices and two physical indices.*/
  class MPO_matrix : public MPX_matrix {
  private:
    void check();
  public:
    MPO_matrix() noexcept : MPX_matrix() {}; 
    MPO_matrix(const Basis& spectrum);
    MPO_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, const SparseMatrix& matrix);
    MPO_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, SparseMatrix&& matrix);

    //MPO_matrix(const MPX_matrix& MPXref);
    MPO_matrix(MPX_matrix&& MPXref) noexcept;
    MPO_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values,bool inverse=0);
    MPO_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<double>& values,bool inverse=0);

    MPO_matrix ExtractMPOBlock(const std::pair<MPXInt,MPXInt>& row_matrix_index_range, const std::pair<MPXInt,MPXInt>& col_matrix_index_range) const;

  };



  struct NamedMPO_matrix {
  public:
    std::string Name;
    ajaj::MPO_matrix Matrix;
    NamedMPO_matrix(const std::string& n,ajaj::MPO_matrix&& m) : Name(n),Matrix(m) {}
  };

  /** Create an identity MPO_matrix. Corresponds to an identity operator.*/
  inline MPO_matrix IdentityMPO_matrix(const Basis& spectrum){
    MPXIndex dummy(1,ajaj::StateArray(1,spectrum[0]-spectrum[0]));
    return MPO_matrix(spectrum,dummy,std::vector<complex<double> >(spectrum.size(),1.0));  //create diagonal matrix
  }

  MPO_matrix UnitaryTransformMPO_matrix(const Basis&, const std::vector<MPXIndex>&, const SparseMatrix&, size_t, double);/**< To unitarily transform a local operator*/

  MPO_matrix load_MPO_matrix(const std::string& filename,const Basis& spectrum);

}
#endif

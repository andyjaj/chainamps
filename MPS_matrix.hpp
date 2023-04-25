/** @file MPS_matrix.hpp
 * MPS_matrix class, derived from MPX_matrix, and helpers
 * Essentially a container for a SparseMatrix, with lots of bells and whistles so that they can be blocked by good quantum numbers, and reshaped according
 * to their indices.
 */
#ifndef MPS_MATRIX_H
#define MPS_MATRIX_H

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

  enum class MPS_matrixCanonicalType : unsigned short int {Left, Right, Non};

  
  /** Inherited class from MPX_matrix. An MPX_matrix with two matrix indices and one physical index.*/
  class MPS_matrix : public MPX_matrix {
  private:
    void check();
  public:
    MPS_matrix(const Basis& spectrum);
    MPS_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, const SparseMatrix& matrix);
    MPS_matrix(const Basis& spectrum, const std::vector<MPXIndex>& indices, SparseMatrix&& matrix);
    MPS_matrix(MPX_matrix&& MPXref) noexcept;

    bool is_left_shape() const {
      return m_Indices[0].Physical() && m_Indices[1].Ingoing() && m_Indices[2].Outgoing();
    }

    bool is_right_shape() const {
      return m_Indices[1].Physical() && m_Indices[0].Ingoing() && m_Indices[2].Outgoing();
    }

    const MPXIndex& getInwardMatrixIndex() const;
    const MPXIndex& getOutwardMatrixIndex() const;
    const MPXIndex& getPhysicalIndex() const;
    MPXInt InwardMatrixIndexNumber() const;
    MPXInt OutwardMatrixIndexNumber() const;
    MPXInt PhysicalIndexNumber() const;
    MPS_matrix&& left_shape();
    MPS_matrix&& right_shape();

    friend void MPS_swap(MPS_matrix& A, MPS_matrix& B);
  };

  /** Used for storing the results of decompositions, such as SVD (Schmidt decomposition) that produce two new MPS_matrix objects.*/
  class MPSDecomposition {
  private:
  public:
    double Truncation;
    std::vector<double> Values;
    MPS_matrix LeftMatrix;
    MPS_matrix RightMatrix;
    MPSDecomposition(MPXDecomposition&& X) noexcept : Truncation(X.Truncation), Values(std::move(X.Values)), LeftMatrix(std::move(X.ColumnMatrix)), RightMatrix(std::move(X.RowMatrix)){};
    MPSDecomposition(const Basis& basis) : Truncation(0.0), Values(), LeftMatrix(basis), RightMatrix(basis){};
    const std::vector<double>& SquareRescale(double sqsum) {
      Truncation=Truncation/sqsum;
      SquareSumRescale(Values,sqsum);
      return Values;
    }
    void printValues() const {
      for (std::vector<double>::const_iterator cit=Values.begin();cit!=Values.end();++cit){
	std::cout << *cit << std::endl;
      }
    }

    bool store(const std::string& LeftName, const std::string& RightName) const;
    bool store(const std::string& Name, uMPXInt nl,uMPXInt nr,bool StoreLambda=1) const;
    bool store_left(const std::string& Name, uMPXInt nl) const;
    bool store_right(const std::string& Name,uMPXInt nr) const;

    std::string LeftName(const std::string& s, uMPXInt n) const {
      std::stringstream leftname;
      leftname << s.c_str() << "_Left_" << n << ".MPS_matrix";
      return leftname.str();
    }
    std::string RightName(const std::string& s, uMPXInt n) const {
      std::stringstream rightname;
      rightname << s.c_str() << "_Right_" << n << ".MPS_matrix";
      return rightname.str();
    }

    void OutputPhysicalIndexDensities(std::ofstream& d) const;
    const Basis& basis() const {
      return LeftMatrix.basis();
    }

    void OutputOneVertexDensityMatrix(std::ofstream& d) const;

  };

  //////////
  //Helpers
  //////////

  MPS_matrix load_MPS_matrix(const std::string& filename,const Basis& spectrum);
  MPS_matrix MakeProductState(const Basis& spectrum, const c_specifier_vector& state_index_vec,State leftstate);
  MPS_matrix MakeProductState(const Basis& spectrum, uMPXInt state_index,State leftstate);
  
  c_specifier_array LoadCNumbers(const std::string& filename);


}
#endif

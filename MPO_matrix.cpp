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
#include "MPO_matrix.hpp"

namespace ajaj{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MPO members
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void MPO_matrix::check(){
    if (m_Indices.size()!=4/* || m_NumRowIndices!=2*/){std::cout << "Incorrect number of total indices for MPO_matrix" << std::endl; exit(1);}
    MPXInt numphysical=0;
    MPXInt nummatrix=0;
    for (std::vector<MPXIndex>::const_iterator cit=m_Indices.begin();cit!=m_Indices.end();++cit){
      if (cit->Physical()){++numphysical;}
      else {++nummatrix;}
    }
    if (numphysical!=2){std::cout << "Incorrect number of physical indices for MPO_matrix: " << numphysical << std::endl; exit(1);}
    if (nummatrix!=2){std::cout << "Incorrect number of matrix indices for MPO_matrix: " << numphysical << std::endl; exit(1);}
  }

  MPO_matrix::MPO_matrix(const Basis& spectrum) : MPX_matrix(spectrum){};

  MPO_matrix::MPO_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, const SparseMatrix& matrix) : MPX_matrix(spectrum,indices,2,matrix)
  {
    check();
  }

  MPO_matrix::MPO_matrix(const Basis& spectrum,const std::vector<MPXIndex>& indices, SparseMatrix&& matrix) : MPX_matrix(spectrum,indices,2,std::move(matrix))
  {
    check();
  }

  MPO_matrix::MPO_matrix(MPX_matrix&& MPXref) noexcept : MPX_matrix(std::move(MPXref))
  {
    check();
  }

  MPO_matrix::MPO_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<complex<double> >& values, bool inverse): MPX_matrix(spectrum) {
    m_Indices.reserve(4);
    m_Indices.emplace_back(1,spectrum);
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,spectrum);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=2;
    m_Matrix=SparseMatrix(values,inverse);
  };

  MPO_matrix::MPO_matrix(const Basis& spectrum, const MPXIndex& index, const std::vector<double>& values, bool inverse): MPX_matrix(spectrum) {
    m_Indices.reserve(4);
    m_Indices.emplace_back(1,spectrum);
    m_Indices.emplace_back(1,index);
    m_Indices.emplace_back(0,spectrum);
    m_Indices.emplace_back(0,index);
    m_NumRowIndices=2;
    m_Matrix=SparseMatrix(values,inverse);
  }

  MPO_matrix MPO_matrix::ExtractMPOBlock(const std::pair<MPXInt,MPXInt>& matrix_index_row_range, const std::pair<MPXInt,MPXInt>& matrix_index_col_range) const {
    //get row range and col range
    std::pair<MPXInt,MPXInt> array_row_range(matrix_index_row_range.first*basis().size(),(matrix_index_row_range.second+1)*basis().size()-1);
    std::pair<MPXInt,MPXInt> array_col_range(matrix_index_col_range.first*basis().size(),(matrix_index_col_range.second+1)*basis().size()-1);

   StateArray row_matrix_index;
   for (MPXInt r=matrix_index_row_range.first;r<=matrix_index_row_range.second;++r){
      row_matrix_index.emplace_back(Index(1)[r]);
    }

   StateArray col_matrix_index;
   for (MPXInt c=matrix_index_col_range.first;c<=matrix_index_col_range.second;++c){
      col_matrix_index.emplace_back(Index(3)[c]);
    }
    
    return MPO_matrix(basis(),std::vector<MPXIndex>({{1,basis()},{1,std::move(row_matrix_index)},{0,basis()},{0,std::move(col_matrix_index)}}),m_Matrix.ExtractSubMatrix(array_row_range,array_col_range));

  }

  MPO_matrix UnitaryTransformMPO_matrix(const Basis& b, const std::vector<MPXIndex>& idxs, const SparseMatrix& s, size_t c, double f){
    if (c>=b.getChargeRules().size()){std::cout <<"Error, requested quantum number outised bounds!" <<std::endl; return MPO_matrix();}

    SparseMatrix T(copy(s)); //make a copy (default constructor is private)
    if (f!=0.0){
      for (MPXInt col=0;col<T.cols();++col){
	std::complex<double> Udagger(cos(b[col][c]*f),sin(-b[col][c]*f));
	for (MPXInt p=T.get_p(col);p<T.get_p(col+1);++p){
	  MPXInt row(T.get_i(p));
	  std::complex<double> U(cos(b[row][c]*f),sin(b[row][c]*f));
	  T.put_x(p)=U*T.get_x(p)*Udagger;
	}
      }
    }
    return MPO_matrix(b,idxs,std::move(T));
  }

  MPO_matrix load_MPO_matrix(const std::string& filename,const Basis& spectrum){
    return MPO_matrix(std::move(load_MPX_matrix_binary(filename,spectrum)));
  }

}

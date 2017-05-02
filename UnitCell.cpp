#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
//#include <utility>
#include <algorithm>

#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp" //specialisations, MPS
#include "UnitCell.hpp" //specialisations, MPS

namespace ajaj {

  void UnitCell::OutputPhysicalIndexDensities(std::ofstream& d) const {
    if (d.is_open() && size()){
      std::vector<std::complex<double> > densities;
      MPX_matrix accumulator(contract(Matrices.at(0),1,Matrices.at(0),0,contract0011));
      for (uMPXInt i=1;i<Matrices.size()-1;++i){
	accumulator=std::move(contract(contract(accumulator,0,Matrices.at(i),1,contract01),0,Matrices.at(i),0,contract0110));
      }
      {
	const MPS_matrix& A=Matrices.back();
	const std::vector<double>& L=Lambdas.at(0);
	MPX_matrix temp(contract(A,0,MPX_matrix(A.GetPhysicalSpectrum(),A.Index(2),L),0,contract20));
	densities=contract_to_sparse(contract(temp,1,accumulator,0,contract10),0,temp,0,contract1221).diagonal();
      }
      for (auto &&i : densities){
	d << real(i) << " ";
      }
      d << std::endl;
    }
  }
  
  void UnitCell::OutputOneVertexDensityMatrix(std::ofstream& d) const {
    if (d.is_open() && size()){
      MPX_matrix accumulator(contract(Matrices.at(0),1,Matrices.at(0),0,contract0011));
      for (uMPXInt i=1;i<Matrices.size()-1;++i){
	accumulator=std::move(contract(contract(accumulator,0,Matrices.at(i),1,contract01),0,Matrices.at(i),0,contract0110));
      }
      {
	const MPS_matrix& A=Matrices.back();
	const std::vector<double>& L=Lambdas.at(0);
	MPX_matrix temp(contract(A,0,MPX_matrix(A.GetPhysicalSpectrum(),A.Index(2),L),0,contract20));
	SparseMatrix S(contract_to_sparse(contract(temp,1,accumulator,0,contract10),0,temp,0,contract1221));
	S.fprint(d);
      }
    }
  }

  void UnitCell::OutputOneVertexDensityMatrix(const std::string& name, uMPXInt l) const {
    std::stringstream DensityMatrixNameStream;
    DensityMatrixNameStream << name << "_" << l << ".dat";
    std::ofstream DensityMatrixFileStream;
    DensityMatrixFileStream.open(DensityMatrixNameStream.str().c_str(),ios::out | ios::trunc);
    OutputOneVertexDensityMatrix(DensityMatrixFileStream);
    DensityMatrixFileStream.close();
  }

  void UnitCell::store(const std::string& filename, uMPXInt l) const {
    std::stringstream FileNameStream;
    FileNameStream << filename << "_" << l;
    store(FileNameStream.str());
  }

  void UnitCell::store(const std::string& filename) const {
    //first do a consistency check
    if (Matrices.size()!=Lambdas.size()){
      std::cout << "Invalid Unit Cell, skipping store." << std::endl; 
    }
    else {
      std::stringstream FileNameStream;
      FileNameStream << filename << ".UNITCELL"; 
      std::ofstream outfile;
      outfile.open(FileNameStream.str().c_str(),ios::out | ios::trunc | ios::binary);
      //first output basis
      basis_ptr_->fprint_binary(outfile);   
      size_t num_in_unit_cell(Matrices.size());
      outfile.write(reinterpret_cast<const char*>(&num_in_unit_cell),sizeof(size_t));
      for (auto&& m : Matrices){
	m.fprint_binary(outfile);
      }
      for (auto&& l : Lambdas){
	size_t lambda_size(l.size());
	outfile.write(reinterpret_cast<const char*>(&lambda_size),sizeof(size_t));	
	for (auto v : l){
	  outfile.write(reinterpret_cast<const char*>(&v),sizeof(double));
	}
      }
    }
  }


  double UnitCell::Entropy() const {
    if (Lambdas.size()){
      return entropy(Lambdas[0]);
    }
    else {
      std::cout << "UnitCell illformed. Can't compute entanglement entropy, returning -1" <<std::endl;
      return -1.0; //indicates failure
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  UnitCell load_UnitCell_binary(std::ifstream& infile, QNVector& charge_rules, Basis& basis){
    UnitCell ans(basis);
    if (infile.is_open()){
      //load in basis
      if (load_Basis_binary(infile,charge_rules,basis)) return ans;
      //read in size of unitcell
      size_t num_in_cell(0);
      infile.read(reinterpret_cast<char*>(&num_in_cell),sizeof(size_t));
      for (size_t n=0;n<num_in_cell;++n){
	ans.Matrices.emplace_back(MPS_matrix(load_MPX_matrix_binary(infile,basis)));
      }
      for (size_t n=0;n<num_in_cell;++n){
	std::vector<double> lambda;
	size_t num_in_lambda(0);
	infile.read(reinterpret_cast<char*>(&num_in_lambda),sizeof(size_t));
	for (size_t l=0;l<num_in_lambda;++l){
	  double val(0.0);
	  infile.read(reinterpret_cast<char*>(&val),sizeof(double));
	  lambda.push_back(val);
	}
	ans.Lambdas.emplace_back(std::move(lambda));
      }
      return ans;
    }
    else {
      return ans;
    }
  }

  UnitCell load_UnitCell_binary(std::ifstream& infile, const QNVector& charge_rules, const Basis& basis){
    UnitCell ans(basis);
    if (infile.is_open()){
      //load in basis form file and see if it is consistent...
      Basis l_basis;
      QNVector l_charge_rules;
      load_Basis_binary(infile,l_charge_rules,l_basis);
      //test Basis
      bool basis_good=1;
      if (l_basis.size()==basis.size()){
	for (size_t b=0; b< l_basis.size();++b){
	  if (l_basis[b]!=basis[b]) {basis_good=0; break;}
	}
      }
      else {
	basis_good=0;
      }

      if (!basis_good){
	std::cout << "Loaded basis doesn't match model!" <<std::endl;
	return ans;							
      }

      //read in size of unitcell
      size_t num_in_cell(0);
      infile.read(reinterpret_cast<char*>(&num_in_cell),sizeof(size_t));
      for (size_t n=0;n<num_in_cell;++n){
	ans.Matrices.emplace_back(MPS_matrix(load_MPX_matrix_binary(infile,basis)));
      }
      for (size_t n=0;n<num_in_cell;++n){
	std::vector<double> lambda;
	size_t num_in_lambda(0);
	infile.read(reinterpret_cast<char*>(&num_in_lambda),sizeof(size_t));
	for (size_t l=0;l<num_in_lambda;++l){
	  double val(0.0);
	  infile.read(reinterpret_cast<char*>(&val),sizeof(double));
	  lambda.push_back(val);
	}
	ans.Lambdas.emplace_back(std::move(lambda));
      }
      return ans;
    }
    else {
      return ans;
    }
  }

}

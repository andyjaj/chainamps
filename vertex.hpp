/** @file vertex.hpp
 * Vertex means site or integrable chain depending on the model studied. Essentially a vertex definition describes the local Hilbert space.
 */
#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <string>
#include <iostream>

#include "ajaj_common.hpp" //defines basic types
#include "states.hpp" //Spectrum etc.
#include "sparse_interface.hpp" //for VertexOperators
#include "dense_interface.hpp"
#include "MPX.hpp"


namespace ajaj{

  typedef EigenStateArray Basis;

/** Vertex Operator class contains the matrix elements and the name of an operator of the local (site or integrable chain) Hilbert space.
 *
 */
  class VertexOperator {
  public:
    std::string Name;
    SparseMatrix MatrixElements;
    VertexOperator(const string& name) : Name(name){};
    VertexOperator(const string& name, const Sparseint& lineardim) : Name(name), MatrixElements(SparseMatrix(lineardim,lineardim)) {};
    void print() const {std::cout << Name << std::endl; MatrixElements.print();}

  };
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
/** A parameter for a Vertex. For example for a site, the local field or for an integrable chain the mass.*/
  class VertexParameter {
  public:
    std::string Name;
    double Value;

    VertexParameter(const string& name,const double& value) : Name(name), Value(value){};

    void print() const {std::cout << Name << " = " << Value << std::endl;}
  };
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
/** Container for vertex operators. Uses std::vector interface.*/
  typedef std::vector<VertexOperator> VertexOperatorArray;
/** Container for vertex parameters. Uses std::vector interface.*/
  typedef std::vector<VertexParameter> VertexParameterArray;
  //////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
/** Vertex defines the local Hilbert space, the local Hamiltonian, any operators etc.*/
  class Vertex {
  public:
    std::string Name;
    VertexParameterArray Parameters; /**< Params used to create this vertex. Could be empty if predefined basis.*/
    Basis Spectrum; /**< Describes basis states*/
    VertexOperatorArray Operators; /**< Collection of all necessary matrix elements. */
    QNVector ChargeRules; //vertex should know about its good quantum numbers

    Vertex(){};
    Vertex(const VertexParameterArray& inputs) : Parameters(inputs){};

    const SparseMatrix& get_operator_matrix(const std::string& Name) const {
      for (auto&& o : Operators){
	if (Name==o.Name) return o.MatrixElements;
      }
      //force failure in non elegant manner
      std::cout << "Cannot find operator named " << Name <<std::endl;
      return Operators.at(Operators.size()).MatrixElements; //use at() to raise exception, by looking past end of vector
    }

    MPO_matrix make_one_site_operator(uMPXInt i) const{
      std::vector<MPXIndex> indices;
      StateArray dummy(1,State(basis().getChargeRules()));
      indices.emplace_back(1,basis());
      indices.emplace_back(1,dummy);
      indices.emplace_back(0,basis());
      indices.emplace_back(0,dummy);
      return MPO_matrix(basis(),indices,Operators.at(i).MatrixElements);
    };

    MPO_matrix make_one_site_operator(const std::string& Name) const{
      std::vector<MPXIndex> indices;
      StateArray dummy(1,State(basis().getChargeRules()));
      indices.emplace_back(1,basis());
      indices.emplace_back(1,dummy);
      indices.emplace_back(0,basis());
      indices.emplace_back(0,dummy);
      return MPO_matrix(basis(),indices,get_operator_matrix(Name));
    }

    const Basis& basis() const {return Spectrum;}

    const QNVector& getChargeRules() const {
      return ChargeRules;
    }

  };

}
#endif

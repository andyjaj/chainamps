/** @file vertex.hpp
 * A model is a combination of the vertex definition (e.g. a single spin or perhaps an exactly solvable chain) and definitions for the couplings between vertices (plus possible local terms applied to a vertex, like a magnetic field). Hence the model describes the whole system.
 */
#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <complex>

//#include "ajaj_common.hpp" //defines basic types
#include "vertex.hpp"
#include "MPX.hpp"
#include "states.hpp"


namespace ajaj{

  struct Model{
  public:    
    Vertex vertex;
    MPO_matrix H_MPO;

    Model() {};
    Model(const VertexParameterArray& vp, Vertex (*generator) (const VertexParameterArray&), const VertexParameterArray& cp,  MPO_matrix (*makeH) (const Vertex&, const VertexParameterArray&)) : vertex(generator(vp)),H_MPO(makeH(vertex,cp)) {
      std::cout << "MODEL'S LOCAL BASIS" <<std::endl;
      basis().print();
      std::cout << "MPO matrix info" <<std::endl;
      H_MPO.print_indices();
    }

    const Basis& basis() const {return vertex.basis();}

  };

  struct Coupling{
  public:
    std::string Name;
    std::string Op1Name;
    std::string Op2Name;
    std::complex<double> Value;
  };

  typedef std::vector<Coupling> CouplingArray;

  std::istream& operator>>(std::istream& is, Coupling& c)
  {
    char check;

    // read c from stream
    //format is CouplingName,(Op1Name,Op2Name):Value

    //drop initial whitespace
    is >> std::ws;

    bool failure(0);
    Coupling temp;


    if (!getline(is,temp.Name,',')) failure=1;
    if (!is.get(check) || check!='(') failure=1;
    if (!getline(is,temp.Op1Name,',')) failure=1;
    if (!getline(is,temp.Op2Name,')')) failure=1;
    if (!is.get(check) || check!=':') failure=1;
    if (!(is >>temp.Value)) failure=1;

    if( failure )
      is.setstate(std::ios::failbit);
    else
      c=std::move(temp);

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const Coupling& c){

    os << c.Name <<",("<<c.Op1Name<<","<<c.Op2Name<<"):"<<c.Value;

    return os;
  }

}

#endif

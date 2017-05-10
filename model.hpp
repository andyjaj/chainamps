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
#include <utility>

#include "vertex.hpp"
#include "MPX.hpp"
#include "states.hpp"


namespace ajaj{

  struct Coupling{
  public:
    std::string Name;
    std::complex<double> Value;
    std::string Op1Name;
    std::string Op2Name;    

    Coupling() {};
    Coupling(const std::string& N, std::complex<double> v) : Name(N),Value(v){};

    bool single_vertex() const {
      return Op2Name.empty() && !Op1Name.empty();
    }

    bool two_vertex() const {
      return !Op2Name.empty() && !Op1Name.empty();
    }
  };

  typedef std::vector<Coupling> CouplingArray;

  struct Model{
  private:
    std::vector<double> times_;
    std::vector<CouplingArray> CAs_;
    MPO_matrix (*makeH_) (const Vertex&, const CouplingArray&);

  public:    
    Vertex vertex;
    MPO_matrix H_MPO;

    Model() : makeH_(nullptr) {};//default, used by user defined
    Model(const VertexParameterArray& vp, Vertex (*generator) (const VertexParameterArray&), const std::vector<CouplingArray>& CAs, MPO_matrix (*makeH) (const Vertex&, const CouplingArray&), const std::vector<double>& t=std::vector<double>()) : CAs_(CAs), times_(t), makeH_(makeH), vertex(generator(vp)),H_MPO(makeH(vertex,CAs_.size() ? CAs_[0] : CouplingArray())) { //used by built in models
      std::cout << "MODEL'S LOCAL BASIS" <<std::endl;
      basis().print();
      std::cout << "MPO MATRIX INFO" <<std::endl;
      H_MPO.print_indices();
      H_MPO.print_sparse_info();
    }
    Model(const std::vector<CouplingArray>& CAs,  MPO_matrix (*makeH) (const Vertex&, const CouplingArray&), const std::vector<double>& t) : CAs_(CAs), times_(t), makeH_(makeH){};//default, used by user defined

    const Basis& basis() const {return vertex.basis();}

    State make_target(const QNVector& vec) const {
      std::cout <<"Setting target quantum numbers" <<std::endl;
      return State(basis().getChargeRules(),vec);
    }

    const std::vector<double>& times() const {return times_;}
    const std::vector<CouplingArray>& coupling_arrays() const {return CAs_;}

    const MPO_matrix& change_H_MPO(const CouplingArray& ca){
      H_MPO=makeH_(vertex,ca);
      return H_MPO;
    }

  };

  std::string trim(const std::string& s){
    size_t start = s.find_first_not_of(' ');

    if (start == string::npos) return std::string();

    size_t end = s.find_last_not_of(' ');

    return s.substr(start, (end-start+1));
  }


  std::istream& operator>>(std::istream& is, Coupling& c)
  {
    char check;

    // read c from stream
    //format is CouplingName,(Op1Name,Op2Name):Value

    //drop initial whitespace
    is >> std::ws;

    bool failure(0);
    Coupling temp;

    if(is.peek()=='\n' ){
      is.setstate(std::ios::failbit);
      char dump;
      is >> dump;
      return is;
    }

    if (!getline(is,temp.Name,',')) failure=1;
    if (!is.get(check) || check!='(') failure=1;
    //Now inside parenthesis
    //need to check if there are two names (comma separated) or one
    std::string opnames;
    if (!getline(is,opnames,')')) failure=1;

    size_t comma_pos=opnames.find(",");
    if (comma_pos==std::string::npos){
      //only one name
      temp.Op1Name=trim(opnames);

    }
    else { //two names

      //check for more commas...
      size_t more=opnames.find(",",comma_pos+1);
      if (more==std::string::npos){
	temp.Op1Name=trim(opnames.substr(0,comma_pos));
	temp.Op2Name=trim(opnames.substr(comma_pos+1,opnames.length()));
      }
      else {
	std::cout << "ERROR: Too many comma separated operator names specified. " << opnames << std::endl;
	std::cout << "Leaving couplings undefined." <<std::endl;
	failure=1;
      }
    }

    is >> std::ws;
    if (!is.get(check) || check!=':') failure=1;

    if (!(is >>temp.Value)) failure=1;

    if( failure )
      is.setstate(std::ios::failbit);
    else
      c=std::move(temp);

    return is;
  }

  std::ostream& operator<<(std::ostream& os, const Coupling& c){
    if (c.Op2Name.empty()){
      os << c.Name <<",("<<c.Op1Name<<"):"<<c.Value;
    }
    else {
      os << c.Name <<",("<<c.Op1Name<<","<<c.Op2Name<<"):"<<c.Value;
    }
    return os;
  }

}

#endif

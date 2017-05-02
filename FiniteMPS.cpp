#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp"
#include "FiniteMPS.hpp"

namespace ajaj{

  FiniteMPS::FiniteMPS(const Basis& model_basis, std::string name, uMPXInt num,c_num_specifier coeffs) : Basis_(model_basis),MPSName_(name),NumVertices_(num),current_(model_basis)
  {

    //check coeffs for legality first?

    if (coeffs.size()==0){
      std::cout << "Error: Empty c number specification for state!" <<std::endl; exit(1);
    }
    for (auto&& c_vect : coeffs){
      if (c_vect.size()==0){
	std::cout << "Error: Empty c number specification for vertex!" <<std::endl; exit(1);
      }
      State test_state(Basis_[c_vect.begin()->first]);
      for (auto&& c_pair : c_vect){
	if (Basis_[c_pair.first]!=test_state){
	  std::cout << "Error: Inconsistent c numbers for vertex state: " << c_vect.begin()->first << " " << c_pair.first  << std::endl; exit(1);
	}
	else {std::cout << "(" << c_pair.first << " : " << c_pair.second << ") ";}
      }
      std::cout << std::endl;
    }

    State LeftState(Basis_.getChargeRules()); //dummy state (all zeros)
    size_t counter=0;
    for (uMPXInt n=1;n<=NumVertices_;++n){
      auto& c_vec=coeffs[counter++];
      current_=MakeProductState(Basis_,c_vec,LeftState);
      LeftState+=Basis_[c_vec.begin()->first];
      if(counter==coeffs.size()){
	counter=0;//back to beginning of c number defns
      }
      std::stringstream NameStream;
      NameStream << MPSName_ << "_Left_" << n << ".MPS_matrix";
      current_.store(NameStream.str());
    }
  }


}

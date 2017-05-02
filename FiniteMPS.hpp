/** @file FiniteMPS.hpp
 * Class to manage disk stored MPS matrices, for finite systems.
 */
#ifndef FMPS_H
#define FMPS_H

#include <vector>
#include <string>
#include <utility>

//#include <sstream>
//#include <iostream>
//#include <fstream>

#include "common_defs.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPXIndex.hpp" //tensor index class
#include "MPX_matrix.hpp" //general tensor class
#include "MPS_matrix.hpp"

namespace ajaj{

  typedef std::vector<std::vector<std::pair<uMPXInt,double> > > c_num_specifier;

  class FiniteMPS{
  private:
    const Basis& Basis_;
    std::string MPSName_;
    uMPXInt NumVertices_;
    MPS_matrix current_;

  public:

    FiniteMPS(const Basis& model_basis, std::string name, uMPXInt num) : Basis_(model_basis),MPSName_(name),NumVertices_(num),current_(model_basis) {}
    FiniteMPS(const Basis& model_basis, std::string name, uMPXInt num,c_num_specifier coeffs);

  };

}

#endif

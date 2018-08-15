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

  enum class CanonicalType : unsigned short int {Left, Right, Mixed, Non, Error};

  class FiniteMPS{
  private:
    const Basis& Basis_;
    std::string MPSName_;
    uMPXInt NumVertices_;
    std::pair<uMPXInt,MPS_matrix> Current_;
    bool Canonical_;
    uMPXInt MixPoint_; 

    void fetch_matrix(uMPXInt i,bool Left=1); /**<Get a specific matrix*/
    std::string filename(uMPXInt i,bool Left=1,const std::string& name=std::string()) const;
    void store_current();
    CanonicalType CheckFilesExist();

  public:

    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num) : Basis_(model_basis),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(0) {} /**< Create a non canonical finite MPS, with no data */
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt mix_idx); /**< Finite MPS, with mixpoint=mix_idx */
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs); /**< Specify a finite MPS product state, makes it left canonical*/
    
    uMPXInt position() const {return Current_.first;}

    const MPS_matrix& matrix() const {return Current_.second;} //const function to get const ref to current buffered matrix
    const MPS_matrix& matrix(uMPXInt p,bool Left=1) { //non const function to buffer a particular matrix
      fetch_matrix(p,Left); //buffer matrix
      return matrix(); //return const ref
    }

    uMPXInt size() const {return NumVertices_;}
    const std::string& name() {return MPSName_;}

    std::complex<double> makeLC(const std::string& new_name=std::string()); /**< 'Ensures' left canonical, and makes an optional copy, returns 0.0 if there is an error. */
    std::complex<double> makeRC(const std::string& new_name=std::string()); /**< 'Ensures' right canonical, and makes an optional copy, returns 0.0 if there is an error. */
    bool valid_files() {return CheckFilesExist()==CanonicalType::Error ? 0 : 1;}
  };

}

#endif

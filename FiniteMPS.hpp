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

  enum class MPSCanonicalType : unsigned short int {Left, Right, Mixed, Non, Error};
  class FiniteMPS;
  class ConstFiniteMPS;
  std::complex<double> ApplySingleVertexOperatorToMPS(const MPO_matrix&, FiniteMPS& F, uMPXInt vertex, const MPSCanonicalType& RequestedCanonization=MPSCanonicalType::Non);
  
  class FiniteMPS{
  private:
    const Basis& Basis_;
    std::string MPSName_;
    uMPXInt NumVertices_;
    std::pair<uMPXInt,MPS_matrix> Current_;
    bool Canonical_;
    uMPXInt MixPoint_;
    MPSCanonicalType Canonization_;
    std::complex<double> Weight_;
    std::vector<MPS_matrixCanonicalType> MatrixCanonizations_;

    
    void fetch_matrix(uMPXInt i,bool Left=1); /**<Get a specific matrix*/
    std::string filename(uMPXInt i,bool Left=1,const std::string& name=std::string()) const;
    void store_current();
    MPSCanonicalType CheckFilesExist(const std::string& newname=std::string()); /**<Check files exist, and optionally copy and store with a new name*/
    void set_matrix_canonization(uMPXInt i,const MPS_matrixCanonicalType& c);
    void update_MPS_canonization_status();
    
  public:

    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num) : Basis_(model_basis),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonical_(0),Canonization_(MPSCanonicalType::Non) {} /**< Create a non canonical finite MPS, with no data */
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt mix_idx); /**< Finite MPS, with mixpoint=mix_idx */
    FiniteMPS(const Basis& model_basis, const std::string& oldname, const std::string& newname, uMPXInt num, bool canon, uMPXInt mix_idx); /**< Create stored copy Finite MPS, with mixpoint=mix_idx */
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs); /**< Specify a finite MPS product state, makes it left canonical*/
    
    const MPS_matrix& matrix() const {return Current_.second;} //const function to get const ref to current buffered matrix
    const MPS_matrix& matrix(uMPXInt p,bool Left=1) { //non const function to buffer a particular matrix
      fetch_matrix(p,Left); //buffer matrix
      return matrix(); //return const ref
    }

    uMPXInt position() const {return Current_.first;}
    uMPXInt size() const {return NumVertices_;}
    const std::string& name() {return MPSName_;}
    std::complex<double> makeLC(const std::string& new_name=std::string()); /**< 'Ensures' left canonical, and makes an optional copy, returns final phase times singular val*/
    std::complex<double> makeRC(const std::string& new_name=std::string()); /**< 'Ensures' right canonical, and makes an optional copy, returns final phase times singular val */
    bool valid_files() {return CheckFilesExist()==MPSCanonicalType::Error ? 0 : 1;}
    std::complex<double> weight() const {return Canonical_ ? Weight_ : 0.0 ;}
    void reset_weight(std::complex<double> w) {if (Weight_!=0.0) Weight_=w; }
    
    friend std::complex<double> ApplySingleVertexOperatorToMPS(const MPO_matrix&, FiniteMPS& F, uMPXInt vertex, const MPSCanonicalType& RequestedCanonization);
    friend class ConstFiniteMPS;
    
  };

  class ConstFiniteMPS{
  private:
    const Basis& Basis_;
    const std::string MPSName_;
    const uMPXInt NumVertices_;
    const uMPXInt MixPoint_;
    const MPSCanonicalType Canonization_;
    const std::complex<double> Weight_;
    const std::vector<MPS_matrixCanonicalType> MatrixCanonizations_;
    
    std::string filename(uMPXInt i,bool Left=1,const std::string& name=std::string()) const;

  public:
    ConstFiniteMPS(const FiniteMPS& F) : Basis_(F.Basis_),MPSName_(F.MPSName_),NumVertices_(F.NumVertices_),MixPoint_(F.MixPoint_),Canonization_(F.Canonization_),Weight_(F.Weight_),MatrixCanonizations_(F.MatrixCanonizations_){};
    ConstFiniteMPS(const Basis& B, const std::string& Name, uMPXInt Num, uMPXInt MP, MPSCanonicalType C, const std::complex<double>& W, const std::vector<MPS_matrixCanonicalType>& MC) : Basis_(B), MPSName_(Name), NumVertices_(Num), MixPoint_(MP), Canonization_(C), Weight_(W), MatrixCanonizations_(MC) {};

    ConstFiniteMPS(const Basis& B, const std::string& Name, uMPXInt Num) : Basis_(B), MPSName_(Name), NumVertices_(Num), MixPoint_(Num), Canonization_(MPSCanonicalType::Left), Weight_(1.0), MatrixCanonizations_(std::vector<MPS_matrixCanonicalType>(Num,MPS_matrixCanonicalType::Left)) {};
    
    uMPXInt size() const {return NumVertices_;}
    std::complex<double> weight() const {return Weight_;}
    const std::string& name() {return MPSName_;}
    MPS_matrix matrix(uMPXInt i) const;
  };
}

#endif

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
    State TotalCharges_;
    std::string MPSName_;
    uMPXInt NumVertices_;
    std::pair<uMPXInt,MPS_matrix> Current_; //currently cached matrix
    //bool Canonical_; //is this still used?
    MPSCanonicalType Canonization_;
    bool HasLambda_;
    //uMPXInt MixPoint_; //for mixed canonical, Should be > 0 and < NumVertices_.
                       //In principle a state could have a valid Lambda matrix, but no longer be truly mixed canonical
                       //e.g. because an operator has been applied to it.
    std::complex<double> Weight_; //norm and possible phase
    std::vector<MPS_matrixCanonicalType> MatrixCanonizations_; // 'list' of canonization status for each matrix    
    uMPXInt LambdaPosition_; //replacement for mix point
    
    void fetch_matrix(uMPXInt i,bool Left=1); /**<Get a specific matrix*/
    MPX_matrix fetch_lambda() const; //For mixed canonical states, must be able to incorporate lambda
    std::string filename(uMPXInt i,bool Left=1,const std::string& name=std::string()) const;
    void store_current();
    MPSCanonicalType check_files_exist(bool canon=0, const std::string& newname=std::string()); /**<Check files exist, and optionally copy and store with a new name*/
    void set_matrix_canonization(uMPXInt i,const MPS_matrixCanonicalType& c);
    /*void set_state_nc(bool hl=0, uMPXInt lp=0){Canonization_=MPS_matrixCanonicalType::Non; HasLambda_=hl; LambdaPosition_=lp;}
    void set_state_mc(uMPXInt lp){Canonization_=MPS_matrixCanonicalType::Mixed; HasLambda_=1; LambdaPosition_=lp;}
    void set_state_lc(bool hl=0){Canonization_=MPS_matrixCanonicalType::Left; HasLambda_=hl; LambdaPosition_=NumVertices_ ;}
    void set_state_rc(bool hl=0){Canonization_=MPS_matrixCanonicalType::Right; HasLambda_=hl; LambdaPosition_=0;}*/
    void update_MPS_canonization_status();
    std::pair<MPX_matrix,MPX_matrix> left_canonize_to(uMPXInt i,const std::string& new_name=std::string());
    std::pair<MPX_matrix,MPX_matrix> right_canonize_to(uMPXInt i,const std::string& new_name=std::string());

  public:

    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num) : Basis_(model_basis),TotalCharges_(Basis_.Identity()),MPSName_(name),NumVertices_(num),Current_(std::pair<uMPXInt,MPS_matrix>(0,MPS_matrix(model_basis))),Canonization_(MPSCanonicalType::Non),HasLambda_(0),Weight_(1.0) {} /**< Create a non canonical finite MPS, with no data */
    FiniteMPS(const Basis& model_basis, const std::string& oldname, const std::string& newname, uMPXInt num, bool canon, uMPXInt lp); /**< Create Finite MPS (a new copy unless newname is an empty string), using stored MPS with oldname */
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num, bool canon, uMPXInt lp); /**< Finite MPS, with LambdaPosition_=lp */
    
    FiniteMPS(const Basis& model_basis, const std::string& name, uMPXInt num,const c_specifier_array& coeffs); /**< Specify a finite MPS product state*/
    
    const MPS_matrix& matrix() const {return Current_.second;} //const function to get const ref to current buffered matrix
    const MPS_matrix& matrix(uMPXInt p,bool Left=1) { //non const function to buffer a particular matrix
      fetch_matrix(p,Left); //buffer matrix
      return matrix(); //return const ref
    }
    MPX_matrix lambda() const {return fetch_lambda();}

    uMPXInt position() const {return Current_.first;}
    uMPXInt size() const {return NumVertices_;}
    const std::string& name() {return MPSName_;}
    std::complex<double> makeLC(const std::string& new_name=std::string()); /**< 'Ensures' left canonical, and makes an optional copy, returns final phase times singular val*/
    std::complex<double> makeRC(const std::string& new_name=std::string()); /**< 'Ensures' right canonical, and makes an optional copy, returns final phase times singular val */
    double mixed(uMPXInt mixposition,const std::string& new_name=std::string());
    State total_charges() {return TotalCharges_;}

    
    bool valid_files() {bool canon = Canonization_!=MPSCanonicalType::Non ? 1 : 0; return check_files_exist(canon)==MPSCanonicalType::Error ? 0 : 1;}
    std::complex<double> weight() const {return (Canonization_!=MPSCanonicalType::Non && Canonization_!=MPSCanonicalType::Error) ? Weight_ : 0.0 ;}
    void reset_weight(std::complex<double> w) {if (Weight_!=0.0) Weight_=w; }
    std::string canon_type_string() const;
    
    friend std::complex<double> ApplySingleVertexOperatorToMPS(const MPO_matrix&, FiniteMPS& F, uMPXInt vertex, const MPSCanonicalType& RequestedCanonization);
    friend class ConstFiniteMPS;
    
  };

  class ConstFiniteMPS{
  private:
    const Basis& Basis_;
    const std::string MPSName_;
    const uMPXInt NumVertices_;
    const uMPXInt LambdaPosition_;
    const MPSCanonicalType Canonization_;
    const std::complex<double> Weight_;
    const std::vector<MPS_matrixCanonicalType> MatrixCanonizations_;
    
    std::string filename(uMPXInt i,bool Left=1,const std::string& name=std::string()) const;

  public:
    ConstFiniteMPS(const FiniteMPS& F) : Basis_(F.Basis_),MPSName_(F.MPSName_),NumVertices_(F.NumVertices_),LambdaPosition_(F.LambdaPosition_),Canonization_(F.Canonization_),Weight_(F.Weight_),MatrixCanonizations_(F.MatrixCanonizations_){};
    ConstFiniteMPS(const Basis& B, const std::string& Name, uMPXInt Num, uMPXInt MP, MPSCanonicalType C, const std::complex<double>& W, const std::vector<MPS_matrixCanonicalType>& MC) : Basis_(B), MPSName_(Name), NumVertices_(Num), LambdaPosition_(MP), Canonization_(C), Weight_(W), MatrixCanonizations_(MC) {};

    ConstFiniteMPS(const Basis& B, const std::string& Name, uMPXInt Num) : Basis_(B), MPSName_(Name), NumVertices_(Num), LambdaPosition_(Num), Canonization_(MPSCanonicalType::Left), Weight_(1.0), MatrixCanonizations_(std::vector<MPS_matrixCanonicalType>(Num,MPS_matrixCanonicalType::Left)) {};
    
    uMPXInt size() const {return NumVertices_;}
    std::complex<double> weight() const {return Weight_;}
    const std::string& name() {return MPSName_;}
    MPS_matrix matrix(uMPXInt i) const;
  };
}

#endif

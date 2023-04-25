#include <vector>
#include <array>

#include "arpack_interface.hpp"
#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"
#include "measurement.hpp"

namespace ajaj {

  static const double IMAGTOL(100.0*std::numeric_limits<double>::epsilon());
  static const int NUMEVALS(1);

  MultiVertexMeasurement::MultiVertexMeasurement(uMPXInt start, const MPO_matrix* Op1Ptr, uMPXInt finish, const MPO_matrix* Op2Ptr) : VertexOperatorPtrs_({{meas_pair(start,Op1Ptr), meas_pair(finish,Op2Ptr)}}), T_(Op1Ptr->GetPhysicalSpectrum()),Result_(0.0) {
    start_=start;
    finish_=finish;
    if (start>finish || start<1){
      std::cout << "Incorrectly defined MultiVertexMeasurement positions, start, finish: " << start << "," << finish <<std::endl;
      exit(1);
    }
  }

  void MultiVertexMeasurement::find_start_and_finish_(){
    if (VertexOperatorPtrs_.size()){
      uMPXInt Ts=VertexOperatorPtrs_.begin()->position();
      uMPXInt Tf=Ts;
      for (auto&& vop : VertexOperatorPtrs_){
	if (vop.position()<Ts) Ts=vop.position();
	else if (vop.position()>Tf) Tf=vop.position();
      }
      start_=Ts;
      finish_=Tf;

      

      if (start_ > finish_) {
	std::cout << "Incorrectly defined MultiVertexMeasurement positions, start, finish: " << start_ << "," << finish_ <<std::endl;
	exit(1);
      }
      
    }
    else {
      start_=0;
      finish_=0;
    }
  }
  
  std::vector<const MPO_matrix*> MultiVertexMeasurement::get_ops(uMPXInt v){
      std::vector<const MPO_matrix*> ans; 
      for (auto&& vop : VertexOperatorPtrs_){
	//if (vop.position()>v) break; //stop searching when we reach a point past v.
	  if (v==vop.position()) {
	    ans.emplace_back(vop.MPO_ptr());
	  }	
	}
      return ans;
    }
  
  void MultiVertexMeasurement::link_(const std::vector<const MPO_matrix*>& ops /*const MPO_matrix* Op*/, const MPX_matrix& A) {
    //if (Op){
    if (ops.size()){
      //make sure to strip 'dummy' indices by contract to sparse and rebuild indices
      //we assume single site operators, so MPO _matrix_ indices are dummies...
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,A.Index(2));
      indices.emplace_back(0,A.Index(2));
      T_=std::move(contract(*ops[0],0,contract(T_,0,A,0,contract11),0,contract21).RemoveDummyIndices(std::vector<MPXInt>({{1,2}})));
      for (size_t o=1;o<ops.size();++o){
	T_=std::move(contract(*ops[o],0,T_,0,contract20).RemoveDummyIndices(std::vector<MPXInt>({{1,2}})));
      }
      T_=std::move(contract(A,1,T_,0,contract0011));
    }
    else {
      T_=std::move(contract(A,1,contract(T_,0,A,0,contract11),0,contract0110));
    }
  }
  
  void MultiVertexMeasurement::start_chain_(const MPX_matrix& A) {
    std::vector<const MPO_matrix*> ops(get_ops(start()));
    if (ops.size()){
      //make sure to strip 'dummy' indices by contract to sparse and rebuild indices
      //we assume single site operators, so MPO _matrix_ indices are dummies...
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,A.Index(2));
      indices.emplace_back(0,A.Index(2));
      //first special
      T_=std::move(contract(*ops[0],0,A,0,contract20).RemoveDummyIndices(std::vector<MPXInt>({{1,2}})));
      for (size_t o=1;o<ops.size();++o){
	T_=std::move(contract(*ops[o],0,T_,0,contract20).RemoveDummyIndices(std::vector<MPXInt>({{1,2}}))); //get rid of the single vertex MPO indices.
      }
      T_=std::move(contract(A,1,T_,0,contract0011));//contract on bra
    }
    else {
      T_=std::move(contract(A,1,A,0,contract0011));//no ops, for canonical matrices this should just give the identity
    }
  }
  
  void MultiVertexMeasurement::finish_chain_(const MPX_matrix& Lambda) {
    std::cout << "Finishing measurement contraction chain" << std::endl;
    //contract density matrix (lambda squared) onto end and trace
    Result_=contract_to_sparse(Lambda,0,contract(T_,0,Lambda,0,contract10),0,contract0110).trace();
    std::cout << "Done" << std::endl;
  }

  void MultiVertexMeasurement::update(uMPXInt v, const MPXDecomposition& D) {
    if (v<finish() && v>start()){
      link_(get_ops(v),D.ColumnMatrix);
    }
    else {
      //depending on v, start measurement chain, or continue it
      if (v==start()){
	//std::vector<const MPO_matrix*> ops(get_ops(v));
	start_chain_(D.ColumnMatrix); //contract first measurement etc.
	
      }
      if (v==finish()){
	//apply and finish off with lambda^2
	if (finish()!=start()) {//one vertex measurement doesn't need this step
	  link_(get_ops(v),D.ColumnMatrix);
	}
	
	finish_chain_(MPX_matrix(D.ColumnMatrix.GetPhysicalSpectrum(),D.ColumnMatrix.Index(2),D.Values));
      }
      else { //v<start() or v>finish()
	//do nothing
      }
    }
  }
  
  UnitCell OrthogonaliseInversionSymmetric(const UnitCell& C){
    std::cout << "Orthogonalising (inversion symmetric version)..." << std::endl;
    const Basis& basis(C.Matrices.front().basis());

    //transfer matrix components with these flags enforces hermiticity of the (reshaped) left eigenvector

    SparseED LeftTdecomp(TransferMatrixComponents(C,1,State(C.basis().getChargeRules())).LeftED(NUMEVALS,LARGESTMAGNITUDE));
    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this
    if (abs(imag(LeftTdecomp.Values.at(0)))>=IMAGTOL*abs(real(LeftTdecomp.Values.at(0)))) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      return UnitCell(basis);
    }

    //decompose into Xdagger X form
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,C.Matrices.front().getInwardMatrixIndex());
    VLIndices.emplace_back(0,C.Matrices.front().getInwardMatrixIndex());
    //decompose it into Xdagger X form
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(basis,VLIndices,1,reshape(LeftTdecomp.EigenVectors.ExtractColumns(std::vector<MPXInt>({0})),C.Matrices.front().getInwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VLDecomp.first.size()) return UnitCell(basis);

    MPX_matrix X(contract(MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    //note that SqrtDR should guarantee that X^dagger is the inverse of X
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));
    
    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,C.Matrices.back().getOutwardMatrixIndex());
    VRIndices.emplace_back(0,C.Matrices.back().getOutwardMatrixIndex());
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(basis,VRIndices,1,reshape(LeftTdecomp.EigenVectors.ExtractColumns(std::vector<MPXInt>({0})),C.Matrices.back().getOutwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VRDecomp.first.size()) return UnitCell(basis);

    MPX_matrix Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));
    MPXDecomposition XLYDecomp(contract(contract(X,0,MPX_matrix(basis,C.Matrices.back().getOutwardMatrixIndex(),C.Lambdas.front()),0,contract10),0,Y,0,contract10).SVD());

    SquareSumRescale(XLYDecomp.Values,1.0);//rescale here to set normalisation
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;

    UnitCell ans(basis);
    std::vector<double> scalevecX(C.Matrices.front().getInwardMatrixIndex().size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));

    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(basis,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,C.Matrices.front(),0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(contract(C.Matrices.back(),0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20));
    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(C.Lambdas.back());
    return ans;
  }

  UnitCell Orthogonalise(const UnitCell& C){
    std::cout << "Orthogonalising..." << std::endl;
    const Basis& basis(C.Matrices.front().basis());

    //transfer matrix components with these flags enforces hermiticity of the (reshaped) left eigenvector
    SparseED LeftTdecomp(TransferMatrixComponents(C,1,C.basis().Identity()/*State(C.basis().getChargeRules())*/).LeftED(NUMEVALS,LARGESTMAGNITUDE));
    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this
    if (abs(imag(LeftTdecomp.Values.at(0)))>=IMAGTOL*abs(real(LeftTdecomp.Values.at(0)))) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      return UnitCell(basis);
    }

    //decompose into Xdagger X form
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,C.Matrices.front().getInwardMatrixIndex());
    VLIndices.emplace_back(0,C.Matrices.front().getInwardMatrixIndex());    
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(basis,VLIndices,1,reshape(LeftTdecomp.EigenVectors.ExtractColumns(std::vector<MPXInt>({0})),C.Matrices.front().getInwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VLDecomp.first.size()) {
      std::cout << "SqrtDR failure" <<std::endl;
      return UnitCell(basis);
    }
    std::cout << "Left vector decomposed" << std::endl;

    MPX_matrix X(contract(MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    //note that SqrtDR should guarantee that X^dagger is the inverse of X

    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));

    SparseED RightTdecompSpecial(TransferMatrixComponents(C,1,C.basis().Identity()/*State(C.basis().getChargeRules())*/).RightED(NUMEVALS,LARGESTMAGNITUDE));

    std::cout << "Leading right eigenvalue of unshifted right transfer matrix: " << RightTdecompSpecial.Values.at(0) << std::endl; //will need to rescale by this
    if (abs(imag(RightTdecompSpecial.Values.at(0)))>=IMAGTOL*abs(real(RightTdecompSpecial.Values.at(0)))) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      return UnitCell(basis);
    }
    
    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,C.Matrices.back().getOutwardMatrixIndex());
    VRIndices.emplace_back(0,C.Matrices.back().getOutwardMatrixIndex());
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(basis,VRIndices,1,reshape(RightTdecompSpecial.EigenVectors.ExtractColumns(std::vector<MPXInt>({0})),C.Matrices.back().getOutwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VRDecomp.first.size()) {
      std::cout << "VRDecomp failure" <<std::endl;
      return UnitCell(basis);
    }

    //std::cout << "Forming lambda * y" <<std::endl;
    
    MPX_matrix LAMBDA_Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));

    //std::cout << "Forming XLY" <<std::endl;
    MPXDecomposition XLYDecomp(contract(X,0,LAMBDA_Y,0,contract10).SVD());
    SquareSumRescale(XLYDecomp.Values,1.0);//rescale here to set normalisation
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;

    UnitCell ans(basis);
    //if norm isn't 1.0, we need to rescale X
    //std::vector<double> scalevecX(C.Matrices.front().getInwardMatrixIndex().size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));
    std::vector<double> scalevecX(X.Index(0).size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));
    
    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(basis,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,C.Matrices.front(),0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(contract(C.Matrices.back(),0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20));
    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(C.Lambdas.back());
    return ans;
  }

  //Only use this as a last resort, because it includes an explicit inverse of singular values
  UnitCell OrthogonaliseInversionSymmetric(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda){
    //makes a copy, in order to preserve the A Lambda B decomposition
    return OrthogonaliseInversionSymmetric(UnitCell(MPSD,PreviousLambda));
  }

  UnitCell Orthogonalise(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda){
    //makes a copy, in order to preserve the A Lambda B decomposition
    return Orthogonalise(UnitCell(MPSD,PreviousLambda));
  }

  MPX_matrix MakeMeasurementTransferMatrix(const MPO_matrix& W1, const MPO_matrix& W2, const MPS_matrix& A1, const MPS_matrix& A2){
    //assumption is the mother of all...
    //MPX_matrix V(contract(contract(A1,1,W1,0,contract00),0,A1,0,contract30));
    MPX_matrix ans(contract(contract(contract(contract(contract(A1,1,W1,0,contract00),0,A1,0,contract30),0,A2,1,contract11),0,W2,0,contract2150),0,A2,0,contract3150));
    ans.ShiftNumRowIndices(2);
    return ans;
  }

  //needs to be cleverer to deal with B type matrices
  MPX_matrix MakeTransferMatrix(const MPS_matrix& A1bra, const MPS_matrix& A2bra, const MPS_matrix& A1ket, const MPS_matrix& A2ket){
    MPX_matrix ans(contract(contract(contract(A1bra,1,A1ket,0,contract00),0,A2bra,1,contract11),0,A2ket,0,contract2130));
    ans.ShiftNumRowIndices(2);
    return ans;
  }

  MPX_matrix MakeLTransferMatrix(const UnitCell& bra, const UnitCell& ket){
    //check lengths are the same
    if (bra.Matrices.size()!=ket.Matrices.size()){
      std::cout << "UnitCell lengths don't match!" << std::endl; exit(1);
    }
    MPX_matrix accumulator(reorder(contract(bra.Matrices.at(0),1,ket.Matrices.at(0),0,contract00),0,reorder0213,2));
    for (uMPXInt i=1;i<bra.Matrices.size();++i){
      accumulator=std::move(contract(contract(accumulator,0,bra.Matrices.at(i),1,contract21),0,ket.Matrices.at(i),0,contract2130));
    }
    accumulator.ShiftNumRowIndices(2);
    return accumulator;
  }

  //inaccurate?
  void ShiftLTransferMatrixToR(MPX_matrix& TransferMatrix, const std::vector<double>& lambda0){
    MPX_matrix LambdaInverse(TransferMatrix.GetPhysicalSpectrum(),TransferMatrix.Index(0),lambda0,1);
    MPX_matrix Lambda(TransferMatrix.GetPhysicalSpectrum(),TransferMatrix.Index(0),lambda0,0);
    TransferMatrix=std::move(contract(contract(LambdaInverse,0,contract(contract(LambdaInverse,0,TransferMatrix,0,contract11),0,Lambda,0,contract21),0,contract01),0,Lambda,0,contract20).ShiftNumRowIndices(2));
  }
  //assumes inversion symmetry!!!!
  void ShiftLTransferMatrixToR(MPX_matrix& TransferMatrix){
    TransferMatrix.Transpose();
  }

  std::complex<double> Overlap(const UnitCell& bra, const UnitCell& ket){
    return Overlap(bra,ket,1).at(0);
  }

  std::vector<std::complex<double> > Overlap(const UnitCell& bra, const UnitCell& ket, uMPXInt nev){
    return TransferMatrixComponents(bra,ket,0,State(ket.basis().getChargeRules())).LeftED(nev,LARGESTMAGNITUDE).Values;
    //return MakeLTransferMatrix(bra,ket).LeftEigs(nev,LARGESTMAGNITUDE).Values;
  }

  std::vector<std::complex<double> > TransferMatrixEigs(const UnitCell& ket, uMPXInt nev,const State& TargetState){
    //make the required bits
    /*std::vector<const MPS_matrix*> ket_ptrs;
    for (auto&& k: ket.Matrices){
      ket_ptrs.emplace_back(&k);
    }*/
    return TransferMatrixComponents(ket,1,TargetState).LeftED(nev,LARGESTMAGNITUDE).Values;
    //return TransferMatrixComponents(ket_ptrs,1,TargetState).LeftED(nev,LARGESTMAGNITUDE).Values;
  }

  std::complex<double> TwoVertexMeasurement(const MPO_matrix& W1, const MPO_matrix& W2, const MPS_matrix& A1, const MPS_matrix& A2, const MPX_matrix& Lambda){
    //make rho matrix
    //MPX_matrix rho(contract(Lambda,0,Lambda,0,contract10));
    std::complex<double> ans(contract_to_sparse(contract(contract(contract(contract(contract(A1,1,A1,0,contract11),0,W1,0,contract0022),0,A2,1,contract01),0,W2,0,contract2130),0,A2,0,contract0130),0,contract(Lambda,0,Lambda,0,contract10),0,contract1130).trace());
    return ans;
  }

  std::complex<double> OneVertexMeasurement(const MPO_matrix& W, const UnitCell& U){ //measures on 2nd site
    MPX_matrix LAMBDA0(U.Matrices.at(0).GetPhysicalSpectrum(),U.Matrices.at(0).Index(1),U.Lambdas.at(0));
    //THIS IS INEFFICIENT BECAUSE MATRICES AREN'T CANONICAL
    //SHOULD BE FIXED IN ORTHO ROUTINE
    std::complex<double> ans(contract_to_sparse(contract(U.Matrices.at(1),1,contract(contract(contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract0011),0,U.Matrices.at(1),0,contract11),0,W,0,contract12),0,contract0210),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0110).trace());
    return ans;
  }

  std::complex<double> TwoVertexMeasurement(const MPO_matrix& W1,const MPO_matrix& W2,const UnitCell& U,uMPXInt separation){
    MPX_matrix LAMBDA0(U.Matrices.at(0).basis(),U.Matrices.at(0).Index(1),U.Lambdas.at(0));
    if (separation==0){
      return contract_to_sparse(contract(U.Matrices.at(1),1,contract(W2,0,contract(W1,0,contract(contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract0011),0,U.Matrices.at(1),0,contract11),0,contract21),0,contract20),0,contract0015),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0150).trace();
      //return contract_to_sparse(contract(U.Matrices.at(1),1,contract(W2,0,contract(W1,0,contract(contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract0011),0,U.Matrices.at(1),0,contract11),0,contract21),0,contract20).RemoveDummyIndices(std::vector<MPXInt>({{1,2,3,4}})),0,contract0011),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0110).trace();
      //return 0.0;
    }
    if (separation==1){
      return TwoVertexMeasurement(W1,W2,U.Matrices.at(0),U.Matrices.at(1),LAMBDA0);
    }
    else {
      MPX_matrix accumulator;
      if (separation % 2) //odd > 1
	//accumulator=contract(U.Matrices.at(1),1,contract(contract(W1,0,contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract11),0,contract0022).RemoveDummyIndices(std::vector<MPXInt>({{0,1}})),0,U.Matrices.at(1),0,contract11),0,contract0110);
	accumulator=contract(U.Matrices.at(1),1,contract(contract(W1,0,contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract11),0,contract0022).RemoveDummyIndices(std::vector<MPXInt>({1})),0,U.Matrices.at(1),0,contract21).RemoveDummyIndices(std::vector<MPXInt>({0})),0,contract0110);
      else //even > 0
	accumulator=contract(U.Matrices.at(1),1,contract(contract(contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract0011),0,U.Matrices.at(1),0,contract11),0,W1,0,contract12),0,contract0210).RemoveDummyIndices(std::vector<MPXInt>({{2,3}}));

      for (size_t i=1;i<separation/2;++i){
	accumulator=std::move(contract(U.Matrices.at(1),1,contract(contract(U.Matrices.at(0),1,contract(accumulator,0,U.Matrices.at(0),0,contract11),0,contract0110),0,U.Matrices.at(1),0,contract11),0,contract0110));
      }
      accumulator=std::move(contract(U.Matrices.at(1),1,contract(contract(contract(U.Matrices.at(0),1,contract(accumulator,0,U.Matrices.at(0),0,contract11),0,contract0110),0,U.Matrices.at(1),0,contract11),0,W2,0,contract12),0,contract0210));
      return contract_to_sparse(accumulator,0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0110).trace();
      //last, use lambda
    }
  }

  std::complex<double> SimpleEnergy(const MPO_matrix& LeftH,const MPO_matrix& RightH,const MPO_matrix& H1,const MPO_matrix& I,const UnitCell& Ortho){
    if (Ortho.Matrices.size()>2){std::cout << "Only two vertex basis accepted for now" <<std::endl; return 0.0;}
    //Need to make an MPX with the lambda
    MPX_matrix LAMBDA0(LeftH.GetPhysicalSpectrum(),Ortho.Matrices.at(0).Index(1),Ortho.Lambdas.at(0));
    //make the onsite/on chain hamiltonian MPO using the onsite (or on chain) part of the Hamiltonian
    //next step is probably a memory hog?
    std::complex<double> BondEnergy(TwoVertexMeasurement(LeftH,RightH,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0));
    //complex<double> ExcessEnergy(0.5*(TwoVertexMeasurement(H1,I,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0)+TwoVertexMeasurement(I,H1,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0)));
    std::complex<double> ExcessEnergy(0.5*(TwoVertexMeasurement(H1,I,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0)+OneVertexMeasurement(H1,Ortho)));
    return BondEnergy-ExcessEnergy;
  }

  std::complex<double> iTwoVertexEnergy(const MPO_matrix& H_MPO,const UnitCell& Ortho){
    const ajaj::MPO_matrix H1(H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(H_MPO.Index(1).size()-1,H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0))); //form the on-vertex part of a Hamiltonian
    const ajaj::MPO_matrix ColX(H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,H_MPO.Index(1).size()-2),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0)));
    const ajaj::MPO_matrix RowX(H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(H_MPO.Index(1).size()-1,H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(1,H_MPO.Index(3).size()-2)));
    return iTwoVertexEnergy(ColX,RowX,H1,Ortho);
  }
  
  std::complex<double> iTwoVertexEnergy(const MPO_matrix& ColX,const MPO_matrix& RowX,const MPO_matrix& H1,const UnitCell& Ortho){
    if (Ortho.Matrices.size()!=2){std::cout << "Only two vertex basis accepted for now" <<std::endl; return 0.0;}

    MPX_matrix LAMBDA0(H1.basis(),Ortho.Matrices.at(0).Index(1),Ortho.Lambdas.at(0));
    std::cout << "Calculating energy" <<std::endl;
    std::complex<double> LocalEnergy=0.0;
    if (!H1.empty()){
      LocalEnergy+=contract_to_sparse(contract(Ortho.Matrices.at(1),1,contract(contract(Ortho.Matrices.at(0),1,contract(H1,0,Ortho.Matrices.at(0),0,contract20),0,contract0013),0,Ortho.Matrices.at(1),0,contract31),0,contract0310),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0130).trace(); //contribution from 1st site
      LocalEnergy+=OneVertexMeasurement(H1,Ortho); //add contribution from 2nd site
    }
    std::cout << "On vertex energy =" << LocalEnergy <<std::endl;
    std::complex<double> Bond1Energy(TwoVertexMeasurement(RowX,ColX,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0));
    std::cout << "Bond 1 energy =" << Bond1Energy <<std::endl;
    //and now a more complicated thing...
    std::complex<double> Bond2Energy(contract_to_sparse(contract(Ortho.Matrices.at(1),1,contract(contract(Ortho.Matrices.at(0),1,contract(contract(Ortho.Matrices.at(0),0,contract(Ortho.Matrices.at(1),1,contract(contract(contract(Ortho.Matrices.at(0),1,Ortho.Matrices.at(0),0,contract0011),0,Ortho.Matrices.at(1),0,contract11),0,RowX,0,contract12),0,contract0210),0,contract11),0,ColX,0,contract0241),0,contract0311),0,Ortho.Matrices.at(1),0,contract11),0,contract0310),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0130).trace());
    std::cout << "Bond 2 energy = " << Bond2Energy <<std::endl;

    return 0.5*(LocalEnergy+Bond1Energy+Bond2Energy);

  }
  
  TransferMatrixComponents::TransferMatrixComponents(const UnitCell& KetCell, bool HV, const State S) : TransferMatrixComponents(KetCell,KetCell,HV,S) {}

  TransferMatrixComponents::TransferMatrixComponents(const UnitCell& BraCell, const UnitCell& KetCell, bool HV, const State S) : CellSize_(BraCell.size()), Hermitian_answer_(HV),TargetState_(S),BraCell_(BraCell), KetCell_(KetCell) {
#ifndef NDEBUG
    std::cout << "Constructing TransferMatrixComponents" << std::endl;
#endif
    if (BraCell.size()!=KetCell.size()){std::cout << "Bra and Ket fragments have different sizes! " << BraCell.size() << " " << KetCell.size() <<std::endl; exit(1);}
    else {
      for (std::pair<std::vector<MPS_matrix>::const_iterator,std::vector<MPS_matrix>::const_iterator> cits={BraCell.Matrices.begin(),KetCell.Matrices.begin()};cits.first!= BraCell.Matrices.end() && cits.second!= KetCell.Matrices.end(); ++cits.first, ++cits.second) {
	BraKetMatrixPtrs_.push_back({std::addressof(*(cits.first)),std::addressof(*(cits.second))});
      }
    }

    left_indices_.emplace_back(1,BraMatrix(0).getInwardMatrixIndex());
    left_indices_.emplace_back(0,KetMatrix(0).getInwardMatrixIndex());
    right_indices_.emplace_back(0,BraMatrix(last()).getOutwardMatrixIndex());
    right_indices_.emplace_back(1,KetMatrix(last()).getOutwardMatrixIndex());

    length_=left_indices_[0].size()*left_indices_[1].size();
    if (length_!=right_indices_[0].size()*right_indices_[1].size()){std::cout << "Index dimensions mismatch between ends of translation unit!" <<std::endl; exit(1);}

    //work out allowed indices
    const MPXIndex& mprimed(left_indices_[0]); 
    const MPXIndex& m(left_indices_[1]);
    vrows_=mprimed.size();
    vcols_=m.size();

    for (MPXInt c=0;c<m.size();++c){
      /*if (TargetStatePtr_){ //have we set a target state here? If so, do we have a match, are we in the target state sector?
	if (m[c] != *TargetStatePtr_){continue;}
	}*/
      for (MPXInt r=0;r<mprimed.size();++r){
	if (m[c]==mprimed[r]+TargetState_){ //check for consistency of charges	      
	  allowed_indices_.push_back(r+vrows_*c);
	  allowed_indices_dagger_.push_back(c+vcols_*r);
	  rows_and_cols_.push_back(std::array<MPXInt,2> {{r,c}});
	}
      }
    }
#ifndef NDEBUG
    std::cout << "End constructing TransferMatrixComponents" << std::endl;
#endif
  }

  MPX_matrix TransferMatrixComponents::accumulate_() const {
    std::vector<MPXPair> ContractPP({{MPXPair(BraMatrix(0).PhysicalIndexNumber(),KetMatrix(0).PhysicalIndexNumber())}});
    //contract will leave in form a'_i-1 a'_i a_i-1 a_i so reorder
    MPX_matrix accumulator(reorder(contract(BraMatrix(0),1,KetMatrix(0),0,ContractPP),0,reorder0213,2));
    for (uMPXInt i=1; i<=last();++i){
      std::vector<MPXPair> BraContract({{MPXPair(2,BraMatrix(i).InwardMatrixIndexNumber())}});
      std::vector<MPXPair> KetContract({{MPXPair(2,KetMatrix(i).InwardMatrixIndexNumber()),MPXPair(3,KetMatrix(i).PhysicalIndexNumber())}});
      accumulator=contract(contract(accumulator,0,BraMatrix(i),1,BraContract),0,KetMatrix(i),0,KetContract);
    }
    accumulator.ShiftNumRowIndices(2);
    return accumulator;
  }

  SparseED TransferMatrixComponents::LeftED(Sparseint numevals, char which[3],SparseMatrix* initial) const {
    std::cout <<"Eigensolver for matrix of length " << length() << std::endl;
    std::cout <<"Using reduced subspace of length " << allowed_indices_.size() <<std::endl;

    if (length()<=SPARSE_THRESHOLD){
      //use simple method
      std::cout << "Simple method" <<std::endl;
      //accumulate
      MPX_matrix TM(accumulate_());

      return TM.LeftEigs(TM.basis().Identity(),numevals,which,initial);
   
      //#ifndef NDEBUG
      //std::cout << "TM is consistent? " << TM.isConsistent() << std::endl;

      //MPX_matrix CheckTMSym=reorder(TM,1,reorder1032,2);
      //std::cout << "Checking Transfer matrix symmetry" << std::endl;
      //(TM.matrix()-CheckTMSym.matrix()).print();
      //#endif

      /*if (Hermitian_answer_){
	//massage answer into Hermitian form
	//collect old index dims
	std::vector<MPXInt> olddims;
	for (const auto& I : left_indices()){
	  olddims.push_back(I.size());
	}
	olddims.push_back(numevals);
	ans.EigenVectors.print();
	//ans.EigenVectors+=reshape(ans.EigenVectors,2,2,olddims,reorder102,1);		
	}*/
      //return ans;

    }
    else {
      SparseED ans(length(),numevals);
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices_.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];

      SparseVectorWithRestriction guess_struct(initial,&allowed_indices_);
      arpack::arpack_eigs<TransferMatrixComponents,SparseVectorWithRestriction> eigensystem(this,&LeftComponentsMultiply,allowed_indices_.size(),initial ? &guess_struct : nullptr, &ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);

      if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}

      if (Hermitian_answer_ && which[0]=='L' && (abs(imag(Evals[0]))>=IMAGTOL*abs(real(Evals[0])))){
	//we want the largest eval to be real, and it isn't quite
	//try a better guess...
	std::cout << "Dominant eigenvalue has non vanishing imaginary part, trying a new guess." <<std::endl;

	SparseMatrix gM(length(),1);
	for (MPXInt i=0;i<allowed_indices_.size();++i){
	  //here is where to drop tiny values...
	  if (abs(Evecs[i])>SPARSETOL){
	    gM.entry(allowed_indices_[i],0,Evecs[i]);
	    gM.entry(allowed_indices_dagger_[i],0,conj(Evecs[i]));
	  }
	}
	gM.finalise();
	SparseVectorWithRestriction newguess(&gM,&allowed_indices_);
	arpack::arpack_eigs<TransferMatrixComponents,SparseVectorWithRestriction> eigensystem2(this,&LeftComponentsMultiply,allowed_indices_.size(),&newguess,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
	if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}
      }

      for (size_t v=0;v<static_cast<size_t>(numevals);++v){
	std::cout << "Eval: " << Evals[v] <<std::endl;
	ans.Values.push_back(Evals[v]);
	
	for (MPXInt i=0;i<allowed_indices_.size();++i){
	  //here is where to drop tiny values...
	  if (abs(Evecs[i+v*allowed_indices_.size()])>SPARSETOL){
	    ans.EigenVectors.entry(allowed_indices_[i],v,Evecs[i+v*allowed_indices_.size()]);
	  }
	  if (v==0 && Hermitian_answer_) //wrecks normalisation, but makes reshape Hermitian to high precision
	    ans.EigenVectors.entry(allowed_indices_dagger_[i],v,conj(Evecs[i+v*allowed_indices_.size()]));
	}
	
      }
      ans.EigenVectors.finalise();
      delete[] Evecs;
      delete[] Evals;
      return ans;
    }

  }

  SparseED TransferMatrixComponents::RightED(Sparseint numevals, char which[3],SparseMatrix* initial) const {
    std::cout <<"Eigensolver for matrix of length " << length() << std::endl;
    std::cout <<"Using reduced subspace of length " << allowed_indices_.size() <<std::endl;

    if (length()<SPARSE_THRESHOLD){
      //use simple method
      std::cout << "Simple method" <<std::endl;
      //accumulate
      MPX_matrix TM(accumulate_());
      return TM.RightEigs(TM.basis().Identity(),numevals,which,initial);
      //std::cout << "TM is consistent? " << TM.isConsistent() << std::endl;
    }
    else {
      SparseED ans(length(),numevals);
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices_.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      
      SparseVectorWithRestriction guess_struct(initial,&allowed_indices_);
      arpack::arpack_eigs<TransferMatrixComponents,SparseVectorWithRestriction> eigensystem(this,&RightComponentsMultiply,allowed_indices_.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
      if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}
      for (size_t v=0;v<static_cast<size_t>(numevals);++v){
	std::cout << "Eval: " << Evals[v] <<std::endl;
	ans.Values.push_back(Evals[v]);
	
	for (MPXInt i=0;i<allowed_indices_.size();++i){
	  if(abs(Evecs[i+v*allowed_indices_.size()])>SPARSETOL){
	    ans.EigenVectors.entry(allowed_indices_[i],v,Evecs[i+v*allowed_indices_.size()]);
	  }
	  if (Hermitian_answer_)
	    ans.EigenVectors.entry(allowed_indices_dagger_[i],v,conj(Evecs[i+v*allowed_indices_.size()]));
	}
	
      }
      ans.EigenVectors.finalise();
      delete[] Evecs;
      delete[] Evals;
      return ans;
    }

  }

  void LeftComponentsMultiply(const TransferMatrixComponents* stuff, std::complex<double> *in, std::complex<double> *out){

    //convert *in to sparse
    SparseMatrix V(stuff->vcols(),stuff->vrows(),stuff->vrows());//undocumented behaviour of SparseMatrix, using transposed rows and cols to avoid unnecessary row ordering...
    for (struct {std::vector<std::array<Sparseint,2> >::const_iterator cit; uMPXInt idx;} itstruct ={stuff->rows_and_cols().begin(), 0};itstruct.cit!=stuff->rows_and_cols().end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    //loop over all matrices in the unitcells, doing contraction
    MPX_matrix accumulator(stuff->KetMatrix(0).basis(),stuff->left_indices(),1,V.cheap_no_transpose_finalise());

    for (uMPXInt s=0;s<stuff->last();++s){
      std::vector<MPXPair> contractK({{1,stuff->KetMatrix(s).InwardMatrixIndexNumber()}});
      std::vector<MPXPair> contractB({{stuff->BraMatrix(s).PhysicalIndexNumber(),1},{stuff->BraMatrix(s).InwardMatrixIndexNumber(),0}});
      accumulator=std::move(contract(stuff->BraMatrix(s),1,contract(accumulator,0,stuff->KetMatrix(s),0,contractK),0,contractB));	
    }

    std::vector<MPXPair> contractK({{1,stuff->KetMatrix(stuff->last()).InwardMatrixIndexNumber()}});
    std::vector<MPXPair> contractB({{stuff->BraMatrix(stuff->last()).PhysicalIndexNumber(),1},{stuff->BraMatrix(stuff->last()).InwardMatrixIndexNumber(),0}});
    DumbExtractWithZeros(reshape(contract_to_sparse(stuff->BraMatrix(stuff->last()),1,contract(accumulator,0,stuff->KetMatrix(stuff->last()),0,contractK),0,contractB),stuff->length()),stuff->allowed_indices(),out);

  }
  void RightComponentsMultiply(const TransferMatrixComponents* stuff, std::complex<double> *in, std::complex<double> *out){
//convert *in to sparse
    SparseMatrix V(stuff->vcols(),stuff->vrows(),stuff->vrows());//undocumented behaviour of SparseMatrix, using transposed rows and cols to avoid unnecessary row ordering...
    for (struct {std::vector<std::array<Sparseint,2> >::const_iterator cit; uMPXInt idx;} itstruct ={stuff->rows_and_cols().begin(), 0};itstruct.cit!=stuff->rows_and_cols().end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    //loop over all matrices in the unitcells, doing contraction
    MPX_matrix accumulator(stuff->KetMatrix(stuff->last()).basis(),stuff->right_indices(),1,V.cheap_no_transpose_finalise());
    for (uMPXInt s=stuff->last();s>0;--s){
      MPXInt KetPhys(stuff->KetMatrix(s).PhysicalIndexNumber());
      MPXInt KetOut(stuff->KetMatrix(s).OutwardMatrixIndexNumber());
      MPXInt BraPhys(stuff->BraMatrix(s).PhysicalIndexNumber());
      MPXInt BraOut(stuff->BraMatrix(s).OutwardMatrixIndexNumber());

      std::vector<MPXPair> contractK({{KetOut,1}});
      std::vector<MPXPair> contractB({{BraPhys,KetPhys},{BraOut,KetOut}});
      accumulator=std::move(contract(stuff->BraMatrix(s),1,contract(stuff->KetMatrix(s),0,accumulator,0,contractK),0,contractB));	
    }
    MPXInt KetPhys(stuff->KetMatrix(0).PhysicalIndexNumber());
    MPXInt KetOut(stuff->KetMatrix(0).OutwardMatrixIndexNumber());
    MPXInt BraPhys(stuff->BraMatrix(0).PhysicalIndexNumber());
    MPXInt BraOut(stuff->BraMatrix(0).OutwardMatrixIndexNumber());
    std::vector<MPXPair> contractK({{KetOut,1}});
    std::vector<MPXPair> contractB({{BraPhys,KetPhys},{BraOut,KetOut}});
    DumbExtractWithZeros(reshape(contract_to_sparse(stuff->BraMatrix(0),1,contract(stuff->KetMatrix(0),0,accumulator,0,contractK),0,contractB),stuff->length()),stuff->allowed_indices(),out);

  }

  double MultiVertexEntropy(const UnitCell& U,uMPXInt v /* number of vertices */){
    //contract together and include the correct lambda at end
    if (!U.Matrices.size() || !U.Lambdas.size() || v<1 || v>2) return -1.0; //error!
    MPX_matrix acc(U.Matrices.back().basis(),U.Matrices.back().Index(2),U.Lambdas.at(0));

    if (v==1){
      MPXDecomposition U1D(U.Matrices.at(0).SVD());
      return entropy(contract(contract(U.Matrices.at(1),0,contract(MPX_matrix(U1D.RowMatrix.basis(),U1D.RowMatrix.Index(0),U1D.Values),0,U1D.RowMatrix,0,contract10),0,contract11),0,acc,0,contract10).ShiftNumRowIndices(1).SVD().Values);
    }

    acc=std::move(contract(U.Matrices.back(),0,acc,0,contract20));

    for (size_t m=1; m<v; ++m){
      //normal contraction
      acc=std::move(contract(U.Matrices.at(U.Matrices.size()-m-1),0,acc,0,contract21));

    }
    return entropy(reorder(acc,0,reorder0213,2).SVD().Values);

    //reshape so that physical indices are rows, matrix indices are cols
    
    //svd

    //get entropy
  }

  bool check_for_H_MPO_file(const std::string& name, uMPXInt sliceindex){
    std::stringstream SliceHMPOnamestream;
    SliceHMPOnamestream << name << "_" << sliceindex << "_H.MPO_matrix";
    
    //check file exists and can be opened
    
    std::ifstream fs;
    
    fs.open(SliceHMPOnamestream.str().c_str(),ios::in);
    
    if (fs.is_open()){
      return 1;
    }
    
    return 0;
  }

  std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Bra, const ConstFiniteMPS& Ket, const std::vector<meas_pair>& ops){
    if (Bra.size()!=Ket.size()){
      std::cout << "States defined for different length systems! Orthogonal but probably not intentional, aborting." <<std::endl<<std::endl;
      exit(1);
    }

    std::vector<std::vector<const MPO_matrix*> > ordered_ops(Ket.size(),std::vector<const MPO_matrix*>());

    for (auto& op_v : ops){
      if (op_v.MPO_ptr()){
	ordered_ops[op_v.position()-1].emplace_back(op_v.MPO_ptr());
      }
    }

    for (auto& o :ordered_ops) { //reverse order as user will expect the rightmost arg to be applied first
      std::reverse(o.begin(),o.end());
    }

    MPX_matrix accumulator(Ket.matrix(1).left_shape());
    for (auto& op1 : ordered_ops.at(0)){
      accumulator=std::move(contract(*op1,0,accumulator,0,contract20).RemoveDummyIndices({{1,2}}).ShiftNumRowIndices(2)); //puts it back into left shaped MPS_matrix form
    }
    accumulator=std::move(contract(Bra.matrix(1).left_shape(),1,accumulator,0,contract0011));
    for (uMPXInt i=2;i<=Ket.size();++i){
      accumulator=std::move(contract(accumulator,0,Ket.matrix(i).left_shape(),0,contract11));
      if (ordered_ops.at(i-1).size()) {
	bool first=1;
	for (auto& opi : ordered_ops.at(i-1)){
	  accumulator=std::move(contract(*opi,0,accumulator,0,first ? contract21 : contract20)).RemoveDummyIndices({{1,2}});
	  first=0;
	}
	accumulator=std::move(contract(Bra.matrix(i).left_shape(),1,accumulator,0,contract0011));
      }
      else {
	accumulator=std::move(contract(Bra.matrix(i).left_shape(),1,accumulator,0,contract0110));
      }
    }
    return accumulator.Trace();
  }

    std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Ket, const std::vector<meas_pair>& ops){

    std::vector<std::vector<const MPO_matrix*> > ordered_ops(Ket.size(),std::vector<const MPO_matrix*>());

    for (auto& op_v : ops){
      if (op_v.MPO_ptr()){
	ordered_ops[op_v.position()-1].emplace_back(op_v.MPO_ptr());
      }
    }

    for (auto& o :ordered_ops) { //reverse order as user will expect the rightmost arg to be applied first
      std::reverse(o.begin(),o.end());
    }

    MPX_matrix accumulator(Ket.matrix(1).left_shape());
    for (auto& op1 : ordered_ops.at(0)){
      accumulator=std::move(contract(*op1,0,accumulator,0,contract20).RemoveDummyIndices({{1,2}}).ShiftNumRowIndices(2)); //puts it back into left shaped MPS_matrix form
    }
    accumulator=std::move(contract(Ket.matrix(1).left_shape(),1,accumulator,0,contract0011));
    for (uMPXInt i=2;i<=Ket.size();++i){
      accumulator=std::move(contract(accumulator,0,Ket.matrix(i).left_shape(),0,contract11));
      if (ordered_ops.at(i-1).size()) {
	bool first=1;
	for (auto& opi : ordered_ops.at(i-1)){
	  accumulator=std::move(contract(*opi,0,accumulator,0,first ? contract21 : contract20)).RemoveDummyIndices({{1,2}});
	  first=0;
	}
	accumulator=std::move(contract(Ket.matrix(i).left_shape(),1,accumulator,0,contract0011));
      }
      else {
	accumulator=std::move(contract(Ket.matrix(i).left_shape(),1,accumulator,0,contract0110));
      }
    }
    return accumulator.Trace();
  }

  
  std::complex<double> GeneralisedOverlap(const ConstFiniteMPS& Bra, const ConstFiniteMPS& Ket, const MPO_matrix& Op, uMPXInt v){
    //zip through the states and contract
    //if Op is defined, apply it at site v.

    //check states have same length!

    if (Bra.size()!=Ket.size()){
      std::cout << "States defined for different length systems! Orthogonal but probably not intentional, aborting." <<std::endl<<std::endl;
      exit(1);
    }
    
    if (!Op.empty()){
      if (v<1 || v>Ket.size()){
	std::cout <<"Specified vertex for application of operator is out of bounds!"<<std::endl<<std::endl;
	exit(1);
      }
    }

    //treat first as special case
    MPX_matrix accumulator(Ket.matrix(1).left_shape());
    if (v==1 && !Op.empty()){
      accumulator=std::move(contract(Op,0,accumulator,0,contract20).RemoveDummyIndices({{1,2}}).ShiftNumRowIndices(2)); //puts it back into left shaped MPS_matrix form
    }
    accumulator=std::move(contract(Bra.matrix(1).left_shape(),1,accumulator,0,contract0011));
    for (uMPXInt i=2;i<=Ket.size();++i){
      accumulator=std::move(contract(accumulator,0,Ket.matrix(i).left_shape(),0,contract11));
      if (i==v){
	accumulator=std::move(contract(Op,0,accumulator,0,contract21)).RemoveDummyIndices({{1,2}});
	accumulator=std::move(contract(Bra.matrix(i).left_shape(),1,accumulator,0,contract0011));
      }
      else {
	accumulator=std::move(contract(Bra.matrix(i).left_shape(),1,accumulator,0,contract0110));
      }
    }

    //finish off. Final index should be a dummy, so really this is a trace of a 1x1.
    
    return accumulator.Trace();
    
  }

  
}

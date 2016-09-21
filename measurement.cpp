#include <vector>

#include "ajaj_common.hpp"
#include "sparse_interface.hpp"
#include "arpack_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"
#include "measurement.hpp"

namespace ajaj {

  TransferMatrixParts::TransferMatrixParts(const UnitCell& B,const UnitCell& K,const State* T) : BraCell(B),KetCell(K),TargetStatePtr(T){m_init();}
  TransferMatrixParts::TransferMatrixParts(const UnitCell& C,const State* T) : BraCell(C),KetCell(C),TargetStatePtr(T){m_init();}

  void TransferMatrixParts::m_init(){
    //check unit cell lengths
    if (BraCell.size() != KetCell.size()) {std::cout << "Malformed unit cells" << std::endl; exit(1);}
    if (BraCell.Matrices.at(0).Index(1).size()!= BraCell.Matrices.back().Index(2).size()){std::cout << "Malformed Bra unit cell" << std::endl;exit(1);}
    if (KetCell.Matrices.at(0).Index(1).size()!= KetCell.Matrices.back().Index(2).size()){std::cout << "Malformed Ket unit cell" << std::endl;exit(1);}
    left_indices.emplace_back(1,BraCell.Matrices.at(0).Index(1));
    left_indices.emplace_back(0,KetCell.Matrices.at(0).Index(1));
    right_indices.emplace_back(0,BraCell.Matrices.back().Index(2));
    right_indices.emplace_back(1,KetCell.Matrices.back().Index(2));

    m_length=left_indices[0].size()*left_indices[1].size();

    //work out allowed indices
    const MPXIndex& mprimed(left_indices[0]); 
    const MPXIndex& m(left_indices[1]);
    vrows=mprimed.size();
    vcols=m.size();

    for (Sparseint c=0;c<m.size();++c){
      if (TargetStatePtr){ //have we set a target state here? If so, do we have a match, are we in the target state sector?
	if (m[c] != *TargetStatePtr){continue;}
      }
      for (Sparseint r=0;r<mprimed.size();++r){
	if (m[c]==mprimed[r]){ //check for consistency of charges	      
	  allowed_indices.push_back(r+vrows*c);
	  allowed_indices_dagger.push_back(c+vrows*r);
	  rows_and_cols.push_back(std::array<Sparseint,2> {{r,c}});
	}
      }
    }
  }

  SparseED TransferMatrixParts::LeftED(Sparseint numevals, char which[3],SparseMatrix* initial){
    //set up workspace (Evals, Evecs) SparseHED
    //MPXInt fulldim(m_length);
    
    if (m_length<500){
      //use dense method
      MPX_matrix TM(MakeLTransferMatrix(BraCell,KetCell));
      return TM.LeftEigs(numevals,which,initial);
    }
    else {
      SparseED ans(m_length,numevals);
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      
      SparseVectorWithRestriction guess_struct(initial,&allowed_indices);
      arpack::arpack_eigs<TransferMatrixParts,SparseVectorWithRestriction> eigensystem(this,&LeftTransferMatrixMultiply,allowed_indices.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
      if (eigensystem.error_status()) {std::cout << "error with tensor arpack" << std::endl;exit(1);}
      for (size_t v=0;v<static_cast<size_t>(numevals);++v){
	std::cout << "Eval: " << Evals[v] <<std::endl;
	ans.Values.push_back(Evals[v]);
	
	for (Sparseint i=0;i<allowed_indices.size();++i){
	  ans.EigenVectors.entry(allowed_indices[i],v,Evecs[i+v*allowed_indices.size()]);
	  ans.EigenVectors.entry(allowed_indices_dagger[i],v,conj(Evecs[i+v*allowed_indices.size()]));
	}
	
      }
      ans.EigenVectors.finalise();
      delete[] Evecs;
      delete[] Evals;
      return ans;
    }
  }

  SparseED TransferMatrixParts::RightED(Sparseint numevals, char which[3],SparseMatrix* initial){
    //set up workspace (Evals, Evecs) SparseHED
    //MPXInt fulldim(m_length);
    SparseED ans(m_length,numevals);
    std::complex<double>* Evecs = new std::complex<double>[allowed_indices.size()*numevals];
    std::complex<double>* Evals = new std::complex<double>[numevals];
    SparseVectorWithRestriction guess_struct(initial,&allowed_indices);
    arpack::arpack_eigs<TransferMatrixParts,SparseVectorWithRestriction> eigensystem(this,&RightTransferMatrixMultiply,allowed_indices.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
    if (eigensystem.error_status()) {std::cout << "error with tensor arpack" << std::endl;exit(1);}
    for (size_t v=0;v<static_cast<size_t>(numevals);++v){
      std::cout << "Eval: " << Evals[v] <<std::endl;
      ans.Values.push_back(Evals[v]);
     
      for (Sparseint i=0;i<allowed_indices.size();++i){
	ans.EigenVectors.entry(allowed_indices[i],v,Evecs[i+v*allowed_indices.size()]);
	ans.EigenVectors.entry(allowed_indices_dagger[i],v,conj(Evecs[i+v*allowed_indices.size()]));
      }
    }
    ans.EigenVectors.finalise();
    delete[] Evecs;
    delete[] Evals;
    return ans;
  }

  void LeftTransferMatrixMultiply(TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out){
    //convert *in to sparse
    SparseMatrix V(stuff->vcols,stuff->vrows,stuff->vrows);//undocumented behaviour of SparseMatrix, using transposed rows and cols to avoid unnecessary row ordering...
    for (struct {std::vector<std::array<Sparseint,2> >::const_iterator cit; uMPXInt idx;} itstruct ={stuff->rows_and_cols.begin(), 0};itstruct.cit!=stuff->rows_and_cols.end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    //loop over all matrices in the unitcells, doing contraction
    MPX_matrix accumulator(stuff->KetCell.Matrices.at(0).GetPhysicalSpectrum(),stuff->left_indices,1,V.cheap_no_transpose_finalise());

    for (uMPXInt s=0;s<stuff->BraCell.size()-1;++s){
      accumulator=std::move(contract(stuff->BraCell.Matrices.at(s),1,contract(accumulator,0,stuff->KetCell.Matrices.at(s),0,contract11),0,contract0110));	
    }

    DumbExtractWithZeros(reshape(contract_to_sparse(stuff->BraCell.Matrices.at(stuff->BraCell.size()-1),1,contract(accumulator,0,stuff->KetCell.Matrices.at(stuff->KetCell.size()-1),0,contract11),0,contract0110),stuff->length()),stuff->allowed_indices,out);
  };

  void RightTransferMatrixMultiply(TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out){
    //convert *in to sparse
    SparseMatrix V(stuff->vcols,stuff->vrows,stuff->vrows);//undocumented behaviour of SparseMatrix, using transposed rows and cols to avoid unnecessary row ordering...
    for (struct {std::vector<std::array<Sparseint,2> >::const_iterator cit; uMPXInt idx;} itstruct ={stuff->rows_and_cols.begin(), 0};itstruct.cit!=stuff->rows_and_cols.end();++itstruct.cit,++itstruct.idx){
      V.entry((*(itstruct.cit))[1],(*(itstruct.cit))[0],in[itstruct.idx]);
    }
    //loop over all matrices in the unitcells, doing contraction
    MPX_matrix accumulator(stuff->KetCell.Matrices.back().GetPhysicalSpectrum(),stuff->right_indices,1,V.cheap_no_transpose_finalise());
    for (uMPXInt s=stuff->BraCell.size()-1;s>0;--s){
      accumulator=std::move(contract(stuff->BraCell.Matrices.at(s),1,contract(accumulator,0,stuff->KetCell.Matrices.at(s),0,contract12),0,contract0120));	
    }
    DumbExtractWithZeros(reshape(contract_to_sparse(stuff->BraCell.Matrices.at(0),1,contract(accumulator,0,stuff->KetCell.Matrices.at(0),0,contract12),0,contract0120),stuff->length()),stuff->allowed_indices,out);
  };

  UnitCell Orthogonalise(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda){
    std::cout << "Orthogonalising..." << std::endl;
    const EigenStateArray& spectrum(MPSD.LeftMatrix.GetPhysicalSpectrum());
    MPS_matrix NLeft(contract(contract(MPX_matrix(spectrum,MPSD.RightMatrix.Index(0),MPSD.Values),0,MPSD.RightMatrix,0,contract10),0,MPX_matrix(spectrum,MPSD.RightMatrix.Index(2),PreviousLambda,1),0,contract20));
    MPX_matrix TransferMatrix(std::move(contract(contract(contract(MPSD.LeftMatrix,1,MPSD.LeftMatrix,0,contract00),0,NLeft,1,contract10),0,NLeft,0,contract2031).ShiftNumRowIndices(2)));
    SparseED LeftTdecomp(contract(contract(contract(MPSD.LeftMatrix,1,MPSD.LeftMatrix,0,contract00),0,NLeft,1,contract10),0,NLeft,0,contract2031).ShiftNumRowIndices(2).LeftEigs(State(spectrum[0].getChargeRules()),1,LARGEST));
    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this
    //make it into an MPX_matrix
    std::vector<double> scalevecX(MPSD.LeftMatrix.Index(1).size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,MPSD.LeftMatrix.Index(1));
    VLIndices.emplace_back(0,MPSD.LeftMatrix.Index(1));
    //decompose it into Xdagger X form
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(spectrum,VLIndices,1,reshape(LeftTdecomp.EigenVectors,MPSD.LeftMatrix.Index(1).size()))));
    MPX_matrix X(contract(MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));
    //record scale factor for later
    //now need to make Y
    ShiftLTransferMatrixToR(TransferMatrix);
    SparseED RightTdecomp(TransferMatrix.RightEigs(State(spectrum[0].getChargeRules()),1,LARGEST));

    std::cout << "Leading right eigenvalue of right transfer matrix: " << RightTdecomp.Values.at(0) << std::endl; //will need to rescale by this if used
    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,MPSD.RightMatrix.Index(2));
    VRIndices.emplace_back(0,MPSD.RightMatrix.Index(2));
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(spectrum,VRIndices,1,reshape(RightTdecomp.EigenVectors,MPSD.RightMatrix.Index(2).size()))));
    //to make Y we need to reorder and conjugate VRDecomp.second()
    MPX_matrix Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));
    //need to form and svd  new lambda
    MPXDecomposition XLYDecomp(contract(contract(X,0,MPX_matrix(spectrum,MPSD.RightMatrix.Index(2),PreviousLambda),0,contract10),0,Y,0,contract10).SVD());
    //Udagger X A
    //NLeft Xinv U
    //new lambda will need to be scaled to correct normalisation
    SquareSumRescale(XLYDecomp.Values,1.0);
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;

    UnitCell ans(spectrum);
    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(spectrum,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,MPSD.LeftMatrix,0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(reorder(contract(NLeft,0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20),0,reorder102,2));
    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(MPSD.Values);
    return ans;
  }

  UnitCell Orthogonalise(const UnitCell& C){
    std::cout << "Orthogonalising..." << std::endl;
    const EigenStateArray& spectrum(C.Matrices.at(0).GetPhysicalSpectrum());
    std::cout << "Forming left transfer matrix" << std::endl;
    MPX_matrix TransferMatrix(MakeLTransferMatrix(C,C));
    SparseED LeftTdecomp(TransferMatrix.LeftEigs(State(spectrum[0].getChargeRules()),1,LARGEST));
    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this
    //make it into an MPX_matrix
    std::vector<double> scalevecX(C.Matrices.at(0).Index(1).size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,C.Matrices.at(0).Index(1));
    VLIndices.emplace_back(0,C.Matrices.at(0).Index(1));
    //decompose it into Xdagger X form
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(spectrum,VLIndices,1,reshape(LeftTdecomp.EigenVectors,C.Matrices.at(0).Index(1).size()).massage())));
    MPX_matrix X(contract(MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));
    //record scale factor for later
    //now need to make Y
    std::cout << "Forming right transfer matrix" << std::endl;
    ShiftLTransferMatrixToR(TransferMatrix); //assumes inversion symmetry!!!!!
    SparseED RightTdecomp(TransferMatrix.RightEigs(State(spectrum[0].getChargeRules()),1,LARGEST));
    std::cout << "Leading right eigenvalue of right transfer matrix: " << RightTdecomp.Values.at(0) << std::endl; //will need to rescale by this if used
    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,C.Matrices.back().Index(2));
    VRIndices.emplace_back(0,C.Matrices.back().Index(2));
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(spectrum,VRIndices,1,reshape(RightTdecomp.EigenVectors,C.Matrices.back().Index(2).size()).massage())));
    //to make Y we need to reorder and conjugate VRDecomp.second()
    MPX_matrix Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));
    //need to form and svd  new lambda
    MPXDecomposition XLYDecomp(contract(contract(X,0,MPX_matrix(spectrum,C.Matrices.back().Index(2),C.Lambdas.at(0)),0,contract10),0,Y,0,contract10).SVD());
    //Udagger X A
    //new lambda will need to be scaled to correct normalisation
    SquareSumRescale(XLYDecomp.Values,1.0);
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;
    UnitCell ans(spectrum);
    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(spectrum,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,C.Matrices.at(0),0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(contract(C.Matrices.back(),0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20));
    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(C.Lambdas.back());
    return ans;
  }

  UnitCell OrthogonaliseInversionSymmetric(const UnitCell& C){
    std::cout << "Orthogonalising..." << std::endl;
    const EigenStateArray& spectrum(C.Matrices.at(0).GetPhysicalSpectrum());
    SparseED LeftTdecomp(TransferMatrixParts(C).LeftED(1,LARGEST));
    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this

    if (abs(imag(LeftTdecomp.Values.at(0)))>=1.0e-14) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      exit(1);
    }
    //make it into an MPX_matrix
    std::vector<double> scalevecX(C.Matrices.at(0).Index(1).size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,C.Matrices.at(0).Index(1));
    VLIndices.emplace_back(0,C.Matrices.at(0).Index(1));
    //decompose it into Xdagger X form
    //std::cout << "XDagger X"<< std::endl; //will need to rescale by this
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(spectrum,VLIndices,1,reshape(LeftTdecomp.EigenVectors,C.Matrices.at(0).Index(1).size()).massage())));
    MPX_matrix X(contract(MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));
    //record scale factor for later
    //now need to make Y
    //if we have inversion symmetry, then we shouldn't need to worry about solving again
    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,C.Matrices.back().Index(2));
    VRIndices.emplace_back(0,C.Matrices.back().Index(2));
    //std::cout << "Y YDagger"<< std::endl; //will need to rescale by this
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(spectrum,VRIndices,1,reshape(LeftTdecomp.EigenVectors,C.Matrices.back().Index(2).size()).massage())));
    //to make Y we need to reorder and conjugate VRDecomp.second()
    MPX_matrix Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));
    //need to form and svd  new lambda
    MPXDecomposition XLYDecomp(contract(contract(X,0,MPX_matrix(spectrum,C.Matrices.back().Index(2),C.Lambdas.at(0)),0,contract10),0,Y,0,contract10).SVD());
    //Udagger X A
    //new lambda will need to be scaled to correct normalisation
    SquareSumRescale(XLYDecomp.Values,1.0);
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;
    UnitCell ans(spectrum);
    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(spectrum,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,C.Matrices.at(0),0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(contract(C.Matrices.back(),0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20));
    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(C.Lambdas.back());
    return ans;
  }

  UnitCell OrthogonaliseInversionSymmetric(const MPSDecomposition& MPSD,const std::vector<double>& PreviousLambda){
    //makes a copy, in order to preserve the A Lambda B decomposition
    return OrthogonaliseInversionSymmetric(UnitCell(MPSD,PreviousLambda));
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
    MPX_matrix accumulator(contract(bra.Matrices.at(0),1,ket.Matrices.at(0),0,contract00));
    for (uMPXInt i=1;i<bra.Matrices.size();++i){
      accumulator=std::move(contract(contract(accumulator,0,bra.Matrices.at(i),1,contract11),0,ket.Matrices.at(i),0,contract2130));
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
    return MakeLTransferMatrix(bra,ket).LeftEigs(1,LARGEST).Values.at(0);
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

  std::complex<double> SophisticatedEnergy(const MPO_matrix& ColX,const MPO_matrix& RowX,const MPO_matrix& H1,const UnitCell& Ortho){
    if (Ortho.Matrices.size()>2){std::cout << "Only two vertex basis accepted for now" <<std::endl; return 0.0;}
    MPX_matrix LAMBDA0(H1.basis(),Ortho.Matrices.at(0).Index(1),Ortho.Lambdas.at(0));

    std::complex<double> LocalEnergy(contract_to_sparse(contract(Ortho.Matrices.at(1),1,contract(contract(Ortho.Matrices.at(0),1,contract(H1,0,Ortho.Matrices.at(0),0,contract20),0,contract0013),0,Ortho.Matrices.at(1),0,contract31),0,contract0310),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0130).trace());
    LocalEnergy+=OneVertexMeasurement(H1,Ortho);
    std::complex<double> Bond1Energy(TwoVertexMeasurement(RowX,ColX,Ortho.Matrices.at(0),Ortho.Matrices.at(1),LAMBDA0));
    //and now a more complicated thing...
    std::complex<double> Bond2Energy(contract_to_sparse(contract(Ortho.Matrices.at(1),1,contract(contract(Ortho.Matrices.at(0),1,contract(contract(Ortho.Matrices.at(0),0,contract(Ortho.Matrices.at(1),1,contract(contract(contract(Ortho.Matrices.at(0),1,Ortho.Matrices.at(0),0,contract0011),0,Ortho.Matrices.at(1),0,contract11),0,RowX,0,contract12),0,contract0210),0,contract11),0,ColX,0,contract0241),0,contract0311),0,Ortho.Matrices.at(1),0,contract11),0,contract0310),0,contract(LAMBDA0,0,LAMBDA0,0,contract10),0,contract0130).trace());

    return 0.5*(LocalEnergy+Bond1Energy+Bond2Energy);

  }


  void MultiVertexMeasurement::link_(const MPO_matrix* Op, const MPX_matrix& A) {
    if (Op){
      //make sure to strip 'dummy' indices by contract to sparse and rebuild indices
      //we assume single site operators, so MPO _matrix_ indices are dummies...
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,A.Index(2));
      indices.emplace_back(0,A.Index(2));
      T_=std::move(MPX_matrix(A.GetPhysicalSpectrum(),indices,1,contract_to_sparse(A,1,contract(*Op,0,contract(T_,0,A,0,contract11),0,contract21),0,contract0013)));
    }
    else {
      T_=std::move(contract(A,1,contract(T_,0,A,0,contract11),0,contract0110));
    }
  }
  
  void MultiVertexMeasurement::start_chain_(const MPO_matrix* Op, const MPX_matrix& A) {
    if (Op){
      //make sure to strip 'dummy' indices by contract to sparse and rebuild indices
      //we assume single site operators, so MPO _matrix_ indices are dummies...
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,A.Index(2));
      //indices.emplace_back(0,Op->Index(3));
      indices.emplace_back(0,A.Index(2));
      T_=std::move(MPX_matrix(A.GetPhysicalSpectrum(),indices,1,contract_to_sparse(A,1,contract(*Op,0,A,0,contract20),0,contract0013)));
    }
    else {
      T_=std::move(contract(A,1,A,0,contract0011));
    }
  }
  
  void MultiVertexMeasurement::finish_chain_(const MPX_matrix& Lambda) {
    std::cout << "Finishing measurement contraction chain" << std::endl;
    Result_=contract_to_sparse(Lambda,0,contract(T_,0,Lambda,0,contract10),0,contract0110).trace();
    std::cout << "Done" << std::endl;

  }

}

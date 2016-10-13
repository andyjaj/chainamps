#include <vector>

#include "ajaj_common.hpp"
#include "sparse_interface.hpp"
#include "arpack_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"
#include "measurement.hpp"

namespace ajaj {

  static const double IMAGTOL(100.0*std::numeric_limits<double>::epsilon());

  TransferMatrixParts::TransferMatrixParts(const UnitCell& B,const UnitCell& K,const State* T) : BraCell(B),KetCell(K),TargetStatePtr(T){
    if (BraCell.size() != KetCell.size()) {
      std::cout << "Malformed unit cells" << std::endl;
      exit(1);
    }
    m_init();
  }
  
  TransferMatrixParts::TransferMatrixParts(const UnitCell& C,const State* T) : BraCell(C),KetCell(C),TargetStatePtr(T){m_init();}

  void TransferMatrixParts::m_init(){
    //check unit cell lengths
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

  SparseED TransferMatrixParts::LeftED(Sparseint numevals, char which[3],SparseMatrix* initial) const {    
    std::cout <<"Eigensolver for matrix of length " << this->length() << std::endl;
    std::cout <<"Using reduced subspace of length " << allowed_indices.size() <<std::endl;
    if (m_length<1000){
      //use simple method
      std::cout << "Simple method" <<std::endl;
      MPX_matrix TM(MakeLTransferMatrix(BraCell,KetCell));
      return TM.LeftEigs(numevals,which,initial);
    }
    else {
      SparseED ans(m_length,numevals);
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      
      SparseVectorWithRestriction guess_struct(initial,&allowed_indices);
      arpack::arpack_eigs<TransferMatrixParts,SparseVectorWithRestriction> eigensystem(this,&LeftTransferMatrixMultiply,allowed_indices.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
      if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}
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

  SparseED TransferMatrixParts::RightED(Sparseint numevals, char which[3],SparseMatrix* initial) const {
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

  void LeftTransferMatrixMultiply(const TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out){
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

  void RightTransferMatrixMultiply(const TransferMatrixParts* stuff, std::complex<double> *in, std::complex<double> *out){
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

    const Basis& basis(MPSD.RightMatrix.basis());

    MPX_matrix PreviousLambdaInverse(basis,MPSD.RightMatrix.Index(2),PreviousLambda,1);
    MPX_matrix CurrentLambda(basis,MPSD.LeftMatrix.Index(2),MPSD.Values);

    //note rightmatrix is stored as a_i-1 sigma_i a_i format...
    //std::vector<MPXPair> contractL({{MPXPair(MPSD.RightMatrix.InwardMatrixIndexNumber(),1)}});
    //MPS_matrix NL(contract(contract(MPSD.RightMatrix,0,CurrentLambda,0,contractL),0,PreviousLambdaInverse,0,contract10));
    //std::vector<MPXPair> contractL({{MPXPair(1,MPSD.RightMatrix.InwardMatrixIndexNumber())}});
    //MPS_matrix NL(MPS_matrix(contract(CurrentLambda,0,contract(MPSD.RightMatrix,0,PreviousLambdaInverse,0,contract20),0,contractL)).left_shape());

    //std::vector<MPXPair> contractR({{MPXPair(MPSD.LeftMatrix.InwardMatrixIndexNumber(),1)}});
    //MPS_matrix NR=(contract(contract(MPSD.LeftMatrix,0,PreviousLambdaInverse,0,contractR),0,CurrentLambda,0,contract10));
    //std::vector<MPXPair> contractR({{MPXPair(1,MPSD.LeftMatrix.InwardMatrixIndexNumber())}});
    //MPS_matrix NR=(contract(PreviousLambdaInverse,0,contract(MPSD.LeftMatrix,0,CurrentLambda,0,contract20),0,std::vector<MPXPair>({{1,MPSD.LeftMatrix.InwardMatrixIndexNumber()}})));

    std::vector<MPXPair> contractL({{MPXPair(1,MPSD.RightMatrix.InwardMatrixIndexNumber())}});
    MPXDecomposition L(MPS_matrix(contract(CurrentLambda,0,MPSD.RightMatrix,0,contractL)).left_shape().SVD());
    //MPS_matrix A2(std::move(L.ColumnMatrix));
    //MPX_matrix P_(contract(MPX_matrix(basis,L.RowMatrix.Index(0),L.Values),0,L.RowMatrix,0,contract10));
    MPS_matrix NL(contract(contract(L.ColumnMatrix,0,contract(MPX_matrix(basis,L.RowMatrix.Index(0),L.Values),0,L.RowMatrix,0,contract10),0,contract20),0,PreviousLambdaInverse,0,contract20));

    SparseED LeftTdecomp(TransferMatrixComponents(std::vector<const MPS_matrix*>({{&MPSD.LeftMatrix,&NL}}),1).LeftED(1,LARGESTMAGNITUDE));

    std::cout << "Leading left eigenvalue of left transfer matrix: " << LeftTdecomp.Values.at(0) << std::endl; //will need to rescale by this
    if (abs(imag(LeftTdecomp.Values.at(0)))>=IMAGTOL*abs(real(LeftTdecomp.Values.at(0)))) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      return UnitCell(basis);
    }

    std::vector<MPXPair> contractR({{MPXPair(MPSD.LeftMatrix.OutwardMatrixIndexNumber(),0)}});
    MPXDecomposition R(MPS_matrix(contract(MPSD.LeftMatrix,0,CurrentLambda,0,contractR)).right_shape().SVD());
    //MPS_matrix B2(std::move(R.RowMatrix)); //rightshaped
    //MPX_matrix Q_(contract(R.ColumnMatrix,0,MPX_matrix(basis,R.ColumnMatrix.Index(1),R.Values),0,contract10));
    MPS_matrix NR(contract(PreviousLambdaInverse,0,contract(contract(R.ColumnMatrix,0,MPX_matrix(basis,R.ColumnMatrix.Index(1),R.Values),0,contract10),0,R.RowMatrix,0,contract10),0,contract10));

    SparseED RightTdecomp(TransferMatrixComponents(std::vector<const MPS_matrix*>({{&NR,&MPSD.RightMatrix}}),1).RightED(1,LARGESTMAGNITUDE));

    std::cout << "Leading right eigenvalue of right transfer matrix: " << RightTdecomp.Values.at(0) << std::endl;
    if (abs(imag(RightTdecomp.Values.at(0)))>=IMAGTOL*abs(real(RightTdecomp.Values.at(0)))) {
      std::cout << "Eigenvalue has non negligible imaginary part, numerical error in transfer matrix contraction?" << std::endl;
      return UnitCell(basis);
    }

    //decompose into Xdagger X form
    std::vector<MPXIndex> VLIndices;
    VLIndices.emplace_back(1,MPSD.LeftMatrix.getInwardMatrixIndex());
    VLIndices.emplace_back(0,MPSD.LeftMatrix.getInwardMatrixIndex());
    //decompose it into Xdagger X form
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(basis,VLIndices,1,reshape(LeftTdecomp.EigenVectors/*.ExtractColumns(std::vector<MPXInt>({0}))*/,MPSD.LeftMatrix.getInwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VLDecomp.first.size()) return UnitCell(basis);

    MPX_matrix X(contract(MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));

    std::vector<MPXIndex> VRIndices;
    VRIndices.emplace_back(1,MPSD.RightMatrix.getOutwardMatrixIndex());
    VRIndices.emplace_back(0,MPSD.RightMatrix.getOutwardMatrixIndex());
    std::pair<std::vector<double>,MPX_matrix> VRDecomp(SqrtDR(MPX_matrix(basis,VRIndices,1,reshape(RightTdecomp.EigenVectors/*.ExtractColumns(std::vector<MPXInt>({0}))*/,MPSD.RightMatrix.getOutwardMatrixIndex().size()))));
    //If failure, return dummy
    if (!VRDecomp.first.size()) return UnitCell(basis);

    MPX_matrix Y(contract(reorder(VRDecomp.second,1,reorder10,1),0,MPX_matrix(basis,VRDecomp.second.Index(0),VRDecomp.first),0,contract10));

    MPXDecomposition XLYDecomp(contract(contract(X,0,MPX_matrix(basis,MPSD.RightMatrix.getOutwardMatrixIndex(),PreviousLambda),0,contract10),0,Y,0,contract10).SVD());//SVD with no args keeps all singular values...
    SquareSumRescale(XLYDecomp.Values,1.0);//rescale here to set normalisation
    std::cout << "Entropy: " << entropy(XLYDecomp.Values) <<", Bond dimension: " << XLYDecomp.Values.size() << std::endl;

    UnitCell ans(basis);
    std::vector<double> scalevecX(MPSD.LeftMatrix.getInwardMatrixIndex().size(),1.0/sqrt(LeftTdecomp.Values.at(0).real()));

    ans.Matrices.emplace_back(reorder(contract(contract(XLYDecomp.ColumnMatrix,1,contract(MPX_matrix(basis,X.Index(0),scalevecX),0,X,0,contract10),0,contract00),0,MPSD.LeftMatrix,0,contract11),0,reorder102,2));
    ans.Matrices.emplace_back(contract(NL,0,contract(Xinv,0,XLYDecomp.ColumnMatrix,0,contract10),0,contract20));

    //ans.Matrices.emplace_back(contract(contract(contract(MPSD.RightMatrix,0,CurrentLambda,0,contractL),0,contract(Y,0,XLYDecomp.RowMatrix,1,contract11),0,contract10),0,MPX_matrix(basis,XLYDecomp.RowMatrix.Index(1),XLYDecomp.Values,1),0,contract20));

    ans.Lambdas.emplace_back(XLYDecomp.Values);
    ans.Lambdas.emplace_back(MPSD.Values);
    return ans;
  }

  UnitCell Orthogonalise(const UnitCell& C){
    std::cout << "Orthogonalising..." << std::endl;
    const EigenStateArray& spectrum(C.Matrices.at(0).GetPhysicalSpectrum());
    std::cout << "Forming left transfer matrix" << std::endl;
    MPX_matrix TransferMatrix(MakeLTransferMatrix(C,C));
    SparseED LeftTdecomp(TransferMatrix.LeftEigs(State(spectrum[0].getChargeRules()),1,LARGESTMAGNITUDE));
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
    SparseED RightTdecomp(TransferMatrix.RightEigs(State(spectrum[0].getChargeRules()),1,LARGESTMAGNITUDE));
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
    SparseED LeftTdecomp(TransferMatrixParts(C).LeftED(1,LARGESTMAGNITUDE));
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
    //std::cout << "XDagger X"<< std::endl;
    std::pair<std::vector<double>,MPX_matrix> VLDecomp(SqrtDR(MPX_matrix(spectrum,VLIndices,1,reshape(LeftTdecomp.EigenVectors,C.Matrices.at(0).Index(1).size()).massage())));
    MPX_matrix X(contract(MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first),0,VLDecomp.second,0,contract10));
    MPX_matrix Xinv(contract(reorder(VLDecomp.second,1,reorder10,1),0,MPX_matrix(spectrum,VLDecomp.second.Index(0),VLDecomp.first,1),0,contract10));
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
    return MakeLTransferMatrix(bra,ket).LeftEigs(1,LARGESTMAGNITUDE).Values.at(0);
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
    if (separation==0) return 0.0;

    MPX_matrix LAMBDA0(U.Matrices.at(0).basis(),U.Matrices.at(0).Index(1),U.Lambdas.at(0));
    if (separation==1){
      return TwoVertexMeasurement(W1,W2,U.Matrices.at(0),U.Matrices.at(1),LAMBDA0);
    }
    else {
      MPX_matrix accumulator;
      if (separation % 2) //odd > 1
	accumulator=contract(U.Matrices.at(1),1,contract(contract(W1,0,contract(U.Matrices.at(0),1,U.Matrices.at(0),0,contract11),0,contract0022).RemoveDummyIndices(std::vector<MPXInt>({{0,1}})),0,U.Matrices.at(1),0,contract11),0,contract0110);
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

  std::complex<double> iTwoVertexEnergy(const MPO_matrix& ColX,const MPO_matrix& RowX,const MPO_matrix& H1,const UnitCell& Ortho){
    if (Ortho.Matrices.size()!=2){std::cout << "Only two vertex basis accepted for now" <<std::endl; return 0.0;}
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

  TransferMatrixComponents::TransferMatrixComponents(const std::vector<const MPS_matrix*>& KetPtrs, bool HV, const State* Target) : TransferMatrixComponents(KetPtrs,KetPtrs,HV,Target) {}

  TransferMatrixComponents::TransferMatrixComponents(const std::vector<const MPS_matrix*>& BraPtrs, const std::vector<const MPS_matrix*>& KetPtrs, bool HV, const State* Target) : TargetStatePtr_(Target), CellSize_(BraPtrs.size()), Hermitian_answer_(HV) {
    //init ptr lists

    //range check
    if (BraPtrs.size()!=KetPtrs.size()){std::cout << "Bra and Ket fragments have different sizes! " << BraPtrs.size() << " " << KetPtrs.size() <<std::endl; exit(1);}
    else {
      for (std::pair<std::vector<const MPS_matrix*>::const_iterator,std::vector<const MPS_matrix*>::const_iterator> cits={BraPtrs.begin(),KetPtrs.begin()};cits.first!= BraPtrs.end() && cits.second!= KetPtrs.end(); ++cits.first, ++cits.second) {
	BraKetMatrixPtrs_.push_back({*(cits.first),*(cits.second)});
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
      if (TargetStatePtr_){ //have we set a target state here? If so, do we have a match, are we in the target state sector?
	if (m[c] != *TargetStatePtr_){continue;}
      }
      for (MPXInt r=0;r<mprimed.size();++r){
	if (m[c]==mprimed[r]){ //check for consistency of charges	      
	  allowed_indices_.push_back(r+vrows_*c);
	  allowed_indices_dagger_.push_back(c+vcols_*r);
	  rows_and_cols_.push_back(std::array<MPXInt,2> {{r,c}});
	}
      }
    }
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

    if (length()<500){
      //use simple method
      std::cout << "Simple method" <<std::endl;
      //accumulate
      MPX_matrix TM(accumulate_());
      return TM.LeftEigs(numevals,which,initial);
    }
    else {
      SparseED ans(length(),numevals);
      std::complex<double>* Evecs = new std::complex<double>[allowed_indices_.size()*numevals];
      std::complex<double>* Evals = new std::complex<double>[numevals];
      
      //No initial guess?
      //try equal weighted indices.      

      SparseVectorWithRestriction guess_struct(initial,&allowed_indices_);
      arpack::arpack_eigs<TransferMatrixComponents,SparseVectorWithRestriction> eigensystem(this,&LeftComponentsMultiply,allowed_indices_.size(),initial ? &guess_struct : NULL,&ConvertSparseVectorWithRestriction,numevals,which,Evals,Evecs);
      if (eigensystem.error_status()) {std::cout << "Error with tensor arpack" << std::endl;exit(1);}
      for (size_t v=0;v<static_cast<size_t>(numevals);++v){
	std::cout << "Eval: " << Evals[v] <<std::endl;
	ans.Values.push_back(Evals[v]);
	
	for (MPXInt i=0;i<allowed_indices_.size();++i){
	  //here is where to drop tiny values...
	  if (abs(Evecs[i+v*allowed_indices_.size()])>SPARSETOL){
	    ans.EigenVectors.entry(allowed_indices_[i],v,Evecs[i+v*allowed_indices_.size()]);
	  }
	  if (Hermitian_answer_) //wrecks normalisation, but makes reshape to Hermitian high precision
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

    if (length()<500){
      //use simple method
      std::cout << "Simple method" <<std::endl;
      //accumulate
      MPX_matrix TM(accumulate_());
      return TM.RightEigs(numevals,which,initial);
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


}

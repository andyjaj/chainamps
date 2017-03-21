class LM_TEBD : public TimeBase {
private:
  const std::string MPSName_;
  const EigenStateArray& Basis_;
  const uMPXInt NumVertices_;
  const MPX_matrix SingleVertexOp_; //for open boundary conditions
  const TrotterDecomposition m_EvolutionOperators;
  //MPX_matrix apply_and_decompose(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS); //returns an MPX to apply right
  //void left_canonise_measure(uMPXInt chi,double minS,std::vector<MultiVertexMeasurement>& measurements);
  
public:
  //LM_TEBD(const MPO_matrix& H, const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //load from file
  //LM_TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //use repeating simple MPS_matrix defined by InitialMPS_matrix

  LM_TEBD(const MPO_matrix& H, const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1) : TimeBase(time_step_size,results),MPSName_(MPSName),Basis_(H.basis()),NumVertices_(NumVertices),SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)) {

    //need to right canonize the state

    //store the state with a new name (right canonical parts)
    for (uMPXInt n=1;n<=NumVertices_/2;++n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Right_" << n << ".MPS_matrix";
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Right_" << n << ".MPS_matrix";
      MPS_matrix current(load_MPS_matrix(Initialnamestream.str(),Basis_));
      current.store(Evolvingnamestream.str());
    }
    //do something with lambda
    std::stringstream Lambdanamestream;
    Lambdanamestream << MPSName_ << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
    MPX_matrix U(load_MPX_matrix(Lambdanamestream.str(),Basis_));
    //left canonize the right hand parts
    for (uMPXInt n=NumVertices_/2;n>0;--n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Left_" << n << ".MPS_matrix";
      //contract U onto next left
      //and reshape to right form
      //do SVD
      MPXDecomposition decomp(MPS_matrix(contract(load_MPS_matrix(Initialnamestream.str(),Basis_),0,U,0,contract20)).right_shape().SVD());
      
      //record row vectors (U part) as new U
      U=std::move(contract(decomp.ColumnMatrix,0,MPX_matrix(Basis_,decomp.ColumnMatrix.Index(1),decomp.Values),0,contract10));
      
      //store fully right canonised initial state for calculating overlap...
      std::stringstream Rightnamestream;
      Rightnamestream << MPSName_ << "_Right_" << NumVertices_-n+1 << ".MPS_matrix";
      decomp.RowMatrix.store(Rightnamestream.str());
      //store updated 'current' vertex matrix
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_-n+1 << ".MPS_matrix";
      decomp.RowMatrix.store(Evolvingnamestream.str());
    }
    //at end Vd should be really trivial, and give sqrt(normalisation) factor
    std::cout << "Norm: " << U.Trace() <<std::endl;
  }

  LM_TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1) : TimeBase(time_step_size,results),MPSName_(MPSName),Basis_(H.basis()),NumVertices_(NumVertices),SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)) {
    if (!InitialMPS_matrix.Index(1).Physical() || !InitialMPS_matrix.Index(0).Ingoing() || !InitialMPS_matrix.Index(1).Ingoing() || !InitialMPS_matrix.Index(2).Outgoing()){
      std::cout << "Initial repeating product state matrix needs to be right shaped" << std::endl; exit(1);
    }
    else {
      for (uMPXInt n=1;n<=NumVertices_;++n){
	std::stringstream Initialnamestream;
	Initialnamestream << MPSName_ << "_Right_" << n << ".MPS_matrix";
	std::stringstream Evolvingnamestream;
	Evolvingnamestream << "Evolving_" << MPSName_ << "_Right_" << n << ".MPS_matrix";
	
	InitialMPS_matrix.store(Initialnamestream.str()); //store a copy of initial state for later
	InitialMPS_matrix.store(Evolvingnamestream.str()); //store a copy for evolving
      }
    }
  }

  void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);
  //void left_info();
  //void right_info();
  
};

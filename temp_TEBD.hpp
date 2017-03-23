class LM_TEBD : public TimeBase {
private:
  const std::string MPSName_;
  const EigenStateArray& Basis_;
  const uMPXInt NumVertices_;
  const MPX_matrix SingleVertexOp_; //for open boundary conditions
  const TrotterDecomposition m_EvolutionOperators;
  //void apply_to_odd_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  //void apply_to_even_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS);
  //void left_canonise(uMPXInt chi=0,double minS=0);
  //void left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements,uMPXInt chi=0,double minS=0);


  void apply_to_odd_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS){
    //Start state should be left canonical
    //At end state is right canonical but truncated.
    //can then use left canonize measure if necessary

    //load each and apply
    std::cout << "Odd bonds" << std::endl;
    //need a starter MPS
    std::stringstream StartNameStream;
    StartNameStream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_ << ".MPS_matrix";
    std::string RightName=StartNameStream.str();
    MPS_matrix R(load_MPS_matrix(RightName,Basis_));

    for (uMPXInt v=NumVertices_;v>1;v-=2){
      std::stringstream LNameStream;
      LNameStream << "Evolving_" << MPSName_ << "_Left_" << v-1 << ".MPS_matrix";
      MPSDecomposition decomp(reorder(contract(BondOp,0,contract(load_MPS_matrix(LNameStream.str(),Basis_),0,R,0,contract21),0,contract2032),0,reorder0213,2).SVD(bond_dimension,minS));
      
      //overwrite with new right part of pair
      {
	std::stringstream StoreNameStream;
	StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+1-v << ".MPS_matrix";
	decomp.RightMatrix.store(StoreNameStream.str());
      }
      //need to 'rotate left part'
      {
	std::stringstream StoreNameStream;
	StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+2-v << ".MPS_matrix";
	MPXDecomposition rot(MPS_matrix(contract(decomp.LeftMatrix,0,MPX_matrix(Basis_,decomp.LeftMatrix.Index(2),decomp.Values),0,contract20)).right_shape().SVD());
	rot.RowMatrix.store(StoreNameStream.str()); //now should be right canonical
	if (v>2){
	  std::stringstream NewNameStream;
	  NewNameStream << "Evolving_" << MPSName_ << "_Left_" << v-2 << ".MPS_matrix";
	  RightName=NewNameStream.str();
	  R=std::move(MPS_matrix(contract(load_MPS_matrix(RightName,Basis_),0,contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10),0,contract20)));
	}
	else {
	  //check norm
	  std::cout << "Norm at end of odd bonds, right canonise: " << contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10).Trace() <<std::endl;
	}
      }
    }
  }

  void apply_to_even_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS){
    //Start state should be left canonical
    //At end state is right canonical but truncated.
    //can then use left canonize measure if necessary

    std::cout << "Even bonds" << std::endl;
    //First need to update the final matrix
    std::stringstream SpecialNameStream;
    SpecialNameStream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_ << ".MPS_matrix";
    //MPS_matrix(std::move(contract(SingleVertexOp_,0,load_MPS_matrix(SpecialNameStream.str(),Basis_),0,contract10).CombineSimilarMatrixIndices())).store(SpecialNameStream.str());
    
    std::stringstream StoreNameStream;
    StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << 1 << ".MPS_matrix";
    
    MPXDecomposition rot(MPS_matrix(std::move(contract(SingleVertexOp_,0,load_MPS_matrix(SpecialNameStream.str(),Basis_),0,contract10).CombineSimilarMatrixIndices())).right_shape().SVD());
    rot.RowMatrix.store(StoreNameStream.str()); //now should be right canonical
      
    std::stringstream StartNameStream;
    StartNameStream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_-1 << ".MPS_matrix";
    std::string RightName=StartNameStream.str();
    MPS_matrix R(contract(load_MPS_matrix(RightName,Basis_),0,contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10),0,contract20));      
    
    if (NumVertices_>2){

      for (uMPXInt v=NumVertices_-1;v>1;v-=2){ //need to start at right hand side off last even bond
	std::stringstream LNameStream;
	LNameStream << "Evolving_" << MPSName_ << "_Left_" << v-1 << ".MPS_matrix";
	MPSDecomposition decomp(reorder(contract(BondOp,0,contract(load_MPS_matrix(LNameStream.str(),Basis_),0,R,0,contract21),0,contract2032),0,reorder0213,2).SVD(bond_dimension,minS));
	
	//overwrite with new right part of pair
	{
	  std::stringstream StoreNameStream;
	  StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+1-v << ".MPS_matrix";
	  decomp.RightMatrix.store(StoreNameStream.str());
	}
      //need to 'rotate left part'
      

	std::stringstream StoreNameStream;
	StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+1-v+1 << ".MPS_matrix";
	MPXDecomposition rot(MPS_matrix(contract(decomp.LeftMatrix,0,MPX_matrix(Basis_,decomp.LeftMatrix.Index(2),decomp.Values),0,contract20)).right_shape().SVD());
	rot.RowMatrix.store(StoreNameStream.str()); //now should be right canonical
	if (v>2){ //this should be true as long as N/2 is even
	  std::stringstream NewNameStream;
	  NewNameStream << "Evolving_" << MPSName_ << "_Left_" << v-2 << ".MPS_matrix";
	  RightName=NewNameStream.str();
	  R=std::move(MPS_matrix(contract(load_MPS_matrix(RightName,Basis_),0,contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10),0,contract20)));
	}
      }
    }
    else {
      std::cout <<"Only two vertices, no even bonds" <<std::endl;
    }
    //need to canonize first vertex still
    
    MPXDecomposition decomp(R.right_shape().SVD());
    
    std::stringstream FirstNameStream;
    FirstNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
    decomp.RowMatrix.store(FirstNameStream.str());
    
    std::cout << "Norm at end of even bonds, right canonise: " << contract(decomp.ColumnMatrix,0,MPX_matrix(Basis_,decomp.ColumnMatrix.Index(1),decomp.Values),0,contract10).Trace() <<std::endl;
  }
  void left_canonise(uMPXInt chi=0,double minS=0){
    std::cout << "Begin left canonisation" << std::endl; 
    MPX_matrix Vd(Basis_);
    {
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Basis_).left_shape().SVD(chi,minS));
      RightEnddecomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << RightEnddecomp.Values.size() << std::endl;
      Vd=std::move(contract(MPX_matrix(Basis_,RightEnddecomp.RowMatrix.Index(0),RightEnddecomp.Values),0,RightEnddecomp.RowMatrix,0,contract10));
      std::stringstream LeftStartstream;
      LeftStartstream << "Evolving_" << MPSName_ << "_Left_1" << ".MPS_matrix";
      RightEnddecomp.ColumnMatrix.store(LeftStartstream.str());
    }
    for (uMPXInt n=NumVertices_-1;n>0;--n){
      std::stringstream RightNamestream;
      RightNamestream << "Evolving_" << MPSName_ << "_Right_" << n << ".MPS_matrix";
      //contract U onto current
      //and reshape to left form
      //do SVD
      const uMPXInt v=NumVertices_-n+1;
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Basis_),0,contract10)).left_shape().SVD(chi,minS));
      decomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
      
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      //store updated 'current' vertex matrix
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      decomp.ColumnMatrix.store(LeftNamestream.str());
    }
    //at end Vd should be trivial
    std::cout << "Norm at end of left canonisation: " << abs(Vd.Trace()) << std::endl; 
  }


  void left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements,uMPXInt chi=0,double minS=0){
    std::cout << "Begin left canonisation" << std::endl; 

    //it can be useful to measure as we step through the system
    std::vector<double> real_results(1,current_time());
    //std::vector<MPX_matrix> measure_tensors(measurements.size(),MPX_matrix(Basis_));

    //we'd also like to measure the overlap with the initial state...
    MPX_matrix overlap_matrix(Basis_);
    MPX_matrix Vd(Basis_);
    {
      //treats first vertex as a special case
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Basis_).left_shape().SVD(chi,minS));
      RightEnddecomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << RightEnddecomp.Values.size() << std::endl;
      if (NumVertices_==2){ //special case for two vertices only
	//measure truncation and entropy
	real_results.push_back(RightEnddecomp.Truncation);
	real_results.push_back(entropy(RightEnddecomp.Values));
      }
      for (auto&& m : measurements){
	if (m.finish()>NumVertices_){
	  std::cout << "MultiVertexMeasurements incorrectly defined, final measurement position is > number of vertices: " << m.finish() << " " << NumVertices_ << std::endl;
	  exit(1);
	}
	m.update(1,RightEnddecomp); //first vertex is position 1!
      }
      {
	const MPX_matrix& current=RightEnddecomp.ColumnMatrix;
	std::stringstream Initialstream;
	Initialstream << MPSName_ << "_Left_" << 1 << ".MPS_matrix";
	MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Basis_));
	std::vector<MPXIndex> indices;
	indices.emplace_back(1,initial.Index(2));
	indices.emplace_back(0,current.Index(2));
	overlap_matrix=std::move(MPX_matrix(Basis_,indices,1,contract_to_sparse(initial,1,current,0,contract00)));
      }
      Vd=std::move(contract(MPX_matrix(Basis_,RightEnddecomp.RowMatrix.Index(0),RightEnddecomp.Values),0,RightEnddecomp.RowMatrix,0,contract10));
      std::stringstream LeftStartstream;
      LeftStartstream << "Evolving_" << MPSName_ << "_Left_1" << ".MPS_matrix";
      RightEnddecomp.ColumnMatrix.store(LeftStartstream.str());

    }
    for (uMPXInt n=NumVertices_-1;n>0;--n){
      std::stringstream RightNamestream;
      RightNamestream << "Evolving_" << MPSName_ << "_Right_" << n << ".MPS_matrix";
      //contract U onto current
      //and reshape to left form
      //do SVD
      const uMPXInt v=NumVertices_-n+1;
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Basis_),0,contract10)).left_shape().SVD(chi,minS));
      decomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
      if (v==NumVertices_/2){ //if we only have two chains, this is never obeyed...
	//measure truncation and entropy
	real_results.push_back(decomp.Truncation);
	real_results.push_back(entropy(decomp.Values));
      }
      for (auto&& m : measurements){
	m.update(v,decomp);
      }
      {
	const MPX_matrix& current=decomp.ColumnMatrix;
	std::stringstream Initialstream;
	Initialstream << MPSName_ << "_Left_" << v << ".MPS_matrix";
	std::cout << Initialstream.str() << std::endl;
	MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Basis_));
	overlap_matrix=std::move(contract(initial,1,contract(overlap_matrix,0,current,0,contract11),0,contract0110));
      }
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      //store updated 'current' vertex matrix
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      decomp.ColumnMatrix.store(LeftNamestream.str());
    }

    std::complex<double> overlap=overlap_matrix.Trace();
    std::vector<std::complex<double> > complex_results;
    real_results.push_back(abs(overlap));
    complex_results.push_back(overlap);
    for (auto&& m : measurements){
      complex_results.push_back(m.result());
    }
    m_results.push(Data(real_results,complex_results));
    //at end Vd should be trivial

    std::cout << "Norm at end of left canonisation: " << abs(Vd.Trace()) << std::endl; 
  }


public:
  //LM_TEBD(const MPO_matrix& H, const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //load from file
  //LM_TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1); //use repeating simple MPS_matrix defined by InitialMPS_matrix
  //void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1);

  LM_TEBD(const MPO_matrix& H, const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1) : TimeBase(time_step_size,results),MPSName_(MPSName),Basis_(H.basis()),NumVertices_(NumVertices),SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)) {

    //canonize and store
    for (uMPXInt n=1;n<=NumVertices_/2;++n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Left_" << n << ".MPS_matrix";
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << n << ".MPS_matrix";
      MPS_matrix current(load_MPS_matrix(Initialnamestream.str(),Basis_));
      current.store(Evolvingnamestream.str());
    }
    //do something with lambda
    std::stringstream Lambdanamestream;
    Lambdanamestream << MPSName_ << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
    MPX_matrix Vd(load_MPX_matrix(Lambdanamestream.str(),Basis_));
    //left canonize the right hand parts
    for (uMPXInt n=NumVertices_/2;n>0;--n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Right_" << n << ".MPS_matrix";
      //contract Vd onto current
      //and reshape to left form
      //do SVD
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(Initialnamestream.str(),Basis_),0,contract10)).left_shape().SVD());
     
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));

      //store for calculating overlap...
      std::stringstream Leftnamestream;
      Leftnamestream << MPSName_ << "_Left_" << NumVertices_-n+1 << ".MPS_matrix";
      decomp.ColumnMatrix.store(Leftnamestream.str());
      //store updated 'current' vertex matrix
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_-n+1 << ".MPS_matrix";
      decomp.ColumnMatrix.store(Evolvingnamestream.str());
    }
    //at end Vd should be really trivial, and give sqrt(normalisation) factor
    std::cout << "Norm: " << Vd.Trace() <<std::endl;

  }

  LM_TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order=1) : TimeBase(time_step_size,results),MPSName_(MPSName),Basis_(H.basis()),NumVertices_(NumVertices),SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)) {
    
    if (!InitialMPS_matrix.Index(0).Physical() || !InitialMPS_matrix.Index(0).Ingoing() || !InitialMPS_matrix.Index(1).Ingoing() || !InitialMPS_matrix.Index(2).Outgoing()){
      std::cout << "Initial matrix needs to be left shaped" << std::endl; exit(1);
    }
    else {
      for (uMPXInt n=1;n<=NumVertices_;++n){
	std::stringstream Initialnamestream;
	Initialnamestream << MPSName_ << "_Left_" << n << ".MPS_matrix";
	std::stringstream Evolvingnamestream;
	Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << n << ".MPS_matrix";
	
	InitialMPS_matrix.store(Initialnamestream.str()); //store a copy of initial state for later
	InitialMPS_matrix.store(Evolvingnamestream.str()); //store a copy for evolving
      }
    }

  }

  void evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension=0, double minS=0.0, uMPXInt measurement_interval=1){
    //do the evolution
    if (m_EvolutionOperators.order()==1){
      std::cout <<"1st order time step evolution" <<std::endl;
      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	//this is the first half of the time step....
	apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	left_canonise();
	apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  left_canonise_measure(measurements);
	}
	else {
	  left_canonise();
	}
      }
    }
    else if (m_EvolutionOperators.order()==2){
      std::cout <<"2nd order time step evolution" <<std::endl;
      //bond order grows rapidly, so need to compress after each application of a set of bond operators
      //first step
      std::cout << "Odd bonds" <<std::endl;
      apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS); //half step
      left_canonise();
      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	left_canonise();
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[2]),bond_dimension,minS);
	  std::cout << "Left canonize measure" <<std::endl;
	  left_canonise_measure(measurements);//measurement
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	  std::cout << "Left canonize" <<std::endl;
	  left_canonise();
	}
	else {
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	  std::cout << "Left canonize" <<std::endl;
	  left_canonise();
	}
      }
    }
    else {
      std::cout << "Haven't implemented higher orders yet in TEBD::evolve()" << std::endl; exit(1);
    }
  }
  //void left_info();
  //void right_info();
  
};

#include <vector>
#include <complex>
#include <iostream>
#include <sstream>

#include "states.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"
#include "TEBD_routines.hpp"

namespace ajaj{

  TrotterDecomposition::TrotterDecomposition(const MPO_matrix& H,double time_step_size,uMPXInt order,const State* blockstate_ptr) : m_H_ptr(&H), m_time_step_size(time_step_size), m_order(order){

    double tau=m_time_step_size;
    if (m_order==0){
      //nothing
    }
    else {
      //make bond operator
      if (m_order>4 || m_order==3){std::cout <<"Only 1st, 2nd and 4th order decompositions are supported. Aborting..." << std::endl; exit(1);}
      //make the Hamiltonian for a single bond
    
      //now make the trotter operators

      MPX_matrix OddBondH(MakeOddBondHamiltonian(*m_H_ptr));
      MPX_matrix EvenBondH(MakeEvenBondHamiltonian(*m_H_ptr));
      
      if (m_order==1){
	
	BondOperators.emplace_back(MakeBondEvolutionOperator(OddBondH,tau,blockstate_ptr));
	BondOperators.emplace_back(MakeBondEvolutionOperator(EvenBondH,tau,blockstate_ptr));

	std::cout << "Done" <<std::endl;
	//safeish to have pointer to vector element if the vector no longer grows
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1));

      }
      else if (m_order==2){ //second order is a special case with a nice simplification for timesteps where no measurement is made
	
	BondOperators.emplace_back(MakeBondEvolutionOperator(OddBondH,0.5*tau));
	BondOperators.emplace_back(MakeBondEvolutionOperator(EvenBondH,tau));
	//safeish to have pointer to vector element if the vector no longer grows
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1));
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
      }

      else if (m_order==4){

	//MPX_matrix BondH(MakeBondHamiltonian(*m_H_ptr));
	
	double f1=0.414490771794375737142354062860761495711774604016707133323;
	double f3=-0.657963087177502948569416251443045982847098416066828533293;
	double f1_3=-0.24347231538312721142706218858228448713532381205012139996;

	BondOperators.emplace_back(MakeBondEvolutionOperator(OddBondH,0.5*tau*f1)); //odd half t1 step
	BondOperators.emplace_back(MakeBondEvolutionOperator(EvenBondH,tau*f1)); //even full t1 step
	BondOperators.emplace_back(MakeBondEvolutionOperator(OddBondH,tau*f1)); //odd full t1 step
	BondOperators.emplace_back(MakeBondEvolutionOperator(OddBondH,0.5*tau*f1_3)); //odd t1+t3 half step
	BondOperators.emplace_back(MakeBondEvolutionOperator(EvenBondH,tau*f3)); //even t3 step

	//11 pieces per full step
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(0)); //Odd
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1)); //Even
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(2)); //Odd
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1)); //E
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(3)); //O
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(4)); //E
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(3)); //O
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1)); //E
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(2)); //O
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(1)); //E
	OrderedOperatorPtrs.emplace_back(&BondOperators.at(0)); //Odd
      }
      else {
	std::cout <<"Unsupported Trotter Order!" <<std::endl;
	exit(1);
      }
    }
  }

  iTEBD::iTEBD(const MPO_matrix& H,const UnitCell& C, double time_step_size, DataOutput& results, const std::string& Name, uMPXInt order) : TimeBase(time_step_size,results),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)),m_initial_unit(C),m_unit(C),Name_(Name) {
    std::ofstream DensityFileStream_;
    //std::stringstream dnamestream;
    //dnamestream << "iTEBD_One_Vertex_Densities.dat";
    DensityFileStream_.open("iTEBD_One_Vertex_Densities.dat",ios::out | ios::trunc);
    for (auto&& i : H.getPhysicalSpectrum().Energies()){
      DensityFileStream_ << i << " ";
    }
    DensityFileStream_ << std::endl;
    DensityFileStream_.close();
  }

  void iTEBD::do_measurements(const UnitCell& ortho, const std::vector<MPO_matrix>& measuredMPOs){
    std::cout << "Measuring...." << std::endl;
    std::vector<complex<double> > complex_results(1,Overlap(m_initial_unit,ortho));
    std::vector<double> real_results(1,current_time());
    real_results.push_back(m_truncation);
    real_results.push_back(entropy(ortho.Lambdas.at(0)));
    real_results.push_back(abs(complex_results.at(0)));
    for (std::vector<MPO_matrix>::const_iterator cit=measuredMPOs.begin();cit!=measuredMPOs.end();++cit){
      complex_results.emplace_back(OneVertexMeasurement(*cit,ortho));
    }
    m_results.push(m_current_time_step,Data(real_results,complex_results));
    std::ofstream DensityFileStream_;
    DensityFileStream_.open("iTEBD_One_Vertex_Densities.dat",ios::out | ios::app);
    ortho.OutputPhysicalIndexDensities(DensityFileStream_);
    DensityFileStream_.close();
  }

  const UnitCell& iTEBD::evolve(uMPXInt num_steps, const std::vector<MPO_matrix>& measuredMPOs, uMPXInt bond_dimension, double minS, uMPXInt measurement_interval){
    
    if (order()==1){
      std::cout <<"1st order time step evolution" <<std::endl;
      for (uMPXInt n=0;n<num_steps;++n){
	//++m_current_time_step;
	update_time();
	std::cout << "Time " << current_time() << std::endl;
	//this is the first half of the time step....
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	//this is the second
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	//do we need to take a measurement?
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  UnitCell ortho(OrthogonaliseInversionSymmetric(m_unit));
	  //ortho.OutputOneVertexDensityMatrix("OneVertexRho",m_current_time_step);
	  if (ortho.size()){ //only if unitcell isn't empty
	    ortho.store(Name_,m_current_time_step);
	    this->do_measurements(ortho,measuredMPOs);
	  }
	  m_unit=std::move(ortho);
	}
      }
    }
    else if (order()==2){
      std::cout <<"2nd order time step evolution" <<std::endl;
      //2nd order special start
      apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
      for (uMPXInt n=0;n<num_steps;++n){
	//++m_current_time_step;
	update_time();
	std::cout << "Time " << current_time() << std::endl;
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[2]),bond_dimension,minS);
	  //swap order
	  m_unit.swap(0,1);
	  //measure etc.
	  UnitCell ortho(OrthogonaliseInversionSymmetric(m_unit));
	  if (ortho.size()){ //if generating unitcell works, then measure and use it
	    ortho.store(Name_,m_current_time_step);
	    this->do_measurements(ortho,measuredMPOs);
	    m_unit=std::move(ortho);
	  }

	  m_unit.OutputOneVertexDensityMatrix("OneVertexRho",m_current_time_step);
	  //complete time step
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	}
	else {
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	}
       
      }
    }
    else if (order()==4){
      std::cout <<"4th order time step evolution" <<std::endl;
      std::cout <<"Canonization is performed at every full time step" <<std::endl;
      //special start
      apply_and_decompose(*m_EvolutionOperators.OrderedOperatorPtrs[0],bond_dimension,minS);
      //
      for (uMPXInt n=0;n<num_steps;++n){
	update_time();
	std::cout << "Time " << current_time() << std::endl;
	for (auto i=1; i<10;++i){
	  apply_and_decompose(*m_EvolutionOperators.OrderedOperatorPtrs[i],bond_dimension,minS);
	  if (i%2) m_unit=std::move(OrthogonaliseInversionSymmetric(m_unit));
	}
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  apply_and_decompose(*m_EvolutionOperators.OrderedOperatorPtrs[10],bond_dimension,minS);
	  m_unit.swap(0,1);//need a swap here because an odd operator follows an odd
	  m_unit=std::move(OrthogonaliseInversionSymmetric(m_unit));
	  if (m_unit.size()){ //only if unitcell isn't empty
	    m_unit.store(Name_,m_current_time_step);
	    this->do_measurements(m_unit,measuredMPOs);
	    m_unit.OutputOneVertexDensityMatrix("OneVertexRho",m_current_time_step);
	  }
	  else {
	    std::cout<<"Orthogonalisation Error" <<std::endl; exit(1);
	  }
	  apply_and_decompose(*m_EvolutionOperators.OrderedOperatorPtrs[0],bond_dimension,minS);
	}
	else {
	  apply_and_decompose(*m_EvolutionOperators.OrderedOperatorPtrs[2],bond_dimension,minS);
	}
	
      }
    }
    else {
      std::cout << "Haven't implemented this order yet in iTEBD::evolve()" << std::endl; exit(1);
    }
    return m_unit;
  }

  const UnitCell& iTEBD::apply_and_decompose(const MPX_matrix& BondOp,uMPXInt bond_dimension, double minS){
    //also swaps order of vertices in cell
    //std::vector<double> lambda0(m_unit.Lambdas[0]);
    //one also needs to retain the object Hastings refers to as 'C'
#ifndef NDEBUG
    std::cout << "Forming C" << std::endl;
#endif
    MPX_matrix C(reorder(contract(BondOp,0,contract(m_unit.Matrices.at(0),0,m_unit.Matrices.at(1),0,contract21),0,contract2032),0,reorder0213,2));
#ifndef NDEBUG
    std::cout << "Apply C to unit cell and decompose" << std::endl;
#endif
    //first populate with parts of decomp that we need to keep
    MPSDecomposition decomp(contract(C,0,MPX_matrix(C.GetPhysicalSpectrum(),C.Index(3),m_unit.Lambdas.at(0)),0,contract30).ShiftNumRowIndices(2).SVD(bond_dimension,minS));
    //we can shamelessly steal from the above
    //start by overwriting parts of the old unitcell
#ifndef NDEBUG
    std::cout << "Updating unit cell" << std::endl;
#endif
    m_unit.Matrices[0]=std::move(decomp.LeftMatrix);
    m_unit.Matrices[1]=std::move(reorder(contract(m_unit.Matrices.at(0),1,C,0,contract0011),0,reorder102,2).Rescale(1.0/sqrt(SquareSum(m_unit.Lambdas[0]))));
    m_unit.Lambdas[1]=std::move(decomp.Values);
#ifndef NDEBUG
    std::cout << "Calculate truncation" << std::endl;
#endif
    m_truncation=decomp.Truncation;//1.0-SquareSumRescale(m_unit.Lambdas.at(1),1.0);
#ifndef NDEBUG
    std::cout << "Swap sites" << std::endl;
#endif
    m_unit.swap(0,1);
#ifndef NDEBUG
    std::cout << "Return new unit cell" << std::endl;
#endif
    return m_unit;
  }

  void iTEBD::change_bond_operator(const MPO_matrix& H, double time_step_size){
    m_time_step_size=time_step_size;
    m_EvolutionOperators=TrotterDecomposition(H,time_step_size,m_EvolutionOperators.order());
  }
  
  void TEBD::apply_to_odd_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS){
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
      if (decomp.Truncation>max_truncation_) max_truncation_=decomp.Truncation;
      //overwrite with new right part of pair
      
      std::stringstream RStoreNameStream;
      RStoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+1-v << ".MPS_matrix";
      decomp.RightMatrix.store(RStoreNameStream.str());
      
      //need to 'rotate left part'
      
      std::stringstream LStoreNameStream;
      LStoreNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_+2-v << ".MPS_matrix";
      MPXDecomposition rot(MPS_matrix(contract(decomp.LeftMatrix,0,MPX_matrix(Basis_,decomp.LeftMatrix.Index(2),decomp.Values),0,contract20)).right_shape().SVD());

      if (rot.Truncation>max_truncation_) max_truncation_=rot.Truncation;

      if (v==2){
	std::complex<double> phase=rot.ColumnMatrix.Trace();
	std::cout << "Phase at end of odd bonds is " << phase << std::endl;	    
	rot.RowMatrix.Rescale(phase);
      }
      rot.RowMatrix.store(LStoreNameStream.str()); //now should be right canonical
      if (v>2){ //get ready for next iteration of loop
	std::stringstream NewNameStream;
	NewNameStream << "Evolving_" << MPSName_ << "_Left_" << v-2 << ".MPS_matrix";
	RightName=NewNameStream.str();
	R=std::move(MPS_matrix(contract(load_MPS_matrix(RightName,Basis_),0,contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10),0,contract20)));
      }
    }
  }

  void TEBD::apply_to_even_bonds(const MPX_matrix& BondOp,uMPXInt  bond_dimension, double minS){
    //Start state should be left canonical
    //At end state is right canonical but truncated.
    //can then use left canonize measure if necessary

    std::cout << "Even bonds" << std::endl;
    //First need to update the final matrix
    std::stringstream SpecialNameStream;
    SpecialNameStream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_ << ".MPS_matrix";
    
    std::stringstream StoreNameStream;
    StoreNameStream << "Evolving_" << MPSName_ << "_Right_" << 1 << ".MPS_matrix";

    //apply end operator and rotate
    //when applying single vertex op, shouldn't truncate at all.
    MPXDecomposition rot(load_MPS_matrix(SpecialNameStream.str(),Basis_).right_shape().SVD());
    //MPXDecomposition rot(MPS_matrix(std::move(contract(SingleVertexOp_,0,load_MPS_matrix(SpecialNameStream.str(),Basis_),0,contract20).CombineSimilarMatrixIndices())).right_shape().SVD());
    rot.RowMatrix.store(StoreNameStream.str()); //now should be right canonical
    
    std::stringstream StartNameStream;
    StartNameStream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_-1 << ".MPS_matrix";
    std::string RightName=StartNameStream.str();
    MPS_matrix R(contract(load_MPS_matrix(RightName,Basis_),0,contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10),0,contract20));      
    
    if (NumVertices_>2){

      for (uMPXInt v=NumVertices_-1;v>1;v-=2){ //need to start at right hand side of last even bond
	std::stringstream LNameStream;
	LNameStream << "Evolving_" << MPSName_ << "_Left_" << v-1 << ".MPS_matrix";
	MPSDecomposition decomp(reorder(contract(BondOp,0,contract(load_MPS_matrix(LNameStream.str(),Basis_),0,R,0,contract21),0,contract2032),0,reorder0213,2).SVD(bond_dimension,minS));
	if (decomp.Truncation>max_truncation_) max_truncation_=decomp.Truncation;

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
	if (rot.Truncation>max_truncation_) max_truncation_=rot.Truncation;

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
    //need to canonize first (leftmost) vertex still
    
    MPXDecomposition decomp(R.right_shape().SVD());
    //if (decomp.Truncation>max_truncation_) max_truncation_=decomp.Truncation;

    std::complex<double> phase=decomp.ColumnMatrix.Trace(); //want to keep the phase, singluar value should be 1, but if <1 is indicative of a loss of norm.
    std::cout << "Phase at end of even bonds is " << phase  <<std::endl;

    decomp.RowMatrix.Rescale(phase);
    
    std::stringstream FirstNameStream;
    FirstNameStream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
    
    decomp.RowMatrix.store(FirstNameStream.str());
   
  }
  
  void TEBD::left_canonise(uMPXInt chi,double minS){
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

      //if final vertex, scale phase into matrix
      if (n==1){
	std::complex<double> phase=Vd.Trace(); //want to keep the phase, singluar value should be 1, but if <1 is indicative of a loss of norm.
	std::cout << "Phase at end of left canonisation: " << phase  <<std::endl;
	decomp.ColumnMatrix.Rescale(phase);
      }
      decomp.ColumnMatrix.store(LeftNamestream.str());
    }
  }


  void TEBD::left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements,uMPXInt chi,double minS,bool overlap_requested){
    
    std::cout << "Begin left canonisation" << std::endl; 

    //it can be useful to measure as we step through the system
    std::vector<double> real_results(1,current_time());

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
	//record truncation and entropy
	real_results.push_back(max_truncation_);
	real_results.push_back(entropy(RightEnddecomp.Values));
      }
      for (auto&& m : measurements){
	if (m.finish()>NumVertices_){
	  std::cout << "MultiVertexMeasurements incorrectly defined, final measurement position is > number of vertices: " << m.finish() << " " << NumVertices_ << std::endl;
	  exit(1);
	}
	m.update(1,RightEnddecomp); //first vertex is position 1!
      }
      if (overlap_requested) {
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
      if (SaveAll_){
	std::stringstream SaveAllStream;
	SaveAllStream << SAVEALLNAME << "_" << m_current_time_step << "_" << MPSName_ << "_Left_1" << ".MPS_matrix";
	RightEnddecomp.ColumnMatrix.store(SaveAllStream.str());
      }

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
	real_results.push_back(max_truncation_);
	real_results.push_back(entropy(decomp.Values));
      }

      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      //store updated 'current' vertex matrix
      //if final vertex, scale phase (from Vd) into matrix
      if (n==1){
	std::complex<double> phase=Vd.Trace(); //want to keep the phase, singluar value should be 1, but if <1 is indicative of a loss of norm.
	std::cout << "Phase at end of left canonisation: " << phase  <<std::endl;
	decomp.ColumnMatrix.Rescale(phase);
      }
      
      for (auto&& m : measurements){
	m.update(v,decomp);
      }
      if (overlap_requested) {
	const MPX_matrix& current=decomp.ColumnMatrix;
	std::stringstream Initialstream;
	Initialstream << MPSName_ << "_Left_" << v << ".MPS_matrix";
	std::cout << Initialstream.str() << std::endl;
	MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Basis_));
	overlap_matrix=std::move(contract(initial,1,contract(overlap_matrix,0,current,0,contract11),0,contract0110));
      }
     
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      decomp.ColumnMatrix.store(LeftNamestream.str());
      if (SaveAll_){
	std::stringstream SaveAllStream;
	SaveAllStream << SAVEALLNAME <<"_" << m_current_time_step << "_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
	decomp.ColumnMatrix.store(SaveAllStream.str());
      }
    }

    std::vector<std::complex<double> > complex_results;

    if (overlap_requested) {
      std::complex<double> overlap=overlap_matrix.Trace();
      real_results.push_back(abs(overlap));
      complex_results.push_back(overlap);
    }
    
    for (auto&& m : measurements){
      complex_results.push_back(m.result());
    }
    m_results.push(m_current_time_step,Data(real_results,complex_results));
  }

  void TEBD::left_canonise_measure_special(std::vector<MultiVertexMeasurement>& measurements, uMPXInt Index){
    std::cout << "Begin left canonisation" << std::endl; 

    //it can be useful to measure as we step through the system
    std::vector<double> real_results;

    MPX_matrix Vd(Basis_);
    {
      //treats first vertex as a special case
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Basis_).left_shape().SVD());
      RightEnddecomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << RightEnddecomp.Values.size() << std::endl;
      /*if (NumVertices_==2){ //special case for two vertices only
	real_results.push_back(entropy(RightEnddecomp.Values));
	}*/
      real_results.push_back(entropy(RightEnddecomp.Values));
      for (auto&& m : measurements){
	if (m.finish()>NumVertices_){
	  std::cout << "MultiVertexMeasurements incorrectly defined, final measurement position is > number of vertices: " << m.finish() << " " << NumVertices_ << std::endl;
	  exit(1);
	}
	m.update(1,RightEnddecomp); //first vertex is position 1!
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
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Basis_),0,contract10)).left_shape().SVD());
      decomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
      
      if (v<NumVertices_) real_results.push_back(entropy(decomp.Values));

      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      //store updated 'current' vertex matrix
      
      //if final vertex, scale phase (from Vd) into matrix
      if (n==1){
	std::complex<double> phase=Vd.Trace(); //want to keep the phase, singluar value should be 1, but if <1 is indicative of a loss of norm.
	std::cout << "Phase at end of left canonisation: " << phase  <<std::endl;
	decomp.ColumnMatrix.Rescale(phase);
      }
      
      for (auto&& m : measurements){
	m.update(v,decomp);
      }
      
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      decomp.ColumnMatrix.store(LeftNamestream.str());
    }

    std::vector<std::complex<double> > complex_results;
    
    for (auto&& m : measurements){
      complex_results.push_back(m.result());
    }
    m_results.push(Index,Data(real_results,complex_results));
  }
  
  TEBD::TEBD(const MPO_matrix& H, FiniteMPS& F, DataOutput& results) : TimeBase(0.0,results),F_(F),MPSName_(F.name()),Basis_(H.basis()),NumVertices_(F.size()),SingleVertexOp_(MPO_matrix()),m_EvolutionOperators(TrotterDecomposition(H,0.0,0)),GoodInitial_(0),SaveAll_(0) {
    std::stringstream Evolvingnamestream;
    Evolvingnamestream << "Evolving_" << MPSName_;
    initial_weight_=F.makeRC(Evolvingnamestream.str()); //combination of initial state norm and overall phase
    //it is useful to do the above for measurements, just to ensure normalisation, and that all files exist.
    
    std::cout << "Initial state weight was " << initial_weight_ <<std::endl;

    if (initial_weight_!=0.0)
      GoodInitial_=1;
  }
  
  TEBD::TEBD(const MPO_matrix& H, FiniteMPS& F, double time_step_size, DataOutput& results, uMPXInt order,const State* blockstate_ptr, bool save_all_flag) : TimeBase(time_step_size,results),F_(F),MPSName_(F.name()),Basis_(H.basis()),NumVertices_(F.size()),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order,blockstate_ptr)),GoodInitial_(0),SaveAll_(save_all_flag) {

    SingleVertexOp_=std::move(MakeSingleSiteEvolutionOperatorFromLowTriMPO(H,time_step_size));
    
    std::stringstream Evolvingnamestream;
    Evolvingnamestream << "Evolving_" << MPSName_;
    initial_weight_=F.makeLC(Evolvingnamestream.str()); //combination of initial state norm and overall phase

    std::cout << "Initial state weight was " << initial_weight_ <<std::endl;

    if (initial_weight_!=0.0)
      GoodInitial_=1;

     if (SaveAll_){
	std::stringstream H_MPO_Filename;
	H_MPO_Filename << SAVEALLNAME << "_"<< current_time_step() <<"_H.MPO_matrix";
	H.store(H_MPO_Filename.str());
	std::stringstream State_Filename;
	State_Filename << SAVEALLNAME << "_"<< current_time_step() << "_" << MPSName_;
	F.makeLC(State_Filename.str());

      }
    
  }

  void TEBD::change_bond_operator(const MPO_matrix& H, double time_step_size, uMPXInt order, const State* blockstate_ptr){
    m_time_step_size=time_step_size;
    SingleVertexOp_=MakeSingleSiteEvolutionOperatorFromLowTriMPO(H,time_step_size);
    
    m_EvolutionOperators=TrotterDecomposition(H,time_step_size,order,blockstate_ptr);
     if (SaveAll_){
	std::stringstream H_MPO_Filename;
	H_MPO_Filename << SAVEALLNAME << "_"<< current_time_step() <<"_H.MPO_matrix";
	H.store(H_MPO_Filename.str());
      }
    
  }

  void TEBD::evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension, double minS, uMPXInt measurement_interval){
    //do the evolution
    if (GoodInitial_){
     
      if (m_EvolutionOperators.order()==0){
	left_canonise_measure_special(measurements, num_steps);
      }
      else if (m_EvolutionOperators.order()==1){
	std::cout <<"1st order time step evolution" <<std::endl;
	for (uMPXInt n=0;n<num_steps;++n){
	  //++m_current_time_step;
	  update_time();
	  std::cout << "Time " << current_time() << std::endl;
	  //this is the first half of the time step....
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	  if (NumVertices_>2) {//if we only have two vertices, then we should skip the even bond part altogether
	    left_canonise();
	    apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	  }
	  if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	    left_canonise_measure(measurements);
	  }
	  else {
	    left_canonise();
	  }
	  max_truncation_=0.0; //reset
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
	  update_time();
	  //	  ++m_current_time_step;
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
	  max_truncation_=0.0; //reset
	}
      }
      else if (m_EvolutionOperators.order()==4){
	std::cout <<"4th order time step evolution" <<std::endl;
	for (uMPXInt n=0;n<num_steps;++n){
	  //++m_current_time_step;
	  update_time();
	  std::cout << "Time " << current_time() << std::endl;
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	  left_canonise();
	  apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	  left_canonise();
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[2]),bond_dimension,minS);
	  left_canonise();
	  apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[3]),bond_dimension,minS);
	  left_canonise();
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[4]),bond_dimension,minS);
	  left_canonise();
	  apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[5]),bond_dimension,minS);
	  left_canonise();
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[6]),bond_dimension,minS);
	  left_canonise();
	  apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[7]),bond_dimension,minS);
	  left_canonise();
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[8]),bond_dimension,minS);
	  left_canonise();
	  apply_to_even_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[9]),bond_dimension,minS);
	  left_canonise();
	  apply_to_odd_bonds(*(m_EvolutionOperators.OrderedOperatorPtrs[10]),bond_dimension,minS);

	  if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	    left_canonise_measure(measurements);
	  }
	  else {
	    left_canonise();
	  }
	  max_truncation_=0.0; //reset
	}
      }
      else {
	std::cout << "Haven't implemented higher orders yet in TEBD::evolve(), doing nothing." << std::endl;
      }
    }
    else {
      std::cout << "Bad initial state, doing nothing. " <<std::endl;
    }
  }

  std::string TEBD::evolution_name() const {
    std::stringstream Evolvingnamestream;
    Evolvingnamestream << "Evolving_" << MPSName_;
    return Evolvingnamestream.str();
  }

  void TEBD::left_info(){
    for (uMPXInt v=1;v<=NumVertices_;++v){
      std::cout << "L" <<v<<std::endl;
      std::stringstream Leftstream;
      Leftstream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      load_MPS_matrix(Leftstream.str(),Basis_).print_indices();
    }
  }

  void TEBD::right_info(){
    for (uMPXInt v=1;v<=NumVertices_;++v){
      std::cout <<"R"<<v<<std::endl;
      std::stringstream Rightstream;
      Rightstream << "Evolving_" << MPSName_ << "_Right_" << v << ".MPS_matrix";
      load_MPS_matrix(Rightstream.str(),Basis_).print_indices();
    }
  }

  /* MPX_matrix MakeBondHamiltonian(const MPO_matrix& H, const std::string& SaveName) {
    std::cout << "Forming bond Hamiltonian" << std::endl;
    MPX_matrix LeftHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(1,H.dimsvector()[1]-1))).ZeroLastBlock());
    MPX_matrix RightHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(3,0))));
    MPX_matrix ans(reorder(contract(LeftHalf,0,RightHalf,0,std::vector<MPXPair>(1,MPXPair(3,1))),0,reorder032415,2));
    if (!SaveName.empty()) ans.store(SaveName);
    return ans;
  }*/
  
  MPX_matrix MakeOddBondHamiltonian(const MPO_matrix& H, const std::string& SaveName) {
    std::cout << "Forming odd bond Hamiltonian" << std::endl;
    //MPX_matrix LeftHalf(H.ExtractSubMPX(std::vector<MPXPair>{{MPXPair(1,H.dimsvector()[1]-1)}}));
    //MPX_matrix RightHalf(H.ExtractSubMPX(std::vector<MPXPair>{{MPXPair(3,0)}}));

    MPX_matrix LeftHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({H.dimsvector().at(1)-1,H.dimsvector().at(1)-1}),std::pair<MPXInt,MPXInt>({0,H.dimsvector().at(3)-1})));
    MPX_matrix RightHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({0,H.dimsvector().at(1)-1}),std::pair<MPXInt,MPXInt>({0,0})));
    
    MPX_matrix ans(reorder(contract(LeftHalf,0,RightHalf,0,std::vector<MPXPair>{{MPXPair(3,1)}}),0,reorder032415,2));
    if (!SaveName.empty()) ans.store(SaveName);
    return ans;
  }

  MPX_matrix MakeEvenBondHamiltonian(const MPO_matrix& H, const std::string& SaveName) {
    std::cout << "Forming even bond Hamiltonian" << std::endl;
    //Need to check here if there is no coupling part in the MPO (i.e. vertices are uncoupled).
    if (H.dimsvector()[1]-1==1 || H.dimsvector()[3]-1==1){
      std::cout << "Even bond is empty" << std::endl;
      //exit(1);
    }
    
    /*std::vector<MPXPair> LeftKeepers;
    LeftKeepers.emplace_back(1,H.dimsvector()[1]-1);
    std::cout <<  LeftKeepers.size() << std::endl;
    for (uMPXInt l=1;l<H.dimsvector()[3]-1;++l){
      LeftKeepers.emplace_back(3,l);
    }
    std::vector<MPXPair> RightKeepers;
    RightKeepers.emplace_back(3,0);
    for (uMPXInt r=1;r<H.dimsvector()[1]-1;++r){
      RightKeepers.emplace_back(1,r);
    }
    MPX_matrix LeftHalf(H.ExtractSubMPX(LeftKeepers));
    MPX_matrix RightHalf(H.ExtractSubMPX(RightKeepers));*/

    MPX_matrix LeftHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({H.dimsvector()[1]-1,H.dimsvector()[1]-1}),std::pair<MPXInt,MPXInt>({1,H.dimsvector()[3]-2})));
    MPX_matrix RightHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({1,H.dimsvector()[1]-2}),std::pair<MPXInt,MPXInt>({0,0})));
    
    MPX_matrix ans(reorder(contract(LeftHalf,0,RightHalf,0,std::vector<MPXPair>{{MPXPair(3,1)}}),0,reorder032415,2));
    if (!SaveName.empty()) ans.store(SaveName);
    return ans;
  }  

  MPX_matrix MakeBondEvolutionOperator(const MPX_matrix& BondH, double timestep,const State* blockstate_ptr){
    std::cout << "Forming bond evolution operator" << std::endl;
    std::vector<MPXIndex> indices;
    indices.emplace_back(BondH.Index(0));
    indices.emplace_back(BondH.Index(1));
    indices.emplace_back(BondH.Index(2));
    indices.emplace_back(BondH.Index(3));
    if (blockstate_ptr!=nullptr){
      return MPX_matrix(BondH.basis(),indices,2,Exponentiate(BondH.Eigs(*blockstate_ptr),std::complex<double>(0.0,-timestep)));
    }
    return MPX_matrix(BondH.basis(),indices,2,Exponentiate(BondH.Eigs(),std::complex<double>(0.0,-timestep)));
  }

  //uses lower triangular MPO
  MPO_matrix MakeSingleSiteEvolutionOperatorFromLowTriMPO(const MPO_matrix& H_MPO, double timestep){
    
    std::vector<MPXIndex> indices;
    indices.emplace_back(H_MPO.Index(0));
    indices.emplace_back(MPXIndex(H_MPO.Index(1),H_MPO.Index(1).size()-1));
    indices.emplace_back(H_MPO.Index(2));
    indices.emplace_back(MPXIndex(H_MPO.Index(3),0));
    return MPO_matrix(H_MPO.basis(),indices,Exponentiate(H_MPO.ExtractMPOBlock(std::pair<ajaj::MPXInt,ajaj::MPXInt>(H_MPO.Index(1).size()-1,H_MPO.Index(1).size()-1),std::pair<ajaj::MPXInt,ajaj::MPXInt>(0,0)).Eigs(),std::complex<double>(0.0,-timestep)));
   
  }

}

#include <vector>
#include <complex>
#include <iostream>

#include "ajaj_common.hpp"
//#include "sparse_interface.hpp"
#include "states.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"
#include "TEBD_routines.hpp"

namespace ajaj{


  UnitCell MakeProductStateUnitCell(const Basis& b,uMPXInt state_index, uMPXInt length){
    return UnitCell(b,MakeProductState(b,state_index),std::vector<double>(1,1.0),2);
  }

  TrotterDecomposition::TrotterDecomposition(const MPO_matrix& H,double time_step_size,uMPXInt order) : m_H_ptr(&H), m_time_step_size(time_step_size), m_order(order){
//make bond operator
    if (m_order>2){std::cout <<"Higher than 2nd order Trotter decomposition not curtrently implemented. Aborting..." << std::endl; exit(1);}
    if (m_order==0){std::cout <<"Zeroth order requested! Assuming 1st order instead" << std::endl; m_order=1;}
    //make the Hamiltonian for a single bond
    MPX_matrix BondH(MakeBondHamiltonian(*m_H_ptr));
    //now make the trotter operators  
    if (m_order==1){
      BondOperators.emplace_back(MakeBondEvolutionOperator(BondH,m_time_step_size));
      //safeish to have pointer to vector element if the vector no longer grows
      OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
    }
    if (m_order==2){ //second order is a special case with a nice simplification for timesteps where no measurement is made
      BondOperators.emplace_back(MakeBondEvolutionOperator(BondH,m_time_step_size/2.0));
      BondOperators.emplace_back(MakeBondEvolutionOperator(BondH,m_time_step_size));
      //safeish to have pointer to vector element if the vector no longer grows
      OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
      OrderedOperatorPtrs.emplace_back(&BondOperators.at(1));
      OrderedOperatorPtrs.emplace_back(&BondOperators.at(0));
    }
  }

  SeparatedTrotterDecomposition::SeparatedTrotterDecomposition(const MPO_matrix& H,double time_step_size,uMPXInt order) : H_(H),TimeStepSize_(time_step_size),Order_(order){
    //first generate a normal decomp
    TrotterDecomposition TD(H_,TimeStepSize_,Order_);
    //decompose bond operators and push_back
    const EigenStateArray& Spectrum=H.GetPhysicalSpectrum();
    for (uMPXInt p=0;p<TD.OrderedOperatorPtrs.size();++p){
      const MPX_matrix& bondop(*TD.OrderedOperatorPtrs[p]);
      bool flag=0;
      uMPXInt order=p;
      for (uMPXInt i=0;i<p;++i){
	if (&bondop==TD.OrderedOperatorPtrs[i]){ //check for pointer equality for bond ops
	  flag=1;
	  order=i;
	  break;
	}
      }
      OperatorOrder_.push_back(order);
      if (!flag) {
	MPXDecomposition decomp(reorder(bondop,0,reorder0213,2).SVD(0,BONDDECOMPTOL*H.GetPhysicalSpectrum().size()));
	for (auto&& s : decomp.Values){
	  s=sqrt(s);
	}
	UOperators_.emplace_back(contract(decomp.ColumnMatrix,0,MPX_matrix(Spectrum,decomp.RowMatrix.Index(0),decomp.Values),0,contract20));
	UbarOperators_.emplace_back(contract(MPX_matrix(Spectrum,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      }
    }
  }

  iTEBD::iTEBD(const MPO_matrix& H,const UnitCell& C, double time_step_size, DataOutput& results, uMPXInt order) : TimeBase(time_step_size,results),m_EvolutionOperators(TrotterDecomposition(H,time_step_size,order)),m_initial_unit(C),m_unit(C) {
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
    m_results.push(Data(real_results,complex_results));
    std::ofstream DensityFileStream_;
    DensityFileStream_.open("iTEBD_One_Vertex_Densities.dat",ios::out | ios::app);
    ortho.OutputPhysicalIndexDensities(DensityFileStream_);
    DensityFileStream_.close();
  }

  const UnitCell& iTEBD::evolve(uMPXInt num_steps, const std::vector<MPO_matrix>& measuredMPOs, uMPXInt bond_dimension, double minS, uMPXInt measurement_interval){
    //depends on order, as 2nd order is special
    if (order()==1){
      std::cout <<"1st order time step evolution" <<std::endl;
      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	//this is the first half of the time step....
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	//this is the second
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	//do we need to take a measurement?
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  UnitCell ortho(OrthogonaliseInversionSymmetric(m_unit));
	  this->do_measurements(ortho,measuredMPOs);
	}
      }
    }
    else if (order()==2){
      std::cout <<"2nd order time step evolution" <<std::endl;
      //2nd order special start
      apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[2]),bond_dimension,minS);
	  //swap order
	  m_unit.swap(0,1);
	  //measure etc.
	  m_unit=std::move(OrthogonaliseInversionSymmetric(m_unit));
	  this->do_measurements(m_unit,measuredMPOs);

	  //complete time step
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[0]),bond_dimension,minS);
	}
	else {
	  apply_and_decompose(*(m_EvolutionOperators.OrderedOperatorPtrs[1]),bond_dimension,minS);
	}
       
      }
    }
    else {
      std::cout << "Haven't implemented higher orders yet in iTEBD::evolve()" << std::endl; exit(1);
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
    m_truncation=1.0-SquareSumRescale(m_unit.Lambdas.at(1),1.0);
#ifndef NDEBUG
    std::cout << "Swap sites" << std::endl;
#endif
    m_unit.swap(0,1);
#ifndef NDEBUG
    std::cout << "Return new unit cell" << std::endl;
#endif
    return m_unit;
  }

  TEBD::TEBD(const MPO_matrix& H,const std::string& MPSName, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order) : TimeBase(time_step_size,results), Spectrum_(H.GetPhysicalSpectrum()), MPSName_(MPSName), NumVertices_(NumVertices), SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),EvolutionOps_(SeparatedTrotterDecomposition(H,time_step_size,order)) {
    //canonize and store
    for (uMPXInt n=1;n<=NumVertices_/2;++n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Left_" << n << ".MPS_matrix";
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << n << ".MPS_matrix";
      MPS_matrix current(load_MPS_matrix(Initialnamestream.str(),Spectrum_));
      current.store(Evolvingnamestream.str());
    }
    //do something with lambda
    std::stringstream Lambdanamestream;
    Lambdanamestream << MPSName_ << "_Lambda_" << NumVertices_/2 << "_" << NumVertices_/2 << ".MPX_matrix";
    MPX_matrix Vd(load_MPX_matrix(Lambdanamestream.str(),Spectrum_));
    //left canonize the right hand parts
    for (uMPXInt n=NumVertices_/2;n>0;--n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName_ << "_Right_" << n << ".MPS_matrix";
      //contract Vd onto current
      //and reshape to left form
      //do SVD
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(Initialnamestream.str(),Spectrum_),0,contract10)).left_shape().SVD());
     
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Spectrum_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));

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

  TEBD::TEBD(const MPO_matrix& H, const std::string& MPSName, const MPS_matrix& InitialMPS_matrix, uMPXInt NumVertices, double time_step_size, DataOutput& results, uMPXInt order) : TimeBase(time_step_size,results), Spectrum_(InitialMPS_matrix.GetPhysicalSpectrum()), MPSName_(MPSName), NumVertices_(NumVertices), SingleVertexOp_(MakeSingleSiteEvolutionOperator(H,time_step_size)),EvolutionOps_(SeparatedTrotterDecomposition(H,time_step_size,order)) {
    //save initial MPS matrices
    //should all be left type!!!

    //SingleVertexOp_.print_indices();
    //SingleVertexOp_.print_sparse_info();


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

  void TEBD::left_info(){
    for (uMPXInt v=1;v<=NumVertices_;++v){
      std::cout << "L" <<v<<std::endl;
      std::stringstream Leftstream;
      Leftstream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      load_MPS_matrix(Leftstream.str(),Spectrum_).print_indices();
    }
  }

  void TEBD::right_info(){
    for (uMPXInt v=1;v<=NumVertices_;++v){
      std::cout <<"R"<<v<<std::endl;
      std::stringstream Rightstream;
      Rightstream << "Evolving_" << MPSName_ << "_Right_" << v << ".MPS_matrix";
      load_MPS_matrix(Rightstream.str(),Spectrum_).print_indices();
    }
  }

  void TEBD::left_canonise(uMPXInt chi,double minS){
    std::cout << "Begin left canonisation" << std::endl; 
    MPX_matrix Vd(Spectrum_);
    {
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Spectrum_).left_shape().SVD(chi,minS));
      RightEnddecomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << RightEnddecomp.Values.size() << std::endl;
      Vd=std::move(contract(MPX_matrix(Spectrum_,RightEnddecomp.RowMatrix.Index(0),RightEnddecomp.Values),0,RightEnddecomp.RowMatrix,0,contract10));
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
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Spectrum_),0,contract10)).left_shape().SVD(chi,minS));
      decomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
      
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Spectrum_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
      //store updated 'current' vertex matrix
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      decomp.ColumnMatrix.store(LeftNamestream.str());
    }
    //at end Vd should be trivial
    std::cout << "Norm at end of left canonisation: " << abs(Vd.Trace()) << std::endl; 
  }

  void TEBD::left_canonise_measure(uMPXInt chi,double minS,std::vector<MultiVertexMeasurement>& measurements){
    std::cout << "Begin left canonisation" << std::endl; 

    //it can be useful to measure as we step through the system
    std::vector<double> real_results(1,current_time());
    //std::vector<MPX_matrix> measure_tensors(measurements.size(),MPX_matrix(Spectrum_));

    //we'd also like to measure the overlap with the initial state...
    MPX_matrix overlap_matrix(Spectrum_);
    MPX_matrix Vd(Spectrum_);

    {
      //treats first vertex as a special case
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Spectrum_).left_shape().SVD(chi,minS));
      RightEnddecomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << RightEnddecomp.Values.size() << std::endl;
      if (NumVertices_==2){ //special case for two vertices only
	//measure truncation and entropy
	real_results.push_back(RightEnddecomp.getRescaleDifference());
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
	MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Spectrum_));
	std::vector<MPXIndex> indices;
	indices.emplace_back(1,initial.Index(2));
	indices.emplace_back(0,current.Index(2));
	overlap_matrix=std::move(MPX_matrix(Spectrum_,indices,1,contract_to_sparse(initial,1,current,0,contract00)));
      }
      Vd=std::move(contract(MPX_matrix(Spectrum_,RightEnddecomp.RowMatrix.Index(0),RightEnddecomp.Values),0,RightEnddecomp.RowMatrix,0,contract10));
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
      MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Spectrum_),0,contract10)).left_shape().SVD(chi,minS));
      decomp.SquareRescale(1.0);
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
      if (v==NumVertices_/2){ //if we only have two chains, this is never obeyed...
	//measure truncation and entropy
	real_results.push_back(decomp.getRescaleDifference());
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
	MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Spectrum_));
	overlap_matrix=std::move(contract(initial,1,contract(overlap_matrix,0,current,0,contract11),0,contract0110));
      }
      //record row vectors (Vdagger part) as new Vd
      Vd=std::move(contract(MPX_matrix(Spectrum_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
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

  void TEBD::right_canonise(double minS){
    //1st is a special case
    std::cout << "Begin right canonisation" << std::endl; 

    MPX_matrix U(Spectrum_);
    {
      std::stringstream LeftEndstream;
      LeftEndstream << "Evolving_" << MPSName_ << "_Left_" << NumVertices_ << ".MPS_matrix";
      MPXDecomposition LeftEnddecomp(load_MPS_matrix(LeftEndstream.str(),Spectrum_).right_shape().SVD(0,minS));
      //LeftEnddecomp.SquareRescale(1.0); //don't need this unless truncating
      std::cout << "Bond dimension: " << LeftEnddecomp.Values.size() << std::endl;

      //LeftEnddecomp.printValues();
      U=std::move(contract(LeftEnddecomp.ColumnMatrix,0,MPX_matrix(Spectrum_,LeftEnddecomp.RowMatrix.Index(0),LeftEnddecomp.Values),0,contract10));
      std::stringstream RightStartstream;
      RightStartstream << "Evolving_" << MPSName_ << "_Right_1" << ".MPS_matrix";
      LeftEnddecomp.RowMatrix.store(RightStartstream.str());
    }

    for (uMPXInt n=NumVertices_-1;n>0;--n){
      std::stringstream LeftNamestream;
      LeftNamestream << "Evolving_" << MPSName_ << "_Left_" << n << ".MPS_matrix";
      //contract U onto current
      //and reshape to left form
      //do SVD
      MPXDecomposition decomp(MPS_matrix(contract(load_MPS_matrix(LeftNamestream.str(),Spectrum_),0,U,0,contract20)).right_shape().SVD(0,minS));
      std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;

      //decomp.SquareRescale(1.0); //don't need this unless truncating
      //record row vectors (Vdagger part) as new Vd
      U=std::move(contract(decomp.ColumnMatrix,0,MPX_matrix(Spectrum_,decomp.ColumnMatrix.Index(1),decomp.Values),0,contract10));
      //store updated 'current' vertex matrix
      std::stringstream RightNamestream;
      RightNamestream << "Evolving_" << MPSName_ << "_Right_" << NumVertices_-n+1 << ".MPS_matrix";
      decomp.RowMatrix.store(RightNamestream.str());
    }
    //at end U should be trivial
    std::cout << "Norm at end of right canonisation: " << abs(U.Trace()) << std::endl; 
  }

  void TEBD::apply_to_odd_bonds(const MPX_matrix& U, const MPX_matrix&  Ubar){
    //load each and apply
    std::cout << "Odd bonds" << std::endl;
    for (uMPXInt v=1;v<=NumVertices_;++v){
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      MPS_matrix current(load_MPS_matrix(Evolvingnamestream.str(),Spectrum_));
      if (v % 2){
	MPS_matrix(std::move(contract(U,0,current,0,contract10).CombineSimilarMatrixIndices())).store(Evolvingnamestream.str());
      }
      else {
	MPS_matrix(std::move(contract(Ubar,0,current,0,contract20).CombineSimilarMatrixIndices())).store(Evolvingnamestream.str());
      }
    }
  }

  void TEBD::apply_to_even_bonds(const MPX_matrix& U, const MPX_matrix&  Ubar){
    //skip first vertex, apply special to last
    std::cout << "Even bonds" << std::endl;
    for (uMPXInt v=2;v<=NumVertices_;++v){
      std::stringstream Evolvingnamestream;
      Evolvingnamestream << "Evolving_" << MPSName_ << "_Left_" << v << ".MPS_matrix";
      MPS_matrix current(load_MPS_matrix(Evolvingnamestream.str(),Spectrum_));
      //current.print_indices();
      if (v==NumVertices_){ //last site
	MPS_matrix(std::move(contract(SingleVertexOp_,0,current,0,contract10).CombineSimilarMatrixIndices())).store(Evolvingnamestream.str());
      }
      else if (v % 2){ //if odd
	MPS_matrix(std::move(contract(Ubar,0,current,0,contract20).CombineSimilarMatrixIndices())).store(Evolvingnamestream.str());
      }
      else { //if even
	MPS_matrix(std::move(contract(U,0,current,0,contract10).CombineSimilarMatrixIndices())).store(Evolvingnamestream.str());
      }     
    }
  }

  void TEBD::evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements, uMPXInt bond_dimension, double minS, uMPXInt measurement_interval){
    //do the evolution
    if (EvolutionOps_.order()==1){
      std::cout <<"1st order time step evolution" <<std::endl;
      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	//this is the first half of the time step....
	apply_to_odd_bonds(EvolutionOps_.getU(0),EvolutionOps_.getUbar(0));
	apply_to_even_bonds(EvolutionOps_.getU(0),EvolutionOps_.getUbar(0));
	right_canonise(minS);
	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  left_canonise_measure(bond_dimension,minS,measurements);
	}
	else {
	  left_canonise(bond_dimension,minS);
	}
      }
    }
    else if (EvolutionOps_.order()==2){
      std::cout <<"2nd order time step evolution" <<std::endl;
      //bond order grows rapidly, so need to compress after each application of a set of bond operators
      //first step
      apply_to_odd_bonds(EvolutionOps_.getU(0),EvolutionOps_.getUbar(0));
      right_canonise(minS);
      left_canonise(bond_dimension,minS);

      for (uMPXInt n=0;n<num_steps;++n){
	++m_current_time_step;
	std::cout << "Time " << current_time() << std::endl;
	apply_to_even_bonds(EvolutionOps_.getU(1),EvolutionOps_.getUbar(1));
	right_canonise(minS);
	left_canonise(bond_dimension,minS);

	if (m_current_time_step % measurement_interval==0) /*make measurement*/ {
	  apply_to_odd_bonds(EvolutionOps_.getU(2),EvolutionOps_.getUbar(2));
	  right_canonise(minS);
	  left_canonise_measure(bond_dimension,minS,measurements);//measurement
	  apply_to_odd_bonds(EvolutionOps_.getU(0),EvolutionOps_.getUbar(0));
	  right_canonise(minS);
	  left_canonise(bond_dimension,minS);
	}
	else {
	  apply_to_odd_bonds(EvolutionOps_.getU(1),EvolutionOps_.getUbar(1));
	  right_canonise(minS);
	  left_canonise(bond_dimension,minS);
	}
      }
    }
    else {
      std::cout << "Haven't implemented higher orders yet in TEBD::evolve()" << std::endl; exit(1);
    }
  }


  MPX_matrix MakeBondHamiltonian(const MPO_matrix& H) {
    std::cout << "Forming bond Hamiltonian" << std::endl;
    MPX_matrix LeftHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(1,H.dimsvector()[1]-1))).RestrictColumnIndex());
    MPX_matrix RightHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(3,0))));
    return reorder(contract(LeftHalf,0,RightHalf,0,std::vector<MPXPair>(1,MPXPair(3,1))),0,reorder032415,2);
  }

  MPX_matrix MakeBondEvolutionOperator(const MPX_matrix& BondH, double timestep){
    std::cout << "Forming bond evolution operator" << std::endl;
    std::vector<MPXIndex> indices;
    indices.emplace_back(BondH.Index(0));
    indices.emplace_back(BondH.Index(1));
    indices.emplace_back(BondH.Index(2));
    indices.emplace_back(BondH.Index(3));
    return MPX_matrix(BondH.GetPhysicalSpectrum(),indices,2,Exponentiate(BondH.Eigs(),std::complex<double>(0.0,-timestep)));
  }

  MPX_matrix MakeSingleSiteEvolutionOperator(const MPO_matrix& H, double timestep){
    std::cout << "Forming single site evolution operator (diagonal)" << std::endl;
    std::vector<std::complex<double> > phase;
    for (auto i : H.GetPhysicalSpectrum().Energies()){
      phase.push_back(std::complex<double>(cos(timestep*i),-sin(timestep*i)));
    }
    return MPX_matrix(H.GetPhysicalSpectrum(),H.Index(0),phase);
  }

  bool CheckMPSFilesExist(const std::string& MPSName,uMPXInt NumVertices){
    //open 
    
    for (auto n=1;n<=NumVertices/2;++n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName << "_Left_" << n << ".MPS_matrix";
      std::ofstream infile;
      infile.open(Initialnamestream.str().c_str(),ios::in | ios::binary); 
      if (!infile.is_open()) {
	std::cout << "Couldn't open file: " << Initialnamestream.str() <<std::endl;
	return 0;
      }
      infile.close();
    }
    //open check opens

    std::stringstream Lambdanamestream;
    Lambdanamestream << MPSName << "_Lambda_" << NumVertices/2 << "_" << NumVertices/2 << ".MPX_matrix";
    std::ofstream Lfile;
    Lfile.open(Lambdanamestream.str().c_str(),ios::in | ios::binary); 
    if (!Lfile.is_open()) {
      std::cout << "Couldn't open file: " << Lambdanamestream.str() <<std::endl;
      return 0;
    }
    Lfile.close();

    for (auto n=1;n<=NumVertices/2;++n){
      std::stringstream Initialnamestream;
      Initialnamestream << MPSName << "_Right_" << n << ".MPS_matrix";
      std::ofstream infile;
      infile.open(Initialnamestream.str().c_str(),ios::in | ios::binary); 
      if (!infile.is_open()) {
	std::cout << "Couldn't open file: " << Initialnamestream.str() <<std::endl;
	return 0;
      }
      infile.close();
    }
    return 1;
  }


}

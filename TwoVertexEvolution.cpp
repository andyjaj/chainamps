#include <sstream>

#include "states.hpp"
#include "MPX.hpp"
#include "data.hpp"
#include "measurement.hpp"
#include "TwoVertexEvolution.hpp"

namespace ajaj{

  MPX_matrix Make2VBondHamiltonian(const MPO_matrix& H) {
    std::cout << "Forming bond Hamiltonian" << std::endl;
    MPX_matrix LeftHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({H.dimsvector()[1]-1,H.dimsvector()[1]-1}),std::pair<MPXInt,MPXInt>({0,H.dimsvector()[3]-1})));
    MPX_matrix RightHalf(H.ExtractMPOBlock(std::pair<MPXInt,MPXInt>({0,H.dimsvector()[1]-1}),std::pair<MPXInt,MPXInt>({0,0})));
    //MPX_matrix LeftHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(1,H.dimsvector()[1]-1))));
    //MPX_matrix RightHalf(H.ExtractSubMPX(std::vector<MPXPair>(1,MPXPair(3,0))));
    return reorder(contract(LeftHalf,0,RightHalf,0,std::vector<MPXPair>(1,MPXPair(3,1))),0,reorder032415,2);
  }

  MPX_matrix Make2VEvolutionOperator(const MPX_matrix& BondH, double timestep,const State& blockstate){    
    std::cout << "Forming evolution operator" << std::endl;
    std::vector<MPXIndex> indices;
    indices.emplace_back(BondH.Index(0));
    indices.emplace_back(BondH.Index(1));
    indices.emplace_back(BondH.Index(2));
    indices.emplace_back(BondH.Index(3));
    return MPX_matrix(BondH.basis(),indices,2,Exponentiate(BondH.Eigs(blockstate),std::complex<double>(0.0,-timestep)));
  }

  TwoVE::TwoVE(const MPO_matrix& HMPO, FiniteMPS& F, double time_step_size, DataOutput& results, DataOutput& fresults ,const State& blockstate) : TimeBase(time_step_size,results,fresults),MPSName_(F.name()),Basis_(HMPO.basis()),EvolutionOperator_(Make2VEvolutionOperator(Make2VBondHamiltonian(HMPO),time_step_size,blockstate)),GoodInitial_(0),BlockState_(blockstate) {
   
    std::stringstream Evolvingnamestream;
    Evolvingnamestream << "Evolving_" << MPSName_;
    std::complex<double> initial_weight(F.makeLC(Evolvingnamestream.str())); //copy state and canonise just in case
    
    std::cout << "Initial state weight was " << initial_weight <<std::endl;
    if (initial_weight!=0.0)
      GoodInitial_=1;
  }

  void TwoVE::change_bond_operator(const MPO_matrix& HMPO, double time_step_size){
    m_time_step_size=time_step_size;
    EvolutionOperator_=Make2VEvolutionOperator(HMPO,time_step_size,BlockState_);
  }

  void TwoVE::evolve(uMPXInt num_steps, std::vector<MultiVertexMeasurement>& measurements){
    //do the evolution
    if (GoodInitial_){
      for (uMPXInt n=0;n<num_steps;++n){
	//++m_current_time_step;
	update_time();
	std::cout << "Time " << current_time() << std::endl;
	//this is the first half of the time step....
	apply();
	left_canonise_measure(measurements);
      }
    }
    else {std::cout << "Bad initial state!"<<std::endl;}
  }

  void TwoVE::apply(){
    std::cout << "Applying evolution operator" << std::endl;
    //Start state should be left canonical
    //At end state is right canonical but truncated.
    //can then use left canonize measure if necessary
    //load each and apply
    //need a starter MPS
    std::stringstream StartNameStream;
    StartNameStream << "Evolving_" << MPSName_ << "_Left_2.MPS_matrix";
    std::string RightName=StartNameStream.str();
    MPS_matrix R(load_MPS_matrix(RightName,Basis_));
    std::stringstream LNameStream;
    LNameStream << "Evolving_" << MPSName_ << "_Left_1.MPS_matrix";
    MPSDecomposition decomp(reorder(contract(EvolutionOperator_,0,contract(load_MPS_matrix(LNameStream.str(),Basis_),0,R,0,contract21),0,contract2032),0,reorder0213,2).SVD());
    
    std::stringstream Store1NameStream;
    Store1NameStream << "Evolving_" << MPSName_ << "_Right_1.MPS_matrix";
    decomp.RightMatrix.store(Store1NameStream.str());
    //need to 'rotate left part'
    
    std::stringstream Store2NameStream;
    Store2NameStream << "Evolving_" << MPSName_ << "_Right_2.MPS_matrix";
    MPXDecomposition rot(MPS_matrix(contract(decomp.LeftMatrix,0,MPX_matrix(Basis_,decomp.LeftMatrix.Index(2),decomp.Values),0,contract20)).right_shape().SVD());
    rot.RowMatrix.store(Store2NameStream.str()); //now should be right canonical
    
    //check norm
    std::cout << "Norm: " << contract(rot.ColumnMatrix,0,MPX_matrix(Basis_,rot.ColumnMatrix.Index(1),rot.Values),0,contract10).Trace() <<std::endl;
  }

  void TwoVE::left_canonise_measure(std::vector<MultiVertexMeasurement>& measurements){
    std::cout << "Begin left canonisation" << std::endl; 
    
    //it can be useful to measure as we step through the system
    std::vector<double> real_results(1,current_time());
    //we'd also like to measure the overlap with the initial state...
    MPX_matrix overlap_matrix(Basis_);
    MPX_matrix Vd(Basis_);
    {
      //treats first vertex as a special case
      std::stringstream RightEndstream;
      RightEndstream << "Evolving_" << MPSName_ << "_Right_2.MPS_matrix";
      MPXDecomposition RightEnddecomp(load_MPS_matrix(RightEndstream.str(),Basis_).left_shape().SVD());
      RightEnddecomp.SquareRescale(1.0);
      //record truncation and entropy
      real_results.push_back(0.0);
      real_results.push_back(entropy(RightEnddecomp.Values));
      
      for (auto&& m : measurements){
	if (m.finish()>2){
	  std::cout << "MultiVertexMeasurements incorrectly defined, final measurement position is > 2: " << m.finish() << std::endl;
	  exit(1);
	}
	m.update(1,RightEnddecomp); //first vertex is position 1!
      }
      
      const MPX_matrix& current=RightEnddecomp.ColumnMatrix;
      std::stringstream Initialstream;
      Initialstream << MPSName_ << "_Left_1.MPS_matrix";
      MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Basis_));
      std::vector<MPXIndex> indices;
      indices.emplace_back(1,initial.Index(2));
      indices.emplace_back(0,current.Index(2));
      overlap_matrix=std::move(MPX_matrix(Basis_,indices,1,contract_to_sparse(initial,1,current,0,contract00)));
      
      Vd=std::move(contract(MPX_matrix(Basis_,RightEnddecomp.RowMatrix.Index(0),RightEnddecomp.Values),0,RightEnddecomp.RowMatrix,0,contract10));
      std::stringstream LeftStartstream;
      LeftStartstream << "Evolving_" << MPSName_ << "_Left_1.MPS_matrix";
      RightEnddecomp.ColumnMatrix.store(LeftStartstream.str());
      
    }
    std::stringstream RightNamestream;
    RightNamestream << "Evolving_" << MPSName_ << "_Right_1.MPS_matrix";
    //contract U onto current
    //and reshape to left form
    //do SVD
    MPXDecomposition decomp(MPS_matrix(contract(Vd,0,load_MPS_matrix(RightNamestream.str(),Basis_),0,contract10)).left_shape().SVD());
    decomp.SquareRescale(1.0);
    std::cout << "Bond dimension: " << decomp.Values.size() << std::endl;
    
    for (auto&& m : measurements){
      m.update(2,decomp);
    }
    const MPX_matrix& current=decomp.ColumnMatrix;
    std::stringstream Initialstream;
    Initialstream << MPSName_ << "_Left_2.MPS_matrix";
    std::cout << Initialstream.str() << std::endl;
    MPS_matrix initial(load_MPS_matrix(Initialstream.str(),Basis_));
    overlap_matrix=std::move(contract(initial,1,contract(overlap_matrix,0,current,0,contract11),0,contract0110));
    
    //record row vectors (Vdagger part) as new Vd
    Vd=std::move(contract(MPX_matrix(Basis_,decomp.RowMatrix.Index(0),decomp.Values),0,decomp.RowMatrix,0,contract10));
    //store updated 'current' vertex matrix
    std::stringstream LeftNamestream;
    LeftNamestream << "Evolving_" << MPSName_ << "_Left_2.MPS_matrix";
    decomp.ColumnMatrix.store(LeftNamestream.str());

    std::vector<std::complex<double> > complex_results;

    std::complex<double> overlap=overlap_matrix.Trace();
    real_results.push_back(abs(overlap));
    complex_results.push_back(overlap);
    
    for (auto&& m : measurements){
      complex_results.push_back(m.result());
    }
    m_results.push(m_current_time_step,Data(real_results,complex_results));
    //at end Vd should be trivial

    std::cout << "Norm at end of measurement process: " << abs(Vd.Trace()) << std::endl; 
  }
  
}

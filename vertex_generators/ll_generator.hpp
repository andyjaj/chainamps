#ifndef THEORY_LL_H
#define THEORY_LL_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

////////////////////////
#ifndef RMK
#define RMK
namespace rmk { //namespace for Robert's ancilliary arrays and definitions
  /*definitions for construction of full states*/
  static const ajaj::uMPXInt MaxStates=2000; /*maximum number of states allowed in the Hilbert space of a single boson*/
  static const ajaj::uMPXInt MaxN=80; /*largest allowed value of N in constructing highest weight states |N,M>*/
  static const ajaj::uMPXInt MaxM=3; /*largest allowed value of M in constructing highest weight states |N,M>*/

  static const ajaj::uMPXInt MaxPart=300; /*maximum number of partitions allowed at a given level; used in denumerating number of chiral states*/
  static const ajaj::uMPXInt MaxChiralStates=600; /*largest allowed number of chiral states*/
  static const ajaj::uMPXInt UL=14;  /*the level up to (but not including) which we keep chiral states*/
  static const ajaj::uMPXInt LL=0;
  
  static const ajaj::uMPXInt MaxNum_a=20; /*the maximum number of chiral a's that we allow in constructing a state, i.e.  
					    for a state, a_{-n_1} ... a_{-n_{k_L}}\bar a_{-n_1} ... \bar a_{-n_{k_R}}|N,M>,
					    both k_L and k_R must be less than MaxNum_a; we must have MaxNum_a <= MaxLev */

  static const ajaj::uMPXInt MaxLev=20; /*the largest chiral level that we can possibly consider in forming states*/

  static ajaj::uMPXInt nstates_lev[MaxLev] = {1,2,4,7,12,19,30,45,67,97,139,195,272,373,508,684,915,1212,1597,2087}; 
  /*a vector containing the number of different partitions at a given level (for levels 1 through MaxLev)*/

  
  double state_en[MaxStates];

  //ajaj::MPXInt nstate; /*index keeping track of how many non-chiral states there are*/
  ajaj::MPXInt state_chiral[3][MaxStates], /*chiral info: left chiral state, right chiral state, total number of left and right a's of state*/
    state_Z[3][MaxStates], /*Z quantum numbers: N, M, mom*/
    state_metaZ[MaxStates], /*Z quantum number incorporating all 3 Z quantum numbers*/
    state_metaZ_ind[MaxStates], /*index of Z quantum number incorporating all 3 Z quantum numbers*/
    state_parity[2][MaxStates]; /*parity quantum numbers: \pm 1 for parity of state; and sign of application of parity operator upon positive component of state*/

  /*routine to compute non-chiral matrix elements*/
  double non_chir_me(int i, int j);

  /*routine to sort in order of energy non-chiral states*/
  void shell_en(int n);

  /*routines to compute non-chiral matrix elements*/
  double non_chir_psi_dagger(int i, int j);
  double non_chir_psi(int i, int j);

  ajaj::MPXInt num_metaZ; /*range of metaZ quantum number: from 0 to num_metaZ-1*/
 
#include "./ll/chiral_states_bosons_fA.hpp"
#include "./ll/element_computation_bosons.hpp"
#include "./ll/me_chir_bosons.hpp"
#include "./ll/non_chiral_states_bosons_no_parity_fA.hpp"
#include "./ll/non_chir_me_no_parity.hpp"
}
#endif

namespace ll{
  ajaj::Vertex VertexGenerator(const ajaj::VertexParameterArray& inputs)
  {
    //initialise vertex object
    ajaj::Vertex ModelVertex(inputs);
    ModelVertex.ChargeRules=ajaj::QNVector({{0,0,0}}); //three Z like quantum numbers, take the first as LL momentum, next as N the M

    const ajaj::QNVector& ChargeRules=ModelVertex.ChargeRules;

    const double R = inputs[0].Value; /*chain length, NOT RADIUS!*/
    const double Beta=inputs[1].Value; /*inverse compactification radius*/
    const double spectrum_cutoff=inputs[2].Value;  /*energy cutoff on states*/
    double tpi_R=2.0*M_PI/R;
    /*compute chiral states: this routine computes the possible ways to form
      chiral states up to a certain level, i.e. states of the form a_{n_1} \cdots a_{n_N}|N,M>*/ 
    rmk::chiral_states();
    /*compute chiral matrix elements*/
    rmk::me_chir(Beta);
    /*compute non-chiral states*/  
    ajaj::uMPXInt nstates=rmk::non_chir_bosons(tpi_R,Beta,spectrum_cutoff,0);

    //populate spectrum
    for (ajaj::uMPXInt n=0;n<nstates;++n){
      //Robert's charges are not in the order I prefer
      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(rmk::state_Z[2][n]),ajaj::QuantumNumberInt(rmk::state_Z[0][n]),ajaj::QuantumNumberInt(rmk::state_Z[1][n])}}),rmk::state_en[n]));
    }

    if (ModelVertex.Spectrum.size()!=nstates){
      std::cout << "Mismatch error " << ModelVertex.Spectrum.size() << " " << nstates << std::endl; exit(1);
    }

    std::cout << "End generating chain spectrum, generated " << nstates << " states" << std::endl;
    std::cout << "Starting Matrix Elements" << std::endl;

    ModelVertex.Operators.push_back(ajaj::VertexOperator("Psi",ModelVertex.Spectrum.size()));      
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Psi_dagger",ModelVertex.Spectrum.size()));      
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Density",ModelVertex.Spectrum.size()));   
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Density_integrated",ModelVertex.Spectrum.size()));   
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Phase",ModelVertex.Spectrum.size()));

    double factor=pow(tpi_R,Beta*Beta);

    for (ajaj::MPXInt col=0;col<ModelVertex.Spectrum.size();++col){
      for (ajaj::MPXInt row=0;row<ModelVertex.Spectrum.size();++row){
	std::complex<double> M_E_=rmk::non_chir_psi(row,col)*factor;
	std::complex<double> M_E_dagger_=rmk::non_chir_psi_dagger(row,col)*factor;
	std::complex<double> density=-1.*rmk::density_me(Beta,row,col)/(2.0 *R*Beta);
	std::complex<double> density_int=-1.*R*rmk::density_int_me(Beta,row,col)/(2.0 *R*Beta);
	std::complex<double> phase=rmk::phase_me(Beta,row,col);

	if (abs(M_E_)>SPARSETOL){
	  ModelVertex.Operators.at(0).MatrixElements.entry(row,col,M_E_);
	}
	if (abs(M_E_dagger_)>SPARSETOL){
	  ModelVertex.Operators.at(1).MatrixElements.entry(row,col,M_E_dagger_);
	}
	if (abs(density)>SPARSETOL){
	  ModelVertex.Operators.at(2).MatrixElements.entry(row,col,density);
	}
	if (abs(density_int)>SPARSETOL){
	  ModelVertex.Operators.at(3).MatrixElements.entry(row,col,density_int);
	}
	if (abs(phase)>SPARSETOL){
	  ModelVertex.Operators.at(4).MatrixElements.entry(row,col,phase);
	}
      }
    }

    for (ajaj::uMPXInt nmode=1;nmode<rmk::UL;++nmode) {
      std::ostringstream mn;
      mn << "Mode_" << nmode;
      ModelVertex.Operators.push_back(ajaj::VertexOperator(mn.str(),ModelVertex.Spectrum.size()));
      for (ajaj::uMPXInt i=0;i<ModelVertex.Spectrum.size();++i){
	double occ=rmk::n_a(i,i,nmode);
	if (occ!=0.0){
	  ModelVertex.Operators.back().MatrixElements.entry(i,i,occ);
	}
	for (ajaj::uMPXInt j=i+1;j<ModelVertex.Spectrum.size();++j){
	  double occ=rmk::n_a(i,j,nmode);
	  if (occ!=0.0){
	    ModelVertex.Operators.back().MatrixElements.entry(i,j,occ);
	    ModelVertex.Operators.back().MatrixElements.entry(j,i,occ);
	  }
	}
      }
    }

    for (ajaj::VertexOperatorArray::iterator it=ModelVertex.Operators.begin();it!=ModelVertex.Operators.end();++it){
      it->MatrixElements.finalise();
#ifndef DNDEBUG
      if (ModelVertex.basis().size()<=5){
	it->print();
      }
#endif
    }

    ModelVertex.Operators.push_back(ajaj::VertexOperator("Vertex_Hamiltonian"));
    ModelVertex.Operators.back().MatrixElements=ajaj::SparseMatrix(ModelVertex.basis().Energies());

    return ModelVertex;
  }

  ajaj::MPO_matrix MakeHamiltonian(const ajaj::Vertex& modelvertex, const ajaj::CouplingArray& couplingparams){
    //Lower triangular MPO
    // I           0    0         0
    // JPsidagger  0    0         0
    // JPsi
    // HV          Psi' Psidagger I        //note the charges Q, are such that Q[Psidagger[i]]+Q[Psi'[i]]=0

    double R(modelvertex.Parameters[0].Value);
    ajaj::QNCombinations differencecombinations(modelvertex.Spectrum,1); //1 means use difference

    ajaj::MPXInt lineardim(modelvertex.Spectrum.size()*2*(1+differencecombinations.size()));
    ajaj::MPXInt offset_to_last_block=modelvertex.Spectrum.size()*(1+2*differencecombinations.size()); //offset to get to the last row of operators
    ajaj::SparseMatrix M(lineardim,lineardim,lineardim);

    //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
    for (size_t i=0;i<modelvertex.Spectrum.size();++i){
      M.entry(i,i,1.0);
      M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
      M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }
    if (couplingparams[0].Value.imag()!=0.0){
      std::cout << "ERROR: tunnelling must have a real value. Imag part=" << couplingparams[0].Value.imag() <<std::endl;
    }
    double tunnelling=couplingparams[0].Value.real();

    //quickest to use psi and its Hermitian conjugate
    ajaj::Sparseint psi_col_offset=1; //+1 for identity matrix in first block
    ajaj::Sparseint psi_row_offset=differencecombinations.size()+1; //reversed order
    ajaj::Sparseint psidagger_col_offset=differencecombinations.size()+1; //+1 for identity matrix in first block
    ajaj::Sparseint psidagger_row_offset=1; //reversed order
    for (ajaj::Sparseint col=0;col<modelvertex.Operators[0].MatrixElements.cols();++col){ //use psi matrix elements
      for (ajaj::Sparseint p=modelvertex.Operators[0].MatrixElements.get_p(col);p<modelvertex.Operators[0].MatrixElements.get_p(col+1);++p){
	ajaj::Sparseint i=modelvertex.Operators[0].MatrixElements.get_i(p);
	ajaj::State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	//now look through charge groups to find MPO indices
	ajaj::Sparseint MPO_subcol=0;
	ajaj::Sparseint MPO_subrow=0;
	ajaj::Sparseint MPO_subcolpsidagger=0;
	ajaj::Sparseint MPO_subrowpsidagger=0;
	if (diffstate==-diffstate){ //involution block
	  for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	    if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
	      MPO_subcol=l;
	      MPO_subrow=l; //involution pairs don't get reversed
	      MPO_subcolpsidagger=l;
	      MPO_subrowpsidagger=l;
	      break;
	    }
	  }
	}
	else {
	  for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	    if (diffstate==differencecombinations.OrderedPairs[l].PairState){
	      MPO_subcol=differencecombinations.InvolutionPairs.size()+l;
	      MPO_subrow=differencecombinations.size()-l-1;
	      MPO_subcolpsidagger=MPO_subrow;
	      MPO_subrowpsidagger=MPO_subcol;
	      break;
	    }
	  }
	}
	std::complex<double> x=modelvertex.Operators[0].MatrixElements.get_x(p); //uses psi
	//do psi first
	//lowest block row
	M.entry(offset_to_last_block+i,modelvertex.Spectrum.size()*(psi_col_offset+MPO_subcol)+col,x);
	//first block col
	M.entry(modelvertex.Spectrum.size()*(psi_row_offset+MPO_subrow)+i,col,R*x*tunnelling);
	//now psi dagger so swap i and col
	//lowest block row
	M.entry(offset_to_last_block+col,modelvertex.Spectrum.size()*(psidagger_col_offset+MPO_subcolpsidagger)+i,conj(x));
	//first block col
	M.entry(modelvertex.Spectrum.size()*(psidagger_row_offset+MPO_subrowpsidagger)+col,i,R*conj(x)*tunnelling);
      }
    }

    //done forming array, now do indices
    ajaj::StateArray b;
    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    for (size_t c=0;c<2;++c){ //do for each operator
      for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	b.push_back(differencecombinations.InvolutionPairs[l].PairState);
      }
      for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	b.push_back(differencecombinations.OrderedPairs[l].PairState);
      }
    }
    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    
    std::vector<ajaj::MPXIndex> indices;
    indices.push_back(ajaj::MPXIndex(1,modelvertex.Spectrum)); //sigma primed
    indices.push_back(ajaj::MPXIndex(1,b)); //b_left
    indices.push_back(ajaj::MPXIndex(0,modelvertex.Spectrum)); //sigma
    indices.push_back(ajaj::MPXIndex(0,b)); //b_right

    return  ajaj::MPO_matrix(modelvertex.Spectrum, indices, M.finalise());
  }
}

#endif

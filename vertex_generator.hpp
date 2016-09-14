#ifndef THEORY_ISING_H
#define THEORY_ISING_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <algorithm>

#include "sparse_interface.hpp"

//static const unsigned short int MAX_NUM_CHAIN_STATES=500;//useful buffer size for forming matrix elements for the occupation number as we generate the spectrum 
inline double NSvacuum(const double Delta, const double R);
inline double Rvacuum(const double Delta, const double R);
inline double Kap(const double t, const double DR);
inline double S(const double DR);
complex<double> calculate_chainmatrixelement(const EigenStateArray& m_states, const MPXInt a, const MPXInt b, const MPXInt modenum_a, const MPXInt modenum_b, const vector<double> kappa_a, const vector<double> kappa_b, const vector<double> theta_a, const vector<double> theta_b, const double Delta,const double DR);

inline void check_mode(const int double_mode_num,const unsigned short int entry,const int measured_occupations, bool** occupations){
  for (int i=0;i<measured_occupations;++i){
    if (double_mode_num==i){occupations[i][entry]=1;break;}//we have a particle in this mode, and it can't be in any others
  }
};

inline void check_mode_ordered_R(const int mode_num,const unsigned short int entry,const int measured_occupations, bool** occupations){
  for (int i=0;i<3;++i){
    if (mode_num==i){occupations[i][entry]=1;break;}//we have a particle in this mode, and it can't be in any others
  }
};

inline void check_mode_ordered_NS(const int mode_num,const unsigned short int entry,const int measured_occupations, bool** occupations){
  for (int i=0;i<3;++i){
    if (mode_num==i){occupations[i+3][entry]=1;break;}//we have a particle in this mode, and it can't be in any others
  }
};

MPO_matrix MakeHamiltonian(const Vertex& modelvertex, const VertexParameterArray& couplingparams){
  //Lower triangular MPO
  // I   0   0
  // JK  0   0
  // HV  K'  I  //note the charges are such that K' and K differ in that Q[K[i]]+Q[K'[i]]=0
  QNCombinations differencecombinations(modelvertex.Spectrum,1); //1 means use difference

  MPXInt lineardim=modelvertex.Spectrum.size()*(2+couplingparams.size()*differencecombinations.size()); //the actual length of the sparse matrix needed
  MPXInt offset_to_last_block=modelvertex.Spectrum.size()*(1+couplingparams.size()*differencecombinations.size()); //offset to get to the last row of operators
  SparseMatrix M(lineardim,lineardim,lineardim);

  //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
  for (size_t i=0;i<modelvertex.Spectrum.size();++i){
    M.entry(i,i,1.0);
    M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
    M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
  }
  //now the more annoying pieces
  //for each coupling operator
  for (size_t c=0;c<couplingparams.size();++c){ //assume each of the coupling params refers to an operator
    Sparseint operator_col_offset=c*differencecombinations.size()+1; //+1 for identity matrix in first block
    Sparseint operator_row_offset=(couplingparams.size()-c-1)*differencecombinations.size()+1; //reversal, +1 for identity

    double operatorparam=couplingparams[c].Value;
    for (Sparseint col=0;col<modelvertex.Operators[c].MatrixElements.cols();++col){
      for (Sparseint p=modelvertex.Operators[c].MatrixElements.get_p(col);p<modelvertex.Operators[c].MatrixElements.get_p(col+1);++p){
	Sparseint i=modelvertex.Operators[c].MatrixElements.get_i(p);
	//
	//this could easily be a lookup function for QNCombinations
	State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	//now look through charge groups to find MPO indices
	Sparseint MPO_subcol=0;
	Sparseint MPO_subrow=0;
	if (diffstate==-diffstate){ //involution block
	  for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	    if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
	      MPO_subcol=l;
	      MPO_subrow=l; //involution pairs don't get reversed
	      break;
	    }
	  }
	}
	else {
	  for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	    if (diffstate==differencecombinations.OrderedPairs[l].PairState){
	      MPO_subcol=l+differencecombinations.InvolutionPairs.size();
	      MPO_subrow=differencecombinations.size()-l-1; //reversed
	      break;
	    }
	  }
	}
	//end of possible lookup function

	std::complex<double> x=modelvertex.Operators[c].MatrixElements.get_x(p);
	//lowest block row
	M.entry(offset_to_last_block+i,modelvertex.Spectrum.size()*(operator_col_offset+MPO_subcol)+col,x);
	//first block col
	M.entry(modelvertex.Spectrum.size()*(operator_row_offset+MPO_subrow)+i,col,x*operatorparam*modelvertex.Parameters[0].Value);
      }
    }
  }
  StateArray b;

  b.push_back(State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
  for (size_t c=0;c<couplingparams.size();++c){ //do for each operator
    for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
      b.push_back(differencecombinations.InvolutionPairs[l].PairState);
    }
    for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
      b.push_back(differencecombinations.OrderedPairs[l].PairState);
    }
  }
  b.push_back(State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    
  std::vector<MPXIndex> indices;
  indices.push_back(MPXIndex(1,modelvertex.Spectrum)); //sigma primed
  indices.push_back(MPXIndex(1,b)); //b_left
  indices.push_back(MPXIndex(0,modelvertex.Spectrum)); //sigma
  indices.push_back(MPXIndex(0,b)); //b_right

  return  MPO_matrix(modelvertex.Spectrum, indices, M.finalise());
}

//////////////////////////////////////////////////////////////
Vertex VertexGenerator(const VertexParameterArray& inputs)
{
  static const unsigned short int MAX_NUM_CHAIN_STATES=500;//useful buffer size for forming matrix elements for the occupation number as we generate the spectrum 
  //initialise vertex object
  Vertex ModelVertex(inputs);
  ModelVertex.ChargeRules.push_back(0); //a momentum like quantum number
  ModelVertex.ChargeRules.push_back(2); //Z_2 quantum number

  const QNVector& ChargeRules=ModelVertex.ChargeRules;

  const double R = inputs[0].Value;
  const double mass=inputs[1].Value;
  const double spectrum_cutoff=inputs[2].Value;

  const double Delta = abs(mass); //use the magnitude of delta
  const double DR=Delta*R;
  
  const double rvac = Rvacuum(Delta,R);
  const double nvac = NSvacuum(Delta,R);
  cout << "NS vacuum energy: " << nvac << ", R vacuum energy: " << rvac  << endl;    
  const double onechain_gs_energy=nvac;
  MPXInt num_chain_states = 0;
  
  //useful containers
  vector<vector<double> >* onechain_theta_ptr=new vector<vector<double> >();
  vector<vector<double> > &onechain_theta=*onechain_theta_ptr;
  vector<vector<double> >* onechain_kappa_ptr=new vector<vector<double> >();
  vector<vector<double> > &onechain_kappa=*onechain_kappa_ptr;  
  vector<int>* onechain_modenumber_ptr=new vector<int>();
  vector<int> &onechain_modenumber=*onechain_modenumber_ptr;
  
  //block for setting up occupation numbers
  int measured_occupations=6; 
  bool** occupations=new bool*[measured_occupations]; //we use double the mode integer, because of the half integer cases
  for (int i=0;i<measured_occupations;++i){
    occupations[i]=new bool[MAX_NUM_CHAIN_STATES];
    for (int j=0;j<MAX_NUM_CHAIN_STATES;++j){
      occupations[i][j]=0;
    }
  }

  vector<int> total_occupation;

  if (mass<0){
    const QuantumNumberInt num1R = 200;
    const QuantumNumberInt num2NS = 100;
    const QuantumNumberInt num3R = 50;
    const QuantumNumberInt num4NS = 50;
    const QuantumNumberInt num5R = 20;
    const QuantumNumberInt num6NS = 20;
    /*NS sector 0 particle mode, continuum chain ground state */    
    QuantumNumberInt gsarray[2] ={0,0};

    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(gsarray,gsarray+2),onechain_gs_energy));

    onechain_modenumber.push_back(0);
    onechain_kappa.push_back(vector<double>());
    onechain_kappa[0].push_back(0.0); //just dummy values for ground state

    onechain_theta.push_back(vector<double>());
    onechain_theta[0].push_back(0.0); //just dummy values for ground state

    cout << "NS0: " << num_chain_states << " * " << ModelVertex.Spectrum[0].en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[0][0]<< " " << onechain_theta[0][0] << endl;
    cout << "mode 0 particles: 1" << endl << endl;
    num_chain_states++;
    total_occupation.push_back(0);

    /*Ramond Sector -- one particle modes*/
    MPXInt mode1part_R =0;
    if (Delta +rvac <= spectrum_cutoff){
      for(QuantumNumberInt i=-num1R;i<=num1R;++i) {
	double en = Delta*sqrt(1.+pow(2.*M_PI*i/DR,2.)) +rvac;
	if (en < spectrum_cutoff) {       
	  mode1part_R++;
	  QuantumNumberInt R1values[2]={i,1};
	  ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R1values,R1values+2),en));
	  check_mode(2*i,num_chain_states,measured_occupations,occupations);
	  
	  onechain_modenumber.push_back(1);
	  onechain_theta.push_back(vector<double>());
	  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI/DR*i));
	  
	  onechain_kappa.push_back(vector<double>());
	  onechain_kappa[num_chain_states].push_back(Kap(asinh(2.*M_PI/DR*i),DR));
	  cout << "R1: " << num_chain_states << " " << i << " " << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << onechain_theta[num_chain_states][0] << endl;
	  num_chain_states++;
	  total_occupation.push_back(1);
	}
      }
    }
    cout << "mode 1 particles: " << mode1part_R << endl << endl;
    
    /*Neveu-Schwarz Sector -- two particle modes*/
    MPXInt mode2part_NS=0;
    if (2*Delta*(sqrt(1.+pow(2.*M_PI*(.5)/DR,2.))) +nvac <= spectrum_cutoff){
       
      for(QuantumNumberInt i=-num2NS-1;i<num2NS;++i) {
	if (Delta +rvac > spectrum_cutoff){break;}

	for(QuantumNumberInt j=i+1;j<=num2NS;++j) {
	  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+nvac;
	  if (en< spectrum_cutoff) {
	    mode2part_NS++;
	    QuantumNumberInt NS2values[2]={static_cast<QuantumNumberInt>(i+j+1),0};
	    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS2values,NS2values+2),en));
	    check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
	    check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
	     
	    onechain_modenumber.push_back(2);
	    onechain_theta.push_back(vector<double>());
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
	     
	    onechain_kappa.push_back(vector<double>());
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
	     
	    cout << "NS2: " << num_chain_states << " " << i << " " << j << " ";
	    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
	    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << endl;
	    num_chain_states++;
	    total_occupation.push_back(2);
	  }
	}
      }
    }
    cout << "mode 2 particles: " << mode2part_NS << endl << endl;

    /*Ramond Sector -- three particle modes*/
    MPXInt mode3part_R = 0;
    if (Delta+2.0*Delta*(sqrt(1.+pow(2.*M_PI/DR,2.))) +rvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num3R;i<num3R-1;++i) {
	for(QuantumNumberInt j=i+1;j<num3R;++j) {
	  for(QuantumNumberInt k=j+1;k<=num3R;++k) {
	    double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*j/DR,2.)) + sqrt(1.+pow(2.*M_PI*k/DR,2.))) + rvac;
	    if (en < spectrum_cutoff) {
	      mode3part_R++;
	      QuantumNumberInt R3values[2]={static_cast<QuantumNumberInt>(i+j+k),1};
	      ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R3values,R3values+2),en));
	      check_mode(2*i,num_chain_states,measured_occupations,occupations);
	      check_mode(2*j,num_chain_states,measured_occupations,occupations);
	      check_mode(2*k,num_chain_states,measured_occupations,occupations);
		  
	      onechain_modenumber.push_back(3);
	      onechain_theta.push_back(vector<double>());
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*i/DR));
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*j/DR));
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*k/DR));
		  
	      onechain_kappa.push_back(vector<double>());
	      onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
	      onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
	      onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		  
	      cout << "R3: " << num_chain_states << " " << i << " " << j << " " << k << " ";
	      cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
	      cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << endl;
	      num_chain_states++;	
	      total_occupation.push_back(3);
	    }
	  }
	}
      }
    }
    cout << "mode 3 particles: " << mode3part_R << endl << endl;

    /*Neveu-Schwarz Sector -- four particle modes*/
    MPXInt mode4part_NS = 0.;
    if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))) +nvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num4NS-1;i<num4NS-2;++i) {
	for(QuantumNumberInt j=i+1;j<num4NS-1;++j) {
	  for(QuantumNumberInt k=j+1;k<num4NS;++k) {
	    for(QuantumNumberInt l=k+1;l<=num4NS;++l) {
	      double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+nvac;
	      if (en < spectrum_cutoff) {
		mode4part_NS++;
		QuantumNumberInt NS4values[2]={static_cast<QuantumNumberInt>(i+j+k+l+2),0};	     
		ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS4values,NS4values+2),en));
		check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
		check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
		check_mode(2*k+1,num_chain_states,measured_occupations,occupations);
		check_mode(2*l+1,num_chain_states,measured_occupations,occupations);

		onechain_modenumber.push_back(4);	    
		onechain_theta.push_back(vector<double>());
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
	    
		onechain_kappa.push_back(vector<double>());
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));	    
		cout << "NS4: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " ";
		cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << endl;
		num_chain_states++;
		total_occupation.push_back(4);
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 4 particles: " << mode4part_NS << endl << endl;

    /*Ramond Sector -- five particle modes*/
    MPXInt mode5part_R = 0;
    if (Delta+2.0*Delta*(sqrt(1.+pow(2.*M_PI/DR,2.))+sqrt(1.+pow(2.*M_PI*2.0/DR,2.))) +rvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num5R;i<num5R-3;++i) {
	for(QuantumNumberInt j=i+1;j<num5R-2;++j) {
	  for(QuantumNumberInt k=j+1;k<num5R-1;++k) {
	    for(QuantumNumberInt l=k+1;l<num5R;++l) {
	      for(QuantumNumberInt m=k+1;m<=num5R;++m) {
		double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*(j)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(m)/DR,2.))) + rvac;//-onechain_gs_energy;
		if (en < spectrum_cutoff) {
		  mode5part_R++;
		  QuantumNumberInt R5values[2]={static_cast<QuantumNumberInt>(i+j+k+l+m),1};
		  ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R5values,R5values+2),en));
		  check_mode(2*i,num_chain_states,measured_occupations,occupations);
		  check_mode(2*j,num_chain_states,measured_occupations,occupations);
		  check_mode(2*k,num_chain_states,measured_occupations,occupations);
		  check_mode(2*l,num_chain_states,measured_occupations,occupations);
		  check_mode(2*m,num_chain_states,measured_occupations,occupations);

		  onechain_modenumber.push_back(5);
		  onechain_theta.push_back(vector<double>());
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*i/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*j/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*k/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*l/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*m/DR));
	      
		  onechain_kappa.push_back(vector<double>());
		  onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		  onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		  onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		  onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));
		  onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][4],DR));
	  
		  cout << "R5: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " ";
		  cout << en << " " << (2.*M_PI)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		  cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << " " << onechain_theta[num_chain_states][4] << endl;
		  num_chain_states++;
		  total_occupation.push_back(5);
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 5 particles: " << mode5part_R << endl << endl;

    /*Neveu-Schwarz Sector -- six particle modes*/
    MPXInt mode6part_NS = 0;
    if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))+sqrt(1.+pow(2.*M_PI*2.5/DR,2.))) +nvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num6NS-1;i<num6NS-4;++i) {
	for(QuantumNumberInt j=i+1;j<num6NS-3;++j) {
	  for(QuantumNumberInt k=j+1;k<num6NS-2;++k) {
	    for(QuantumNumberInt l=k+1;l<num6NS-1;++l) {
	      for(QuantumNumberInt m=l+1;m<num6NS;++m) {
		for(QuantumNumberInt n=m+1;n<=num6NS;++n) {
		  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n+0.5)/DR,2.))+nvac;
		  if (en < spectrum_cutoff) {
		    mode6part_NS++;
		    QuantumNumberInt NS6values[2]={static_cast<QuantumNumberInt>(i+j+k+l+m+n+3),0};
		    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS6values,NS6values+2),en));
		    check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
		    check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
		    check_mode(2*k+1,num_chain_states,measured_occupations,occupations);
		    check_mode(2*l+1,num_chain_states,measured_occupations,occupations);
		    check_mode(2*m+1,num_chain_states,measured_occupations,occupations);
		    check_mode(2*n+1,num_chain_states,measured_occupations,occupations);

		    onechain_modenumber.push_back(6);
		    onechain_theta.push_back(vector<double>());
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n+0.5)/DR));		
		    onechain_kappa.push_back(vector<double>());
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][4],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][5],DR));
 
		    cout << "NS6: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << " ";
		    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << " " << onechain_theta[num_chain_states][4] << " " << onechain_theta[num_chain_states][5] << endl;
		    num_chain_states++;
		    total_occupation.push_back(6);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 6 particles: " << mode6part_NS << endl;
  }
  else if (mass>0){ //ordered or critical chains...

    const QuantumNumberInt num2R = 100;
    const QuantumNumberInt num2NS = 100;
    const QuantumNumberInt num4R = 50;
    const QuantumNumberInt num4NS = 50;
    const QuantumNumberInt num6R = 20;
    const QuantumNumberInt num6NS = 20;
    
    //push back the two vacuua
    QuantumNumberInt nvacvalues[2]={0,0};		    
    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(nvacvalues,nvacvalues+2),nvac));

    onechain_modenumber.push_back(0);
    onechain_theta.push_back(vector<double>());
    onechain_kappa.push_back(vector<double>());
    onechain_kappa[0].push_back(0.0); //just dummy values
    onechain_theta[0].push_back(0.0); //just dummy values
    num_chain_states++;
    total_occupation.push_back(0);

    QuantumNumberInt rvacvalues[2] ={0,1}; //1 is ramond sector
    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(rvacvalues,rvacvalues+2),rvac));

    onechain_modenumber.push_back(0);
    onechain_theta.push_back(vector<double>());
    onechain_kappa.push_back(vector<double>());
    onechain_kappa[0].push_back(0.0); //just dummy values
    onechain_theta[0].push_back(0.0); //just dummy values
    num_chain_states++;
    total_occupation.push_back(0);

    std::cout << "mode 0 particles: " << total_occupation.size() << std::endl << std::endl;
  
    /*Now do 2 particle modes*/	    
    /* R sector*/
    MPXInt mode2part_R=0;
    if (Delta*(1.0+sqrt(1.+pow(2.*M_PI/DR,2.))) +rvac <= spectrum_cutoff){ //check lowest state
      for(QuantumNumberInt i=-num2R-1;i<num2R;++i) {
	if (Delta +rvac > spectrum_cutoff){break;}
	for(QuantumNumberInt j=i+1;j<=num2R;++j) {
	  double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.))+sqrt(1.+pow(2.*M_PI*j/DR,2.)))+nvac;
	  if (en< spectrum_cutoff) {
	    mode2part_R++;
	    QuantumNumberInt R2values[2]={static_cast<QuantumNumberInt>(i+j),1};
	    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R2values,R2values+2),en));

	    check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
	    check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
	     
	    onechain_modenumber.push_back(2);
	    onechain_theta.push_back(vector<double>());
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
	     
	    onechain_kappa.push_back(vector<double>());
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
	     
	    cout << "R2: " << num_chain_states << " " << i << " " << j << " ";
	    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
	    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << endl;
	    num_chain_states++;
	    total_occupation.push_back(2);
	  }
	}
      }
    }
    cout << "mode 2 R particles: " << mode2part_R << endl << endl;
    /* NS sector*/
    MPXInt mode2part_NS=0;
    if (2*Delta*(sqrt(1.+pow(2.*M_PI*(.5)/DR,2.))) +nvac <= spectrum_cutoff){
      for(QuantumNumberInt i=-num2NS-1;i<num2NS;++i) {
	if (Delta +rvac > spectrum_cutoff){break;}

	for(QuantumNumberInt j=i+1;j<=num2NS;++j) {
	  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+nvac;
	  if (en< spectrum_cutoff) {
	    mode2part_NS++;
	    QuantumNumberInt NS2values[2]={static_cast<QuantumNumberInt>(i+j+1),0};
	    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS2values,NS2values+2),en));

	    check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
	    check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
	     
	    onechain_modenumber.push_back(2);
	    onechain_theta.push_back(vector<double>());
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
	     
	    onechain_kappa.push_back(vector<double>());
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
	    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
	     
	    cout << "NS2: " << num_chain_states << " " << i << " " << j << " ";
	    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
	    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << endl;
	    num_chain_states++;
	    total_occupation.push_back(2);
	  }
	}
      }
    }
    cout << "mode 2 NS particles: " << mode2part_NS << endl << endl;

    /*four particle modes*/
    /*Ramond sector -- four particle modes*/
    MPXInt mode4part_R = 0.;
    if (Delta*(1.0+2.0*sqrt(1.+pow(2.*M_PI/DR,2.))+sqrt(1.+pow(2.*M_PI*2.0/DR,2.))) +rvac <= spectrum_cutoff){
      for(QuantumNumberInt i=-num4R-1;i<num4R-2;++i) {
	for(QuantumNumberInt j=i+1;j<num4R-1;++j) {
	  for(QuantumNumberInt k=j+1;k<num4R;++k) {
	    for(QuantumNumberInt l=k+1;l<=num4R;++l) {
	      double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+rvac;
	      if (en < spectrum_cutoff) {
		mode4part_R++;
		QuantumNumberInt R4values[2]={static_cast<QuantumNumberInt>(i+j+k+l),1};
		ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R4values,R4values+2),en));

		check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_R(k,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_R(l,num_chain_states,measured_occupations,occupations);

		onechain_modenumber.push_back(4);	    
		onechain_theta.push_back(vector<double>());
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l)/DR));
	    
		onechain_kappa.push_back(vector<double>());
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));	    
		cout << "R4: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " ";
		cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << endl;
		num_chain_states++;
		total_occupation.push_back(4);
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 4 R particles: " << mode4part_R << endl << endl;

    /*Neveu-Schwarz Sector -- four particle modes*/
    MPXInt mode4part_NS = 0.;
    if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))) +nvac <= spectrum_cutoff){
      for(QuantumNumberInt i=-num4NS-1;i<num4NS-2;++i) {
	for(QuantumNumberInt j=i+1;j<num4NS-1;++j) {
	  for(QuantumNumberInt k=j+1;k<num4NS;++k) {
	    for(QuantumNumberInt l=k+1;l<=num4NS;++l) {
	      double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+nvac;
	      if (en < spectrum_cutoff) {
		mode4part_NS++;
		QuantumNumberInt NS4values[2]={static_cast<QuantumNumberInt>(i+j+k+l+2),0};
		ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS4values,NS4values+2),en));

		check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_NS(k,num_chain_states,measured_occupations,occupations);
		check_mode_ordered_NS(l,num_chain_states,measured_occupations,occupations);

		onechain_modenumber.push_back(4);	    
		onechain_theta.push_back(vector<double>());
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
	    
		onechain_kappa.push_back(vector<double>());
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));	    
		cout << "NS4: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " ";
		cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << endl;
		num_chain_states++;
		total_occupation.push_back(4);
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 4 NS particles: " << mode4part_NS << endl << endl;

    /*6 particle modes */

    /*Ramond Sector -- six particle modes*/
    MPXInt mode6part_R = 0;
    if (Delta*(1.0+2.0*sqrt(1.+pow(2.*M_PI/DR,2.))+2.0*sqrt(1.+pow(2.*M_PI*2.0/DR,2.))+sqrt(1.+pow(2.*M_PI*3.0/DR,2.))) +rvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num6R-1;i<num6R-4;++i) {
	for(QuantumNumberInt j=i+1;j<num6R-3;++j) {
	  for(QuantumNumberInt k=j+1;k<num6R-2;++k) {
	    for(QuantumNumberInt l=k+1;l<num6R-1;++l) {
	      for(QuantumNumberInt m=l+1;m<num6R;++m) {
		for(QuantumNumberInt n=m+1;n<=num6R;++n) {
		  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n)/DR,2.))+nvac;
		  if (en < spectrum_cutoff) {
		    mode6part_R++;
		    QuantumNumberInt R6values[2]={static_cast<QuantumNumberInt>(i+j+k+l+m+n),1};
		    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(R6values,R6values+2),en));

		    check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_R(k,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_R(l,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_R(m,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_R(n,num_chain_states,measured_occupations,occupations);

		    onechain_modenumber.push_back(6);
		    onechain_theta.push_back(vector<double>());
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n)/DR));		
		    onechain_kappa.push_back(vector<double>());
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][4],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][5],DR));
		
		    cout << "R6: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << " ";
		    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << " " << onechain_theta[num_chain_states][4] << " " << onechain_theta[num_chain_states][5] << endl;
		    num_chain_states++;
		    total_occupation.push_back(6);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 6 R particles: " << mode6part_R << endl;

    /*Neveu-Schwarz Sector -- six particle modes*/
    MPXInt mode6part_NS = 0;
    if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))+sqrt(1.+pow(2.*M_PI*2.5/DR,2.))) +nvac <= spectrum_cutoff){

      for(QuantumNumberInt i=-num6NS-1;i<num6NS-4;++i) {
	for(QuantumNumberInt j=i+1;j<num6NS-3;++j) {
	  for(QuantumNumberInt k=j+1;k<num6NS-2;++k) {
	    for(QuantumNumberInt l=k+1;l<num6NS-1;++l) {
	      for(QuantumNumberInt m=l+1;m<num6NS;++m) {
		for(QuantumNumberInt n=m+1;n<=num6NS;++n) {
		  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n+0.5)/DR,2.))+nvac;
		  if (en < spectrum_cutoff) {
		    mode6part_NS++;
		    QuantumNumberInt NS6values[2]={static_cast<QuantumNumberInt>(i+j+k+l+m+n+3),0};
		    ModelVertex.Spectrum.push_back(EigenState(ChargeRules,QNVector(NS6values,NS6values+2),en));

		    check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_NS(k,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_NS(l,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_NS(m,num_chain_states,measured_occupations,occupations);
		    check_mode_ordered_NS(n,num_chain_states,measured_occupations,occupations);

		    onechain_modenumber.push_back(6);
		    onechain_theta.push_back(vector<double>());
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m+0.5)/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n+0.5)/DR));		
		    onechain_kappa.push_back(vector<double>());
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][0],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][1],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][2],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][3],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][4],DR));
		    onechain_kappa[num_chain_states].push_back(Kap(onechain_theta[num_chain_states][5],DR));
		
		    cout << "NS6: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << " ";
		    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " ";
		    cout << onechain_theta[num_chain_states][0] << " " << onechain_theta[num_chain_states][1] << " " << onechain_theta[num_chain_states][2] << " " << onechain_theta[num_chain_states][3] << " " << onechain_theta[num_chain_states][4] << " " << onechain_theta[num_chain_states][5] << endl;
		    num_chain_states++;
		    total_occupation.push_back(6);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "mode 6 NS particles: " << mode6part_NS << endl;
  }
  else {
    cout << "Critical chains not supported yet" << endl; exit(1);
  }

  //copy spectrum into pairs, find key listing
  std::vector<std::pair<size_t, double> > sort_keys;

  for (auto&& ES : ModelVertex.Spectrum.Energies()){
    sort_keys.push_back(std::pair<size_t, double>(sort_keys.size(),ES));
  }

  std::sort(sort_keys.begin(),sort_keys.end(),[] (std::pair<size_t, double> a, std::pair<size_t, double> b) {
      return b.second > a.second;
    });

  EigenStateArray TempSpectrum;
  for (size_t s=0;s<ModelVertex.Spectrum.size();++s){
    TempSpectrum.push_back(ModelVertex.Spectrum[sort_keys[s].first]); 
  }
  if (ModelVertex.Spectrum.size()!=static_cast<size_t>(num_chain_states)){
    cout << "Mismatch error" << endl; exit(1);
  }

    cout << "End generating chain spectrum, generated " << num_chain_states << " states" << endl;
    cout << "Starting Matrix Elements" << endl;
#if defined NDEBUG
    double sr = S(DR); //really slow right now
    cout << "S(R): " << sr << endl;
#else
    double sr=1.0;
    cout << "WARNING APPROXIMATING S(R) AS " << sr <<  endl;
#endif
    ModelVertex.Operators.push_back(VertexOperator("Spin Operator",ModelVertex.Spectrum.size()));      
    //off diagonal only
    for (size_t col=0;col<ModelVertex.Spectrum.size();++col){
      size_t s_col=sort_keys[col].first;
      for (size_t row=col+1;row<ModelVertex.Spectrum.size();++row){
	size_t s_row=sort_keys[row].first;
	if(TempSpectrum[row][1]!=TempSpectrum[col][1]){ //off diagonal in Ising 'sector' charge.
	  complex<double> M_E_=sr*calculate_chainmatrixelement(TempSpectrum,static_cast<MPXInt>(s_row),static_cast<MPXInt>(s_col),static_cast<MPXInt>(onechain_modenumber[s_row]),static_cast<MPXInt>(onechain_modenumber[s_col]),onechain_kappa[s_row],onechain_kappa[s_col],onechain_theta[s_row],onechain_theta[s_col],Delta,DR);
	  ModelVertex.Operators.back().MatrixElements.entry(row,col,M_E_);
	  ModelVertex.Operators.back().MatrixElements.entry(col,row,conj(M_E_));
	}
      }
    }
    ModelVertex.Operators.back().MatrixElements.finalise();

    //make total number operator
    ModelVertex.Operators.push_back(VertexOperator("Total Number",ModelVertex.Spectrum.size()));
    for (size_t j=0;j<ModelVertex.Spectrum.size();++j){
      ModelVertex.Operators.back().MatrixElements.entry(j,j,static_cast<complex<double> >(total_occupation[sort_keys[j].first]));
    }
    ModelVertex.Operators.back().MatrixElements.finalise();

    //make number operators
    for (MPXInt i=0;i<measured_occupations;++i){
      ModelVertex.Operators.push_back(VertexOperator("Mode Occupation",ModelVertex.Spectrum.size()));
      SparseMatrix occ_number(ModelVertex.Spectrum.size(),ModelVertex.Spectrum.size(),ModelVertex.Spectrum.size());
      for (size_t j=0;j<ModelVertex.Spectrum.size();++j){
	if(occupations[i][sort_keys[j].first]==1){ModelVertex.Operators.back().MatrixElements.entry(j,j,1.0);}
      }
      ModelVertex.Operators.back().MatrixElements.finalise();
    }

    ModelVertex.Spectrum=std::move(TempSpectrum);

    //clean up workspace
    delete onechain_theta_ptr;
    delete onechain_kappa_ptr;
    delete onechain_modenumber_ptr;
    for (MPXInt i=0;i<measured_occupations;++i){
      delete[] occupations[i];
    }
    delete[] occupations;

  //Operators[i]
  //0 spin
  //1 total chain fermion occupation number
  //2+ individual chain mode occupations

  return ModelVertex;
}


//user functions for making the spectrum of the chain
inline double NSvacuum(const double Delta, const double R)
{
  double term, den, th,DR;
  int i, num;
  DR=Delta*R;
  term = 0.;
  num = 100000.;
  den = 5./num;
  for(i=0;i<num;++i) {
    th = i*den;
    term += cosh(th)*log(1.+exp(-1.*DR*cosh(th)));
  }
  return(-1.*term*den*Delta/(M_PI));
}

inline double Rvacuum(const double Delta, const double R)
{
  double term, den, th,DR;
  int i, num;
  DR=Delta*R;
  term = 0.;
  num = 100000.;
  den = 5./num;
  for(i=0;i<num;++i) { //even function so just do positive range and compensate with factor of 2
    th = i*den;
    term += cosh(th)*log(1.-exp(-1.*DR*cosh(th)));
  }
  return(-1.*term*den*Delta/(M_PI));
}

inline double Kap(const double t, const double DR)
{
  double term, den, th;
  int i, num;
  term = 0.;
  num = 10000.;
  den = 10./num;
  for(i=-num;i<num;++i) {
    th = i*den;
    term += (1./cosh(t-th)*log((1.-exp(-1.*DR*cosh(th)))/(1.+exp(-1.*DR*cosh(th)))));
  }
  return(term*den/(2.0*M_PI));
}

inline double S(const double DR)
{
  double term, den, t1, t2, err,jterm;
  double ct1,st1,ct2,st2,ft12,sdc1;
  int i, j, num;

  err = .0000001;
  term = 0.;
  num = 6000.;
  den = 4./num;
    
  for(i=-num;i<num;++i) {
    t1 = i*den;
    ct1=cosh(t1);
    st1=sinh(t1);
    sdc1=st1/sinh(DR*ct1);
    jterm=0.;
    for(j=-num;j<num;++j) {
      t2 = j*den;
      ct2=cosh(t2);
      st2=sinh(t2);
      ft12=abs(t1-t2);
      jterm -= (st2/sinh(DR*ct2))*log(tanh(err+.5*ft12));
    }
    term +=sdc1*jterm;
  }
  return(exp(pow(DR*den/M_PI,2.)*.5*.25*term));
}

inline complex<double> calculate_chainmatrixelement(const EigenStateArray& onechain_states, const MPXInt a, const MPXInt b, const MPXInt modenum_a, const MPXInt modenum_b, const vector<double> kappa_a, const vector<double> kappa_b, const vector<double> theta_a, const vector<double> theta_b, const double Delta,const double DR)
{
  double term = 1.0;
  double sign;
  (onechain_states[a][1] == 1) ? sign=1.0 : sign=-1.0;
  for(MPXInt i=0;i<modenum_a;++i) {
    term *= (exp(sign*kappa_a[i])/sqrt(DR*cosh(theta_a[i])));
  }
  (onechain_states[b][1] == 1) ? sign=1.0 : sign=-1.0;
  for(MPXInt i=0;i<modenum_b;++i) {
    term *= (exp(sign*kappa_b[i])/sqrt(DR*cosh(theta_b[i])));
  }
  if (modenum_a > 1) {
    for(MPXInt i=0;i<modenum_a;++i) {
      for(MPXInt j=i+1;j<modenum_a;++j) {
	term *= (tanh((theta_a[i]-theta_a[j])*.5));
      }
    }
  }
  if (modenum_b > 1) {
    for(MPXInt i=0;i<modenum_b;++i) {
      for(MPXInt j=i+1;j<modenum_b;++j) {
	term *= (tanh((theta_b[i]-theta_b[j])*.5));
      }
    }
  }
  if ((modenum_a > 0) && (modenum_b > 0)) {
    for(MPXInt i=0;i<modenum_a;++i) {
      for(MPXInt j=0;j<modenum_b;++j) {
	term /= tanh((theta_a[i]-theta_b[j])*.5);
      }
    }
  }
  return(complex<double>(term*1.35783834*pow(Delta,.125),0.0));
}
#endif

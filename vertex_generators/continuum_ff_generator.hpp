#ifndef THEORY_FREE_FERMION_H
#define THEORY_FREE_FERMION_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>
#include <iostream>
#include <sstream>

#include "../ajaj_common.hpp"
#include "../sparse_interface.hpp"

namespace continuumff{

  unsigned long int makebit(ajaj::QuantumNumberInt i){

    return i == 0 ? 1 : 1 << int(2*abs(i)-1 +(i>0));
  }

  unsigned int Hamming(unsigned long int b){
    unsigned int count ;
    for (count=0; b; count++)
      b&=b-1;
    return count;
  }

  int signfactor(unsigned long int b /* state bit string */, int m /* the mode we have*/){
    b >>= (m==0)+(2*abs(m) +(m>0)); //shift away bits for lower mode numbers
    return (Hamming(b) % 2) ? -1 : 1;
  }  

  inline void check_mode(const int i,const unsigned short int entry,const int measured_occupations, bool** occupations){
    occupations[measured_occupations+i-1][entry]=1;
    //for (int i=-measured_occupations+1;i<measured_occupations;++i){
    //if (mode_num==i){occupations[measured_occupations+i-1][entry]=1;break;}//we have a particle in this mode, and it can't be in any others
    //}
  };

  ajaj::MPO_matrix MakeHamiltonian(const ajaj::Vertex& modelvertex, const ajaj::VertexParameterArray& couplingparams){
    //Lower triangular MPO
    // I           0    0         0
    // JPsidagger  0    0         0
    // JPsi
    // HV          Psi' Psidagger I        //note the charges Q, are such that Q[Psidagger[i]]+Q[Psi'[i]]=0

    ajaj::QNCombinations differencecombinations(modelvertex.Spectrum,1); //1 means use difference

    //Unlike the Ising case (where the coupling operator is spin, which is Hermitian) the coupling here is a tunneling: psi^dagger psi, where psi is not Hermitian.
    //Therefore there is an extra factor of 2, as A_n and A_n^dagger are different
    //Additionally there is a different tunnelling for each fermionic mode, n, i.e. a tunnelling term looks like A_{n,i}^dagger A_{n,i+1} where i is the vertex index

    ajaj::MPXInt number_of_measured_modes=modelvertex.Operators.size()-2; //first operator is the total occupation and last is vertex hamiltonian.
    //This is a flaky step, and should be checked first if there are any errors
    //would be better to count up the number of 'weighted annihilation operators'

    ajaj::MPXInt lineardim=modelvertex.Spectrum.size()*(2+2*number_of_measured_modes*differencecombinations.size()); //the actual length of the sparse matrix needed
    ajaj::MPXInt offset_to_last_block=modelvertex.Spectrum.size()*(1+2*number_of_measured_modes*differencecombinations.size()); //offset to get to the last row of operators
    ajaj::SparseMatrix M(lineardim,lineardim,lineardim);

    //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
    for (size_t i=0;i<modelvertex.Spectrum.size();++i){
      M.entry(i,i,1.0);
      M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
      M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }

    //now the more annoying pieces

    double tunnelling=couplingparams[0].Value;
    double Delta=modelvertex.Parameters[1].Value;

    ajaj::Sparseint A_col_offset=1; //+1 for identity matrix in first block
    ajaj::Sparseint A_row_offset=differencecombinations.size()*number_of_measured_modes+1; //reversed order
    ajaj::Sparseint Adagger_col_offset=differencecombinations.size()*number_of_measured_modes+1; //+1 for identity matrix in first block
    ajaj::Sparseint Adagger_row_offset=1; //reversed order
    for (size_t m=0;m<number_of_measured_modes;++m){
      for (ajaj::Sparseint col=0;col<modelvertex.Operators[1+m].MatrixElements.cols();++col){
	for (ajaj::Sparseint p=modelvertex.Operators[1+m].MatrixElements.get_p(col);p<modelvertex.Operators[1+m].MatrixElements.get_p(col+1);++p){
	  ajaj::Sparseint i=modelvertex.Operators[1+m].MatrixElements.get_i(p);
	  ajaj::State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	  //now look through charge groups to find MPO indices
	  ajaj::Sparseint MPO_subcol=0;
	  ajaj::Sparseint MPO_subrow=0;
	  ajaj::Sparseint MPO_subcolAdagger=0;
	  ajaj::Sparseint MPO_subrowAdagger=0;
	  if (diffstate==-diffstate){ //involution block
	    for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	      if (diffstate==differencecombinations.InvolutionPairs[l].PairState){
		MPO_subcol=l;
		MPO_subrow=l; //involution pairs don't get reversed
		MPO_subcolAdagger=l;
		MPO_subrowAdagger=l;
		break;
	      }
	    }
	  }
	  else {
	    for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	      if (diffstate==differencecombinations.OrderedPairs[l].PairState){
		MPO_subcol=differencecombinations.InvolutionPairs.size()+l;
		MPO_subrow=differencecombinations.size()-l-1; //reversed
		MPO_subcolAdagger=MPO_subrow;
		MPO_subrowAdagger=MPO_subcol;
		break;
	      }
	    }
	  }
	  std::complex<double> x=modelvertex.Operators[1+m].MatrixElements.get_x(p);
	  //do A first
	  //lowest block row
	  M.entry(offset_to_last_block+i,modelvertex.Spectrum.size()*(A_col_offset+differencecombinations.size()*m+MPO_subcol)+col,1.0*x);
	  //first block col
	  M.entry(modelvertex.Spectrum.size()*(A_row_offset+differencecombinations.size()*m+MPO_subrow)+i,col,x*tunnelling); //no factor of R
	  //now A dagger so swap i and col
	  //lowest block row
	  M.entry(offset_to_last_block+col,modelvertex.Spectrum.size()*(Adagger_col_offset+differencecombinations.size()*m+MPO_subcolAdagger)+i,conj(x));
	  //first block col
	  M.entry(modelvertex.Spectrum.size()*(Adagger_row_offset+differencecombinations.size()*m+MPO_subrowAdagger)+col,i,conj(x)*tunnelling); //no factor of R
	}
      }
    }

    ajaj::StateArray b;

    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    for (size_t c=0;c<2;++c){ //do for each operator
      for (size_t m=0;m<number_of_measured_modes;++m){ //for each mode
	for (size_t l=0;l<differencecombinations.InvolutionPairs.size();++l){
	  b.push_back(differencecombinations.InvolutionPairs[l].PairState);
	}
	for (size_t l=0;l<differencecombinations.OrderedPairs.size();++l){
	  b.push_back(differencecombinations.OrderedPairs[l].PairState);
	}
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

  ajaj::Vertex VertexGenerator(const ajaj::VertexParameterArray& inputs)
  {
    static const unsigned short int MAX_NUM_CHAIN_STATES=500;//useful buffer size for forming matrix elements for the occupation number as we generate the spectrum 
    //initialise vertex object
    ajaj::Vertex ModelVertex(inputs);
    ModelVertex.ChargeRules.push_back(0); //just a momentum quantum number
    ModelVertex.ChargeRules.push_back(0); //number

    const ajaj::QNVector& ChargeRules=ModelVertex.ChargeRules; //useful reference, saves typing ModelVertex

    const double R = inputs[0].Value;
    const double Delta=inputs[1].Value;
    const double spectrum_cutoff=inputs[2].Value;

    if (Delta<=0.0) {
      std::cout << "Negative mass parameter! Aborting. " <<std::endl; exit(1);
    }
    const double DR=Delta*R;
    const double vac = 0.0;
    const ajaj::MPXInt num1R = 100;
    const ajaj::MPXInt num2NS = 100;
    const ajaj::MPXInt num3R = 100;
    const ajaj::MPXInt num4NS = 10;
    const ajaj::MPXInt num5R = 2;
    const ajaj::MPXInt num6NS = 2;
    
    const double onechain_gs_energy=vac;
    ajaj::MPXInt num_chain_states = 0;

    //useful containers
    //std::vector<int>* onechain_modenumber_ptr=new std::vector<int>();
    //std::vector<int> &onechain_modenumber=*onechain_modenumber_ptr;
    //block for setting up occupation numbers
    //need to guess largest k
    int measured_occupations=0;
    for (int i=0; i<100;++i){
      if (sqrt(Delta*Delta+(2.0*2.0*M_PI*M_PI*i*i/R/R))+vac<spectrum_cutoff){
	measured_occupations=i;
      }
      else {break;}
    }
    ++measured_occupations; //add on 1 to make array sizes right

    cout << "Number of measured occupations including zero mode: " << 2*measured_occupations-1 << endl;

    bool** occupations=new bool*[2*measured_occupations-1];
    for (int i=0;i<2*measured_occupations-1;++i){
      occupations[i]=new bool[MAX_NUM_CHAIN_STATES];
      for (int j=0;j<MAX_NUM_CHAIN_STATES;++j){
	occupations[i][j]=0;
      }
    }

    std::vector<int> total_occupation;
    std::vector<unsigned long int> occ_bitstring;

    /*NS sector 0 particle mode, continuum chain ground state */
    ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(0),0}}),onechain_gs_energy));

    //onechain_modenumber.push_back(0);
    total_occupation.push_back(0);
    occ_bitstring.push_back(0);
    cout << "NS0: " << num_chain_states << " * " << ModelVertex.Spectrum[0].en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[0][0] <<  " " << occ_bitstring.back() << endl;
    cout << "vacuum states: 1" << endl << endl;
    num_chain_states++;

    /*Ramond Sector -- one particle modes*/

    ajaj::MPXInt mode1part_R =0;
    if (Delta +vac <= spectrum_cutoff){
      for(ajaj::MPXInt i=-num1R;i<=num1R;++i) {
	double en = Delta*sqrt(1.+pow(2.*M_PI*i/DR,2.)) +vac;
	if (en < spectrum_cutoff) {       
	  mode1part_R++;
	  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i),1}}),en));
	  check_mode(i,num_chain_states,measured_occupations,occupations);
  
	  total_occupation.push_back(1);
	  occ_bitstring.push_back(makebit(i));
	  cout << "R1: " << num_chain_states << " " << i << " " << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back() << endl;
	  num_chain_states++;
	      
	}
      }
    }
    cout << "1 particle states: " << mode1part_R << endl << endl;
    
    /*Neveu-Schwarz Sector -- two particle modes*/
    ajaj::MPXInt mode2part_NS=0;
    if (2.0*Delta +vac <= spectrum_cutoff){
       
      for(ajaj::MPXInt i=-num2NS;i<num2NS;++i) {

	for(ajaj::MPXInt j=i+1;j<=num2NS;++j) {
	  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+vac;
	  if (en< spectrum_cutoff) {
	    mode2part_NS++;
	    ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i+j),2}}),en));
	    check_mode(i,num_chain_states,measured_occupations,occupations);
	    check_mode(j,num_chain_states,measured_occupations,occupations);
	     
	    //onechain_modenumber.push_back(2);

	    total_occupation.push_back(2);
	    occ_bitstring.push_back(makebit(i)+makebit(j)); 
	    cout << "NS2: " << num_chain_states << " " << i << " " << j << " ";
	    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back()<< endl;
	    num_chain_states++;

	  }
	}
      }
    }
    cout << "2 particle states: " << mode2part_NS << endl << endl;

    /*Ramond Sector -- three particle modes*/
    ajaj::MPXInt mode3part_R = 0;
    if (3.0*Delta+vac <= spectrum_cutoff){

      for(ajaj::MPXInt i=-num3R;i<num3R-1;++i) {
	for(ajaj::MPXInt j=i+1;j<num3R;++j) {
	  for(ajaj::MPXInt k=j+1;k<=num3R;++k) {
	    double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*j/DR,2.)) + sqrt(1.+pow(2.*M_PI*k/DR,2.))) + vac;
	    if (en < spectrum_cutoff) {
	      mode3part_R++;
	      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i+j+k),3}}),en));
	      check_mode(i,num_chain_states,measured_occupations,occupations);
	      check_mode(j,num_chain_states,measured_occupations,occupations);
	      check_mode(k,num_chain_states,measured_occupations,occupations);
		  
	      //onechain_modenumber.push_back(3);
	      total_occupation.push_back(3);
	      occ_bitstring.push_back(makebit(i)+makebit(j)+makebit(k));  
	      cout << "R3: " << num_chain_states << " " << i << " " << j << " " << k << " ";
	      cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back()<<endl;
	      num_chain_states++;	
	        
	    }
	  }
	}
      }
    }
    cout << "3 particle states: " << mode3part_R << endl << endl;

    /*Neveu-Schwarz Sector -- four particle modes*/
    ajaj::MPXInt mode4part_NS = 0.;
    if (4.0*Delta +vac <= spectrum_cutoff){

      for(ajaj::MPXInt i=-num4NS-1;i<num4NS-2;++i) {
	for(ajaj::MPXInt j=i+1;j<num4NS-1;++j) {
	  for(ajaj::MPXInt k=j+1;k<num4NS;++k) {
	    for(ajaj::MPXInt l=k+1;l<=num4NS;++l) {
	      double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+vac;
	      if (en < spectrum_cutoff) {
		mode4part_NS++;
		ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i+j+k+l),4}}),en));
		check_mode(i,num_chain_states,measured_occupations,occupations);
		check_mode(j,num_chain_states,measured_occupations,occupations);
		check_mode(k,num_chain_states,measured_occupations,occupations);
		check_mode(l,num_chain_states,measured_occupations,occupations);

		//onechain_modenumber.push_back(4);	    
	    
		total_occupation.push_back(4);
		occ_bitstring.push_back(makebit(i)+makebit(j)+makebit(k)+makebit(l));    
		cout << "NS4: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " ";
		cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back()<< endl;
		num_chain_states++;
		
	      
	      }
	    }
	  }
	}
      }
    }
    cout << "4 particle states: " << mode4part_NS << endl << endl;

    /*Ramond Sector -- five particle modes*/
    ajaj::MPXInt mode5part_R = 0;
    if (5.0*Delta +vac <= spectrum_cutoff){

      for(ajaj::MPXInt i=-num5R;i<num5R-3;++i) {
	for(ajaj::MPXInt j=i+1;j<num5R-2;++j) {
	  for(ajaj::MPXInt k=j+1;k<num5R-1;++k) {
	    for(ajaj::MPXInt l=k+1;l<num5R;++l) {
	      for(ajaj::MPXInt m=k+1;m<=num5R;++m) {
		double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*(j)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(m)/DR,2.))) + vac;//-onechain_gs_energy;
		if (en < spectrum_cutoff) {
		  mode5part_R++;
		  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i+j+k+l+m),5}}),en));
		  check_mode(i,num_chain_states,measured_occupations,occupations);
		  check_mode(j,num_chain_states,measured_occupations,occupations);
		  check_mode(k,num_chain_states,measured_occupations,occupations);
		  check_mode(l,num_chain_states,measured_occupations,occupations);
		  check_mode(m,num_chain_states,measured_occupations,occupations);

		  //onechain_modenumber.push_back(5);
		  total_occupation.push_back(5);
		  occ_bitstring.push_back(makebit(i)+makebit(j)+makebit(k)+makebit(l)+makebit(m)); 
		  cout << "R5: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " ";
		  cout << en << " " << (2.*M_PI)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back()<< endl;
		  num_chain_states++;
		     

		}
	      }
	    }
	  }
	}
      }
    }
    cout << "5 particle states: " << mode5part_R << endl << endl;

    /*Neveu-Schwarz Sector -- six particle modes*/
    ajaj::MPXInt mode6part_NS = 0;
    if (6.0*Delta +vac <= spectrum_cutoff){

      for(ajaj::MPXInt i=-num6NS-1;i<num6NS-4;++i) {
	for(ajaj::MPXInt j=i+1;j<num6NS-3;++j) {
	  for(ajaj::MPXInt k=j+1;k<num6NS-2;++k) {
	    for(ajaj::MPXInt l=k+1;l<num6NS-1;++l) {
	      for(ajaj::MPXInt m=l+1;m<num6NS;++m) {
		for(ajaj::MPXInt n=m+1;n<=num6NS;++n) {
		  double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n)/DR,2.))+vac;
		  if (en < spectrum_cutoff) {
		    mode6part_NS++;
		    ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector({{ajaj::QuantumNumberInt(i+j+k+l+m+n),6}}),en));
		    check_mode(i,num_chain_states,measured_occupations,occupations);
		    check_mode(j,num_chain_states,measured_occupations,occupations);
		    check_mode(k,num_chain_states,measured_occupations,occupations);
		    check_mode(l,num_chain_states,measured_occupations,occupations);
		    check_mode(m,num_chain_states,measured_occupations,occupations);
		    check_mode(n,num_chain_states,measured_occupations,occupations);

		    //onechain_modenumber.push_back(6);

		    total_occupation.push_back(6);
		    occ_bitstring.push_back(makebit(i)+makebit(j)+makebit(k)+makebit(l)+makebit(m)+makebit(n)); 

		    cout << "NS6: " << num_chain_states << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << " ";
		    cout << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << occ_bitstring.back() << endl;
		    num_chain_states++;
		       

		  }
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "6 particle states: " << mode6part_NS << endl;

    if (ModelVertex.Spectrum.size()!=static_cast<size_t>(num_chain_states)){
      cout << "Mismatch error" << endl; exit(1);
    }

    cout << "End generating spectrum, generated " << num_chain_states << " states" << endl;

    /*for (int m=0;m<2*measured_occupations-1;++m){
      for (int i=0;i<ModelVertex.Spectrum.size();++i){
	std::cout << occupations[m][i] << " ";
      }
      std::cout << std::endl;
    }*/

    cout << "Starting Matrix Elements" << endl;

    ModelVertex.Operators.push_back(ajaj::VertexOperator("Total_Number",ModelVertex.Spectrum.size()));
    for (int j=0;j<ModelVertex.Spectrum.size();++j){
      ModelVertex.Operators.back().MatrixElements.entry(j,j,static_cast<complex<double> >(total_occupation[j]));
    }
    ModelVertex.Operators.back().MatrixElements.finalise();

    // add some more matrices
    for (int m=0;m<2*measured_occupations-1;++m){
      std::stringstream name;
      name << "Weighted_Annihilation_Operator_" << m+1;
      ModelVertex.Operators.push_back(ajaj::VertexOperator(name.str(),ModelVertex.Spectrum.size())); //weighted by sqrt(mass/energy)
    }
  
    //populate annihilation operators
    for (int j=0;j<ModelVertex.Spectrum.size();++j){
      for (int i=0;i<ModelVertex.Spectrum.size();++i){//inefficient
	/*if (total_occupation[i]==total_occupation[j]-1){ //bra has one less particle than ket
	  for (int m=0;m<2*measured_occupations-1;++m){
	    int mode=-measured_occupations+1+m;
	    if (occupations[m][j] == 1 && occupations[m][i]==0){ //search through modes to find which one is missing a particle
	      //if (ModelVertex.Spectrum[j].en<=0.0) {cout << "Incorrect occupation operator" << endl; exit(1);}
	      ModelVertex.Operators[m+1].MatrixElements.entry(i,j,signfactor(occ_bitstring[j],mode)*sqrt(Delta/ModelVertex.Spectrum[j].en)); //note the offset
	      break; //if this was a match, we can skip to next
	    }
	  }
	  }*/

	if (total_occupation[i]==total_occupation[j]-1){ //bra has one less particle than ket
	  unsigned long int where=occ_bitstring[i]^occ_bitstring[j];
	  if (Hamming(where)==1){ //is this the right one?	  
	    //convert where back to mode
	    int mode=0;
	    while (where!=1) {
	      if (mode>=0) mode++;
	      mode*=-1;
	      where >>= 1;
	    }
	    ModelVertex.Operators[mode+measured_occupations].MatrixElements.entry(i,j,signfactor(occ_bitstring[j],mode)/sqrt(sqrt(1.+4.0*M_PI*M_PI*(mode*mode)/DR/DR))); //note the offset	    
	  }
	}
      }
    }
  
    //finish off
    for (int m=0;m<2*measured_occupations-1;++m){
      ModelVertex.Operators[m+1].MatrixElements.finalise();
      //ModelVertex.Operators[m+1].MatrixElements.print();
    }

    //clean up workspace
    //delete onechain_modenumber_ptr;
    for (int i=0;i<2*measured_occupations-1;++i){
      delete[] occupations[i];
    }
    delete[] occupations;

    //push back Vertex Hamiltonain
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Vertex_Hamiltonian"));
    ModelVertex.Operators.back().MatrixElements=ajaj::SparseMatrix(ModelVertex.basis().Energies());

    //Operators[i]
    //0 total chain fermion occupation number
    //1+ annihilation operators
    //last is vertex hamiltonian

    return ModelVertex;

  };
}
#endif

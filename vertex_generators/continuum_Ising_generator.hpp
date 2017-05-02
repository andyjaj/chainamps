#ifndef THEORY_ISING_H
#define THEORY_ISING_H

#include <cmath>
#include <cstdlib>
#include <complex>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>

#include "../ajaj_common.hpp"
#include "../sparse_interface.hpp"

namespace continuumIsing {

  //static const unsigned short int MAX_NUM_CHAIN_STATES=500;//useful buffer size for forming matrix elements for the occupation number as we generate the spectrum 
  inline double NSvacuum(const double Delta, const double R);
  inline double Rvacuum(const double Delta, const double R);
  inline double Kap(const double t, const double DR);
  inline double S(const double DR);
  complex<double> calculate_chainmatrixelement(const ajaj::EigenStateArray& m_states, const ajaj::MPXInt a, const ajaj::MPXInt b, const ajaj::MPXInt modenum_a, const ajaj::MPXInt modenum_b, const std::vector<double> kappa_a, const std::vector<double> kappa_b, const std::vector<double> theta_a, const std::vector<double> theta_b, const double Delta,const double DR);

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
      term +=log((1.-exp(-1.*DR*cosh(th)))/(1.+exp(-1.*DR*cosh(th))))/cosh(t-th);
    }
    return(term*den/(2.0*M_PI));
  }

  inline double S(const double DR)
  {
    double term, den, t1, t2, err,jterm;
    double ct1,st1,drct2,st2,hft12,sdc1;
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
	drct2=DR*cosh(t2);
	st2=sinh(t2);
	hft12=0.5*abs(t1-t2)+err;
	jterm -= (st2)*log(tanh(hft12))/sinh(drct2);
      }
      term +=sdc1*jterm;
    }
    return(exp(pow(DR*den/M_PI,2.)*.5*.25*term));
  }

  inline complex<double> calculate_chainmatrixelement(const ajaj::EigenStateArray& onechain_states, const ajaj::MPXInt a, const ajaj::MPXInt b, const ajaj::MPXInt modenum_a, const ajaj::MPXInt modenum_b, const std::vector<double> kappa_a, const std::vector<double> kappa_b, const std::vector<double> theta_a, const std::vector<double> theta_b, const double Delta,const double DR)
  {
    double term = 1.0;
    double sign;
    (onechain_states[a][1] == 1) ? sign=1.0 : sign=-1.0;
    for(ajaj::MPXInt i=0;i<modenum_a;++i) {
      term *= (exp(sign*kappa_a[i])/sqrt(DR*cosh(theta_a[i])));
    }
    (onechain_states[b][1] == 1) ? sign=1.0 : sign=-1.0;
    for(ajaj::MPXInt i=0;i<modenum_b;++i) {
      term *= (exp(sign*kappa_b[i])/sqrt(DR*cosh(theta_b[i])));
    }
    if (modenum_a > 1) {
      for(ajaj::MPXInt i=0;i<modenum_a;++i) {
	for(ajaj::MPXInt j=i+1;j<modenum_a;++j) {
	  term *= (tanh((theta_a[i]-theta_a[j])*.5));
	}
      }
    }
    if (modenum_b > 1) {
      for(ajaj::MPXInt i=0;i<modenum_b;++i) {
	for(ajaj::MPXInt j=i+1;j<modenum_b;++j) {
	  term *= (tanh((theta_b[i]-theta_b[j])*.5));
	}
      }
    }
    if ((modenum_a > 0) && (modenum_b > 0)) {
      for(ajaj::MPXInt i=0;i<modenum_a;++i) {
	for(ajaj::MPXInt j=0;j<modenum_b;++j) {
	  term /= tanh((theta_a[i]-theta_b[j])*.5);
	}
      }
    }
    return(complex<double>(term*1.35783834*pow(Delta,.125),0.0));
  }

  ajaj::Vertex VertexGenerator(const ajaj::VertexParameterArray& inputs)
  {
    static const unsigned short int MAX_NUM_CHAIN_STATES=500;//useful buffer size for forming matrix elements for the occupation number as we generate the spectrum 
    //initialise vertex object
    ajaj::Vertex ModelVertex(inputs);
    ModelVertex.ChargeRules.push_back(0); //a momentum like quantum number
    ModelVertex.ChargeRules.push_back(2); //Z_2 quantum number
   

    const ajaj::QNVector& ChargeRules=ModelVertex.ChargeRules;

    const double R = inputs[0].Value;
    const double mass=inputs[1].Value;
    const double spectrum_cutoff=inputs[2].Value;

    bool use_sector(1); //default true
    double transverse_field(0.0);
    if (inputs.size()>3){
      use_sector=static_cast<int>(inputs[3].Value);
    }

    const double Delta = abs(mass); //use the magnitude of delta
    const double DR=Delta*R;
  
    const double rvac = Rvacuum(Delta,R);
    const double nvac = NSvacuum(Delta,R);
    cout << "NS vacuum energy: " << nvac << ", R vacuum energy: " << rvac  << endl;    
    const double onechain_gs_energy=nvac;
    ajaj::MPXInt num_chain_states = 0;
  
    //useful containers
    std::vector<std::vector<double> >* onechain_theta_ptr=new std::vector<std::vector<double> >();
    std::vector<std::vector<double> > &onechain_theta=*onechain_theta_ptr;
    std::vector<std::vector<double> >* onechain_kappa_ptr=new std::vector<std::vector<double> >();
    std::vector<std::vector<double> > &onechain_kappa=*onechain_kappa_ptr;  
    std::vector<int>* onechain_modenumber_ptr=new std::vector<int>();
    std::vector<int> &onechain_modenumber=*onechain_modenumber_ptr;
  
    //block for setting up occupation numbers
    int measured_occupations=6; 
    bool** occupations=new bool*[measured_occupations]; //we use double the mode integer, because of the half integer cases
    for (int i=0;i<measured_occupations;++i){
      occupations[i]=new bool[MAX_NUM_CHAIN_STATES];
      for (int j=0;j<MAX_NUM_CHAIN_STATES;++j){
	occupations[i][j]=0;
      }
    }

    std::vector<int> total_occupation;

    if (mass<0){
      const ajaj::QuantumNumberInt num1R = 200;
      const ajaj::QuantumNumberInt num2NS = 100;
      const ajaj::QuantumNumberInt num3R = 50;
      const ajaj::QuantumNumberInt num4NS = 50;
      const ajaj::QuantumNumberInt num5R = 20;
      const ajaj::QuantumNumberInt num6NS = 20;
      /*NS sector 0 particle mode, continuum chain ground state */    
      ajaj::QuantumNumberInt gsarray[2] ={0,0};

      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(gsarray,gsarray+2),onechain_gs_energy));

      onechain_modenumber.push_back(0);
      onechain_kappa.push_back(std::vector<double>());
      onechain_kappa[0].push_back(0.0); //just dummy values for ground state

      onechain_theta.push_back(std::vector<double>());
      onechain_theta[0].push_back(0.0); //just dummy values for ground state

      cout << "NS0: " << num_chain_states << " * " << ModelVertex.Spectrum[0].en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[0][0]<< " " << onechain_theta[0][0] << endl;
      cout << "mode 0 particles: 1" << endl << endl;
      num_chain_states++;
      total_occupation.push_back(0);

      /*Ramond Sector -- one particle modes*/
      ajaj::MPXInt mode1part_R =0;
      if (Delta +rvac <= spectrum_cutoff){
	for(ajaj::QuantumNumberInt i=-num1R;i<=num1R;++i) {
	  double en = Delta*sqrt(1.+pow(2.*M_PI*i/DR,2.)) +rvac;
	  if (en < spectrum_cutoff) {       
	    mode1part_R++;
	    ajaj::QuantumNumberInt R1values[2]={i,1};
	    ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R1values,R1values+2),en));
	    check_mode(2*i,num_chain_states,measured_occupations,occupations);
	  
	    onechain_modenumber.push_back(1);
	    onechain_theta.push_back(std::vector<double>());
	    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI/DR*i));
	  
	    onechain_kappa.push_back(std::vector<double>());
	    onechain_kappa[num_chain_states].push_back(Kap(asinh(2.*M_PI/DR*i),DR));
	    cout << "R1: " << num_chain_states << " " << i << " " << en << " " << (2.*M_PI/R)*ModelVertex.Spectrum[num_chain_states][0] << " " << onechain_theta[num_chain_states][0] << endl;
	    num_chain_states++;
	    total_occupation.push_back(1);
	  }
	}
      }
      cout << "mode 1 particles: " << mode1part_R << endl << endl;
    
      /*Neveu-Schwarz Sector -- two particle modes*/
      ajaj::MPXInt mode2part_NS=0;
      if (2*Delta*(sqrt(1.+pow(2.*M_PI*(.5)/DR,2.))) +nvac <= spectrum_cutoff){
       
	for(ajaj::QuantumNumberInt i=-num2NS-1;i<num2NS;++i) {
	  if (Delta +rvac > spectrum_cutoff){break;}

	  for(ajaj::QuantumNumberInt j=i+1;j<=num2NS;++j) {
	    double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+nvac;
	    if (en< spectrum_cutoff) {
	      mode2part_NS++;
	      ajaj::QuantumNumberInt NS2values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+1),0};
	      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS2values,NS2values+2),en));
	      check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
	      check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
	     
	      onechain_modenumber.push_back(2);
	      onechain_theta.push_back(std::vector<double>());
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
	     
	      onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode3part_R = 0;
      if (Delta+2.0*Delta*(sqrt(1.+pow(2.*M_PI/DR,2.))) +rvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num3R;i<num3R-1;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num3R;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<=num3R;++k) {
	      double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*j/DR,2.)) + sqrt(1.+pow(2.*M_PI*k/DR,2.))) + rvac;
	      if (en < spectrum_cutoff) {
		mode3part_R++;
		ajaj::QuantumNumberInt R3values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k),1};
		ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R3values,R3values+2),en));
		check_mode(2*i,num_chain_states,measured_occupations,occupations);
		check_mode(2*j,num_chain_states,measured_occupations,occupations);
		check_mode(2*k,num_chain_states,measured_occupations,occupations);
		  
		onechain_modenumber.push_back(3);
		onechain_theta.push_back(std::vector<double>());
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*i/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*j/DR));
		onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*k/DR));
		  
		onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode4part_NS = 0.;
      if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))) +nvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num4NS-1;i<num4NS-2;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num4NS-1;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num4NS;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<=num4NS;++l) {
		double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+nvac;
		if (en < spectrum_cutoff) {
		  mode4part_NS++;
		  ajaj::QuantumNumberInt NS4values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+2),0};	     
		  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS4values,NS4values+2),en));
		  check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
		  check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
		  check_mode(2*k+1,num_chain_states,measured_occupations,occupations);
		  check_mode(2*l+1,num_chain_states,measured_occupations,occupations);

		  onechain_modenumber.push_back(4);	    
		  onechain_theta.push_back(std::vector<double>());
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
	    
		  onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode5part_R = 0;
      if (Delta+2.0*Delta*(sqrt(1.+pow(2.*M_PI/DR,2.))+sqrt(1.+pow(2.*M_PI*2.0/DR,2.))) +rvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num5R;i<num5R-3;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num5R-2;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num5R-1;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<num5R;++l) {
		for(ajaj::QuantumNumberInt m=k+1;m<=num5R;++m) {
		  double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.)) + sqrt(1.+pow(2.*M_PI*(j)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(m)/DR,2.))) + rvac;//-onechain_gs_energy;
		  if (en < spectrum_cutoff) {
		    mode5part_R++;
		    ajaj::QuantumNumberInt R5values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+m),1};
		    ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R5values,R5values+2),en));
		    check_mode(2*i,num_chain_states,measured_occupations,occupations);
		    check_mode(2*j,num_chain_states,measured_occupations,occupations);
		    check_mode(2*k,num_chain_states,measured_occupations,occupations);
		    check_mode(2*l,num_chain_states,measured_occupations,occupations);
		    check_mode(2*m,num_chain_states,measured_occupations,occupations);

		    onechain_modenumber.push_back(5);
		    onechain_theta.push_back(std::vector<double>());
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*i/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*j/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*k/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*l/DR));
		    onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*m/DR));
	      
		    onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode6part_NS = 0;
      if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))+sqrt(1.+pow(2.*M_PI*2.5/DR,2.))) +nvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num6NS-1;i<num6NS-4;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num6NS-3;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num6NS-2;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<num6NS-1;++l) {
		for(ajaj::QuantumNumberInt m=l+1;m<num6NS;++m) {
		  for(ajaj::QuantumNumberInt n=m+1;n<=num6NS;++n) {
		    double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n+0.5)/DR,2.))+nvac;
		    if (en < spectrum_cutoff) {
		      mode6part_NS++;
		      ajaj::QuantumNumberInt NS6values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+m+n+3),0};
		      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS6values,NS6values+2),en));
		      check_mode(2*i+1,num_chain_states,measured_occupations,occupations);
		      check_mode(2*j+1,num_chain_states,measured_occupations,occupations);
		      check_mode(2*k+1,num_chain_states,measured_occupations,occupations);
		      check_mode(2*l+1,num_chain_states,measured_occupations,occupations);
		      check_mode(2*m+1,num_chain_states,measured_occupations,occupations);
		      check_mode(2*n+1,num_chain_states,measured_occupations,occupations);

		      onechain_modenumber.push_back(6);
		      onechain_theta.push_back(std::vector<double>());
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n+0.5)/DR));		
		      onechain_kappa.push_back(std::vector<double>());
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

      const ajaj::QuantumNumberInt num2R = 100;
      const ajaj::QuantumNumberInt num2NS = 100;
      const ajaj::QuantumNumberInt num4R = 50;
      const ajaj::QuantumNumberInt num4NS = 50;
      const ajaj::QuantumNumberInt num6R = 20;
      const ajaj::QuantumNumberInt num6NS = 20;
    
      //push back the two vacuua
      ajaj::QuantumNumberInt nvacvalues[2]={0,0};		    
      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(nvacvalues,nvacvalues+2),nvac));

      onechain_modenumber.push_back(0);
      onechain_theta.push_back(std::vector<double>());
      onechain_kappa.push_back(std::vector<double>());
      onechain_kappa[0].push_back(0.0); //just dummy values
      onechain_theta[0].push_back(0.0); //just dummy values
      num_chain_states++;
      total_occupation.push_back(0);

      ajaj::QuantumNumberInt rvacvalues[2] ={0,1}; //1 is ramond sector
      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(rvacvalues,rvacvalues+2),rvac));

      onechain_modenumber.push_back(0);
      onechain_theta.push_back(std::vector<double>());
      onechain_kappa.push_back(std::vector<double>());
      onechain_kappa[0].push_back(0.0); //just dummy values
      onechain_theta[0].push_back(0.0); //just dummy values
      num_chain_states++;
      total_occupation.push_back(0);

      std::cout << "mode 0 particles: " << total_occupation.size() << std::endl << std::endl;
  
      /*Now do 2 particle modes*/	    
      /* R sector*/
      ajaj::MPXInt mode2part_R=0;
      if (Delta*(1.0+sqrt(1.+pow(2.*M_PI/DR,2.))) +rvac <= spectrum_cutoff){ //check lowest state
	for(ajaj::QuantumNumberInt i=-num2R-1;i<num2R;++i) {
	  if (Delta +rvac > spectrum_cutoff){break;}
	  for(ajaj::QuantumNumberInt j=i+1;j<=num2R;++j) {
	    double en = Delta*(sqrt(1.+pow(2.*M_PI*i/DR,2.))+sqrt(1.+pow(2.*M_PI*j/DR,2.)))+nvac;
	    if (en< spectrum_cutoff) {
	      mode2part_R++;
	      ajaj::QuantumNumberInt R2values[2]={static_cast<ajaj::QuantumNumberInt>(i+j),1};
	      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R2values,R2values+2),en));

	      check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
	      check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
	     
	      onechain_modenumber.push_back(2);
	      onechain_theta.push_back(std::vector<double>());
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
	     
	      onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode2part_NS=0;
      if (2*Delta*(sqrt(1.+pow(2.*M_PI*(.5)/DR,2.))) +nvac <= spectrum_cutoff){
	for(ajaj::QuantumNumberInt i=-num2NS-1;i<num2NS;++i) {
	  if (Delta +rvac > spectrum_cutoff){break;}

	  for(ajaj::QuantumNumberInt j=i+1;j<=num2NS;++j) {
	    double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+nvac;
	    if (en< spectrum_cutoff) {
	      mode2part_NS++;
	      ajaj::QuantumNumberInt NS2values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+1),0};
	      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS2values,NS2values+2),en));

	      check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
	      check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
	     
	      onechain_modenumber.push_back(2);
	      onechain_theta.push_back(std::vector<double>());
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
	      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
	     
	      onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode4part_R = 0.;
      if (Delta*(1.0+2.0*sqrt(1.+pow(2.*M_PI/DR,2.))+sqrt(1.+pow(2.*M_PI*2.0/DR,2.))) +rvac <= spectrum_cutoff){
	for(ajaj::QuantumNumberInt i=-num4R-1;i<num4R-2;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num4R-1;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num4R;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<=num4R;++l) {
		double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+rvac;
		if (en < spectrum_cutoff) {
		  mode4part_R++;
		  ajaj::QuantumNumberInt R4values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l),1};
		  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R4values,R4values+2),en));

		  check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_R(k,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_R(l,num_chain_states,measured_occupations,occupations);

		  onechain_modenumber.push_back(4);	    
		  onechain_theta.push_back(std::vector<double>());
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l)/DR));
	    
		  onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode4part_NS = 0.;
      if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))) +nvac <= spectrum_cutoff){
	for(ajaj::QuantumNumberInt i=-num4NS-1;i<num4NS-2;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num4NS-1;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num4NS;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<=num4NS;++l) {
		double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+nvac;
		if (en < spectrum_cutoff) {
		  mode4part_NS++;
		  ajaj::QuantumNumberInt NS4values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+2),0};
		  ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS4values,NS4values+2),en));

		  check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_NS(k,num_chain_states,measured_occupations,occupations);
		  check_mode_ordered_NS(l,num_chain_states,measured_occupations,occupations);

		  onechain_modenumber.push_back(4);	    
		  onechain_theta.push_back(std::vector<double>());
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		  onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
	    
		  onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode6part_R = 0;
      if (Delta*(1.0+2.0*sqrt(1.+pow(2.*M_PI/DR,2.))+2.0*sqrt(1.+pow(2.*M_PI*2.0/DR,2.))+sqrt(1.+pow(2.*M_PI*3.0/DR,2.))) +rvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num6R-1;i<num6R-4;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num6R-3;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num6R-2;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<num6R-1;++l) {
		for(ajaj::QuantumNumberInt m=l+1;m<num6R;++m) {
		  for(ajaj::QuantumNumberInt n=m+1;n<=num6R;++n) {
		    double en = Delta*(sqrt(1.+pow(2.*M_PI*(i)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n)/DR,2.))+nvac;
		    if (en < spectrum_cutoff) {
		      mode6part_R++;
		      ajaj::QuantumNumberInt R6values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+m+n),1};
		      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(R6values,R6values+2),en));

		      check_mode_ordered_R(i,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_R(j,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_R(k,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_R(l,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_R(m,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_R(n,num_chain_states,measured_occupations,occupations);

		      onechain_modenumber.push_back(6);
		      onechain_theta.push_back(std::vector<double>());
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n)/DR));		
		      onechain_kappa.push_back(std::vector<double>());
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
      ajaj::MPXInt mode6part_NS = 0;
      if (2.0*Delta*(sqrt(1.+pow(2.*M_PI*0.5/DR,2.))+sqrt(1.+pow(2.*M_PI*1.5/DR,2.))+sqrt(1.+pow(2.*M_PI*2.5/DR,2.))) +nvac <= spectrum_cutoff){

	for(ajaj::QuantumNumberInt i=-num6NS-1;i<num6NS-4;++i) {
	  for(ajaj::QuantumNumberInt j=i+1;j<num6NS-3;++j) {
	    for(ajaj::QuantumNumberInt k=j+1;k<num6NS-2;++k) {
	      for(ajaj::QuantumNumberInt l=k+1;l<num6NS-1;++l) {
		for(ajaj::QuantumNumberInt m=l+1;m<num6NS;++m) {
		  for(ajaj::QuantumNumberInt n=m+1;n<=num6NS;++n) {
		    double en = Delta*(sqrt(1.+pow(2.*M_PI*(i+.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(j+.5)/DR,2.)))+ sqrt(1.+pow(2.*M_PI*(k+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(l+0.5)/DR,2.))+sqrt(1.+pow(2.*M_PI*(m+0.5)/DR,2.)) + sqrt(1.+pow(2.*M_PI*(n+0.5)/DR,2.))+nvac;
		    if (en < spectrum_cutoff) {
		      mode6part_NS++;
		      ajaj::QuantumNumberInt NS6values[2]={static_cast<ajaj::QuantumNumberInt>(i+j+k+l+m+n+3),0};
		      ModelVertex.Spectrum.push_back(ajaj::EigenState(ChargeRules,ajaj::QNVector(NS6values,NS6values+2),en));

		      check_mode_ordered_NS(i,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_NS(j,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_NS(k,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_NS(l,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_NS(m,num_chain_states,measured_occupations,occupations);
		      check_mode_ordered_NS(n,num_chain_states,measured_occupations,occupations);

		      onechain_modenumber.push_back(6);
		      onechain_theta.push_back(std::vector<double>());
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(i+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(j+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(k+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(l+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(m+0.5)/DR));
		      onechain_theta[num_chain_states].push_back(asinh(2.*M_PI*(n+0.5)/DR));		
		      onechain_kappa.push_back(std::vector<double>());
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

    ajaj::EigenStateArray TempSpectrum;
    for (size_t s=0;s<ModelVertex.Spectrum.size();++s){
      TempSpectrum.push_back(ModelVertex.Spectrum[sort_keys[s].first]); 
    }
    if (ModelVertex.Spectrum.size()!=static_cast<size_t>(num_chain_states)){
      cout << "Mismatch error" << endl; exit(1);
    }
    //TempSpectrum is now ordered by energy

    cout << "End generating chain spectrum, generated " << num_chain_states << " states" << endl;
    cout << "Starting Matrix Elements" << endl;
#if defined NDEBUG
    double sr = S(DR); //really slow right now
    cout << "S(R): " << sr << endl;
#else
    double sr=1.0;
    cout << "WARNING APPROXIMATING S(R) AS " << sr <<  endl;
#endif
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Spin_Operator",ModelVertex.Spectrum.size()));      
    //off diagonal only
    for (size_t col=0;col<TempSpectrum.size();++col){
      size_t s_col=sort_keys[col].first;
      for (size_t row=col+1;row<TempSpectrum.size();++row){
	size_t s_row=sort_keys[row].first;
	if(TempSpectrum[row][1]!=TempSpectrum[col][1]){ //off diagonal in Ising 'sector' charge.
	  complex<double> M_E_=sr*calculate_chainmatrixelement(ModelVertex.Spectrum,static_cast<ajaj::MPXInt>(s_row),static_cast<ajaj::MPXInt>(s_col),static_cast<ajaj::MPXInt>(onechain_modenumber[s_row]),static_cast<ajaj::MPXInt>(onechain_modenumber[s_col]),onechain_kappa[s_row],onechain_kappa[s_col],onechain_theta[s_row],onechain_theta[s_col],Delta,DR);
	  ModelVertex.Operators.back().MatrixElements.entry(row,col,M_E_);
	  ModelVertex.Operators.back().MatrixElements.entry(col,row,conj(M_E_));
	}
      }
    }
    ModelVertex.Operators.back().MatrixElements.finalise();

    //make total number operator
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Total_Number",ModelVertex.Spectrum.size()));
    for (size_t j=0;j<ModelVertex.Spectrum.size();++j){
      ModelVertex.Operators.back().MatrixElements.entry(j,j,static_cast<complex<double> >(total_occupation[sort_keys[j].first]));
    }
    ModelVertex.Operators.back().MatrixElements.finalise();

    //make number operators
    for (ajaj::MPXInt i=0;i<measured_occupations;++i){
      std::stringstream name;
      name << "Mode_Occupation_" << i+1;
      ModelVertex.Operators.push_back(ajaj::VertexOperator(name.str(),ModelVertex.Spectrum.size()));
      //ajaj::SparseMatrix occ_number(ModelVertex.Spectrum.size(),ModelVertex.Spectrum.size(),ModelVertex.Spectrum.size());
      for (size_t j=0;j<ModelVertex.Spectrum.size();++j){
	if(occupations[i][sort_keys[j].first]==1){ModelVertex.Operators.back().MatrixElements.entry(j,j,1.0);}
      }
      ModelVertex.Operators.back().MatrixElements.finalise();
    }

    if (!use_sector){
      //if we don't want to use
      ModelVertex.ChargeRules.resize(1);
    }

    ModelVertex.Spectrum=std::move(TempSpectrum);

    //clean up workspace
    delete onechain_theta_ptr;
    delete onechain_kappa_ptr;
    delete onechain_modenumber_ptr;
    for (ajaj::MPXInt i=0;i<measured_occupations;++i){
      delete[] occupations[i];
    }
    delete[] occupations;

    //push back Vertex Hamiltonain
    ModelVertex.Operators.push_back(ajaj::VertexOperator("Vertex_Hamiltonian"));
    ModelVertex.Operators.back().MatrixElements=ajaj::SparseMatrix(ModelVertex.basis().Energies());

    //Operators[i]
    //0 spin
    //1 total chain fermion occupation number
    //2+ individual chain mode occupations
    //last is vertex hamiltonian

    return ModelVertex;
  }

  ajaj::MPO_matrix MakeHamiltonian(const ajaj::Vertex& modelvertex, const ajaj::VertexParameterArray& couplingparams){
    //Lower triangular MPO
    // I   0   0
    // JK  0   0
    // HV  K'  I  //note the charges are such that K' and K differ in that Q[K[i]]+Q[K'[i]]=0
    ajaj::QNCombinations differencecombinations(modelvertex.Spectrum,1); //1 means use difference

    ajaj::MPXInt lineardim=modelvertex.Spectrum.size()*(2+couplingparams.size()*differencecombinations.size()); //the actual length of the sparse matrix needed
    ajaj::MPXInt offset_to_last_block=modelvertex.Spectrum.size()*(1+couplingparams.size()*differencecombinations.size()); //offset to get to the last row of operators
    ajaj::SparseMatrix M(lineardim,lineardim,lineardim);

    //start with the really easy bits, the Identities, I, and the vertex Hamiltonian HV
    for (size_t i=0;i<modelvertex.Spectrum.size();++i){
      M.entry(i,i,1.0);
      M.entry(i+offset_to_last_block,i,modelvertex.Spectrum[i].en);
      M.entry(i+offset_to_last_block,i+offset_to_last_block,1.0);
    }
    //now the more annoying pieces
    if (couplingparams[1].Value!=0.0){
      //must check that use sector is turned off, or quit
      if (modelvertex.Spectrum.getChargeRules().size()!=1){
	std::cout << "Error: if a local longitudinal field is applied the CONSERVE_SECTOR flag must be 0!" <<std::endl;
	exit(1);
      }

      //also put in a longitudinal field
      const double long_field(couplingparams[1].Value);
      //loop over Operators[0] and multiply by param
      const ajaj::SparseMatrix& spin=modelvertex.Operators[0].MatrixElements;
      for (ajaj::Sparseint col=0;col<spin.cols();++col){
	for (ajaj::Sparseint p=spin.get_p(col);p<spin.get_p(col+1);++p){
	  ajaj::MPXInt row=spin.get_i(p);
	  if (modelvertex.Spectrum[row]==modelvertex.Spectrum[col]){ //check momenta are equal
	    M.entry(row+offset_to_last_block,col,long_field*spin.get_x(p));
	  }
	}
      }
    }


    //for each coupling operator
    for (size_t c=0;c<1;++c){ //assume each of the coupling params refers to an operator
      ajaj::Sparseint operator_col_offset=c*differencecombinations.size()+1; //+1 for identity matrix in first block
      ajaj::Sparseint operator_row_offset=(couplingparams.size()-c-1)*differencecombinations.size()+1; //reversal, +1 for identity

      double operatorparam=couplingparams[c].Value;
      if (operatorparam!=0.0){
	for (ajaj::Sparseint col=0;col<modelvertex.Operators[c].MatrixElements.cols();++col){
	  for (ajaj::Sparseint p=modelvertex.Operators[c].MatrixElements.get_p(col);p<modelvertex.Operators[c].MatrixElements.get_p(col+1);++p){
	    ajaj::Sparseint i=modelvertex.Operators[c].MatrixElements.get_i(p);
	    //
	    //this could easily be a lookup function for QNCombinations
	    ajaj::State diffstate=modelvertex.Spectrum[i]-modelvertex.Spectrum[col];
	    //now look through charge groups to find MPO indices
	    ajaj::Sparseint MPO_subcol=0;
	    ajaj::Sparseint MPO_subrow=0;
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
    }

    ajaj::StateArray b;
    b.push_back(ajaj::State(modelvertex.Spectrum[0].getChargeRules())); //push back 'zero state'
    for (size_t c=0;c<couplingparams.size();++c){ //do for each operator
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
//////////////////////////////////////////////////////////////


#endif

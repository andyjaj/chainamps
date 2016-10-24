#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>

/*definitions for construction of full states*/
#define MaxStates 2000 /*maximum number of states allowed in the Hilbert space of a single boson*/
#define MaxN 20 /*largest allowed value of N in constructing highest weight states |N,M>*/
#define MaxM 20 /*largest allowed value of M in constructing highest weight states |N,M>*/

double state_en[MaxStates];

int nstate, /*index keeping track of how many non-chiral states there are*/
  state_chiral[3][MaxStates], /*chiral info: left chiral state, right chiral state, total number of left and right a's of state*/
  state_Z[3][MaxStates], /*Z quantum numbers: N, M, mom*/
  state_metaZ[MaxStates], /*Z quantum number incorporating all 3 Z quantum numbers*/
  state_metaZ_ind[MaxStates], /*index of Z quantum number incorporating all 3 Z quantum numbers*/
  state_parity[2][MaxStates]; /*parity quantum numbers: \pm 1 for parity of state; and sign of application of parity operator upon positive component of state*/

/*routine to compute non-chiral matrix elements*/
double non_chir_me(int i, int j);

/*routine to sort in order of energy non-chiral states*/
void shell_en(int n);

/*routines to compute non-chiral matrix elements*/
double non_chir_psi_dagger(int i, int j),
  non_chir_psi(int i, int j),
  density_me(double Beta, int i, int j),
  density_int_me(double Beta, int i, int j),
  a_n(int k,int l,int n),
  bar_a_n(int k,int l,int n),
  n_a(int k,int l,int n),
  bar_n_a(int k,int l,int n);


int num_metaZ; /*range of metaZ quantum number: from 0 to num_metaZ-1*/
 
#include "chiral_states_bosons_fA.hpp"
#include "element_computation_bosons.hpp"
#include "me_chir_bosons.hpp"

//#include "non_chiral_states_bosons_fA.hpp"
#include "non_chiral_states_bosons_no_parity_fA.hpp"
//#include "non_chir_me.hpp"
#include "non_chir_me_no_parity.hpp"

int main()
{

  double R=4.0;
  double tpi_R=2.0*3.1415926/R;
  double pi=3.1415926;

  double Beta= 1.; /*inverse of compactification radius; Beta^2 = 1/(4K) where K is Luttinger parameter of Eqn. 114 in
			     arXiv:1101.5337*/
  double en_cutoff = 8.; /*energy cutoff on states*/

  /*compute chiral states: this routine computes the possible ways to form
   chiral states up to a certain level, i.e. states of the form a_{n_1} \cdots a_{n_N}|N,M>*/ 

  chiral_states();

  /*compute chiral matrix elements*/
  me_chir(Beta);
  

  /*compute non-chiral states*/  
  non_chir_bosons(tpi_R,Beta,en_cutoff);

  double *psi_dagger = (double*)malloc(nstate*nstate*sizeof(double));
  double *psi = (double*)malloc(nstate*nstate*sizeof(double));
  double *density = (double*)malloc(nstate*nstate*sizeof(double));
  double *density_int = (double*)malloc(nstate*nstate*sizeof(double));
  double *phase = (double*)malloc(nstate*nstate*sizeof(double));
  double **a_left = (double**)malloc((2*UL-1)*sizeof(double*));
  double **a_right = (double**)malloc((2*UL-1)*sizeof(double*));
  double **na_left = (double**)malloc((2*UL-1)*sizeof(double*));
  double **na_right = (double**)malloc((2*UL-1)*sizeof(double*));
  for(int i=0;i<(2*UL-1);++i) {
    a_left[i] = (double*)malloc(nstate*nstate*sizeof(double));
    a_right[i] = (double*)malloc(nstate*nstate*sizeof(double));
    na_left[i] = (double*)malloc(nstate*nstate*sizeof(double));
    na_right[i] = (double*)malloc(nstate*nstate*sizeof(double));
  }

  /*compute non-chiral matrix elements*/

  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      // <i|psi^\dagger|j>
      psi_dagger[i*nstate+j] = non_chir_psi_dagger(i,j)*pow(tpi_R,(Beta*Beta));
      // <i|psi|j>
      psi[i*nstate+j] = non_chir_psi(i,j)*pow(tpi_R,(Beta*Beta));;
    }
  }

  /*checking for hermiticity*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      if (fabs(psi_dagger[i*nstate+j]-psi[j*nstate+i])>0.000001)
	printf("Warning: Hermiticity broken: %d %d %3.8f %3.8f\n",i,j,psi_dagger[i*nstate+j],psi[j*nstate+i]);
    }
  }

  /*compute density matrix elements*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      // <i|rho|j> 
      density[i*nstate+j] = -1.*density_me(Beta,i,j)*tpi_R/(4.*pi*Beta);
      printf("%3.3f ",density[i*nstate+j]);
    }
  }

  /*compute integrated density matrix elements*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      // <i|rho|j> 
      density_int[i*nstate+j] = -1.*R*density_int_me(Beta,i,j)*tpi_R/(4.*pi*Beta);
      printf("%3.3f ",density_int[i*nstate+j]);
    }
  }

  /*compute phase matrix elements*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      // <i|phase|j> 
      phase[i*nstate+j] = phase_me(Beta,i,j);
      printf("%3.3f ",phase[i*nstate+j]);
    }
  }

  //for(int i=0;i<nstate;++i) {printf("%d %3.8f %3.8f\n",i,psi_dagger[i*nstate+i],psi[i*nstate+i]);}

  /*compute a_n for left movers*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      for(int n=-UL+1;n<UL;++n) {
	a_left[n+UL-1][i*nstate+j] = a_n(i,j,n);
	if (fabs(a_n(i,j,n)-a_n(j,i,-n))>0.001)
	  printf("Warning: Hermiticity broken at %d %d %d\n",n,i,j);
      }
    }
  }

  /*compute \bar a_n for right movers*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      for(int n=-UL+1;n<UL;++n) {
	a_right[n+UL-1][i*nstate+j] = bar_a_n(i,j,n);
	if (fabs(bar_a_n(i,j,n)-bar_a_n(j,i,-n))>0.001)
	  printf("Warning: Hermiticity broken at %d %d %d\n",n,i,j);
      }
    }
  }

  /*compute mode occupancy na for left movers*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      for(int n=1;n<UL;++n) {
	na_left[n][i*nstate+j] = n_a(i,j,n);
	if (fabs(n_a(i,j,n)-n_a(j,i,n))>0.001)
	  printf("Warning: Hermiticity broken for n_a left movers at %d %d %d\n",n,i,j);
      }
    }
  }

  /*compute mode occupancy bar_na for right movers*/
  for(int i=0;i<nstate;++i) {
    for(int j=0;j<nstate;++j) {
      for(int n=1;n<UL;++n) {
	na_right[n][i*nstate+j] = bar_n_a(i,j,n);
	if (fabs(bar_n_a(i,j,n)-bar_n_a(j,i,n))>0.001)
	  printf("Warning: Hermiticity broken for bar_n_a right movers at %d %d %d\n",n,i,j);
      }
    }
  }

}






/** @file arpack_interface.hpp 
 * Interface to arpack fortran calls
 * Can do sparse and dense eigensolving for a few eigenvalues and vectors
 */
#ifndef arpack_interface_H
#define arpack_interface_H

#include <utility>
#include <algorithm>
#include <complex> //std::complex<double>
#include <iostream>
#include <cs.h> //cs_cl* type

namespace arpack {

#if defined(LAPACK_ILP64)
  /** If LAPACK_ILP64 defined, use long integers */
  typedef long arpack_int; 
#elif defined(MKL_ILP64)
  typedef long long arpack_int; 
#else
  /** Else use normal int type */
  typedef int arpack_int;
#endif

  /** external arpack function */
  extern "C" void znaupd_(arpack_int *ido, const char *bmat, arpack_int *n, const char *which, arpack_int *nev, double *tol, std::complex<double> *resid, arpack_int *ncv, std::complex<double> *v, arpack_int *ldv, arpack_int *iparam, arpack_int *ipntr, std::complex<double> *workd, std::complex<double> *workl, arpack_int *lworkl, double *rwork, arpack_int *info);
  /** external arpack function */
  extern "C" void zneupd_(arpack_int *rvec, const char *All, arpack_int *select, std::complex<double> *d, std::complex<double> *v, arpack_int *ldv, double *sigma, std::complex<double> *workev, const char* bmat, arpack_int *n, const char *which, arpack_int *nev, double *tol, std::complex<double> *resid, arpack_int *ncv, std::complex<double> *v2, arpack_int *ldv2, arpack_int *iparam, arpack_int *ipntr, std::complex<double> *workd, std::complex<double> *workl, arpack_int *lworkl, double *rwork, arpack_int *ierr);

  /** sparse matrix * vector routine */
  void av(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out);
  /** dense matrix * vector routine */
  void av(std::complex<double>* dense,arpack_int n, std::complex<double> *in, std::complex<double> *out);
  /** transpose matrix * vector */
  void aTv(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out);

  class arpack_workspace{
  public:
    static const char* all;
    static const char* bmat;

    arpack_int n;
    arpack_int nev;
    arpack_int rvec;
    arpack_int info;

    bool converge_flag;
    
    double tol;

    arpack_int maxiter;
    arpack_int ido;
    arpack_int ncv;
    arpack_int ldv;
    arpack_int ierr;
    arpack_int lworkl;

    
    double sigma;

    const char* which;
    std::complex<double>* Evals;
    std::complex<double>* Evecs;
    std::complex<double>* resid;

    std::complex<double> *v;
    arpack_int *iparam;
    arpack_int *ipntr;
    std::complex<double> *workd;
    std::complex<double> *workl;
    double *rwork;
    arpack_int *select;
    std::complex<double> *d;
    std::complex<double> *workev;

    arpack_workspace(arpack_int length, arpack_int num_e_vals, char which_e_vals[3], std::complex<double>* Evals_ptr, arpack_int need_e_vectors, std::complex<double>* Evecs_ptr, std::complex<double>* resid_ptr, bool converge=1, double use_tolerance=-0.0);
    ~arpack_workspace();
    void init();
    void reset(arpack_int ncv_delta=0);    
  };

  /** class to find eigenvalues and maybe vectors using arpack */
  template <typename ArrayType,typename GuessType>
  class arpack_eigs{
  public:
    const ArrayType* m_array_stuff; //input
    arpack_int m_length;
    GuessType* m_initial_guess; //possible input
    std::complex<double>* m_evals; //output
    std::complex<double>* m_evecs; //output
    arpack_workspace m_workspace; //workspace storage objects and params

    //need to set up workspace and get pointers to array (or components to make array) and ptrs to output containers
    arpack_eigs(const ArrayType* array_stuff, void (*MV)(const ArrayType*,std::complex<double>*,std::complex<double>*), arpack_int length, GuessType* initial_guess, void (*converter)(GuessType*,std::complex<double>*), arpack_int num_e_vals, char which_e_vals[3], std::complex<double> *Evals, std::complex<double> *Evecs=nullptr, bool converge=1, double tol=-0.0);
    ~arpack_eigs(){ delete[] resid; }

    void do_znaupd();
    void do_zneupd();

    arpack_int error_status(){return m_workspace.info+m_workspace.ierr;}
    arpack_int iterations() const;

  private:
    std::complex<double>* resid;
    arpack_int m_cumulative_iterations;
    void (*m_MV)(const ArrayType*,std::complex<double>*,std::complex<double>*);
    void (*m_converter)(GuessType*,std::complex<double>*);
    void sort();//always check the order...

  };
  
  template <typename ArrayType, typename GuessType>
  inline void arpack_eigs<ArrayType,GuessType>::arpack_eigs::sort(){
    //first create an indexed list
    std::vector<std::pair<arpack_int,std::complex<double> > > indexed_e_vals;
    for (arpack_int i=0; i<m_workspace.nev; ++i){
      indexed_e_vals.emplace_back(i,m_workspace.d[i]);
    }
    //next sort, needs to know if mode is SR or LM
    if (m_workspace.which[0]=='S' && m_workspace.which[1]=='R'){
      std::sort(indexed_e_vals.begin(),indexed_e_vals.end(),
		[](const std::pair<arpack_int,std::complex<double> >& a, const std::pair<arpack_int,std::complex<double> >& b)
		{ return real(a.second) < real(b.second);}
		);
    }
    else if (m_workspace.which[0]=='L' && m_workspace.which[1]=='M'){
      std::sort(indexed_e_vals.begin(),indexed_e_vals.end(),
		[](const std::pair<arpack_int,std::complex<double> >& a, const std::pair<arpack_int,std::complex<double> >& b)
		{ return abs(a.second) > abs(b.second);}
		);
    }
    else {
      std::cout <<"Unsupported arpack mode!" <<std::endl;
      exit(1);
    }

    //now check if anything changed
    arpack_int count=0;
    bool firstpass=1;
    for (auto& i : indexed_e_vals){
      if (i.first!=count) {
	std::cout <<"ARPACK has returned eigenvalues in the wrong order!"<<std::endl;
	if (firstpass){
	  for (auto i=0;i<m_workspace.nev;++i){
	    std::cout << m_workspace.d[i] <<std::endl;
	  }
	  firstpass=0;
	}
	/*for(auto& i : indexed_e_vals){
	  std::cout << i.first <<": " << i.second <<std::endl; 
	  }*/
	std::cout << "Reordering..." <<std::endl;
	//swap the data, and then update the indices in the list.
	for (auto j=count+1;j<m_workspace.nev;++j){
	  if (indexed_e_vals[j].first==count){
	    std::swap_ranges(&m_workspace.d[i.first],&m_workspace.d[i.first]+1,&m_workspace.d[count]);
	    std::swap_ranges(&m_evecs[i.first*m_length],&m_evecs[i.first*m_length]+m_length,&(m_evecs[count*m_length]));
	    indexed_e_vals[j].first=i.first; //update moved idx
	    i.first=count;//not necessary
	    i.second=m_workspace.d[count];
	    break;
	  }
	}
      }
      ++count;
    }

    /*std::cout <<"Ordered eigenvalues from ARPACK"<<std::endl;
    for (auto i=0;i<m_workspace.nev;++i){
      std::cout << m_workspace.d[i] <<std::endl;
      }*/
    
  }
  
  template <typename ArrayType, typename GuessType>
  arpack_eigs<ArrayType,GuessType>::arpack_eigs(const ArrayType* array_stuff, void (*MV)(const ArrayType*,std::complex<double>*,std::complex<double>*), arpack_int length, GuessType* initial_guess, void (*converter)(GuessType*,std::complex<double>*), arpack_int num_e_vals, char which_e_vals[3], std::complex<double> *Evals, std::complex<double> *Evecs, bool converge, double tol) : m_array_stuff(array_stuff), m_MV(MV), m_length(length), m_initial_guess(initial_guess),m_converter(converter),m_evals(Evals),m_evecs(Evecs),m_workspace(m_length,num_e_vals,which_e_vals,Evals,Evecs ? 1 : 0, Evecs,nullptr,converge, tol), m_cumulative_iterations(0){
    //form resid
    resid = new std::complex<double>[m_length];
    m_workspace.resid=resid;
    std::cout << "Arpack on matrix of length " << m_length << std::endl;
    if (m_initial_guess) {//if there is an initial guess, we need to feed it into resid
      std::cout << "Using initial guess" << std::endl;
      m_converter(m_initial_guess,resid);//convert guess format to dense vector
      m_workspace.info=1; //tell arpack to use initial guess vector
    }
    else { //no initial guess?
      std::cout << "No guess, trying a random starting vector..." << std::endl;
      m_workspace.info=0; //tell arpack to use random initial vector
    }

    arpack_int ncv_diff=0;

    do { //enter loop
      if (m_workspace.info==-9){//if info is -9 on rentering loop then alter things
	m_workspace.reset(); //restart the workspace
	m_workspace.info=0; //tell arpack to use random initial vector
      }
      else if (ncv_diff && m_workspace.info==1){
	m_workspace.reset(ncv_diff); //restart the workspace, changing ncv
      }

      do_znaupd();//call arpack

      if (m_workspace.info==-9)
	std::cout << "Matrix vector product is zero, trying a different starting vector..." << std::endl;

      else if (m_workspace.info==1){//if info is 1 at this point, then we need more iterations
	std::cout << iterations() << " Arnoldi iterations taken" << std::endl;
	if (iterations()<m_length){
	  std::cout << "Poor convergence, altering ncv..." << std::endl;
	  ncv_diff+=3;
	}
	else {
	  std::cout << "Too many arpack iterations, aborting..." << std::endl;
	  exit(1);
	}
      }

    } while (m_workspace.info==1 || m_workspace.info==-9);

    std::cout << "znaupd done, " << iterations() << " Arnoldi iterations taken" << std::endl;
    
    if (m_workspace.info!=0) {std::cout << "Arpack Error in znaupd: " << m_workspace.info << std::endl;}
    else {
      
      for (arpack_int s=0;s<m_workspace.nev;++s){
	m_workspace.select[s]=1;
      }
      do_zneupd();
      if (m_workspace.ierr!=0) {std::cout << "Arpack Error in zneupd: " << m_workspace.ierr << std::endl;}
      //check order

      sort();
 
      std::swap_ranges(m_workspace.d,m_workspace.d+m_workspace.nev,m_evals);
      
    }
  }

  //need to be able to call recursively
  template <typename ArrayType, typename GuessType>
  inline void arpack_eigs<ArrayType,GuessType>::do_znaupd(){
    size_t itnum=1;
    do {
      znaupd_(&m_workspace.ido, m_workspace.bmat, &m_workspace.n, m_workspace.which, &m_workspace.nev, &m_workspace.tol, resid, &m_workspace.ncv, m_workspace.v, &m_workspace.ldv, m_workspace.iparam, m_workspace.ipntr, m_workspace.workd, m_workspace.workl,&m_workspace.lworkl, m_workspace.rwork, &m_workspace.info);
      if ((m_workspace.ido==1)||(m_workspace.ido==-1)) m_MV(m_array_stuff,m_workspace.workd+m_workspace.ipntr[0]-1, m_workspace.workd+m_workspace.ipntr[1]-1);
      std::cout << "." <<std::flush; 
    } while (m_workspace.ido==1 || m_workspace.ido==-1);
    //update total number of iterations
    m_cumulative_iterations+=m_workspace.iparam[2];
    std::cout << std::endl;
  }

  template <typename ArrayType, typename GuessType>
  inline void arpack_eigs<ArrayType,GuessType>::do_zneupd(){
    zneupd_(&m_workspace.rvec, m_workspace.all, m_workspace.select, m_workspace.d, m_evecs, &m_workspace.ldv, &m_workspace.sigma, m_workspace.workev,m_workspace.bmat, &m_workspace.n, m_workspace.which, &m_workspace.nev, &m_workspace.tol, resid, &m_workspace.ncv, m_workspace.v, &m_workspace.ldv,m_workspace.iparam, m_workspace.ipntr, m_workspace.workd, m_workspace.workl, &m_workspace.lworkl, m_workspace.rwork, &m_workspace.ierr);
  }

  template <typename ArrayType,typename GuessType>
  inline arpack_int arpack_eigs<ArrayType,GuessType>::iterations() const{
    return m_cumulative_iterations;
  }

  inline void sparse_matrix_vector_mult(const cs_cl* array, std::complex<double> *in, std::complex<double> *out){av(const_cast<cs_cl*>(array),in,out);}//annoying const_cast required because cxsparse doesn't use const in its sparse matrix vector routines.
  void cs_cl_to_dense(cs_cl* sparsevec, std::complex<double>* densevec);

  /** sparse matrix * vector routine */
  void av(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out);
  /** dense matrix * vector routine */
  void av(std::complex<double>* dense,arpack_int n, std::complex<double> *in, std::complex<double> *out);
  /** transpose matrix * vector */
  void aTv(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out);

  /** Sparse eigensolver for which nev eigenvalues and eigenvectors, optional initial guess vector in sparse form*/
  bool cpparpack(const cs_cl* sparse,arpack_int n, arpack_int nev, std::complex<double> *Evals, std::complex<double> *Evecs, char which[3],cs_cl* initial=NULL);
}
#endif

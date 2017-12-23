#include <cmath>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#if defined(USETBB)
#include <tbb/tbb.h>
#endif
#include <cs.h>
#include "dense_matrix_functions.hpp"
#include "arpack_interface.hpp"

#define ARPACKTOL -0.0

using namespace std;

namespace arpack {

  const char* arpack_workspace::all="A";
  const char* arpack_workspace::bmat="I";

  void arpack_workspace::init(){
      
      ldv = n;
      maxiter = converge_flag ? 50 : 3;
      if (!converge_flag) tol=0.001;
      ido = 0;
      iparam = new arpack_int[11];
      iparam[0] = 1;
      iparam[2] = maxiter;
      iparam[6] = 1;
      ipntr = new arpack_int[14];
      workd = new std::complex<double>[3*n];
      d = new std::complex<double>[nev+1];

      ncv = 20*nev;
      if (ncv>n) ncv = n;
      v = new std::complex<double>[ldv*ncv];

      lworkl = 3*ncv*ncv + 5*ncv;
      workl = new std::complex<double>[lworkl];
      rwork = new double[ncv];
      select = new arpack_int[ncv];
      workev = new std::complex<double>[2*ncv];
      //
      std::cout << "Arpack workspace set, using at least " << double(sizeof(std::complex<double>)*(3*n+lworkl+nev+1+2*ncv+ldv*ncv+nev*n))/(1024.0*1024.0) << " megabytes (including solution space size)" << std::endl; 

  }

  arpack_workspace::arpack_workspace(arpack_int length, arpack_int num_e_vals, char which_e_vals[3], std::complex<double>* Evals_ptr, arpack_int need_e_vectors, std::complex<double>* Evecs_ptr, std::complex<double>* resid_ptr,bool converge, double use_tolerance) : n(length), nev(num_e_vals),which(which_e_vals),Evals(Evals_ptr),rvec(need_e_vectors),Evecs(Evecs_ptr),resid(resid_ptr),info(resid_ptr ? 1 : 0),converge_flag(converge),tol(use_tolerance){
    init();
  }

  arpack_workspace::~arpack_workspace(){
      delete[] v;
      delete[] iparam;
      delete[] ipntr;
      delete[] workd;
      delete[] workl;
      delete[] rwork;
      delete[] select;
      delete[] d;
      delete[] workev;
  }

  void arpack_workspace::reset(arpack_int ncv_delta){
    ido=0;
    if (iparam){//iparam not null
      iparam[0] = 1;
      iparam[2] = maxiter;
      iparam[6] = 1;
    }
    else {
      std::cout << "Null iparam! Allocation error" << std::endl; exit(1);
    }

    if (ncv_delta!=0){
      ncv = ncv+ncv_delta*nev;
      if (ncv>n) ncv = n;
      else if (ncv<2*nev) ncv=2*nev;

      lworkl = 3*ncv*ncv + 5*ncv;

      delete[] v;
      v = new std::complex<double>[ldv*ncv];
      delete[] workl;
      workl = new std::complex<double>[lworkl];
      delete[] rwork;
      rwork = new double[ncv];
      delete[] select;
      select = new arpack_int[ncv];
      delete[] workev;
      workev = new std::complex<double>[2*ncv];
    }

   
  }
  
  bool cpparpack(const cs_cl* sparse,arpack_int n, arpack_int nev, std::complex<double> *Evals, std::complex<double> *Evecs, char which[3],cs_cl* initial){
    return arpack_eigs<cs_cl,cs_cl>(sparse,&sparse_matrix_vector_mult,n,initial,&cs_cl_to_dense,nev,which,Evals,Evecs).error_status();
  }

  void cs_cl_to_dense(cs_cl* sparsevec, std::complex<double>* densevec){
    if (sparsevec->n!=1) {std::cout << "Not a sparse VECTOR!" << std::endl; exit(1);}
    std::fill(densevec,densevec+sparsevec->m,0.0);
    for (SuiteSparse_long p=sparsevec->p[0];p<sparsevec->p[1];++p){
      densevec[sparsevec->i[p]]=sparsevec->x[p];
    }
  }

#if defined(USETBB)
   struct reduce_body {
     std::complex<double>* y_; // accumulating vector
     cs_cl* C_; // reference to a matrix
     std::complex<double>* x_;    // reference to a vector
     reduce_body( cs_cl* C, std::complex<double>* x)  : C_(C), x_(x) {
       y_=new std::complex<double>[C_->n];
      //for (arpack_int row=0; row<C_->m; ++row) y_[row] = 0.0; // prepare for accumulation
      std::fill(y_,y_+C_->m,0.0);
    }
    // splitting constructor required by TBB
    reduce_body( reduce_body& rb, tbb::split ) : C_(rb.C_), x_(rb.x_) {
      y_=new std::complex<double>[C_->n];
      //for (arpack_int row=0; row<C_->m; ++row) y_[row] = 0.0;
      std::fill(y_,y_+C_->m,0.0);
    }
    ~reduce_body(){
      delete[] y_;
    }
    // the main computation method
    void operator()(const tbb::blocked_range<SuiteSparse_long>& r) { 
      // closely resembles the original serial loop
      for (SuiteSparse_long col=r.begin(); col<r.end(); ++col) // iterates over a subrange in [0,N)
	for (SuiteSparse_long p=C_->p[col]; p<C_->p[col+1]; ++p)
	  y_[C_->i[p]] += C_->x[p]*x_[col];
    }
    // the method to reduce computations accumulated in two bodies
    void join( reduce_body& rb ) {
      for (SuiteSparse_long row=0; row<C_->m; ++row) y_[row] += rb.y_[row];
    }
  };
#endif  
  //matrix vector multiplication functions
  //sparse case
  void av(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out) //complex overload
  {
#if defined(USETBB)
    reduce_body body(sparse, in);
    tbb::parallel_reduce(tbb::blocked_range<SuiteSparse_long>(0,sparse->n), body);
    std::copy(body.y_,body.y_+body.C_->m,out);
#else
    std::fill(out,out+sparse->m,0.0);
    for (SuiteSparse_long col=0;col<sparse->n;++col){
      for (SuiteSparse_long p=sparse->p[col];p<sparse->p[col+1];++p){
	out[sparse->i[p]]+=sparse->x[p]*in[col];
      }
    }
#endif
  }
  //dense case
  void av(std::complex<double> *dense,arpack_int n, std::complex<double> *in, std::complex<double> *out) //complex overload
  {
    std::fill(out,out+n,0.0);
    densefuncs::square_gemv(n,dense,in,out);
  }
  void aTv(cs_cl* sparse, std::complex<double> *in, std::complex<double> *out) //complex overload
  {
    std::fill(out,out+sparse->n,0.0);
    cs_cl* trans=cs_cl_transpose(sparse,1);
    if (!cs_cl_gaxpy(trans,in,out)){
      out=NULL; 
      std::cout << "av ERROR" <<std::endl;
      exit(1);
    }
    cs_cl_spfree(trans);
  }
}

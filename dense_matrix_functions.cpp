//file that translates lapack and blas calls
//mostly hideous
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#if defined(__APPLE__) && defined(__MACH__)
#include <Accelerate/Accelerate.h> //MacOSX blas,lapack include
#elif defined(__linux__)
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))
static const int max_num_threads_to_use=10;
#if defined(MKL_ILP64)
#define MKL_INT long long
#define MKL_Complex16 std::complex<double>
#define MKL_SET_INTERFACE_LAYER MKL_INTERFACE_ILP64
#include "mkl.h"
#include "mkl_lapacke.h"
#define set_num_threads(a) (mkl_set_num_threads(a))
#else
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_CPP
#include <cblas.h>
#include <lapacke.h> //linux lapack with c extensions
extern "C" void openblas_set_num_threads(int num_threads);
#define set_num_threads(a) (openblas_set_num_threads(a))
#endif
#else     
#error Platform not supported
#endif

#include "dense_matrix_functions.hpp"

using namespace std;

namespace densefuncs {

void gemm(const dense_int M,  const dense_int K, const dense_int N, const complex<double> alpha, const complex<double> beta, complex<double> *A, complex<double> *B, complex<double> *C)
{
  complex<double> _ALPHA=alpha;
  complex<double> _BETA=beta;
#if defined(__APPLE__) && defined(__MACH__)
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,M, N, K, &_ALPHA, A, M, B, K, &_BETA, C, M);
#elif __linux__
  //openblas annoyingly takes c type __complexdouble for complex
  set_num_threads(max_num_threads_to_use);
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, reinterpret_cast<double*>(&_ALPHA), reinterpret_cast<double*>(A), M, reinterpret_cast<double*>(B), K, reinterpret_cast<double*>(&_BETA), reinterpret_cast<double*>(C), M);
#else     
#error Platform not supported
#endif
}

void square_gemv(dense_int M, complex<double> *A, complex<double> *x, complex<double> *y)
{
  complex<double>* ONE =new complex<double>(1.0,0.0);
  complex<double>* ZERO =new complex<double>(0.0,0.0);
#if defined(__APPLE__) && defined(__MACH__)
  cblas_zgemv(CblasColMajor, CblasNoTrans, M, M, ONE, A, M, x, 1, ZERO, y, 1);
#elif __linux__
  //openblas annoyingly takes c type __complexdouble for complex
  set_num_threads(max_num_threads_to_use);
  cblas_zgemv(CblasColMajor, CblasNoTrans, M, M, reinterpret_cast<double*>(ONE), reinterpret_cast<double*>(A), M, reinterpret_cast<double*>(x), 1, reinterpret_cast<double*>(ZERO), reinterpret_cast<double*>(y), 1);
#else     
#error Platform not supported
#endif
  delete ONE;
  delete ZERO;
}

void matrix_dump(const dense_int rows,  const dense_int cols, const complex<double> *A)
{
  ofstream realoutfile ("realdump.dat");
  realoutfile << std::setprecision(16);
  if (realoutfile.is_open()){
    for (dense_int i=0; i<rows; i++){
      for (dense_int j=0; j<cols; j++){
	realoutfile << real(A[i+j*rows]) << " " ;
      }
      realoutfile << endl;
    }
    realoutfile << endl;
    realoutfile.close();
  }
  ofstream imagoutfile ("imagdump.dat");
  imagoutfile << std::setprecision(16);
  if (imagoutfile.is_open()){
    for (dense_int i=0; i<rows; i++){
      for (dense_int j=0; j<cols; j++){
	imagoutfile << imag(A[i+j*rows]) << " " ;
      }
      imagoutfile << endl;
    }
    imagoutfile << endl;
    imagoutfile.close();
  }
}

//complex overload
  void diagonalise_with_lapack(dense_int lineardim, complex<double> *matrixin, complex<double> *evecs, double *evals, dense_int il,dense_int iu){
  /* Locals */
  dense_int n=lineardim;
  //dense_int il=1;
  //dense_int iu=lineardim;
  dense_int lda=lineardim;
  dense_int ldz=lineardim;
  dense_int info=0;
  dense_int m=0;
  double abstol, vl, vu;
  /* Local arrays */
  vl=0.0;
  vu=0.0;
  //needed because zheevr expects a full array for evals, even if it only calculates some of them
  double* w= new double[lineardim];
  /* Negative abstol means using the default value */
  abstol = -1.0; 

#ifndef NDEBUG
  std::cout << "Calling zheevr on matrix " << lineardim  << "*" << lineardim << std::endl;
#endif

  /* Query and allocate the optimal workspace */
#if defined(__APPLE__) && defined(__MACH__)
  dense_int lwork=-1;
  dense_int liwork=-1;
  dense_int lrwork=-1;
  dense_int iwkopt=0;
  double rwkopt=0.0;
  complex<double> wkopt (0.0,0.0);
  char V[] = {'V','\0'};
  char I[] = {'I','\0'};
  char A[] = {'A','\0'};
  char U[] = {'U','\0'};
  dense_int* isuppz =new dense_int[2*lineardim];
  char* which=0;
  if (iu-il+1==n){which=A;}
  else {which=I;}
  //std::cout << lda << std::endl;
  zheevr_(V,which,U, &n, (__CLPK_doublecomplex*)matrixin, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, (__CLPK_doublecomplex*)evecs, &ldz, isuppz, (__CLPK_doublecomplex*) &wkopt, &lwork, &rwkopt, &lrwork, &iwkopt, &liwork, &info );
  if (info!=0){
    std::cout << "Something went wrong setting up the workspace in zheevr" << std::endl;
  }
  lwork = (dense_int)real(wkopt);
  complex<double>* work = new complex<double>[lwork];
  if (work==NULL) {
    cout << "Can't allocate work, diagonalise_with_lapack" <<endl;
    exit(1);
  }
  lrwork=(dense_int) rwkopt;
  double* rwork = new double[lrwork];
  liwork = iwkopt;
  dense_int* iwork = new dense_int[liwork];
  if (iwork==NULL) {
    cout << "Can't allocate iwork, diagonalise_with_lapack" << endl;
    exit(1);
  }
  zheevr_(V,which,U, &n, (__CLPK_doublecomplex*)matrixin, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, (__CLPK_doublecomplex*)evecs, &ldz, isuppz, (__CLPK_doublecomplex*)work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
  /* Free workspace */
  delete[] isuppz;
  delete[] iwork;
  delete[] work;
  delete[] rwork;
#elif defined(__linux__)
  char V = 'V';
  char I = 'I';
  char U = 'U';
  dense_int* isuppz =new dense_int[2*lineardim];
  set_num_threads(max_num_threads_to_use);
  info=LAPACKE_zheevr(LAPACK_COL_MAJOR,V,I,U,n, reinterpret_cast<lapack_complex_double*>(matrixin), lda, vl, vu, il, iu, abstol, &m, w, reinterpret_cast<lapack_complex_double*>(evecs), ldz, isuppz);
  delete[] isuppz;
#else     
#error Platform not supported
#endif
  std::copy(w,w+iu-il+1,evals);
  delete[] w;
  /* Check for convergence */
  if( info != 0 ) {
    cout << "The algorithm failed to compute eigenvalues." << endl;
    cout << "Info " << info << endl;
    exit( 1 );
  }
}

void diagonalise_with_lapack_nh(dense_int lineardim, complex<double> *matrixin, complex<double> *evals){
  /* Locals */
  dense_int n=lineardim;
  dense_int info=0;
  std::complex<double>* dummy_ptr=0;

#ifndef NDEBUG
  std::cout << "Calling zgeev on matrix " << lineardim  << "*" << lineardim << std::endl;
#endif

  /* Query and allocate the optimal workspace */
#if defined(__APPLE__) && defined(__MACH__)
  char NO[]={'N','\0'};
  dense_int lwork=-1;
  std::complex<double> wkopt (0.0,0.0);
  double* rwork=new double[2*lineardim];
  zgeev_(NO,NO, &n, (__CLPK_doublecomplex*)matrixin, &n, (__CLPK_doublecomplex*)evals, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*) &wkopt, &lwork, rwork, &info);
  lwork = (dense_int)real(wkopt);
  complex<double>* work = new complex<double>[lwork];
  if (work==NULL) {
    cout << "Can't allocate work, diagwithlapack" <<endl;
    exit(1);
  }
  zgeev_(NO,NO, &n, (__CLPK_doublecomplex*)matrixin, &n, (__CLPK_doublecomplex*)evals, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*) work, &lwork, rwork, &info);
  lwork = (dense_int)real(wkopt);
  /* Free workspace */
  delete[] work;
  delete[] rwork;
#elif defined(__linux__)
  set_num_threads(max_num_threads_to_use);
  char NO='N';
  info=LAPACKE_zgeev(LAPACK_COL_MAJOR,NO,NO,n,reinterpret_cast<lapack_complex_double*>(matrixin),n,reinterpret_cast<lapack_complex_double*>(evals),reinterpret_cast<lapack_complex_double*>(dummy_ptr),n, reinterpret_cast<lapack_complex_double*>(dummy_ptr),n);

  if (info!=0){cout << "Error: " << info << endl;exit(1);}

#else     
#error Platform not supported
#endif
  /* Check for convergence */
  if( info != 0 ) {
    cout << "The algorithm failed to compute eigenvalues." << endl;
    cout << "Info " << info << endl;
    exit( 1 );
  }
}

void diagonalise_with_lapack_nh(dense_int lineardim, complex<double> *matrixin, complex<double> *revecs, complex<double> *evals){
  /* Locals */
  dense_int n=lineardim;
  dense_int info=0;
  complex<double>* dummy_ptr=0;

#ifndef NDEBUG
  std::cout << "Calling zgeev on matrix " << lineardim  << "*" << lineardim << std::endl;
#endif

  /* Query and allocate the optimal workspace */
#if defined(__APPLE__) && defined(__MACH__)
  char NO[]={'N','\0'};
  char RV[]={'V','\0'};
  dense_int lwork=-1;
  complex<double> wkopt (0.0,0.0);
  double* rwork=new double[2*lineardim];
  zgeev_(NO,RV, &n, (__CLPK_doublecomplex*)matrixin, &n, (__CLPK_doublecomplex*)evals, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*)revecs, &n, (__CLPK_doublecomplex*) &wkopt, &lwork, rwork, &info);
  lwork = (dense_int)real(wkopt);
  complex<double>* work = new complex<double>[lwork];
  if (work==NULL) {
    cout << "Can't allocate work, diagwithlapack" <<endl;
    exit(1);
  }
  zgeev_(NO,RV, &n, (__CLPK_doublecomplex*)matrixin, &n, (__CLPK_doublecomplex*)evals, (__CLPK_doublecomplex*)dummy_ptr, &n, (__CLPK_doublecomplex*)revecs, &n, (__CLPK_doublecomplex*) work, &lwork, rwork, &info);
  lwork = (dense_int)real(wkopt);
  /* Free workspace */
  delete[] work;
  delete[] rwork;
#elif defined(__linux__)
  set_num_threads(max_num_threads_to_use);
  char NO='N';
  char RV='V';
  info=LAPACKE_zgeev(LAPACK_COL_MAJOR,NO,RV,n,reinterpret_cast<lapack_complex_double*>(matrixin),n,reinterpret_cast<lapack_complex_double*>(evals),reinterpret_cast<lapack_complex_double*>(dummy_ptr),n, reinterpret_cast<lapack_complex_double*>(revecs),n);
  if (info!=0){cout << "Error: " << info << endl;exit(1);}
#else     
#error Platform not supported
#endif
  /* Check for convergence */
  if( info != 0 ) {
    cout << "The algorithm failed to compute eigenvalues." << endl;
    cout << "Info " << info << endl;
    exit( 1 );
  }
}

//complex overload
void svd_with_lapack(const dense_int rows, const dense_int cols, complex<double> *original_matrix, complex<double> *u, complex<double> *vt, double *svals){
  dense_int min = MIN(rows,cols);
  dense_int max = MAX(rows,cols);
  dense_int m = rows, n = cols, lda =rows, ldu = rows, ldvt = min;
  dense_int info=0;
  //cout << "Doing SVD" << endl;
  //figure workspace
  //for absolute safety, we copy the matrix into a new array, because it might be needed again if the routine doesn't converge...
  complex<double>* matrix=new complex<double>[rows*cols];
  for (dense_int ci=0;ci<rows*cols;++ci){
    matrix[ci]=original_matrix[ci];
  }
#if defined(__APPLE__) && defined(__MACH__)
  dense_int lwork;
  complex<double> wkopt = complex<double> (0.0,0.0);
  //char JOB[]={'A','L','L','\0'};
  char JOB[]={'S','\0'};
  
  //double* rwork =new double[MIN(m,n)*MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1)];
  double* rwork =new double[min*MAX(5*min+7,2*max+2*min+1)];

  dense_int* iwork =new dense_int[8*min]; 
  lwork = -1;
  if (rows==0) {std::cout<< "zero rows!" << std::endl; exit(1);}
  zgesdd_( JOB, &m, &n, reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt, reinterpret_cast<__CLPK_doublecomplex*>(&wkopt),&lwork, rwork, iwork, &info );
  if (info==0){
    lwork = (dense_int)real(wkopt);
    complex<double>* work = new complex<double>[lwork];
    zgesdd_( JOB, &m, &n, reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt, reinterpret_cast<__CLPK_doublecomplex*>(work), &lwork, rwork, iwork, &info );
    //delete workspace
    delete[] work;
  }
  delete[] rwork;
  delete[] iwork;
#elif defined(__linux__)
  //  char JOB='A';
  char JOB='S';
  set_num_threads(max_num_threads_to_use);
#ifndef NDEBUG
  std::cout << "Calling LAPACKE_zgesdd" << std::endl;
  std::cout<< m  << " " << n << " " << lda << " " << ldu << std::endl;
#endif
  info=LAPACKE_zgesdd(LAPACK_COL_MAJOR,JOB, m, n, reinterpret_cast<lapack_complex_double*>(matrix), lda, svals, reinterpret_cast<lapack_complex_double*>(u), ldu, reinterpret_cast<lapack_complex_double*>(vt), ldvt);
#else
#error Platform not supported
  exit(1);
#endif
  if (info < 0){std::cout << "Illegal argument to zgesdd, " << abs(info) << std::endl;
    std::cout << rows << "," << cols <<std::endl;
  }
  else if( info > 0 ) {
    cout << "The algorithm (zgesdd) computing SVD failed to converge." << endl;
    cout << "DUMP INFO " << info << endl;
    cout << rows << " " << cols << " " << endl;
    cout << m << " " << n << " " << lda << " " << ldu << " " << ldvt << " " << endl;
    matrix_dump(rows,cols,original_matrix);
    cout << "!!!!!!!!!!!!!!!! Trying another SVD algorithm..." << endl;
    svd_with_lapack_fallback(rows,cols,original_matrix,u,vt,svals);
  }
  delete[] matrix;
  //cout << "Finished SVD" << endl;
}

void svd_with_lapack_fallback(const dense_int rows, const dense_int cols, complex<double> *matrix, complex<double> *u, complex<double> *vt, double *svals){
  dense_int min=MIN(rows,cols);
  dense_int max=MAX(rows,cols);
  dense_int m = rows, n = cols, lda =rows, ldu = rows, ldvt = min;
  dense_int info=0;
  //figure workspace
#if defined(__APPLE__) && defined(__MACH__)
  dense_int lwork;
  complex<double> wkopt = complex<double> (0.0,0.0);
  //char a[]={'A','\0'};
  char JOB[]={'S','\0'};
  double* rwork =new double[5*min];
  lwork = -1;
  zgesvd_(JOB,JOB, &m, &n,reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt,reinterpret_cast<__CLPK_doublecomplex*>(&wkopt),&lwork,rwork,&info);
  if (info==0){
    lwork = (dense_int)real(wkopt);
    complex<double>*  work = new complex<double>[lwork];
    zgesvd_(JOB,JOB, &m, &n,reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt, reinterpret_cast<__CLPK_doublecomplex*>(work),&lwork,rwork,&info);
    delete[] work;
    if (info < 0){cout << "Illegal argument to zgesvd, " << abs(info) << endl; delete[] rwork; exit(1);}
    else if( info > 0 ) {cout << "Still failed... " << info << endl; delete[] rwork; exit(1);}
  }
  delete[] rwork;
#elif defined(__linux__)
  set_num_threads(max_num_threads_to_use);
  //char A='A';
  char JOB='S';
  double superb[5*min];
  info=LAPACKE_zgesvd(LAPACK_COL_MAJOR,JOB,JOB,m,n,reinterpret_cast<lapack_complex_double*>(matrix),lda,svals,reinterpret_cast<lapack_complex_double*>(u),ldu,reinterpret_cast<lapack_complex_double*>(vt),ldvt,superb);
#else     
#error Platform not supported
  exit(1);
#endif
  if (info < 0){cout << "Illegal argument to zgesvd, " << abs(info) << endl; exit(1);}
  else if( info > 0 ) {cout << "Still failed..." << endl; exit(1);}
}

void svd_with_lapack(const dense_int rows, const dense_int cols, complex<double> *original_matrix, double *svals){
  dense_int m = rows, n = cols, lda =rows, ldu = rows, ldvt = cols;
  dense_int info=0;
  //cout << "Doing SVD" << endl;
  //figure workspace
  //for absolute safety, we copy the matrix into a new array, because it might be needed again if the routine doesn't converge...
  complex<double>* matrix=new complex<double>[rows*cols];

  for (dense_int ci=0;ci<rows*cols;++ci){
    matrix[ci]=original_matrix[ci];
  }

  complex<double>* u=NULL;
  complex<double>* vt=NULL;

#if defined(__APPLE__) && defined(__MACH__)
  dense_int lwork;
  complex<double> wkopt = complex<double> (0.0,0.0);
  char none[]={'N','\0'};
  dense_int min,max;
  rows >= cols ? min=cols : min=rows;
  rows >= cols ? max=rows : max=cols; 
  double* rwork =new double[MIN(m,n)*MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1)];
  //double* rwork =new double[5*MIN(rows,cols)*MIN(rows,cols)+7*MIN(rows,cols)+1];
  dense_int* iwork =new dense_int[8*MIN(m,n)]; 
  lwork = -1;
  zgesdd_( none, &m, &n, reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt, reinterpret_cast<__CLPK_doublecomplex*>(&wkopt),&lwork, rwork, iwork, &info );
  if (info==0){
    lwork = (dense_int)real(wkopt);
    complex<double>* work = new complex<double>[lwork];
    zgesdd_( none, &m, &n, reinterpret_cast<__CLPK_doublecomplex*>(matrix), &lda, svals, reinterpret_cast<__CLPK_doublecomplex*>(u), &ldu, reinterpret_cast<__CLPK_doublecomplex*>(vt), &ldvt, reinterpret_cast<__CLPK_doublecomplex*>(work), &lwork, rwork, iwork, &info );
    //delete workspace
    delete[] work;
  }
  delete[] rwork;
  delete[] iwork;
#elif defined(__linux__)
  char None='N';
  set_num_threads(max_num_threads_to_use);
  info=LAPACKE_zgesdd(LAPACK_COL_MAJOR,None, m, n, reinterpret_cast<lapack_complex_double*>(matrix), lda, svals, reinterpret_cast<lapack_complex_double*>(u), ldu, reinterpret_cast<lapack_complex_double*>(vt), ldvt);
#else
#error Platform not supported
  exit(1);
#endif

  if (info < 0){std::cout << "Illegal argument to zgesdd, " << abs(info) << std::endl;
    std::cout << rows << "," << cols <<std::endl;
  }
  else if( info > 0 ) {
    cout << "The algorithm (zgesdd) computing SVD failed to converge." << endl;
    cout << "DUMP INFO " << info << endl;
    cout << rows << " " << cols << " " << endl;
    cout << m << " " << n << " " << lda << " " << ldu << " " << ldvt << " " << endl;
    matrix_dump(rows,cols,original_matrix);
    cout << "!!!!!!!!!!!!!!!! Trying another SVD algorithm..." << endl;
    svd_with_lapack_fallback(rows,cols,original_matrix,u,vt,svals);
  }
  delete[] matrix;
  //cout << "Finished SVD" << endl;
}

void inverse_with_lapack(const dense_int lineardim, complex<double> *matrixin, complex<double> *matrixout){
  dense_int lda=lineardim;
  for (dense_int i=0;i<lineardim*lineardim;++i){
    matrixout[i]=matrixin[i]; //copy first
  }
  dense_int info=0;
  dense_int* ipiv=new dense_int[lineardim];

#ifndef NDEBUG
  std::cout << "Calling zgetri on matrix " << lineardim  << "*" << lineardim << std::endl;
#endif

#if defined(__APPLE__) && defined(__MACH__)
  zgetrf_(&lda,&lda,reinterpret_cast<__CLPK_doublecomplex*>(matrixout),&lda,ipiv,&info);
  if (info < 0){cout << "illegal argument to zgetrf, " << abs(info) << endl; exit(1);}
  if (info > 0) {cout << "singular U... " << abs(info) << endl; exit(1);}
  dense_int lwork=lda;
  info=0; //reset
  complex<double>* work = new complex<double>[lwork];
  zgetri_(&lda,reinterpret_cast<__CLPK_doublecomplex*>(matrixout),&lda,ipiv,reinterpret_cast<__CLPK_doublecomplex*>(work),&lwork,&info);
  delete[] work;
#elif defined(__linux__)
  set_num_threads(max_num_threads_to_use);
  info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,lda,lda,reinterpret_cast<lapack_complex_double*>(matrixout),lda,ipiv);
  if (info < 0){cout << "illegal argument to zgetrf, " << abs(info) << endl; /*square_matrix_print(lineardim,matrixin)*/; exit(1);}
  if (info > 0) {cout << "singular U... " << abs(info) << endl; exit(1);}
  info=LAPACKE_zgetri(LAPACK_COL_MAJOR,lda,reinterpret_cast<lapack_complex_double*>(matrixout),lda,ipiv);
#else  
#error Platform not supported
  exit(1);
#endif
  delete[] ipiv;
  if (info < 0){cout << "illegal argument to zgetri, " << abs(info) << endl; exit(1);}
  if (info > 0) {cout << "singular U... " << abs(info) << endl; exit(1);}
}

double dense_max_element(const dense_int dimA, const dense_int dimB, complex<double>* matrix){
  dense_int a=dimA;
  dense_int b=dimB;
#if defined(__APPLE__) && defined(__MACH__)
  char M[]={'M','\0'};
  double work[1];
  return(zlange_(M,&a,&b,reinterpret_cast<__CLPK_doublecomplex*>(matrix),&a,work));
#elif defined(__linux__)
  set_num_threads(max_num_threads_to_use);
#ifdef MKL_ILP64
  char M[]={'M','\0'};
  double work[1];
  return(zlange_(M,&a,&b,matrix,&a,work));
#else
  char M='M';
  return(LAPACKE_zlange(LAPACK_COL_MAJOR,M,a,b,reinterpret_cast<lapack_complex_double*>(matrix),a));
#endif
#else  
#error Platform not supported
  exit(1);
#endif
}
}

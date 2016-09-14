#ifndef DMAT_FUNC_H
#define DMAT_FUNC_H
#include <cstdint>
#include <complex>

using namespace std;

#if defined(MKL_ILP64) //MKL uses long long
#define dense_int long long
#elif defined(LAPACK_ILP64) //everyone else knows that ILP64 means long long=long
#define dense_int long
#else
#define dense_int int
#endif

namespace densefuncs {
  void gemm(dense_int M, dense_int K, dense_int N, complex<double> alpha, complex<double> beta, complex<double> *A, complex<double> *B, complex<double> *C);
  void square_gemv(dense_int M,complex<double> *A, complex<double> *x, complex<double> *y);
  void diagonalise_with_lapack(dense_int lineardim, complex<double> *matrixin, complex<double> *evecs, double *evals, dense_int il, dense_int iu);
  void diagonalise_with_lapack_nh(dense_int lineardim, complex<double> *matrixin, complex<double> *evals);
  void diagonalise_with_lapack_nh(dense_int lineardim, complex<double> *matrixin,  complex<double> *revecs, complex<double> *evals);
  void svd_with_lapack(dense_int rows, dense_int cols, complex<double> *matrix, complex<double> *u, complex<double> *vt, double *svals);
  void svd_with_lapack(dense_int rows, dense_int cols, complex<double> *matrix, double *svals);
  void svd_with_lapack_fallback(dense_int rows, dense_int cols, complex<double> *matrix, complex<double> *u, complex<double> *vt, double *svals);
  void inverse_with_lapack(dense_int lineardim, complex<double> *matrixin, complex<double> *matrixout);
  double dense_max_element(dense_int dimA, dense_int dimB, complex<double>* matrix);
}
#endif

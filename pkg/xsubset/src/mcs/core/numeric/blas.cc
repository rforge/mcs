/**
 * @file blas.cc
 */
#ifndef MCS_CORE_NUMERIC_BLAS_CC
#define MCS_CORE_NUMERIC_BLAS_CC


#include "../../mcs.hh"

#include "blas.hh"


extern "C" {

void MCS_F77_NAME(dcopy)(const int* n, const double* x, const int* incx,
			 double* y, const int* incy);

void MCS_F77_NAME(dgemv)(const char* trans, const int* m, const int* n,
			 const double* alpha, const double* a, const int* lda,
			 const double* x, const int* incx, const double* beta,
			 double* y, const int* incy);

void MCS_F77_NAME(drot)(const int* n, double* x, const int* incx, double* y,
			const int* incy, const double* c, const double* s);

void MCS_F77_NAME(drotg)(double* x, double* y, double* c, double *s);

}  // extern "C"


namespace mcs     {
namespace core    {
namespace numeric {


template<>
void blas<double>::copy(size_type n, const_reference_type x, size_type incx,
			reference_type y, size_type incy)
{
  MCS_ASSERT(n >= 0, "invalid argument: n (blas::copy)");
  MCS_ASSERT(incx > 0, "invalid argument: incx (blas::copy)");
  MCS_ASSERT(incy > 0, "invalid argument: incy (blas::copy)");

  MCS_F77_CALL(dcopy)(&n, &x, &incx, &y, &incy);
}

template<>
void blas<double>::gemv(const char* const trans, size_type m, size_type n,
			value_type alpha, const_reference_type a, size_type lda,
			const_reference_type x, size_type incx, value_type beta,
			reference_type y, size_type incy)
{
  MCS_ASSERT(m >= 0, "invalid argument: m (blas::copy)");
  MCS_ASSERT(n >= 0, "invalid argument: n (blas::copy)");
  MCS_ASSERT(lda > 0, "invalid argument: lda (blas::copy)");
  MCS_ASSERT(incx > 0, "invalid argument: incx (blas::copy)");
  MCS_ASSERT(incy > 0, "invalid argument: incy (blas::copy)");

  MCS_F77_CALL(dgemv)(trans, &m, &n, &alpha, &a, &lda, &x, &incx, &beta, &y, &incy);
}

template<>
void blas<double>::rotg(reference_type x, reference_type y,
			reference_type c, reference_type s)
{
  MCS_F77_CALL(drotg)(&x, &y, &c, &s);
}

template<>
void blas<double>::rot(size_type n, reference_type x, size_type incx, reference_type y,
		       size_type incy, value_type c, value_type s)
{
  MCS_ASSERT(n >= 0, "invalid argument: n (blas::rot)");
  MCS_ASSERT(incx > 0, "invalid argument: incx (blas::rot)");
  MCS_ASSERT(incy > 0, "invalid argument: incy (blas::rot)");

  MCS_F77_CALL(drot)(&n, &x, &incx, &y, &incy, &c, &s);
}


}  // numeric
}  // core
}  // mcs


#endif

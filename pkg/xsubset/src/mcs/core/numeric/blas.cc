/**
 * @file blas.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_BLAS_CC
#define MCS_CORE_NUMERIC_BLAS_CC


#include "../../mcs.hh"

#include "range.hh"
#include "vector.hh"
#include "matrix.hh"
#include "blas.hh"


#define RANGE range<Size>
#define VECTOR vector<Value, Alloc>
#define MATRIX matrix<Value, Alloc>


extern "C"
{


  void
  scopy_(const int* n,
	 const float* x, const int* incx,
	 float* y, const int* incy);


  void
  dcopy_(const int* n,
	 const double* x, const int* incx,
	 double* y, const int* incy);


  void
  sgemv_(const char* trans,
	 const int* m, const int* n,
	 const float* alpha,
	 const float* a, const int* lda,
	 const float* x, const int* incx,
	 const float* beta,
	 float* y, const int* incy);


  void
  dgemv_(const char* trans,
	 const int* m, const int* n,
	 const double* alpha,
	 const double* a, const int* lda,
	 const double* x, const int* incx,
	 const double* beta,
	 double* y, const int* incy);


  void
  srot_(const int* n,
	float* x, const int* incx,
	float* y, const int* incy,
	const float* c, const float* s);


  void
  drot_(const int* n,
	double* x, const int* incx,
	double* y, const int* incy,
	const double* c, const double* s);


  void
  srotg_(float* x, float* y,
	 float* c, float *s);


  void
  drotg_(double* x, double* y,
	 double* c, double *s);


}


namespace mcs
{

  namespace core
  {

    namespace numeric
    {

      namespace blas
      {


	void
	copy(const int n, const float& x, const int incx,
	     float& y, const int incy)
	{
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(incx > 0);
	  MCS_ASSERT(incy > 0);

	  scopy_(&n, &x, &incx, &y, &incy);
	}


	void
	copy(const int n, const double& x, const int incx,
	     double& y, const int incy)
	{
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(incx > 0);
	  MCS_ASSERT(incy > 0);

	  dcopy_(&n, &x, &incx, &y, &incy);
	}


	template<typename Value,
		 typename Alloc>
	void
	copy(const VECTOR& x, VECTOR&& y)
	{
	  MCS_ASSERT(x.len() == y.len());

	  copy(x.len(), x(0), x.inc(), y(0), y.inc());
	}


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	copy(const RANGE r, const VECTOR& x, VECTOR&& y)
	{
	  MCS_ASSERT(r.len() <= x.len());
	  MCS_ASSERT(r.len() <= y.len());

	  copy(x(r), y(r));
	}


	void
	gemv(const char* const trans, const int m, const int n, const float alpha,
	     const float& a, const int lda, const float& x, const int incx,
	     const float beta, float& y, const int incy)
	{
	  MCS_ASSERT(m >= 0);
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(lda > 0);
	  MCS_ASSERT(incx > 0);
	  MCS_ASSERT(incy > 0);

	  sgemv_(trans, &m, &n, &alpha, &a, &lda,
		 &x, &incx, &beta, &y, &incy);
	}


	void
	gemv(const char* const trans, const int m, const int n, const double alpha,
	     const double& a, const int lda, const double& x, const int incx,
	     const double beta, double& y, const int incy)
	{
	  MCS_ASSERT(m >= 0);
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(lda > 0);
	  MCS_ASSERT(incx > 0);
	  MCS_ASSERT(incy > 0);

	  dgemv_(trans, &m, &n, &alpha, &a, &lda,
		 &x, &incx, &beta, &y, &incy);
	}


	template<typename Value,
		 typename Alloc>
	void
	gemv(const char* const trans, const Value alpha, const MATRIX& a,
	     const VECTOR& x, const Value beta, VECTOR&& y)
	{
	  MCS_ASSERT(a.ncol() == x.len());
	  MCS_ASSERT(a.nrow() == y.len());

	  gemv(trans, a.nrow(), a.ncol(), alpha,
	       a(0, 0), a.ldim(), x(0), x.inc(),
	       beta, y(0), y.inc());
	}


	void
	rot(const int n, float& x, const int incx, float& y, const int incy,
	    const float c, const float s)
	{
	  MCS_ASSERT(n >= 0);

	  srot_(&n, &x, &incx, &y, &incy, &c, &s);
	}

	void
	rot(const int n, double& x, const int incx, double& y, const int incy,
	    const double c, const double s)
	{
	  MCS_ASSERT(n >= 0);

	  drot_(&n, &x, &incx, &y, &incy, &c, &s);
	}


	template<typename Value,
		 typename Alloc>
	void
	rot(VECTOR&& x, VECTOR&& y, const Value c, const Value s)
	{
	  MCS_ASSERT(x.len() == y.len());

	  rot(x.len(), x(0), x.inc(), y(0), y.inc(), c, s);
	}


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	rot(const Size i, VECTOR&& x, VECTOR&& y)
	{
	  MCS_ASSERT(x.len() == y.len());

	  const RANGE r(i + 1);
	  Value c, s;

	  rotg(x(i), y(i), c, s);
	  rot(x(r), y(r), c, s);
	}


	template<typename Value,
		 typename Alloc>
	void
	rot(VECTOR&& x, VECTOR&& y)
	{
	  typedef typename VECTOR::size_type size_type;

	  MCS_ASSERT(x.len() == y.len());

	  rot(size_type(0), x, y);
	}


	void
	rotg(float& x, float& y, float& c, float& s)
	{
	  srotg_(&x, &y, &c, &s);
	}


	void
	rotg(double& x, double& y, double& c, double& s)
	{
	  drotg_(&x, &y, &c, &s);
	}


      }

    }

  }

}


#undef RANGE
#undef VECTOR
#undef MATRIX


#endif

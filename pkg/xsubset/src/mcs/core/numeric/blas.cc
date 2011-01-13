/**
 * @file blas.cc
 */
#ifndef MCS_CORE_NUMERIC_BLAS_CC
#define MCS_CORE_NUMERIC_BLAS_CC


#include "../../mcs.hh"

#include "slice.hh"
#include "vector.hh"
#include "matrix.hh"
#include "blas.hh"


extern "C"
{


  // FIXME: single precision routines fail to link
  //        with RLapack
  // void
  // MCS_F77_NAME(scopy)(const int* n,
  //                     const float* x, const int* incx,
  //                     float* y, const int* incy);


  void
  MCS_F77_NAME(dcopy)(const int* n,
                      const double* x, const int* incx,
                      double* y, const int* incy);


  // FIXME: single precision routines fail to link
  //        with RLapack
  // void
  // MCS_F77_NAME(sgemv)(const char* trans,
  //                     const int* m, const int* n,
  //                     const float* alpha,
  //                     const float* a, const int* lda,
  //                     const float* x, const int* incx,
  //                     const float* beta,
  //                     float* y, const int* incy);


  void
  MCS_F77_NAME(dgemv)(const char* trans,
                      const int* m, const int* n,
                      const double* alpha,
                      const double* a, const int* lda,
                      const double* x, const int* incx,
                      const double* beta,
                      double* y, const int* incy);


  // FIXME: single precision routines fail to link
  //        with RLapack
  // void
  // MCS_F77_NAME(srot)(const int* n,
  //                    float* x, const int* incx,
  //                    float* y, const int* incy,
  //                    const float* c, const float* s);


  void
  MCS_F77_NAME(drot)(const int* n,
                     double* x, const int* incx,
                     double* y, const int* incy,
                     const double* c, const double* s);


  // FIXME: single precision routines fail to link
  //        with RLapack
  // void
  // MCS_F77_NAME(srotg)(float* x, float* y,
  //                     float* c, float *s);


  void
  MCS_F77_NAME(drotg)(double* x, double* y,
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
	copy(const int n,
             const float& x, const int incx,
	     float& y, const int incy)
	{
          MCS_ASSERT(n >= 0, "invalid argument (blas::copy)");
          MCS_ASSERT(incx > 0, "invalid argument (blas::copy)");
          MCS_ASSERT(incy > 0, "invalid argument (blas::copy)");

          // FIXME: single precision routines fail to link
          //        with RLapack
          // MCS_F77_CALL(scopy)(&n, &x, &incx, &y, &incy);
          MCS_ASSERT(false, "routine unavailable (scopy)");
	}


	void
	copy(const int n,
             const double& x, const int incx,
	     double& y, const int incy)
	{
          MCS_ASSERT(n >= 0, "invalid argument (blas::copy)");
          MCS_ASSERT(incx > 0, "invalid argument (blas::copy)");
          MCS_ASSERT(incy > 0, "invalid argument (blas::copy)");

	  MCS_F77_CALL(dcopy)(&n, &x, &incx, &y, &incy);
	}


        

	template<typename Value,
                 typename Size,
		 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
        copy(const vector_base<Value, Size, Derived1>& x,
             vector_base<Value, Size, Derived2>&& y)
	{
	  MCS_ASSERT(x.len() == y.len(), "invalid argument (blas::copy)");

	  copy(x.len(), x(0), x.inc(), y(0), y.inc());
	}


	template<typename Value,
		 typename Size,
		 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	copy(const slice<Size>& s,
             const vector_base<Value, Size, Derived1>& x,
             vector_base<Value, Size, Derived2>&& y)
	{
	  MCS_ASSERT(s.pos() + s.len() <= x.len(),
                     "invalid argument (blas::copy)");
	  MCS_ASSERT(s.pos() + s.len() <= y.len(),
                     "invalid argument (blas::copy)");

	  copy(x(s), y(s));
	}


	void
	gemv(const char* const trans,
             const int m, const int n,
             const float alpha, const float& a, const int lda,
             const float& x, const int incx,
	     const float beta, float& y, const int incy)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (blas::gemv)");
	  MCS_ASSERT(n >= 0, "invalid argument (blas::gemv)");
	  MCS_ASSERT(lda > 0, "invalid argument (blas::gemv)");
	  MCS_ASSERT(incx > 0, "invalid argument (blas::gemv)");
	  MCS_ASSERT(incy > 0, "invalid argument (blas::gemv)");

          // FIXME: single precision routines fail to link
          //        with RLapack
	  // MCS_F77_CALL(sgemv)(trans, &m, &n, &alpha, &a, &lda,
          //                     &x, &incx, &beta, &y, &incy);
          MCS_ASSERT(false, "routine unavailable (sgemv)");
	}


	void
	gemv(const char* const trans,
             const int m, const int n,
             const double alpha, const double& a, const int lda,
             const double& x, const int incx,
	     const double beta, double& y, const int incy)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (blas::copy)");
	  MCS_ASSERT(n >= 0, "invalid argument (blas::copy)");
	  MCS_ASSERT(lda > 0, "invalid argument (blas::copy)");
	  MCS_ASSERT(incx > 0, "invalid argument (blas::copy)");
	  MCS_ASSERT(incy > 0, "invalid argument (blas::copy)");

	  MCS_F77_CALL(dgemv)(trans, &m, &n, &alpha, &a, &lda,
                              &x, &incx, &beta, &y, &incy);
	}


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2,
                 template<typename V,
                          typename S>
                 class Derived3>
	void
	gemv(const char* trans,
             Value alpha, const matrix_base<Value, Size, Derived1>& a,
	     const vector_base<Value, Size, Derived2>& x,
             Value beta,  vector_base<Value, Size, Derived3>&& y)
	{
	  MCS_ASSERT(a.ncol() == x.len(), "invalid argument (blas::copy)");
	  MCS_ASSERT(a.nrow() == y.len(), "invalid argument (blas::copy)");

	  gemv(trans, a.nrow(), a.ncol(), alpha, a(0, 0), a.ldim(),
               x(0), x.inc(), beta, y(0), y.inc());
	}


	void
	rotg(float& x, float& y,
             float& c, float& s)
	{
          // FIXME: single precision routines fail to link
          //        with RLapack
          // MCS_F77_CALL(srotg)(&x, &y, &c, &s);
          MCS_ASSERT(false, "routine unavailable (srotg)");
	}


	void
	rotg(double& x, double& y,
             double& c, double& s)
	{
	  MCS_F77_CALL(drotg)(&x, &y, &c, &s);
	}


	void
	rot(const int n,
            float& x, const int incx,
            float& y, const int incy,
	    const float c, const float s)
	{
	  MCS_ASSERT(n >= 0, "invalid argument (blas::rot)");
	  MCS_ASSERT(incx > 0, "invalid argument (blas::rot)");
	  MCS_ASSERT(incy > 0, "invalid argument (blas::rot)");

          // FIXME: single precision routines fail to link
          //        with RLapack
          // MCS_F77_CALL(srot)(&n, &x, &incx, &y, &incy, &c, &s);
          MCS_ASSERT(false, "routine unavailable (srot)");
	}

	void
	rot(const int n,
            double& x, const int incx,
            double& y, const int incy,
	    const double c, const double s)
	{
	  MCS_ASSERT(n >= 0, "invalid argument (blas::rot)");
	  MCS_ASSERT(incx > 0, "invalid argument (blas::rot)");
	  MCS_ASSERT(incy > 0, "invalid argument (blas::rot)");

	  MCS_F77_CALL(drot)(&n, &x, &incx, &y, &incy, &c, &s);
	}


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	rot(vector_base<Value, Size, Derived1>&& x,
            vector_base<Value, Size, Derived2>&& y,
            Value c, Value s)
	{
	  MCS_ASSERT(x.len() == y.len(), "invalid argument (blas::rot)");

	  rot(x.len(), x(0), x.inc(), y(0), y.inc(), c, s);
	}


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	rot(vector_base<Value, Size, Derived1>&& x,
            vector_base<Value, Size, Derived2>&& y)
	{
	  MCS_ASSERT(x.len() == y.len(), "invalid argument (blas::rot)");

	  rot(0, x, y);
	}


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	rot(const Size i,
            vector_base<Value, Size, Derived1>&& x,
            vector_base<Value, Size, Derived2>&& y)
	{
	  MCS_ASSERT(x.len() == y.len(), "invalid argument (blas::rot)");

          const Size n = x.len();
	  Value c, s;

	  rotg(x(i), y(i), c, s);
	  rot(x({i + 1, n}), y({i + 1, n}), c, s);
	}


      }

    }

  }

}


#endif

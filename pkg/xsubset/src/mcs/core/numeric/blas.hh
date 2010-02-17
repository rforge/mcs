/**
 * @file blas.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_BLAS_HH
#define MCS_CORE_NUMERIC_BLAS_HH


#define RANGE range<Size>
#define VECTOR vector<Value, Alloc>
#define MATRIX matrix<Value, Alloc>


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      template<typename Size>
      class range;


      template<typename Value,
	       typename Alloc>
      class vector;


      template<typename Value,
	       typename Alloc>
      class matrix;


      namespace blas
      {


	void
	copy(int n, const float& x, int incx,
	     float& y, int incy);


	void
	copy(int n, const double& x, int incx,
	     double& y, int incy);


	template<typename Value,
		 typename Alloc>
	void
	copy(const VECTOR& x, VECTOR&& y);


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	copy(RANGE r, const VECTOR& x, VECTOR&& y);


	void
	gemv(const char* trans, int m, int n, float alpha,
	     const float& a, int lda, const float& x, int incx,
	     float beta, float& y, int incy);


	void
	gemv(const char* trans, int m, int n, double alpha,
	     const double& a, int lda, const double& x, int incx,
	     double beta, double& y, int incy);


	template<typename Value,
		 typename Alloc>
	void
	gemv(const char* trans, Value alpha, const MATRIX& a,
	     const VECTOR& x, Value beta, VECTOR&& y);


	void
	rot(int n, float& x, int incx, float& y, int incy,
	    float c, float s);


	void
	rot(int n, double& x, int incx, double& y, int incy,
	    double c, double s);


	template<typename Value,
		 typename Alloc>
	void
	rot(VECTOR&& x, VECTOR&& y, Value c, Value s);


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	rot(Size i, VECTOR&& x, VECTOR&& y);


	template<typename Value,
		 typename Alloc>
	void
	rot(VECTOR&& x, VECTOR&& y);


	void
	rotg(float& x, float& y, float& c, float& s);


	void
	rotg(double& x, double& y, double& c, double& s);


      }

    }

  }

}


#undef RANGE
#undef VECTOR
#undef MATRIX


#include "blas.cc"
#endif

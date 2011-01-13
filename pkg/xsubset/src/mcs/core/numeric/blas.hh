/**
 * @file blas.hh
 */
#ifndef MCS_CORE_NUMERIC_BLAS_HH
#define MCS_CORE_NUMERIC_BLAS_HH


#include "slice.hh"
#include "vector.hh"
#include "matrix.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {


      namespace blas
      {


	void
	copy(int n,
             const float& x, int incx,
	     float& y, int incy);


	void
	copy(int n,
             const double& x, int incx,
	     double& y, int incy);
        

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
             vector_base<Value, Size, Derived2>&& y);


	template<typename Value,
		 typename Size,
		 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	copy(const slice<Size>& r,
             const vector_base<Value, Size, Derived1>& x,
             vector_base<Value, Size, Derived2>&& y);


	void
	gemv(const char* trans,
             int m, int n,
             float alpha, const float& a, int lda,
             const float& x, int incx,
	     float beta, float& y, int incy);


	void
	gemv(const char* trans,
             int m, int n,
             double alpha,  const double& a, int lda,
             const double& x, int incx,
	     double beta, double& y, int incy);


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
             Value beta,  vector_base<Value, Size, Derived3>&& y);


	void
	rotg(float& x, float& y,
             float& c, float& s);


	void
	rotg(double& x, double& y,
             double& c, double& s);


	void
	rot(int n,
            float& x, int incx,
            float& y, int incy,
	    float c, float s);


	void
	rot(int n,
            double& x, int incx,
            double& y, int incy,
	    double c, double s);


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
            Value c, Value s);


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
            vector_base<Value, Size, Derived2>&& y);


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	rot(Size i,
            vector_base<Value, Size, Derived1>&& x,
            vector_base<Value, Size, Derived2>&& y);


      }

    }

  }

}


#include "blas.cc"
#endif

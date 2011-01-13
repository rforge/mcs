/**
 * @file lapack.hh
 */
#ifndef MCS_CORE_NUMERIC_LAPACK_HH
#define MCS_CORE_NUMERIC_LAPACK_HH


#include "slice.hh"
#include "vector.hh"
#include "matrix.hh"


namespace mcs
{

  namespace core
  {

    namespace numeric
    {

      namespace lapack
      {


	void
	lacpy(const char* uplo, int m, int n,
	      const float& a, int lda,
	      float& b, int ldb);


	void
	lacpy(const char* uplo, int m, int n,
	      const double& a, int lda,
	      double& b, int ldb);
        
        
	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
        void
        lacpy(const char* uplo,
              const matrix_base<Value, Size, Derived1>& a,
              matrix_base<Value, Size, Derived2>&& b);


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived1,
                 template<typename V,
                          typename S>
                 class Derived2>
	void
	lacpy(const char* uplo,
              const slice<Size>& s1, const slice<Size>& s2,
	      const matrix_base<Value, Size, Derived1>& a,
              matrix_base<Value, Size, Derived2>&& b);


	void
	geqrf(int m, int n,
              float& a, int lda, float& tau,
	      float& work, int lwork);


	void
	geqrf(int m, int n,
              double& a, int lda, double& tau,
	      double& work, int lwork);


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived>
        void
        geqrf(matrix_base<Value, Size, Derived>&& a);


	template<typename Value,
		 typename Size,
		 template<typename V,
                          typename S>
                 class Derived>
	void
        geqrf(const slice<Size>& s,
              matrix_base<Value, Size, Derived>&& a);


      }

    }

  }

}


#include "lapack.cc"
#endif

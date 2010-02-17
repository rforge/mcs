/**
 * @file lapack.hh
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_LAPACK_HH
#define MCS_CORE_NUMERIC_LAPACK_HH


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
		 typename Alloc>
	void
	lacpy(const char* uplo, const MATRIX& a, MATRIX&& b);


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	lacpy(const char* uplo, RANGE r1, RANGE r2,
	      const MATRIX& a, MATRIX&& b);


	void
	geqrf(int m, int n, float& a, int lda, float& tau,
	      float& work, int lwork);


	void
	geqrf(int m, int n, double& a, int lda, double& tau,
	      double& work, int lwork);


	template<typename Value,
		 typename Alloc>
	void
	geqrf(MATRIX&& a);


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	geqrf(RANGE r, MATRIX&& a);


      }

    }

  }

}


#undef RANGE
#undef VECTOR
#undef MATRIX


#include "lapack.cc"
#endif

/**
 * @file lapack.cc
 *
 * @author Marc Hofmann
 */

#ifndef MCS_CORE_NUMERIC_LAPACK_CC
#define MCS_CORE_NUMERIC_LAPACK_CC


#include <algorithm>

#include "../../mcs.hh"

#include "range.hh"
#include "vector.hh"
#include "matrix.hh"
#include "lapack.hh"


#define RANGE range<Size>
#define VECTOR vector<Value, Alloc>
#define MATRIX matrix<Value, Alloc>


extern "C"
{


  // FIXME: single precision routines fail to link
  //   with RLapack
  // void
  // slacpy_(const char* uplo,
  //         const int* m, const int* n,
  //         const float* a, const int* lda,
  //         float* b, const int* ldb);


  void
  dlacpy_(const char* uplo,
	  const int* m, const int* n,
	  const double* a, const int* lda,
	  double* b, const int* ldb);


  // FIXME: sgeqrf_ fails to link
  // void
  // sgeqrf_(const int* m, const int* n,
  //         float* a, const int* lda,
  //         float* tau,
  //         float* work, const int* lwork,
  //         int* info);


  void
  dgeqrf_(const int* m, const int* n,
	  double* a, const int* lda,
	  double* tau,
	  double* work, const int* lwork,
	  int* info);


}


namespace mcs
{

  namespace core
  {

    namespace numeric
    {

      namespace lapack
      {

	namespace detail
	{


	  int info;


	}


        // FIXME: single precision routines fail to link
        //   with RLapack
	// void
	// lacpy(const char* uplo, const int m, const int n,
	//       const float& a, const int lda,
	//       float& b, const int ldb)
	// {
	//   MCS_ASSERT(m >= 0);
	//   MCS_ASSERT(n >= 0);
	//   MCS_ASSERT(lda >= std::max(1, m));
	//   MCS_ASSERT(ldb >= std::max(1, m));

	//   slacpy_(uplo, &m, &n, &a, &lda, &b, &ldb);
	// }


	void
	lacpy(const char* uplo, const int m, const int n,
	      const double& a, const int lda,
	      double& b, const int ldb)
	{
	  MCS_ASSERT(m >= 0);
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(lda >= std::max(1, m));
	  MCS_ASSERT(ldb >= std::max(1, m));

	  dlacpy_(uplo, &m, &n, &a, &lda, &b, &ldb);
	}


	template<typename Value,
		 typename Alloc>
	void
	lacpy(const char* uplo, const MATRIX& a, MATRIX&& b)
	{
	  MCS_ASSERT(a.nrow() == b.nrow());
	  MCS_ASSERT(a.ncol() == b.ncol());

	  lacpy(uplo, a.nrow(), a.ncol(), a(0, 0), a.ldim(),
		b(0, 0), b.ldim());
	}


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	lacpy(const char* uplo, const RANGE r1, const RANGE r2,
	      const MATRIX& a, MATRIX&& b)
	{
	  lacpy(uplo, a(r1, r2), b(r1, r2));
	}


        // FIXME: sgeqrf_ fails to link
	// void
	// geqrf(const int m, const int n, float& a, const int lda, float& tau,
	//       float& work, const int lwork)
	// {
	//   MCS_ASSERT(m >= 0);
	//   MCS_ASSERT(n >= 0);
	//   MCS_ASSERT(lda >= std::max(1, m));
	//   MCS_ASSERT(lwork >= std::max(1, n));

	//   sgeqrf_(&m, &n, &a, &lda, &tau, &work, &lwork, &detail::info);
	// }


	void
	geqrf(const int m, const int n, double& a, const int lda, double& tau,
	      double& work, const int lwork)
	{
	  MCS_ASSERT(m >= 0);
	  MCS_ASSERT(n >= 0);
	  MCS_ASSERT(lda >= std::max(1, m));
	  MCS_ASSERT(lwork >= std::max(1, n));

	  dgeqrf_(&m, &n, &a, &lda, &tau, &work, &lwork, &detail::info);
	}


	template<typename Value,
		 typename Alloc>
	void
	geqrf(MATRIX&& a)
	{
	  typedef typename MATRIX::size_type size_type;

	  const size_type m = a.nrow();
	  const size_type n = a.ncol();
	  Value tau[std::min(m, n)];
	  size_type lwork = std::max(size_type(1), n);
	  Value work[lwork];

	  geqrf(m, n, a(0, 0), a.ldim(), tau[0], work[0], lwork);
	}


	template<typename Size,
		 typename Value,
		 typename Alloc>
	void
	geqrf(const RANGE r, MATRIX& a)
	{
	  geqrf(a(r, r));
	}


      }

    }

  }

}


#undef RANGE
#undef VECTOR
#undef MATRIX


#endif

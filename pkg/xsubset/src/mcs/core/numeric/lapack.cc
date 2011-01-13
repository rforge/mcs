/**
 * @file lapack.cc
 */
#ifndef MCS_CORE_NUMERIC_LAPACK_CC
#define MCS_CORE_NUMERIC_LAPACK_CC


#include <algorithm>

#include "../../mcs.hh"

#include "slice.hh"
#include "vector.hh"
#include "matrix.hh"
#include "lapack.hh"


extern "C"
{


  // FIXME: single precision routines fail to link
  //        with RLapack
  // void
  // MCS_F77_NAME(slacpy)(const char* uplo,
  //                      const int* m, const int* n,
  //                      const float* a, const int* lda,
  //                      float* b, const int* ldb);


  void
  MCS_F77_NAME(dlacpy)(const char* uplo,
                       const int* m, const int* n,
                       const double* a, const int* lda,
                       double* b, const int* ldb);


  // FIXME: single precision routines fail to link
  //        with RLapack
  // MCS_F77_NAME(sgeqrf(const int* m, const int* n,
  //                     float* a, const int* lda, float* tau,
  //                     float* work, const int* lwork, int* info);


  void
  MCS_F77_NAME(dgeqrf)(const int* m, const int* n,
                       double* a, const int* lda, double* tau,
                       double* work, const int* lwork, int* info);


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


	void
	lacpy(const char* uplo, const int m, const int n,
	      const float& a, const int lda,
	      float& b, const int ldb)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(n >= 0, "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(lda >= std::max(1, m),
                     "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(ldb >= std::max(1, m),
                     "invalid argument (lapack::lacpy)");

          // FIXME: single precision routines fail to link
          //        with RLapack
          // MCS_F77_CALL(slacpy)(uplo, &m, &n, &a, &lda, &b, &ldb);
          MCS_ASSERT(false, "routine unavailable (slacpy)");
	}


	void
	lacpy(const char* uplo, const int m, const int n,
	      const double& a, const int lda,
	      double& b, const int ldb)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(n >= 0, "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(lda >= std::max(1, m),
                     "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(ldb >= std::max(1, m),
                     "invalid argument (lapack::lacpy)");

	  MCS_F77_CALL(dlacpy)(uplo, &m, &n, &a, &lda, &b, &ldb);
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
        lacpy(const char* uplo,
              const matrix_base<Value, Size, Derived1>& a,
              matrix_base<Value, Size, Derived2>&& b)
	{
	  MCS_ASSERT(a.nrow() == b.nrow(), "invalid argument (lapack::lacpy)");
	  MCS_ASSERT(a.ncol() == b.ncol(), "invalid argument (lapack::lacpy)");

	  lacpy(uplo, a.nrow(), a.ncol(), a(0, 0), a.ldim(),
		b(0, 0), b.ldim());
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
	lacpy(const char* uplo,
              const slice<Size>& sr, const slice<Size>& sc,
	      const matrix_base<Value, Size, Derived1>& a,
              matrix_base<Value, Size, Derived2>&& b)
	{
          MCS_ASSERT(sr.pos() + sr.len() <= a.nrow(),
                     "invalid argument (lapack::lacpy)");
          MCS_ASSERT(sc.pos() + sc.len() <= a.ncol(),
                     "invalid argument (lapack::lacpy)");
          MCS_ASSERT(sr.pos() + sr.len() <= b.nrow(),
                     "invalid argument (lapack::lacpy)");
          MCS_ASSERT(sc.pos() + sc.len() <= b.ncol(),
                     "invalid argument (lapack::lacpy)");

	  lacpy(uplo, a(sr, sc), b(sr, sc));
	}


	void
	geqrf(const int m, const int n,
              float& a, const int lda, float& tau,
	      float& work, const int lwork)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (lapack::geqrf)");
	  MCS_ASSERT(n >= 0, "invalid argument (lapack::geqrf)");
	  MCS_ASSERT(lda >= std::max(1, m),
                     "invalid argument(lapack::geqrf)");
	  MCS_ASSERT(lwork >= std::max(1, n),
                     "invalid argument(lapack::geqrf)");

          // FIXME: single precision routines fail to link
          //        with RLapack
          // MCS_F77_CALL(sgeqrf)(&m, &n, &a, &lda, &tau,
          //                      &work, &lwork, &detail::info);
          MCS_ASSERT(false, "routine unavailable (sgeqrf)");
        }


	void
	geqrf(const int m, const int n,
              double& a, const int lda, double& tau,
	      double& work, const int lwork)
	{
	  MCS_ASSERT(m >= 0, "invalid argument (lapack::geqrf)");
	  MCS_ASSERT(n >= 0, "invalid argument (lapack::geqrf)");
	  MCS_ASSERT(lda >= std::max(1, m),
                     "invalid argument(lapack::geqrf)");
	  MCS_ASSERT(lwork >= std::max(1, n),
                     "invalid argument(lapack::geqrf)");

	  MCS_F77_CALL(dgeqrf)(&m, &n, &a, &lda, &tau,
                               &work, &lwork, &detail::info);
	}


	template<typename Value,
		 typename Size,
                 template<typename V,
                          typename S>
                 class Derived>
	void
	geqrf(matrix_base<Value, Size, Derived>&& a)
	{
	  const Size m = a.nrow();
	  const Size n = a.ncol();
	  Value tau[std::min(m, n)];
	  const Size lwork = std::max(Size(1), n);
	  Value work[lwork];

	  geqrf(m, n, a(0, 0), a.ldim(), tau[0], work[0], lwork);
	}


	template<typename Value,
		 typename Size,
		 template<typename V,
                          typename S>
                 class Derived>
	void
        geqrf(const slice<Size>& s,
              matrix_base<Value, Size, Derived>&& a)
	{
          MCS_ASSERT(s.pos() + s.len() <= a.nrow(),
                     "invalid argument (lapack::geqrf)");
          MCS_ASSERT(s.pos() + s.len() <= a.ncol(),
                     "invalid argument (lapack::geqrf)");

	  geqrf(a(s, s));
	}


      }

    }

  }

}


#endif

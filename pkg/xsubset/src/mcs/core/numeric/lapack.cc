/**
 * @file lapack.cc
 */
#ifndef MCS_CORE_NUMERIC_LAPACK_CC
#define MCS_CORE_NUMERIC_LAPACK_CC


#include "../../mcs.hh"

#include "lapack.hh"


extern "C" {

void MCS_F77_NAME(dlacpy)(const char* uplo, const int* m, const int* n,
			  const double* a, const int* lda, double* b, const int* ldb);

void MCS_F77_NAME(dgeqrf)(const int* m, const int* n, double* a, const int* lda, double* tau,
			  double* work, const int* lwork, int* info);


}


namespace mcs     {
namespace core    {
namespace numeric {


template<>
void lapack<double>::lacpy(const char* uplo, const size_type m, const size_type n,
			   const_reference a, const size_type lda, reference b,
			   const size_type ldb)
{
  MCS_ASSERT(m >= 0, "invalid argument: m (lapack::lacpy)");
  MCS_ASSERT(n >= 0, "invalid argument: n (lapack::lacpy)");
  MCS_ASSERT(lda >= std::max(1, m), "invalid argument: lda (lapack::lacpy)");
  MCS_ASSERT(ldb >= std::max(1, m), "invalid argument: ldb (lapack::lacpy)");

  MCS_F77_CALL(dlacpy)(uplo, &m, &n, &a, &lda, &b, &ldb);
}

template<>
void lapack<double>::geqrf(const size_type m, const size_type n, reference a,
			   const size_type lda, reference tau, reference work,
			   const size_type lwork)
{
  MCS_ASSERT(m >= 0, "invalid argument: m (lapack::geqrf)");
  MCS_ASSERT(n >= 0, "invalid argument: n (lapack::geqrf)");
  MCS_ASSERT(lda >= std::max(1, m), "invalid argument: lda (lapack::geqrf)");
  MCS_ASSERT(lwork >= std::max(1, n), "invalid argument: lwork (lapack::geqrf)");

  MCS_F77_CALL(dgeqrf)(&m, &n, &a, &lda, &tau,  &work, &lwork, &info);
}


}  // namespace numeric
}  // namespace core
}  // namespace mcs


#endif

#ifndef MCS_SUBSET_DETAIL_LAPACK_CC
#define MCS_SUBSET_DETAIL_LAPACK_CC


#include "mcs/subset/detail/lapack.hh"

#include "mcs/detail/fortran.hh" // MCS_F77_NAME, MCS_F77_CALL

#include <algorithm>  // std::min, std::max



extern "C" {

  void MCS_F77_NAME(dlacpy)(const char* uplo, const int* m, const int* n,
                            const double* a, const int* lda, double* b, const int* ldb);

  void MCS_F77_NAME(dgeqrf)(const int* m, const int* n, double* a, const int* lda, double* tau,
                            double* work, const int* lwork, int* info);


}  // extern "C"



namespace mcs    {
namespace subset {
namespace detail {



  template<>
  int
  Lapack<double>::info = 0;


  template<>
  void
  Lapack<double>::lacpy(const int m, const int n,
                        const double* const a, const int lda,
                              double* const b, const int ldb)
  {
    MCS_F77_CALL(dlacpy)("ALL", &m, &n, a, &lda, b, &ldb);
  }


  template<>
  void
  Lapack<double>::lacpy(const char* uplo, const int m, const int n,
                        const double* const a, const int lda,
                              double* const b, const int ldb)
  {
    MCS_F77_CALL(dlacpy)(uplo, &m, &n, a, &lda, b, &ldb);
  }


  template<>
  void
  Lapack<double>::geqrf(const int m, const int n,
                        double* const a, const int lda, double* const tau,
                        double* const work, const int lwork)
  {
    MCS_F77_CALL(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
  }


  template<>
  void
  Lapack<double>::geqrf(const int m, const int n, double* const a, const int lda)
  {
    const int lwork = std::max(1, n);
    const int ltau = std::min(m, n);

    double tau[ltau];
    double work[lwork];

    geqrf(m, n, a, lda, tau, work, lwork);
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

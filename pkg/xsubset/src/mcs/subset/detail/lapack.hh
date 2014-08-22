#ifndef MCS_SUBSET_DETAIL_LAPACK_HH
#define MCS_SUBSET_DETAIL_LAPACK_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  struct Lapack
  {

    static
    int info;


    static
    void
    lacpy(int m, int n,
          const TReal* a, int lda,
                TReal* b, int ldb);


    static
    void
    lacpy(const char* uplo, int m, int n,
          const TReal* a, int lda,
                TReal* b, int ldb);


    static
    void
    geqrf(int m, int n, TReal* a, int lda, TReal* tau,
          TReal* work, int lwork);


    static
    void
    geqrf(int m, int n, TReal* a, int lda);

  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/lapack.cc"
#endif

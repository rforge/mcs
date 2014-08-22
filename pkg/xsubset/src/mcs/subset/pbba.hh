#ifndef MCS_SUBSET_PBBA_HH
#define MCS_SUBSET_PBBA_HH


#include "mcs/subset/detail/pbba.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  pbba(const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const index, TReal* const rss, int* const s,
       const TCrit<TReal>& c)
  {
    detail::Preorder::Single<TReal> p(size, pmin);

    return detail::pbba(m, size, mark, nbest, v, ay, lday, index, rss, s, c, p);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  pbba(const int size, const int mark, const int nbest, const int pmin,
      const int* const v, const TReal* const rz, const int ldrz,
      int* const index, TReal* const rss, int* const s,
      const TCrit<TReal>& c)
  {
    detail::Preorder::Single<TReal> p(size, pmin);

    return detail::pbba(size, mark, nbest, v, rz, ldrz, index, rss, s, c, p);
  }


  template<typename TReal>
  int
  pbba(const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const index, TReal* const rss, int* const s)
  {
    detail::Preorder::Single<TReal> p(size, pmin);

    return detail::pbba(m, size, mark, nbest, v, ay, lday, index, rss, s, p);
  }


  template<typename TReal>
  int
  pbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const index, TReal* const rss, int* const s)
  {
    detail::Preorder::Single<TReal> p(size, pmin);

    return detail::pbba(size, mark, nbest, v, rz, ldrz, index, rss, s, p);
  }



}  // namespace subset
}  // namespace mcs



#endif

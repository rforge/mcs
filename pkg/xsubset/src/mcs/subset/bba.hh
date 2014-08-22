#ifndef MCS_SUBSET_BBA_HH
#define MCS_SUBSET_BBA_HH


#include "mcs/subset/detail/bba.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  bba(const int m, const int size, const int mark, const int nbest,
      const int* const v, const TReal* const ay, const int lday,
      int* const index, TReal* const rss, int* const s,
      const TCrit<TReal>& c)
  {
    return detail::bba(m, size, mark, nbest, v, ay, lday, index, rss, s, c);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  bba(const int size, const int mark, const int nbest,
      const int* const v, const TReal* const rz, const int ldrz,
      int* const index, TReal* const rss, int* const s,
      const TCrit<TReal>& c)
  {
    return detail::bba(size, mark, nbest, v, rz, ldrz, index, rss, s, c);
  }


  template<typename TReal>
  int
  bba(const int m, const int size, const int mark, const int nbest,
      const int* const v, const TReal* const ay, const int lday,
      int* const index, TReal* const rss, int* const s)
  {
    return detail::bba(m, size, mark, nbest, v, ay, lday, index, rss, s);
  }


  template<typename TReal>
  int
  bba(const int size, const int mark, const int nbest,
      const int* const v, const TReal* const rz, const int ldrz,
      int* const index, TReal* const rss, int* const s)
  {
    return detail::bba(size, mark, nbest, v, rz, ldrz, index, rss, s);
  }



}  // namespace subset
}  // namespace mcs



#endif

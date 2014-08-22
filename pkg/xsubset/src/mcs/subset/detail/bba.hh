#ifndef MCS_SUBSET_DETAIL_BBA_HH
#define MCS_SUBSET_DETAIL_BBA_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  bba(int m, int size, int mark, int nbest,
      const int* v, const TReal* ay, int lday,
      int* index, TReal* crit, int* s,
      const TCrit<TReal>& c);


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  bba(int size, int mark, int nbest,
      const int* v, const TReal* rz, int ldrz,
      int* index, TReal* crit, int* s,
      const TCrit<TReal>& c);


  template<typename TReal>
  int
  bba(int m, int size, int mark, int nbest,
      const int* v, const TReal* ay, int lday,
      int* index, TReal* rss, int* s);


  template<typename TReal>
  int
  bba(int size, int mark, int nbest,
      const int* v, const TReal* rz, int ldrz,
      int* index, TReal* rss, int* s);



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/bba.cc"
#endif

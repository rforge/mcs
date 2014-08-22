#ifndef MCS_SUBSET_DETAIL_PBBA_HH
#define MCS_SUBSET_DETAIL_PBBA_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCrit,
    template<typename R>
    class TPreorder
  >
  int
  pbba(int m, int size, int mark, int nbest,
       const int* v, const TReal* ay, int lday,
       int* index, TReal* crit, int* s,
       const TCrit<TReal>& c, const TPreorder<TReal>& p);


  template<
    typename TReal,
    template<typename R>
    class TCrit,
    template<typename R>
    class TPreorder
  >
  int
  pbba(int size, int mark, int nbest,
       const int* v, const TReal* rz, int ldrz,
       int* index, TReal* crit, int* s,
       const TCrit<TReal>& c, const TPreorder<TReal>& p);


  template<
    typename TReal,
    template<typename R>
    class TPreorder
  >
  int
  pbba(int m, int size, int mark, int nbest,
       const int* v, const TReal* ay, int lday,
       int* index, TReal* rss, int* s,
       const TPreorder<TReal>& p);


  template<
    typename TReal,
    template<typename R>
    class TPreorder
  >
  int
  pbba(int size, int mark, int nbest,
       const int* v, const TReal* rz, int ldrz,
       int* index, TReal* rss, int* s,
       const TPreorder<TReal>& p);



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/pbba.cc"
#endif

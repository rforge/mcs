#ifndef MCS_SUBSET_DCA_HH
#define MCS_SUBSET_DCA_HH


#include "mcs/subset/detail/dca.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int dca(const int m, const int size, const int mark, const int nbest,
          const int* const v, const TReal* const ay, const int lday,
	  int* const index, TReal* const rss, int* const s,
	  const TCrit<TReal>& c)
  {
    return detail::dca(m, size, mark, nbest, v, ay, lday, index, rss, s, c);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int dca(const int size, const int mark, const int nbest,
          const int* const v, const TReal* const rz, const int ldrz,
	  int* const index, TReal* const rss, int* const s,
	  const TCrit<TReal>& c)
  {
    return detail::dca(size, mark, nbest, v, rz, ldrz, index, rss, s, c);
  }


  template<typename TReal>
  int dca(const int m, const int size, const int mark, const int nbest,
          const int* const v, const TReal* const ay, const int lday,
	  int* const index, TReal* const rss, int* const s)
  {
    return detail::dca(m, size, mark, nbest, v, ay, lday, index, rss, s);
  }


  template<typename TReal>
  int dca(const int size, const int mark, const int nbest,
          const int* const v, const TReal* const rz, const int ldrz,
	  int* const index, TReal* const rss, int* const s)
  {
    return detail::dca(size, mark, nbest, v, rz, ldrz, index, rss, s);
  }



}  // namespace subset
}  // namespace mcs



#endif

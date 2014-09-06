#ifndef MCS_SUBSET_DETAIL_ALGO_CC
#define MCS_SUBSET_DETAIL_ALGO_CC


#include "mcs/subset/detail/algo.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<typename InputIterator,
           typename RandomIterator,
           typename OutputIterator>
  void
  Algo::permute_n(InputIterator first, int n,
                  RandomIterator in,
                  OutputIterator out)
  {
    for (int i = 0; i < n; ++i)
    {
      *(out++) = in[*(first++)];
    }
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

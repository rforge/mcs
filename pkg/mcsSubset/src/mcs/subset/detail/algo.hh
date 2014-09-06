#ifndef MCS_SUBSET_DETAIL_ALGO_HH
#define MCS_SUBSET_DETAIL_ALGO_HH



namespace mcs    {
namespace subset {
namespace detail {



  struct Algo
  {



    template<typename InputIterator,
             typename RandomIterator,
             typename OutputIterator>
    static
    void
    permute_n(InputIterator first, int n,
              RandomIterator in,
              OutputIterator out);


  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/algo.cc"
#endif

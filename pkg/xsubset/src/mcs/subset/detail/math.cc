#ifndef MCS_SUBSET_DETAIL_MATH_CC
#define MCS_SUBSET_DETAIL_MATH_CC


#include "mcs/subset/detail/math.hh"

#include <limits>  // std::numeric_limits



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  const TReal Math<TReal>::MAX = std::numeric_limits<TReal>::max();


  template<typename TReal>
  TReal
  Math<TReal>::sign(TReal dx)
  {
    return (dx > 0) - (dx < 0);
  }


  template<typename TReal>
  TReal
  Math<TReal>::sqr(TReal val)
  {
    return val * val;
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

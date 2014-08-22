#ifndef MCS_SUBSET_DETAIL_MATH_HH
#define MCS_SUBSET_DETAIL_MATH_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  struct Math
  {

    static const TReal MAX;

    static
    TReal
    sign(TReal dx);

    static
    TReal
    sqr(TReal val);

  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/math.cc"
#endif

#ifndef MCS_SUBSET_DETAIL_LOG_LIK_CC
#define MCS_SUBSET_DETAIL_LOG_LIK_CC


#include "mcs/subset/detail/log_lik.cc"

#include <cmath>  // std::log



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  const TReal LogLik<TReal>::LOG_2PI_ = std::log(2 * Math<TReal>::PI);


  template<typename TReal>
  LogLik<TReal>::LogLik(const int nobs) :
    nobs_   {nobs},
    logNobs_{std::log(nobs)}
  {
  }


  template<typename TReal>
  TReal
  LogLik<TReal>::value(const TReal rss) const
  {
    return -0.5 * nobs_ * (LOG_2PI_ - logNobs_ + std::log(rss) + 1);
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs


#endif

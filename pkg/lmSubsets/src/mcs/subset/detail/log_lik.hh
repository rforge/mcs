#ifndef MCS_SUBSET_DETAIL_LOG_LIK_HH
#define MCS_SUBSET_DETAIL_LOG_LIK_HH



namespace mcs    {
namespace subset {
namespace detail {



  template<typename TReal>
  class LogLik
  {

  private:

    static const TReal LOG_2PI_;


  private:

    const int   nobs_;
    const TReal logNobs_;


  public:

    LogLik(int nobs);

    TReal
    value(TReal rss) const;

  };



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#include "mcs/subset/detail/log_lik.cc"
#endif

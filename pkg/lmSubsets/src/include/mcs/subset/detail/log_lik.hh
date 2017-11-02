#ifndef MCS_SUBSET_DETAIL_LOG_LIK_HH
#define MCS_SUBSET_DETAIL_LOG_LIK_HH



#define _USE_MATH_DEFINES  // M_PI
#include <cmath>  // std::log, M_PI



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class log_lik
{

private:

    static constexpr Scalar LOG_2PI_ = std::log(Scalar(2.0) * Scalar(M_PI));



private:

    Scalar nobs_half_;

    Scalar log_nobs_;



public:

    constexpr
    log_lik(const int nobs) noexcept :
        nobs_half_(Scalar(0.5) * Scalar(nobs)),
        log_nobs_(std::log(nobs))
    {
    }



    constexpr Scalar
    operator ()(const Scalar rss) const noexcept
    {
        const Scalar log_rss = std::log(rss);

        return -nobs_half_ * (LOG_2PI_ - log_nobs_ + log_rss + Scalar(1.0));
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif

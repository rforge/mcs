#ifndef MCS_SUBSET_DETAIL_AIC_HH
#define MCS_SUBSET_DETAIL_AIC_HH



#include "mcs/subset/detail/log_lik.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class aic
{

    using log_lik = detail::log_lik<Scalar>;



private:

    Scalar k_;

    log_lik ll_;



public:

    constexpr
    aic(
        const Scalar k,
        const int nobs
    ) noexcept :
        k_(k),
        ll_(nobs)
    {
    }


    constexpr Scalar
    operator ()(
        const int size,
        const Scalar rss
    ) const noexcept
    {
        // size + 1  to account for sd as estimated parameter
        const int npar = size + 1;

        return Scalar(-2.0) * ll_(rss) + k_ * Scalar(npar);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs


#endif

#ifndef MCS_SUBSET_ABBA_HH
#define MCS_SUBSET_ABBA_HH



#include <utility>  // std::pair
#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/core/matrix.hh"

#include "mcs/subset/table.hh"

#include "mcs/subset/detail/abba.hh"
#include "mcs/subset/detail/dca_preo.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {



template<typename Scalar,
         typename CostFunc>
std::pair<table_best<Scalar>, int>
abba_best(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const CostFunc& cost_func,
    const Scalar tau,
    const int nbest,
    const int prad
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = detail::dca_state_best<Scalar, CostFunc, preo_radius>;

    dca_state state(ay_mat, mark, cost_func, nbest, preo_radius(prad));
    const int nodes = detail::abba_best<Scalar, dca_state>(state, tau);

    return { state.table(), nodes };
}



template<typename Scalar>
std::pair<table_all<Scalar>, int>
abba_all(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    gsl::span<const Scalar> tau,
    const int nbest,
    const int prad
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = detail::dca_state_all<Scalar, preo_radius>;

    dca_state state(ay_mat, mark, nbest, preo_radius(prad));
    const int nodes = detail::abba_all<Scalar, dca_state>(
        state, {tau.begin(), tau.end()});

    return { state.table(), nodes };
}



}  // namespace subset
}  // namespace mcs



#endif

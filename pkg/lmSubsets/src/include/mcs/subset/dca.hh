#ifndef MCS_SUBSET_DCA_HH
#define MCS_SUBSET_DCA_HH



#include <utility>  // std::pair
#include <vector>



#include "mcs/core/matrix.hh"

#include "mcs/subset/table.hh"

#include "mcs/subset/detail/dca.hh"
#include "mcs/subset/detail/dca_preo.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {



template<typename Scalar,
         typename CostFunc>
std::pair<table_best<Scalar>, int>
dca_best(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const CostFunc& cost_func,
    const int nbest,
    const int prad = 0
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = dca_state_best<Scalar, CostFunc, preo_radius>;

    dca_state state(ay_mat, mark, cost_func, nbest, preo_radius(prad));
    const int nodes = detail::dca_best<Scalar, dca_state>(state);

    return { state.table(), nodes };
}



template<typename Scalar>
std::pair<table_all<Scalar>, int>
dca_all(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const int nbest,
    const int prad = 0
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = dca_state_all<Scalar, preo_radius>;

    dca_state state(ay_mat, mark, nbest, preo_radius(prad));
    const int nodes = detail::dca_all<Scalar, dca_state>(state);

    return { state.table(), nodes };
}



}  // namespace subset
}  // namespace mcs



#endif

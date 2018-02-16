// Copyright 2018  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs/>).
//
// 'mcs' is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'mcs' is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'mcs'.  If not, see <http://www.gnu.org/licenses/>.



#ifndef MCS_SUBSET_HBBA_HH
#define MCS_SUBSET_HBBA_HH



#include <utility>  // std::pair
#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/core/matrix.hh"

#include "mcs/subset/table.hh"

#include "mcs/subset/detail/hbba.hh"
#include "mcs/subset/detail/dca_preo.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {



template<typename Scalar,
         typename CostFunc>
std::pair<table_best<Scalar>, int>
hbba_best(
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
    const int nodes = detail::hbba_best<Scalar, dca_state>(state, tau);

    return { state.table(), nodes };
}



template<typename Scalar>
std::pair<table_all<Scalar>, int>
hbba_all(
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
    using dca_state = dca_state_all<Scalar, preo_radius>;

    dca_state state(ay_mat, mark, nbest, preo_radius(prad));
    const int nodes = detail::hbba_all<Scalar, dca_state>(
        state, {tau.begin(), tau.end()});

    return { state.table(), nodes };
}



}  // namespace subset
}  // namespace mcs



#endif

#ifndef MCS_SUBSET_DETAIL_SUBSET_HH
#define MCS_SUBSET_DETAIL_SUBSET_HH



#include <vector>



#include "mcs/detail/interrupt.hh"

#include "mcs/util/algo.hh"  // transform



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename DcaState>
void
subset_best(
    DcaState& state,
    const Scalar tau
) noexcept
{
    const Scalar approx = (tau - 1) * state.cost_inf();

    while (!state.is_final())
    {
        MCS_INTERRUPT_BREAK();

        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        const Scalar min_cost = state.min_cost() + approx;

        for (int j = k; j < n - 1; ++j)
        {
            if (tau * state.cost_bound(j) >= min_cost)
            {
                break;
            }

            state.drop_column(j);
        }
    }
}



template<typename Scalar,
         typename DcaState>
void
subset_all(
    DcaState& state,
    const std::vector<Scalar>& tau
) noexcept
{
    const auto approx = mcs::util::transform(
        tau,
        [rss_inf = state.rss_inf()](const Scalar tau) {
            return (tau - 1) * rss_inf;
        }
    );

    while (!state.is_final())
    {
        MCS_INTERRUPT_BREAK();

        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        for (int j_sup = n - 1; j_sup > k; --j_sup)
        {
            if (tau[j_sup - 1] * state.rss_bound() <
                (state.min_rss(j_sup) + approx[j_sup - 1]))
            {
                for (int j = k; j < j_sup; ++j)
                {
                    state.drop_column(j);
                }

                break;
            }
        }
    }
}



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
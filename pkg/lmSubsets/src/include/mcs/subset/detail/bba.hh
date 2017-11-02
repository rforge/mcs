#ifndef MCS_SUBSET_DETAIL_BBA_HH
#define MCS_SUBSET_DETAIL_BBA_HH



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename DcaState>
int
bba_best(DcaState& state) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        const Scalar min_cost = state.min_cost();

        for (int j = k; j < n - 1; ++j)
        {
            if (state.cost_bound(j) >= min_cost)
            {
                break;
            }

            state.drop_column(j);
        }

        ++node_cnt;
    }

    return node_cnt;
}



template<typename Scalar,
         typename DcaState>
int
bba_all(DcaState& state) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        for (int h = n - 1; h > k; --h)
        {
            if (state.rss_bound() < state.min_rss(h))
            {
                for (int j = k; j < h; ++j)
                {
                    state.drop_column(j);
                }

                break;
            }
        }

        ++node_cnt;
    }

    return node_cnt;
}



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

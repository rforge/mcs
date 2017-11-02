#ifndef MCS_SUBSET_DETAIL_DCA_HH
#define MCS_SUBSET_DETAIL_DCA_HH



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename DcaState>
int
dca_impl(DcaState& state) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        for (int j = k; j < n - 1; ++j)
        {
            state.drop_column(j);
        }

        ++node_cnt;
    }

    return node_cnt;
}



template<typename Scalar,
         typename DcaState>
int
dca_best(DcaState& state) noexcept
{
    return dca_impl<Scalar, DcaState>(state);
}



template<typename Scalar,
         typename DcaState>
int
dca_all(DcaState& state) noexcept
{
    return dca_impl<Scalar, DcaState>(state);
}



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

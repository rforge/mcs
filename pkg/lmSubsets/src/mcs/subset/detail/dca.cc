#ifndef MCS_SUBSET_DETAIL_DCA_CC
#define MCS_SUBSET_DETAIL_DCA_CC


#include "mcs/subset/detail/dca.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  dca(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table)
  {
    while (!state.isFinal())
      {
        state.nextNode();
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

        for (int j = k; j < n - 1; ++j)
          {
            state.dropColumn(j);
          }
      }

    return state.nodeCount();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

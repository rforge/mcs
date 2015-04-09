#ifndef MCS_SUBSET_DETAIL_BBA_CC
#define MCS_SUBSET_DETAIL_BBA_CC


#include "mcs/subset/detail/bba.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  bba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
      const TCriterion<TReal>& c)
  {
    while (!state.isDone())
      {
        state.nextNode();
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal max = table.maxBound();

        for (int j = k; j < n - 1; ++j)
          {
	    if (state.minBound(j + 1, c) >= max)
	      {
		break;
	      }

            state.dropColumn(j);
          }
      }

    table.sortSubsets();

    return state.nodeCount();
  }


  template<typename TReal>
  int
  bba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table)
  {
    while (!state.isDone())
      {
        state.nextNode();
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal min = state.minBound();

        for (int j = k; j < n - 1; ++j)
          {
	    if (min >= table.maxBound(j + 1))
	      {
		break;
	      }

            state.dropColumn(j);
          }
      }

    table.sortSubsets();

    return state.nodeCount();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

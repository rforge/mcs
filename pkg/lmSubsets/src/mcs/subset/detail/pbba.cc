#ifndef MCS_SUBSET_DETAIL_PBBA_CC
#define MCS_SUBSET_DETAIL_PBBA_CC


#include "mcs/subset/detail/pbba.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCriterion,
    template<typename R>
    class TPreorder
  >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
       const TCriterion<TReal>& crit, const TPreorder<TReal>& preo)
  {
    while (!state.isFinal())
      {
        state.nextNode(preo);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal best = table.bestValue();

        for (int j = k; j < n - 1; ++j)
          {
	    if (state.lowerBound(j + 1, crit) >= best)
	      {
		break;
	      }

            state.dropColumn(j);
          }
      }

    return state.nodeCount();
  }


  template<
    typename TReal,
    template<typename R>
    class TPreorder
   >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
       const TPreorder<TReal>& preo)
  {
    while (!state.isFinal())
      {
        state.nextNode(preo);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal bound = state.lowerBound();

        for (int j = k; j < n - 1; ++j)
          {
	    if (bound >= table.bestValue(j + 1))
	      {
		break;
	      }

            state.dropColumn(j);
          }
      }

    return state.nodeCount();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

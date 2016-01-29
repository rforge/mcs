#ifndef MCS_SUBSET_DETAIL_HPBBA_CC
#define MCS_SUBSET_DETAIL_HPBBA_CC


#include "mcs/subset/detail/hpbba.hh"



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
  hpbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
	const TCriterion<TReal>& crit, const TPreorder<TReal>& preo,
	const TReal tau)
  {
    while (!state.isFinal())
      {
        state.nextNode(preo);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal best = table.bestValue();

	for (int h = n - 1; h > k; --h)
	  {
	    if (tau * state.lowerBound(h, crit) < best)
	      {
		for (int j = k; j < h; ++j)
		  {
		    state.dropColumn(j);
		  }

		break;
	      }
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
  hpbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
	const TPreorder<TReal>& preo, const TReal* const tau)
  {
    while (!state.isFinal())
      {
        state.nextNode(preo);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	for (int h = n - 1; h > k; --h)
	  {
	    if (state.lowerBound(h, tau) < table.bestValue(h))
	      {
		for (int j = k; j < h; ++j)
		  {
		    state.dropColumn(j);
		  }

		break;
	      }
	  }
      }

    return state.nodeCount();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

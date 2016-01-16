#ifndef MCS_SUBSET_DETAIL_HBBA_CC
#define MCS_SUBSET_DETAIL_HBBA_CC


#include "mcs/subset/detail/hbba.hh"



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
  hbba(DcaState<TReal>& state, DcaTable<TReal,TCriterion>& table,
       const TCriterion<TReal>& crit, const TPreorder<TReal>& preo,
       const TReal tau)
  {
    while (!state.isFinal())
      {
        state.nextNode(preo);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

	const TReal max = table.maxBound();

	for (int h = n - 1; h > k; --h)
	  {
	    if (tau * state.minBound(h, crit) < max)
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
  hbba(DcaState<TReal>& state, DcaTable<TReal,Criteria::None>& table,
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
	    if (state.minBound(h, tau) < table.maxBound(h))
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

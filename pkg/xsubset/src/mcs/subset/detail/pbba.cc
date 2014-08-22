#ifndef MCS_SUBSET_DETAIL_PBBA_CC
#define MCS_SUBSET_DETAIL_PBBA_CC


#include "mcs/subset/detail/pbba.hh"

#include "mcs/subset/detail/algo.hh"
#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/criteria.hh"
#include "mcs/subset/detail/preorder.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCrit,
    template<typename R>
    class TPreorder
  >
  int
  pbba(const int m, const int size, const int mark, const int nbest,
       const int* const v, const TReal* const ay, const int lday,
       int* const index, TReal* const crit, int* const s,
       const TCrit<TReal>& c, const TPreorder<TReal>& p)
  {
    Preorder::Full<TReal> p1(size, p.pmin());

    DcaState<TReal>        state(m, size, mark, v, ay, lday, p1);
    DcaTable<TReal, TCrit> table(size, nbest, index, crit, s, c);

    return pbba(state, table, c);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit,
    template<typename R>
    class TPreorder
  >
  int
  pbba(const int size, const int mark, const int nbest,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const index, TReal* const crit, int* const s,
       const TCrit<TReal>& c, const TPreorder<TReal>& p)
  {
    Preorder::Full<TReal> p1(size, p.pmin());

    DcaState<TReal>        state(size, mark, v, rz, ldrz, p1);
    DcaTable<TReal, TCrit> table(size, nbest, index, crit, s, c);

    return pbba(state, table, c, p);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit,
    template<typename R>
    class TPreorder
  >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal, TCrit>& table,
       const TCrit<TReal>& c, const TPreorder<TReal>& p)
  {
    while (!state.isDone())
      {
        state.nextNode(p);
        state.reportSubleading(table);

        const int n = state.currentSize();
        const int k = state.currentMark();

        for (int j = k; j < n - 1; ++j)
          {
	    if (state.minBound(j + 1, c) >= table.maxBound())
	      {
		break;
	      }

            state.dropColumn(j);
          }
      }

    table.sortSubsets();

    return state.nodeCount();
  }


  template<
    typename TReal,
    template<typename R>
    class TPreorder
  >
  int
  pbba(const int m, const int size, const int mark, const int nbest,
       const int* const v, const TReal* const ay, const int lday,
       int* const index, TReal* const rss, int* const s,
       const TPreorder<TReal>& p)
  {
    Preorder::Full<TReal> p1(size, p.pmin());

    DcaState<TReal>                 state(m, size, mark, v, ay, lday, p1);
    DcaTable<TReal, Criteria::None> table(size, nbest, index, rss, s);

    return pbba(state, table, p);
  }


  template<
    typename TReal,
    template<typename R>
    class TPreorder
  >
  int
  pbba(const int size, const int mark, const int nbest,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const index, TReal* const rss, int* const s,
       const TPreorder<TReal>& p)
  {
    Preorder::Full<TReal> p1(size, p.pmin());

    DcaState<TReal>                 state(size, mark, v, rz, ldrz, p1);
    DcaTable<TReal, Criteria::None> table(size, nbest, index, rss, s);

    return pbba(state, table, p);
  }


  template<
    typename TReal,
    template<typename R>
    class TPreorder
   >
  int
  pbba(DcaState<TReal>& state, DcaTable<TReal, Criteria::None>& table,
       const TPreorder<TReal>& p)
  {
    while (!state.isDone())
      {
        state.nextNode(p);
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

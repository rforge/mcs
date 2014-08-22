#ifndef MCS_SUBSET_DETAIL_DCA_CC
#define MCS_SUBSET_DETAIL_DCA_CC


#include "mcs/subset/detail/dca.hh"

#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/criteria.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {
namespace detail {



  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  dca(const int m, const int size, const int mark, const int nbest,
      const int* const v, const TReal* const ay, const int lday,
      int* const index, TReal* const crit, int* const s,
      const TCrit<TReal>& c)
  {
    DcaState<TReal>        state(m, size, mark, v, ay, lday);
    DcaTable<TReal, TCrit> table(size, nbest, index, crit, s, c);

    return dca(state, table, c);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  dca(const int size, const int mark, const int nbest,
      const int* const v, const TReal* const rz, const int ldrz,
      int* const index, TReal* const crit, int* const s,
      const TCrit<TReal>& c)
  {
    DcaState<TReal>        state(size, mark, v, rz, ldrz);
    DcaTable<TReal, TCrit> table(size, nbest, index, crit, s, c);

    return dca(state, table, c);
  }


  template<
    typename TReal,
    template<typename R>
    class TCrit
  >
  int
  dca(DcaState<TReal>& state, DcaTable<TReal, TCrit>& table,
      const TCrit<TReal>& c)
  {
    while (!state.isDone())
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

    table.sortSubsets();

    return state.nodeCount();
  }


  template<typename TReal>
  int
  dca(const int m, const int size, const int mark, const int nbest,
      const int* const v, const TReal* const ay, const int lday,
      int* const index, TReal* const rss, int* const s)
  {
    DcaState<TReal>                 state(m, size, mark, v, ay, lday);
    DcaTable<TReal, Criteria::None> table(size, nbest, index, rss, s);

    return dca(state, table);
  }


  template<typename TReal>
  int
  dca(const int size, const int mark, const int nbest,
      const int* const v, const TReal* const rz, const int ldrz,
      int* const index, TReal* const rss, int* const s)
  {
    DcaState<TReal>                 state(size, mark, v, rz, ldrz);
    DcaTable<TReal, Criteria::None> table(size, nbest, index, rss, s);

    return dca(state, table);
  }


  template<typename TReal>
  int
  dca(DcaState<TReal>& state, DcaTable<TReal, Criteria::None>& table)
  {
    while (!state.isDone())
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

    table.sortSubsets();

    return state.nodeCount();
  }



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif

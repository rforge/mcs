#ifndef MCS_SUBSET_PBBA_HH
#define MCS_SUBSET_PBBA_HH


#include "mcs/subset/detail/algo.hh"
#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/preorder.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"

#include "mcs/subset/detail/pbba.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  pbba(const int m, const int size, const int mark, const int nbest,
       const int pmin, const int* const v, const TReal* const ay,
       const int lday, int* const sIndex, TReal* const sRss,
       TReal* const sVal, int* const sWhich, const TCriterion<TReal>& crit,
       int& nodes)
  {
    using namespace detail;

    Preorder::Full<TReal> preo(size, pmin);

    DcaState<TReal>            state(m, size, mark, v, ay, lday);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);

    nodes = pbba(state, table, crit, preo);

    table.sortIndex(sIndex);
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  pbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, TReal* const sVal,
       int* const sWhich, const TCriterion<TReal>& crit, int& nodes)
  {
    using namespace detail;

    Preorder::Full<TReal> preo(size, pmin);

    DcaState<TReal>            state(size, mark, v, rz, ldrz);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);

    nodes = pbba(state, table, crit, preo);

    table.sortIndex(sIndex);
  }


  template<typename TReal>
  void
  pbba(const int m, const int size, const int mark, const int nbest,
       const int pmin, const int* const v, const TReal* const ay,
       const int lday, int* const sIndex, TReal* const sRss, int* const sWhich,
       int& nodes)
  {
    using namespace detail;

    Preorder::Full<TReal> preo(size, pmin);

    DcaState<TReal>                 state(m, size, mark, v, ay, lday);
    DcaTable<TReal, Criteria::None> table(size, nbest, sIndex, sRss, sWhich);

    nodes = pbba(state, table, preo);

    table.sortIndex(sIndex);
  }


  template<typename TReal>
  void
  pbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, int* const sWhich, int& nodes)
  {
    using namespace detail;

    Preorder::Full<TReal> preo(size, pmin);

    DcaState<TReal>                 state(size, mark, v, rz, ldrz);
    DcaTable<TReal, Criteria::None> table(size, nbest, sIndex, sRss, sWhich);

    nodes = pbba(state, table, preo);

    table.sortIndex(sIndex);
  }



}  // namespace subset
}  // namespace mcs



#endif

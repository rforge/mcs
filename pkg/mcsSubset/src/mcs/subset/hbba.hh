#ifndef MCS_SUBSET_HBBA_HH
#define MCS_SUBSET_HBBA_HH


#include "mcs/subset/detail/algo.hh"
#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/preorder.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"

#include "mcs/subset/detail/hbba.hh"



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  hbba(const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const sIndex, TReal* const sRss, TReal* const sCrit,
       int* const sSize, int* const s, const TCriterion<TReal>& c, const TReal tau)
  {
    using namespace detail;

    Preorder::Full<TReal> p(size, pmin);

    DcaState<TReal>            state(m, size, mark, v, ay, lday);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);

    return hbba(state, table, c, p, tau);
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  int
  hbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, TReal* const sCrit,
       int* const sSize, int* const s, const TCriterion<TReal>& c, const TReal tau)
  {
    using namespace detail;

    Preorder::Full<TReal> p(size, pmin);

    DcaState<TReal>            state(size, mark, v, rz, ldrz);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss, sCrit, sSize, s, c);

    return hbba(state, table, c, p, tau);
  }


  template<typename TReal>
  int
  hbba(const int m, const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const ay, const int lday,
       int* const sIndex, TReal* const sRss, int* const s, const TReal* tau)
  {
    using namespace detail;

    Preorder::Full<TReal> p(size, pmin);

    DcaState<TReal>                state(m, size, mark, v, ay, lday);
    DcaTable<TReal,Criteria::None> table(size, mark, nbest, sIndex, sRss, s);

    return hbba(state, table, p, tau);
  }


  template<typename TReal>
  int
  hbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, int* const s, const TReal* tau)
  {
    using namespace detail;

    Preorder::Full<TReal> p (size, pmin);

    DcaState<TReal>                state(size, mark, v, rz, ldrz);
    DcaTable<TReal,Criteria::None> table(size, mark, nbest, sIndex, sRss, s);

    return hbba(state, table, p, tau);
  }



}  // namespace subset
}  // namespace mcs



#endif

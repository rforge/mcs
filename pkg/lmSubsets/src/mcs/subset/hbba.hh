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
<<<<<<< .mine
  void
  hbba(const int m, const int size, const int mark, const int nbest,
       const int* const v, const TReal* const ay, const int lday,
       int* const sIndex, TReal* const sRss, TReal* const sVal,
       int* const sWhich, const TCriterion<TReal>& crit, const TReal tau,
       int& nodes)
=======
  void
  hbba(const int m, const int size, const int mark, const int nbest,
       const int pmin, const int* const v, const TReal* const ay,
       const int lday, int* const sIndex, TReal* const sRss, TReal* const sVal,
       int* const sWhich, const TCriterion<TReal>& crit, const TReal tau,
       int& nodes)
>>>>>>> .r101
  {
    using namespace detail;

<<<<<<< .mine
=======
    Preorder::Full<TReal> preo(size, pmin);

>>>>>>> .r101
    DcaState<TReal>            state(m, size, mark, v, ay, lday);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);

<<<<<<< .mine
    nodes = hbba(state, table, crit, tau);

    table.sortIndex(sIndex);
=======
    nodes = hbba(state, table, crit, preo, tau);

    table.sortIndex(sIndex);
>>>>>>> .r101
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
<<<<<<< .mine
  void
  hbba(const int size, const int mark, const int nbest, const int* const v,
       const TReal* const rz, const int ldrz, int* const sIndex,
       TReal* const sRss, TReal* const sVal, int* const sWhich,
       const TCriterion<TReal>& crit, const TReal tau, int& nodes)
=======
  void
  hbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, TReal* const sVal,
       int* const sWhich, const TCriterion<TReal>& crit, const TReal tau,
       int& nodes)
>>>>>>> .r101
  {
    using namespace detail;

<<<<<<< .mine
=======
    Preorder::Full<TReal> preo(size, pmin);

>>>>>>> .r101
    DcaState<TReal>            state(size, mark, v, rz, ldrz);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);

<<<<<<< .mine
    nodes = hbba(state, table, crit, tau);

    table.sortIndex(sIndex);
=======
    nodes = hbba(state, table, crit, preo, tau);

    table.sortIndex(sIndex);
>>>>>>> .r101
  }


  template<typename TReal>
<<<<<<< .mine
  void
  hbba(const int m, const int size, const int mark, const int nbest,
=======
  void
  hbba(const int m, const int size, const int mark, const int nbest, const int pmin,
>>>>>>> .r101
       const int* const v, const TReal* const ay, const int lday,
<<<<<<< .mine
       int* const sIndex, TReal* const sRss, int* const sWhich,
       const TReal* tau, int& nodes)
=======
       int* const sIndex, TReal* const sRss, int* const sWhich, const TReal* tau,
       int& nodes)
>>>>>>> .r101
  {
    using namespace detail;

<<<<<<< .mine
=======
    Preorder::Full<TReal> preo(size, pmin);

>>>>>>> .r101
    DcaState<TReal>                state(m, size, mark, v, ay, lday);
    DcaTable<TReal,Criteria::None> table(size, nbest, sIndex, sRss, sWhich);

<<<<<<< .mine
    nodes = hbba(state, table, tau);

    table.sortIndex(sIndex);
=======
    nodes = hbba(state, table, preo, tau);

    table.sortIndex(sIndex);
>>>>>>> .r101
  }


  template<typename TReal>
<<<<<<< .mine
  void
  hbba(const int size, const int mark, const int nbest, const int* const v,
       const TReal* const rz, const int ldrz, int* const sIndex,
       TReal* const sRss, int* const sWhich, const TReal* tau, int& nodes)
=======
  void
  hbba(const int size, const int mark, const int nbest, const int pmin,
       const int* const v, const TReal* const rz, const int ldrz,
       int* const sIndex, TReal* const sRss, int* const sWhich,
       const TReal* tau, int& nodes)
>>>>>>> .r101
  {
    using namespace detail;

<<<<<<< .mine
=======
    Preorder::Full<TReal> preo(size, pmin);

>>>>>>> .r101
    DcaState<TReal>                state(size, mark, v, rz, ldrz);
    DcaTable<TReal,Criteria::None> table(size, nbest, sIndex, sRss, sWhich);

<<<<<<< .mine
    nodes = hbba(state, table, tau);

    table.sortIndex(sIndex);
=======
    nodes = hbba(state, table, preo, tau);

    table.sortIndex(sIndex);
>>>>>>> .r101
  }



}  // namespace subset
}  // namespace mcs



#endif

#ifndef MCS_SUBSET_XBBA_HH
#define MCS_SUBSET_XBBA_HH


#include "mcs/subset/detail/algo.hh"
#include "mcs/subset/detail/math.hh"
#include "mcs/subset/detail/givens.hh"
#include "mcs/subset/detail/lapack.hh"
#include "mcs/subset/detail/qrz.hh"
#include "mcs/subset/detail/preorder.hh"
#include "mcs/subset/detail/dca_table.hh"
#include "mcs/subset/detail/dca_state.hh"

#include "mcs/subset/detail/hpbba.hh"

#include <cstring>



namespace mcs    {
namespace subset {



  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  xbba(const char* const algo, const int m, const int size, const int mark,
       const int nbest, const int pmin, const int* const v,
       const TReal* const ay, const int lday, int* const sIndex,
       TReal* const sRss, TReal* const sVal, int* const sWhich,
       const TCriterion<TReal>& crit, const TReal tau, int& info, int& nodes)
  {
    using namespace detail;


    Preorder::Full<TReal> preo1(size);

    DcaState<TReal>            state(m, size, mark, v, ay, lday, preo1);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);


    const int pmax = size - mark - 1;

    if (std::strcmp(algo, "xbba1") == 0) {
      Preorder::Single1<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, crit, preo, tau);
    } else if (std::strcmp(algo, "xbba2") == 0) {
      Preorder::Single2<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, crit, preo, tau);
    } else {
      info = -1;
    }
  }


  template<
    typename TReal,
    template<typename R>
    class TCriterion
  >
  void
  xbba(const char* const algo, const int size, const int mark, const int nbest,
       const int pmin, const int* const v, const TReal* const rz,
       const int ldrz, int* const sIndex, TReal* const sRss, TReal* const sVal,
       int* const sWhich, const TCriterion<TReal>& crit, const TReal tau,
       int& info, int& nodes)
  {
    using namespace detail;


    Preorder::Full<TReal> preo1(size);

    DcaState<TReal>            state(size, mark, v, rz, ldrz, preo1);
    DcaTable<TReal,TCriterion> table(size, nbest, sIndex, sRss,
                                     sVal, sWhich, crit);


    const int pmax = size - mark - 1;

    if (std::strcmp(algo, "xbba1") == 0) {
      Preorder::Single1<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, crit, preo, tau);
    } else if (std::strcmp(algo, "xbba2") == 0) {
      Preorder::Single2<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, crit, preo, tau);
    } else {
      info = -1;
    }
  }


  template<typename TReal>
  void
  xbba(const char* const algo, const int m, const int size, const int mark,
       const int nbest, const int pmin, const int* const v,
       const TReal* const ay, const int lday, int* const sIndex,
       TReal* const sRss, int* const sWhich, const TReal* tau, int& info,
       int& nodes)
  {
    using namespace detail;


    Preorder::Full<TReal> preo1(size);

    DcaState<TReal>                 state(m, size, mark, v, ay, lday, preo1);
    DcaTable<TReal, Criteria::None> table(size, nbest, sIndex, sRss, sWhich);


    const int pmax = size - mark - 1;

    if (std::strcmp(algo, "xbba1") == 0) {
      Preorder::Single1<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, preo, tau);
    } else if (std::strcmp(algo, "xbba2") == 0) {
      Preorder::Single2<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, preo, tau);
    } else {
      info = -1;
    }
  }


  template<typename TReal>
  void
  xbba(const char* const algo, const int size, const int mark, const int nbest,
       const int pmin, const int* const v, const TReal* const rz,
       const int ldrz, int* const sIndex, TReal* const sRss, int* const sWhich,
       const TReal* tau, int& info, int& nodes)
  {
    using namespace detail;


    Preorder::Full<TReal> preo1(size);

    DcaState<TReal>                 state(size, mark, v, rz, ldrz, preo1);
    DcaTable<TReal, Criteria::None> table(size, nbest, sIndex, sRss, sWhich);


    const int pmax = size - mark - 1;

    if (std::strcmp(algo, "xbba1") == 0) {
      Preorder::Single1<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, preo, tau);
    } else if (std::strcmp(algo, "xbba2") == 0) {
      Preorder::Single2<TReal> preo(size, pmin, pmax);

      nodes = hpbba(state, table, preo, tau);
    } else {
      info = -1;
    }
  }



}  // namespace subset
}  // namespace mcs



#endif
